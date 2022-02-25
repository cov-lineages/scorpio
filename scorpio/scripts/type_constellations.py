#!/usr/bin/env python3

import csv
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import sys
import json
import re
import logging
import copy
import math
import multiprocessing as mp

if sys.version_info[0] < 3:
    raise Exception("Python 3 or a more recent version is required.")

global_aliases = {"spike": "s", "s": "spike",
                  "envelope": "e", "e": "envelope",
                  "membrane": "m", "m": "membrane",
                  "nucleocapsid": "n", "n": "nucleocapsid"}


def load_feature_coordinates(reference_json):
    """
    Loads a JSON file and extracts a dictionary of coordinates for features to
    be used to translate from amino acid into nucleotide coordinates

    nuc_pos is an integer which is 1-based start pos of codon
    """
    in_json = open(reference_json, 'r')
    json_dict = json.load(in_json, strict=False)

    if "genome" in json_dict:
        refseq = json_dict["genome"]
    else:
        sys.stderr.write("No reference sequence (key \"genome\") provided in JSON %s " % reference_json)
        sys.exit(1)

    features_dict = {}
    for feature in ["genes", "proteins"]:  #, "features"]:
        if feature in json_dict:
            for item in json_dict[feature]:
                name = item.lower()
                if "name" in json_dict[feature][item]:
                    name = json_dict[feature][item]["name"].lower()
                if name in features_dict or name in global_aliases and global_aliases[name] in features_dict:
                    continue

                if "coordinates" in json_dict[feature][item]:
                    if "from" in json_dict[feature][item]["coordinates"]:
                        start = int(json_dict[feature][item]["coordinates"]["from"])
                        end = int(json_dict[feature][item]["coordinates"]["to"])
                    elif "start" in json_dict[feature][item]["coordinates"]:
                        start = int(json_dict[feature][item]["coordinates"]["start"])
                        end = int(json_dict[feature][item]["coordinates"]["end"])

                    if "gene" in json_dict[feature][item]:
                        features_dict[name] = (start, end, json_dict[feature][item]["gene"])
                    else:
                        features_dict[name] = (start, end)
                    logging.info("Found reference feature %s with coordinates %s" % (name, features_dict[name]))
    if len(features_dict) == 0:
        sys.stderr.write("No features (keys \"genes\", \"proteins\" or \"features\" ) provided in JSON %s " %
                         reference_json)
        sys.exit(1)

    in_json.close()
    return refseq, features_dict


def resolve_ambiguous_cds(cds, aa_pos, features_dict):
    cds = cds.lower()

    if cds in features_dict:
        if len(features_dict[cds]) == 3:
            cds, aa_pos = resolve_ambiguous_cds(features_dict[cds][2], features_dict[cds][0]+aa_pos-1, features_dict)
        return cds, aa_pos

    if cds in global_aliases and global_aliases[cds] in features_dict:
        return global_aliases[cds], aa_pos

    if cds[0].isdigit():
        cds = "orf" + cds
        if cds in features_dict:
            return cds, aa_pos

    prefix = cds[:-2]
    potential = [key for key in features_dict if key.startswith(prefix)]
    count = 0
    for key in potential:
        aa_length = (features_dict[key][1] + 1 - features_dict[key][0])/3
        if count <= aa_pos < aa_length:
            return key, aa_pos
        aa_pos -= aa_length
        aa_pos = int(aa_pos)
        if aa_pos < 0:
            return cds, None
    return cds, None


def get_nuc_position_from_aa_description(cds, aa_pos, features_dict):
    """
    given a CDS (eg. S) and the number of an amino acid in it, get the
    1-based start position of that codon in ref coordinates

    nuc_pos is an integer which is 1-based start pos of codon
    """
    if aa_pos is None or cds not in features_dict.keys():
        sys.stderr.write("I don't know about cds: %s \n" % cds)
        sys.stderr.write("please use one of: %s" % ",".join(features_dict.keys()))
        sys.exit(1)

    cds_tuple = features_dict[cds]
    nuc_pos = cds_tuple[0] + ((aa_pos - 1) * 3)

    if nuc_pos > cds_tuple[1]:
        sys.stderr.write("invalid amino acid position for cds %s : %d" % (cds, aa_pos))
        sys.exit(1)

    return int(nuc_pos)


def variant_to_variant_record(l, refseq, features_dict, ignore_fails=False):
    """
    convert a variant in one of the following formats

    snp:T6954C
    nuc:T6954C
    del:11288:9
    aa:orf1ab:T1001I
    aa:orf1ab:T1001del
    aa:orf1ab:T1001 # this is for ambiguous AA change, NOT DELETION

    to a dict
    """
    #print("Parsing variant %s" %l)
    info = {}

    if "#" in l:
        l = l.split("#")[0].strip()
        if l == "":
            return info
    lsplit = l.split(":")

    if "+" in l:
        m = re.match(r'[aa:]*(?P<cds>\w+):(?P<pos>\d+)\+(?P<alt_allele>[a-zA-Z]+)', l)
        if not m:
            sys.stderr.write("Warning: couldn't parse the following string: %s\n" % l)
            if not ignore_fails:
                sys.exit(1)
        info = m.groupdict()
        info["type"] = "ins"
        info["ref_allele"] = ""
        if info["cds"] not in ["snp", "nuc"]:
            cds, pos = resolve_ambiguous_cds(info["cds"], info["pos"], features_dict)
            info["ref_start"] = get_nuc_position_from_aa_description(cds, pos, features_dict)
            info["name"] = "%s:%s%d+%s" % (info["cds"], info["ref_allele"], info["pos"], info["alt_allele"])
        else:
            info["ref_start"] = info["pos"]
            info["name"] = l
        logging.debug("Warning: found variant of type insertion, which will be ignored during typing")
    elif lsplit[0] in ["snp", "nuc"]:
        info = {"name": l, "type": "snp"}
        m = re.match(r'(?P<ref_allele>[ACGTUN]+)(?P<ref_start>\d+)(?P<alt_allele>[AGCTUN]*)', l[4:])
        if not m:
            sys.stderr.write("Warning: couldn't parse the following string: %s\n" % l)
            if not ignore_fails:
                sys.exit(1)
        info.update(m.groupdict())
        info["ref_start"] = int(info["ref_start"])
        ref_allele_check = refseq[info["ref_start"] - 1]

        if info["ref_allele"] != '?' and info["ref_allele"] != ref_allele_check:
            sys.stderr.write(
                "variants file says reference nucleotide at position %d is %s, but reference sequence has %s, "
                "context %s\n" % (info["ref_start"], info["ref_allele"], ref_allele_check, refseq[info["ref_start"] - 4:info["ref_start"] + 3]))
            if not ignore_fails:
                sys.exit(1)

    elif lsplit[0] == "del":
        length = int(lsplit[2])
        info = {"name": l, "type": lsplit[0], "ref_start": int(lsplit[1]), "length": length,
                "ref_allele": refseq[int(lsplit[1]) - 1:int(lsplit[1]) + length - 1], "space": "nuc"}

    else:
        m = re.match(r'[aa:]*(?P<cds>\w+):(?P<ref_allele>[a-zA-Z-*]+)(?P<aa_pos>\d+)(?P<alt_allele>[a-zA-Z-*]*)', l)
        if not m:
            sys.stderr.write("Warning: couldn't parse the following string: %s\n" % l)
            if not ignore_fails:
                sys.exit(1)
            return info

        info = m.groupdict()
        info["aa_pos"] = int(info["aa_pos"])
        info["name"] = "%s:%s%d%s" % (info["cds"], info["ref_allele"], info["aa_pos"], info["alt_allele"])
        info["type"] = "aa"

        cds, aa_pos = resolve_ambiguous_cds(info["cds"], info["aa_pos"], features_dict)
        ref_start = get_nuc_position_from_aa_description(cds, aa_pos, features_dict)

        ref_allele = Seq(refseq[ref_start - 1:ref_start - 1 + 3*len(info["ref_allele"])])
        ref_allele_check = ref_allele.translate()
        if info["ref_allele"] != '?' and info["ref_allele"] != ref_allele_check:
            sys.stderr.write("variants file says reference amino acid in CDS %s at position %d is %s, but "
                             "reference sequence has %s\n" % (cds, aa_pos, info["ref_allele"], ref_allele_check))
            if not ignore_fails:
                sys.exit(1)

        info["cds"] = cds
        info["aa_pos"] = aa_pos
        info["ref_start"] = ref_start
        if info["alt_allele"] in ['del', '-']:
            info["type"] = "del"
            info["space"] = "aa"
            info["length"] = len(info["ref_allele"])
            info["before"] = str(Seq(refseq[ref_start - 4:ref_start - 1]).translate())
            info["after"] = str(Seq(refseq[ref_start - 1 + 3*info["length"]:ref_start - 1 + 3*info["length"]+3]).translate())
        elif info["alt_allele"] == '':
            info["fuzzy"] = True
        else:
            info["fuzzy"] = False

    #print("Found variant %s of type %s" % (info["name"], info["type"]))
    return info


def parse_name_from_file(constellation_file):
    name = constellation_file.split('/')[-1]
    name = name.replace(".json", "").replace(".csv", "").replace(".txt", "")
    return name


def parse_json_in(refseq, features_dict, variants_file, constellation_names=None, include_ancestral=False, label=None, ignore_fails=False):
    """
    returns variant_list name and rules
    """
    variant_list = []
    name = None
    output_name = None
    rules = None
    mrca_lineage = ""
    incompatible_lineage_calls = ""
    parent_lineage = None
    lineage_name = None

    in_json = open(variants_file, 'r')
    json_dict = json.load(in_json, strict=False)

    if "type" in json_dict and json_dict["type"] in json_dict:
        if "mrca_lineage" in json_dict[json_dict["type"]]:
            m = re.match(r'[A-Z0-9.]*', json_dict[json_dict["type"]]["mrca_lineage"])
            if not m:
                sys.stderr.write("Warning: mrca_lineage %s not in acceptable format - ignoring\n" % json_dict[json_dict["type"]]["mrca_lineage"])
            else:
                mrca_lineage = json_dict[json_dict["type"]]["mrca_lineage"]
        if "incompatible_lineage_calls" in json_dict[json_dict["type"]]:
            incompatible_lineage_calls = "|".join(json_dict[json_dict["type"]]["incompatible_lineage_calls"])
        if "parent_lineage" in json_dict[json_dict["type"]]:
            parent_lineage = json_dict[json_dict["type"]]["parent_lineage"]
        if "lineage_name" in json_dict[json_dict["type"]]:
            lineage_name = json_dict[json_dict["type"]]["lineage_name"]

    if "name" in json_dict:
        name = json_dict["name"]
    elif "label" in json_dict:
        name = json_dict["label"]
    else:
        name = parse_name_from_file(variants_file)

    if label:
        if "type" in json_dict and json_dict["type"] in json_dict and label in json_dict[json_dict["type"]]:
            output_name = json_dict[json_dict["type"]][label]
    elif "label" in json_dict:
        output_name = json_dict["label"]
    elif name:
        output_name = name


    if not name:
        return variant_list, name, output_name, rules, mrca_lineage, incompatible_lineage_calls, parent_lineage, lineage_name
    #if constellation_names and name not in constellation_names and output_name not in constellation_names:
    #    return variant_list, name, output_name, rules, mrca_lineage, incompatible_lineage_calls, parent_lineage, lineage_name

    logging.debug("")
    logging.debug("Parsing constellation JSON file %s" % variants_file)

    if "sites" in json_dict:
        for site in json_dict["sites"]:
            record = variant_to_variant_record(site, refseq, features_dict, ignore_fails=ignore_fails)
            if record != {}:
                variant_list.append(record)
    if include_ancestral and "ancestral" in json_dict:
        for site in json_dict["ancestral"]:
            record = variant_to_variant_record(site, refseq, features_dict, ignore_fails=ignore_fails)
            if record != {}:
                variant_list.append(record)

    if "rules" in json_dict:
        if type(json_dict["rules"]) == dict and "default" in json_dict["rules"]:
            rules = json_dict["rules"]
        else:
            rules = {"default": json_dict["rules"]}

    in_json.close()
    sorted_variants = sorted(variant_list, key=lambda x: int(x["ref_start"]))
    #for var in sorted_variants:
    #    print(var)

    return sorted_variants, name, output_name, rules, mrca_lineage, incompatible_lineage_calls, parent_lineage, lineage_name


def parse_csv_in(refseq, features_dict, variants_file, constellation_names=None, ignore_fails=False):
    """
    returns variant_list and name
    """
    variant_list = []
    compulsory = []

    name = parse_name_from_file(variants_file)
    if constellation_names and name not in constellation_names:
        return variant_list, name, compulsory

    logging.debug("\n")
    logging.debug("Parsing constellation CSV file %s" % variants_file)

    csv_in = open("%s" % variants_file, 'r')
    reader = csv.DictReader(csv_in, delimiter=",")

    if "id" not in reader.fieldnames:
        csv_in.close()
        csv_in = open("%s" % variants_file, 'r', encoding="utf-8-sig")
        reader = csv.DictReader(csv_in, delimiter=",")

        if "id" not in reader.fieldnames:
            csv_in.close()
            logging.info("Warning: CSV headerline does not contain 'id': %s - ignoring" % reader.fieldnames)
            return variant_list, name, compulsory

    for row in reader:
        if ":" not in row["id"] and "gene" in reader.fieldnames:
            var = "%s:%s" % (row["gene"], row["id"])
        else:
            var = row["id"]
        record = variant_to_variant_record(var, refseq, features_dict, ignore_fails=ignore_fails)
        if record != {}:
            variant_list.append(record)
        if "compulsory" in reader.fieldnames and row["compulsory"] in ["True", True, "Y", "y", "Yes", "yes", "YES"]:
            compulsory.append(record["name"])

    csv_in.close()
    rules = None
    if len(compulsory) > 0:
        rules = {"default": {}}
        for var in compulsory:
            rules["default"][var] = "alt"
    sorted_variants = sorted(variant_list, key=lambda x: int(x["ref_start"]))
    return sorted_variants, name, rules


def parse_textfile_in(refseq, features_dict, variants_file, constellation_names=None, ignore_fails=False):
    """
    returns variant_list and name
    """
    variant_list = []

    name = parse_name_from_file(variants_file)
    if constellation_names and name not in constellation_names:
        return variant_list, name

    logging.debug("\n")
    logging.debug("Parsing constellation text file %s" % variants_file)

    with open("%s" % variants_file, "r") as f:
        for line in f:
            l = line.split("#")[0].strip()  # remove comments from the line
            if len(l) > 0:  # skip blank lines (or comment only lines)
                record = variant_to_variant_record(l, refseq, features_dict, ignore_fails=ignore_fails)
                if record != {}:
                    variant_list.append(record)
    sorted_variants = sorted(variant_list, key=lambda x: int(x["ref_start"]))

    return sorted_variants, name


def parse_variants_in(refseq, features_dict, variants_file, constellation_names=None, include_ancestral=False, label=None, ignore_fails=False):
    """
    read in a variants file and parse its contents and
    return something sensible.

    format of list_file is:
    snp:T6954C
    del:11288:9
    aa:orf1ab:T1001I
    aa:orf1ab:T1001del

    reference_json requires key "sites":[]

    csv_file requires columns "id" and "gene"

    returns variant_list which is a list of dicts of snps, aas and dels,
    one dict per variant. format of subdict varies by variant type
    """
    variant_list = []
    rule_dict = None
    output_name = None
    mrca_lineage = ""
    incompatible_lineage_calls = ""
    parent_lineage = None
    lineage_name = None

    if variants_file.endswith(".json"):
        variant_list, name, output_name, rule_dict, mrca_lineage, incompatible_lineage_calls, parent_lineage, lineage_name = parse_json_in(refseq, features_dict, variants_file, constellation_names, include_ancestral=include_ancestral,label=label,ignore_fails=ignore_fails)
    elif variants_file.endswith(".csv"):
        variant_list, name, rule_dict = parse_csv_in(refseq, features_dict, variants_file, constellation_names, ignore_fails=ignore_fails)
        output_name = name

    if len(variant_list) == 0 and not variants_file.endswith(".json"):
        variant_list, name = parse_textfile_in(refseq, features_dict, variants_file, constellation_names, ignore_fails=ignore_fails)
        output_name = name

    return name, output_name, variant_list, rule_dict, mrca_lineage, incompatible_lineage_calls, parent_lineage, lineage_name


def parse_mutations_in(mutations_file):
    logging.info("\n")
    logging.info("Parsing mutations file %s" % mutations_file)

    mutations_list = []
    with open("%s" % mutations_file, "r") as f:
        for line in f:
            l = line.split("#")[0].strip().split(',')[0]  # remove comments from the line and assume in first column
            if len(l) > 0:  # skip blank lines (or comment only lines)
                if l.startswith('id'):
                    continue
                mutations_list.append(l)
    logging.info("Found %d mutations" % len(mutations_list))
    return mutations_list


def parse_mutations(refseq, features_dict, mutations_list, ignore_fails=False):
    """
    Parse the mutations specified on command line and make a mutations constellation for them

    returns variant_list which is a list of dicts of snps, aas and dels,
    one dict per variant. format of subdict varies by variant type
    """
    variant_list = []
    problematic = []

    for mutation in mutations_list:
        record = variant_to_variant_record(mutation, refseq, features_dict, ignore_fails=ignore_fails)
        if record != {}:
            variant_list.append(record)
        else:
            problematic.append(mutation)

    if len(problematic) > 0:
        sys.stderr.write("The following mutations were not provided in an acceptable format: %s\n\n"
                         % ",".join(problematic))
        sys.exit(1)
    return variant_list


def call_variant_from_fasta(record_seq, var, ins_char="?", oth_char=None, codon=False):
    #print("Call variant for ", var)
    call = None
    query_allele = None

    if var["type"] == "ins":
        call = "oth"
        query_allele = ins_char
    elif var["type"] == "snp":
        query_allele = record_seq.upper()[var["ref_start"] - 1]
        #query_allele_minus = record_seq.upper()[var["ref_start"] - 2]
        #query_allele_plus = record_seq.upper()[var["ref_start"]]
        #print("Found", query_allele, query_allele_minus, query_allele_plus)
        #print(var["ref_allele"], query_allele == var["ref_allele"], var["alt_allele"], query_allele == var["alt_allele"])
        if query_allele == var["ref_allele"]:
            call = 'ref'
        elif query_allele == "N":
            call = 'ambig'
        elif query_allele == var["alt_allele"] or var["alt_allele"] == "":
            call = 'alt'
        else:
            call = 'oth'
        #print(call, query_allele)

    elif var["type"] == "aa":
        try:
            query = record_seq.upper()[var["ref_start"] - 1:var["ref_start"] - 1 + 3 * len(var["ref_allele"])]
            query_allele = query.translate(gap = "-")
            #query_allele_minus = record_seq.upper()[var["ref_start"] - 2:var["ref_start"] + 1].translate()
            #query_allele_plus = record_seq.upper()[var["ref_start"]:var["ref_start"] + 3].translate()
            #print("Found", query_allele, query_allele_minus, query_allele_plus)
            #print(var["ref_allele"], query_allele == var["ref_allele"], var["alt_allele"],query_allele == var["alt_allele"])
            if query_allele == var["ref_allele"]:
                call = 'ref'
            elif query_allele == "X":
                call = 'ambig'
            elif query_allele == var["alt_allele"] or var["fuzzy"]:
                call = 'alt'
            else:
                call = 'oth'
            #print(call, query_allele)
        except:
            #print("Except")
            call = 'oth'
        #print(call, query_allele)

    elif var["type"] == "del" and var["space"] == "aa":
        query_allele = str(record_seq.upper()[var["ref_start"] - 4:var["ref_start"] + 3*var["length"] + 2])
        query = query_allele.replace("-","")
        if len(query) % 3 != 0:
            query = query.replace("N","")
        if len(query) % 3 != 0:
            query = query_allele.replace("-","N")
        if len(query) % 3 != 0:
            logging.warning("Warning: while typing variant %s (before,ref,after) = (%s,%s,%s) found sequence with query allele %s treated as %s. Handling by adding Ns which will result in ambiguous calls" %(var["name"], var["before"], var["ref_allele"], var["after"], query_allele, query))
        query_allele = query
        while len(query_allele) % 3 != 0:
            query_allele += "N"
        query_allele = Seq(query_allele).translate()
        #print("call for del in aa space with before %s, ref %s, after %s, length %d and query_allele %s" %(var["before"], var["ref_allele"], var["after"], var["length"], query_allele))
        if query_allele == var["before"] + var["ref_allele"] + var["after"]:
            call = 'ref'
            query_allele = 0
        elif query_allele == var["before"] + var["after"]:
            call = 'alt'
            query_allele = var["length"]
        elif "X" in query_allele:
            call = 'ambig'
            query_allele = "X"
        else:
            call = 'oth'
            if not oth_char:
                query_allele = "X"
    elif var["type"] == "del" and var["space"] == "nuc":
        query_allele = record_seq.upper()[var["ref_start"] - 1:var["ref_start"] + var["length"] - 1]
        #print("call for del in nuc space with ref %s, length %d and query_allele %s" %(var["ref_allele"], var["length"], query_allele))

        if query_allele == var["ref_allele"]:
            call = 'ref'
            query_allele = 0
        elif query_allele == "-" * var["length"] or query_allele == "N" * var["length"]:
            call = 'alt'
            query_allele = max(int(var["length"] / 3), 1)
        elif "N" in query_allele:
            call = 'ambig'
            if not oth_char:
                query_allele = "X"
            else:
                query_allele = "N"
        else:
            call = 'oth'
            if not oth_char:
                query_allele = "X"
        #print(call, query_allele)

    if call == "oth" and var["type"] != "ins" and oth_char:
        query_allele = oth_char

    return call, query_allele


def var_follows_rules(call, rule):
    # rules allowed include "ref", "alt", "not ref", "not alt"
    rule_parts = rule.split()
    if len(rule_parts) > 1:
        rule_call = rule_parts[1]
    else:
        rule_call = rule
    if rule_parts[0] == "not":
        return call != rule_call
    else:
        return call == rule_call


def counts_follow_rules(counts, rules, key):
    # rules allowed include "max_ref", "min_alt", "min_snp_alt"
    is_rule_follower = True
    notes = []
    for rule in rules[key]:
        if ":" in rule:
            continue
        elif str(rule).startswith("min") or str(rule).startswith("max"):
            rule_parts = rule.split("_")
            if len(rule_parts) <= 1:
                continue
            elif len(rule_parts) == 2:
                if rule_parts[0] == "min" and counts[rule_parts[1]] < rules[key][rule]:
                    is_rule_follower = False
                elif rule_parts[0] == "max" and counts[rule_parts[1]] > rules[key][rule]:
                    is_rule_follower = False
                else:
                    counts["rules"][key] += 1
            elif len(rule_parts) == 3:
                part = None
                if rule_parts[1] in ["substitution", "snp"]:
                    part = "substitution"
                elif rule_parts[1] in ["indel"]:
                    part = "indel"
                if not part:
                    is_rule_follower = False
                elif rule_parts[0] == "min" and counts[part][rule_parts[2]] < rules[key][rule]:
                    is_rule_follower = False
                    notes.append("%s_%s_count=%i is less than %i" % (part, rule_parts[2], counts[part][rule_parts[2]], rules[key][rule]))
                elif rule_parts[0] == "max" and counts[part][rule_parts[2]] > rules[key][rule]:
                    is_rule_follower = False
                    notes.append("%s_%s_count=%i is more than %i" % (part, rule_parts[2], counts[part][rule_parts[2]], rules[key][rule]))
                else:
                    counts["rules"][key] += 1
        else:
            logging.warning("Warning: Ignoring rule %s:%s" % (rule, str(rules[key][rule])))
    return is_rule_follower, ";".join(notes)


def count_and_classify(record_seq, variant_list, rules):
    assert rules is not None
    counts = {'ref': 0, 'alt': 0, 'ambig': 0, 'oth': 0, 'rules': {},
              'substitution': {'ref': 0, 'alt': 0, 'ambig': 0, 'oth': 0},
              'indel': {'ref': 0, 'alt': 0, 'ambig': 0, 'oth': 0}}
    is_rule_follower_dict = {}
    for key in rules:
        is_rule_follower_dict[key] = True
        counts["rules"][key] = 0

    for var in variant_list:
        call, query_allele = call_variant_from_fasta(record_seq, var)
        #print(var, call, query_allele)
        counts[call] += 1
        if var['type'] in ["aa", "snp"]:
            counts["substitution"][call] += 1
        elif var['type'] in ["ins", "del"]:
            counts["indel"][call] += 1
        for key in rules:
            if var["name"] in rules[key]:
                if var_follows_rules(call, rules[key][var["name"]]):
                    counts['rules'][key] += 1
                elif is_rule_follower_dict[key]:
                    is_rule_follower_dict[key] = False

    counts['support'] = round(counts['alt']/float(counts['alt'] + counts['ref'] + counts['ambig'] + counts['oth']),4)
    counts['conflict'] = round(counts['ref'] /float(counts['alt'] + counts['ref'] + counts['ambig'] + counts['oth']),4)

    for key in rules:
        if not is_rule_follower_dict[key]:
            continue
        else:
            call, note = counts_follow_rules(counts, rules, key)
            if call:
                counts["rules"] = counts["rules"][key]
                call = key
                return counts, call, note
    counts["rules"] = counts["rules"]["default"]
    return counts, False, ""


def generate_barcode(record_seq, variant_list, ref_char=None, ins_char="?", oth_char="X",constellation_count_dict=None):
    barcode_list = []
    counts = {'ref': 0, 'alt': 0, 'ambig': 0, 'oth': 0}
    sorted_alt_sites = []

    for var in variant_list:
        call, query_allele = call_variant_from_fasta(record_seq, var, ins_char, oth_char)
        # print(var, call, query_allele)
        counts[call] += 1
        if constellation_count_dict and "constellations" in var:
            for constellation in var["constellations"]:
                constellation_count_dict[constellation][call] += 1
            if call == "alt":
                sorted_alt_sites.append(var["constellations"])

        if ref_char is not None and call == 'ref':
            barcode_list.append(str(ref_char))
        else:
            barcode_list.append(str(query_allele))

    counts['support'] = round(counts['alt'] / float(counts['alt'] + counts['ref'] + counts['ambig'] + counts['oth']), 4)
    counts['conflict'] = round(counts['ref'] / float(counts['alt'] + counts['ref'] + counts['ambig'] + counts['oth']), 4)

    return barcode_list, counts, constellation_count_dict, sorted_alt_sites


def load_constellations(list_constellation_files, constellation_names, reference_seq, features_dict, label,
                        include_ancestral=True, rules_required=False, ignore_fails=False):
    constellation_dict = {}
    name_dict = {}
    rule_dict = {}
    mrca_lineage_dict = {}
    incompatible_dict = {}
    parent_lineage_dict = {}
    lineage_name_dict = {}
    for constellation_file in list_constellation_files:
        constellation, output_name, variants, rules, mrca_lineage, \
        incompatible_lineage_calls, parent_lineage, lineage_name = \
            parse_variants_in(reference_seq,
                              features_dict,
                              constellation_file,
                              None,
                              include_ancestral=include_ancestral,
                              label=label,
                              ignore_fails=ignore_fails)
        if not constellation:
            continue
        #if constellation_names and constellation not in constellation_names and output_name not in constellation_names:
        #    continue
        #else:
        name_dict[constellation] = output_name
        if rules_required and not rules:
            logging.warning("Warning: No rules provided to classify %s - ignoring" % constellation)
            continue
        else:
            rule_dict[constellation] = rules
        if len(variants) > 0:
            constellation_dict[constellation] = variants
            mrca_lineage_dict[constellation] = mrca_lineage
            incompatible_dict[constellation] = incompatible_lineage_calls
            if parent_lineage:
                parent_lineage_dict[constellation] = parent_lineage
            if lineage_name:
                lineage_name_dict[lineage_name] = constellation
            if constellation_names and constellation not in constellation_names and output_name not in constellation_names:
                continue
            logging.info("\n")
            logging.info("Found file %s for constellation %s containing %i defining mutations" % (
                constellation_file, constellation, len([v["name"] for v in variants])))
            if rules_required:
                logging.info("Rules %s" % rule_dict[constellation])
        else:
            logging.warning("Warning: %s is not a valid constellation file - ignoring" % constellation_file)

    # If we specified only a subset of constellations on the commandline, restrict to those, handling the case where a
    # constellation has a parent definition, but the parent is not included in the list
    # if constellation_names:
    #     for constellation_name in constellation_names:
    #         current_constellation = constellation_name
    #         updated = False
    #         while current_constellation in parent_lineage_dict:
    #             parent = name_dict[lineage_name_dict[parent_lineage_dict[constellation_name]]]
    #             if parent not in constellation_names:
    #                 logging.info("\n")
    #                 logging.info("Add variants for parent %s to constellation %s" % (parent, constellation_name))
    #                 constellation_dict[constellation_name] = add_parent_variants(constellation_dict[constellation_name],
    #                                                                              constellation_dict[parent])
    #                 updated = True
    #                 if rules_required:
    #                     logging.info("Add rules for parent %s to constellation %s" % (parent, constellation_name))
    #                     rule_dict[constellation_name] = add_parent_rules(rule_dict[constellation_name], rule_dict[parent])
    #             current_constellation = parent
    #         if rules_required and updated:
    #             logging.info("Updated rules for %s: %s" % (constellation_name, rule_dict[constellation_name]))
    #
    #     names_to_ignore = [name for name in constellation_dict.keys() if name not in constellation_names and name_dict[name] not in constellation_names]
    #     for name in names_to_ignore:
    #         del constellation_dict[name]
    #     for name in constellation_names:
    #         if (name in parent_lineage_dict and parent_lineage_dict[name] in lineage_name_dict \
    #             and lineage_name_dict[parent_lineage_dict[name]] in names_to_ignore) \
    #                 or (name in parent_lineage_dict and parent_lineage_dict[name] in names_to_ignore):
    #             del parent_lineage_dict[name]

    return constellation_dict, name_dict, rule_dict, mrca_lineage_dict, incompatible_dict, parent_lineage_dict, \
           lineage_name_dict


def parse_mutations_list(mutations_list, reference_seq, features_dict):
    new_mutations_list = []
    for entry in mutations_list:
        if '.' in entry:
            new_mutations_list.extend(parse_mutations_in(entry)) # this is a file
        else:
            new_mutations_list.append(entry)
    mutations_list = new_mutations_list
    mutation_variants = parse_mutations(reference_seq, features_dict, mutations_list)
    return mutations_list, mutation_variants


def combine_constellations(constellation_dict):
    variant_dict = {}
    constellation_count_dict = {}
    for constellation in constellation_dict:
        constellation_count_dict[constellation] = {"total": len(constellation_dict[constellation]), 'ref': 0, 'alt': 0,
                                                   'ambig': 0, 'oth': 0}
        for variant in constellation_dict[constellation]:
            if variant["name"] in variant_dict:
                variant_dict[variant["name"]]["constellations"].append(constellation)
            else:
                variant_dict[variant["name"]] = variant
                variant_dict[variant["name"]]["constellations"] = [constellation]
    sorted_variants = sorted(variant_dict.values(), key=lambda x: int(x["ref_start"]))
    return {"union": sorted_variants}, constellation_count_dict


def combine_constellations_by_name(constellation_dict, name_dict):
    new_constellation_dict = {}
    for constellation in constellation_dict:
        constellation_name = name_dict[constellation]
        if not constellation_name:
            continue
        elif constellation_name == "None":
            constellation_name = constellation
        if constellation_name not in new_constellation_dict:
            new_constellation_dict[constellation_name] = constellation_dict[constellation]
        else:
            for variant in constellation_dict[constellation]:
                if variant not in new_constellation_dict[constellation_name]:
                    new_constellation_dict[constellation_name].append(variant)
    return new_constellation_dict

def add_parent_variants(constellation_variants, parent_variants):
    for variant in parent_variants:
        if variant not in constellation_variants:
              constellation_variants.append(variant)
    return constellation_variants

def add_parent_rules(constellation_rules, parent_rules):
    for key in parent_rules:
        if key in constellation_rules:
            for rule in parent_rules[key]:
                if rule.startswith("min") or rule.startswith("max") and rule in constellation_rules[key]:
                    constellation_rules[key][rule] += parent_rules[key][rule]
                else:
                    constellation_rules[key][rule] = parent_rules[key][rule]
    return constellation_rules


def get_number_switches(barcode_list, ref_char="-", ambig_char = "X"):
    previous = None
    current = None
    number_switches = 0
    len_non_ambig_barcode = 0
    for letter in barcode_list:
        if letter == ref_char:
            current = 1
        elif letter == ambig_char:
            continue
        else:
            current = 0

        len_non_ambig_barcode += 1

        if previous == None:
            previous = current
            continue

        if current != previous:
            number_switches += 1
            previous = current

    if number_switches > 1:
        number_switches -= 1

    return number_switches,len_non_ambig_barcode

#def prob_number_errors_or_fewer(counts):
#    # H0: missing is due to errors
#    # H1: missing is due to recombination
#    prob = counts["substitution"]["alt"]^prob_snp_error + prob_del_error * counts["indel"]["alt"]

#def prob_number_switches_or_fewer(number_switches, len_non_ambig_barcode):
#    # H0: sample is random draws from mixture
#    # H1: sample is a due to recombination
#    return


#def get_interspersion(top_scoring, counts, barcode_list, ref_char="-", ambig_char = "X"):
#    number_switches, len_non_ambig_barcode = get_number_switches(barcode_list, ref_char="-", ambig_char = "X"
#    # H0: missing is due to sequencing errors
#    # H1: sample is a random mixture of
#    return


def file_writer(out_csv, list_of_constellations, q):
    """listens for messages on the q, writes to files. """
    handle_dict = {}
    if len(list_of_constellations) == 1:
        handle_dict[list_of_constellations[0]] = open(out_csv, "w")
    else:
        handle_dict["summary"] = open(out_csv, "w")
        for constellation_name in list_of_constellations:
            clean_name = re.sub("[^a-zA-Z0-9_\-.]", "_", constellation_name)
            handle_dict[constellation_name] = open("%s.%s.csv" % (out_csv.replace(".csv", ""), clean_name), "w")

    while 1:
        handle, message = q.get()
        if handle not in handle_dict:
            print("Failed to write to unopen file %s:\n%s" %(handle, message))
        if message == 'kill':
            break
        handle_dict[handle].write(message)

    for handle in handle_dict:
        handle_dict[handle].close()


def type_record(record, reference_seq, constellation_names, constellation_dict, name_dict, constellation_count_dict, output_counts,
                append_genotypes, combination, ref_char, q):
    if len(record.seq) != len(reference_seq):
        sys.stderr.write("The fasta record for %s isn't as long as the reference, is this fasta file aligned?\n"
                         % record.id)
        sys.exit(1)

    out_list = [record.id]
    for constellation in constellation_dict:
        logging.debug("Consider constellation %s" % constellation)
        barcode_list, counts, sample_constellation_count_dict, sorted_alt_sites = generate_barcode(record.seq,
                                                                                                   constellation_dict[constellation],
                                                                                                   ref_char,
                                                                                                   constellation_count_dict=copy.deepcopy(constellation_count_dict))
        if output_counts or append_genotypes:
            columns = [record.id]
            if output_counts:
                columns.append("%i,%i,%i,%i,%f,%f" % (counts['ref'], counts['alt'],
                                                      counts['ambig'], counts['oth'], counts['support'],
                                                      counts['conflict']))
            if append_genotypes:
                columns.extend(barcode_list)
            if combination:
                scores = {}
                for candidate in sample_constellation_count_dict:
                    if sample_constellation_count_dict[candidate]["alt"] > 0:
                        summary = "%s:%i|%i|%i|%i" % (
                            candidate, sample_constellation_count_dict[candidate]["ref"],
                            sample_constellation_count_dict[candidate]["alt"],
                            sample_constellation_count_dict[candidate]["ambig"],
                            sample_constellation_count_dict[candidate]["oth"])
                        score = float(sample_constellation_count_dict[candidate]["alt"]) / \
                                sample_constellation_count_dict[candidate]["total"]
                        scores[score] = summary
                sorted_scores = sorted(scores, key=lambda x: float(x), reverse=True)
                columns.append("; ".join([scores[score] for score in sorted_scores]))
            res = "%s\n" % ','.join(columns)
            q.put((constellation, res))

        out_list.append(''.join(barcode_list))

    res = "%s\n" % ",".join(out_list)
    if len(constellation_dict) > 1 or not (output_counts or append_genotypes):
        q.put(("summary", res))
    return res

def restrict_to_named_constellations(constellation_names, constellation_dict, name_dict):
    to_remove = []
    if constellation_names:
        for name in constellation_dict:
            if name in constellation_names or (name in name_dict and name_dict[name] in constellation_names):
                continue
            to_remove.append(name)
    for name in to_remove:
        del constellation_dict[name]
    return constellation_dict

def type_constellations(in_fasta, list_constellation_files, constellation_names, out_csv, reference_json, ref_char=None,
                        output_counts=False, label=None, append_genotypes=False, mutations_list=None, dry_run=False,
                        combination=False, interspersion=False, threads=1):
    reference_seq, features_dict = load_feature_coordinates(reference_json)

    constellation_dict, name_dict, rule_dict, mrca_lineage_dict, \
    incompatible_dict, parent_lineage_dict, lineage_name_dict = load_constellations(list_constellation_files,
                                                                                    constellation_names,
                                                                                    reference_seq,
                                                                                    features_dict,
                                                                                    label,
                                                                                    include_ancestral=False,
                                                                                    rules_required=False)
    restrict_to_named_constellations(constellation_names, constellation_dict, name_dict)

    if combination:
        constellation_names = ["union"]
        constellation_dict, constellation_count_dict = combine_constellations(constellation_dict)
        name_dict["union"] = "union"
    else:
        constellation_count_dict = None

    if mutations_list:
        constellation_names.append("mutations")
        mutations_list, mutation_variants = parse_mutations_list(mutations_list, reference_seq, features_dict)
        if len(constellation_dict) == 1 and "mutations" not in constellation_dict:
            constellation = list(constellation_dict)[0]
            new_constellation = "%s+%s" % (constellation, '|'.join(mutations_list))
            constellation_dict[new_constellation] = constellation_dict[constellation] + mutation_variants
            del constellation_dict[constellation]
        else:
            constellation_dict["mutations"] = mutation_variants
            name_dict["mutations"] = "mutations"

    if dry_run:
        return

    logging.info("\n")
    logging.info("Update constellation dict")
    constellation_dict = combine_constellations_by_name(constellation_dict, name_dict)
    logging.debug(constellation_dict)
    logging.info("Have %i constellations to type: %s" % (len(constellation_dict), list(constellation_dict.keys())))
    if constellation_count_dict:
        logging.info("Have %i candidate constellations to collect counts: %s" % (len(constellation_count_dict),
                                                                                 list(constellation_count_dict.keys())))

    # setup manager
    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(threads)

    # put listener to work first
    list_constellation_names = constellation_dict.keys()
    if not output_counts and not append_genotypes:
        list_constellation_names = []
    if threads > 1:
        pool.apply_async(file_writer, (out_csv, list_constellation_names, q))

    if len(constellation_dict) > 1 or not (output_counts or append_genotypes):
        header = "query,%s\n" % ",".join(list_constellation_names)
        q.put(("summary", header))

    if output_counts or append_genotypes:
        for constellation in constellation_dict:
            columns = ["query"]
            if output_counts:
                columns.append("ref_count,alt_count,ambig_count,other_count,support,conflict")
            if append_genotypes:
                columns.extend([var["name"] for var in constellation_dict[constellation]])
            if combination:
                columns.append("notes")
            header = "%s\n" % ','.join(columns)
            q.put((constellation, header))

    # fire off workers
    jobs = []
    with open(in_fasta, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            job = pool.apply_async(type_record, (record, reference_seq, constellation_names, constellation_dict,
                                                 name_dict, constellation_count_dict, output_counts, append_genotypes,
                                                 combination, ref_char, q))
            jobs.append(job)

    # collect results from the workers through the pool result queue
    for job in jobs:
        job.get()

    # now we are done, kill the listener
    if len(list_constellation_names) > 1 or not (append_genotypes or output_counts):
        q.put(("summary", 'kill'))

    for constellation_name in list_constellation_names:
        q.put((constellation_name, "kill"))

    if threads == 1:
        pool.apply_async(file_writer, (out_csv, list_constellation_names, q))

    pool.close()
    pool.join()


def combine_counts_call_notes(counts1, call1, note1, counts2, call2, note2):
    counts = {'ref': 0, 'alt': 0, 'ambig': 0, 'oth': 0, 'rules': 0,
              'substitution': {'ref': 0, 'alt': 0, 'ambig': 0, 'oth': 0},
              'indel': {'ref': 0, 'alt': 0, 'ambig': 0, 'oth': 0},
              'support': 0, 'conflict': 0}
    for key in counts:
        if key in ["substitution", "indel"]:
            for subkey in counts[key]:
                counts[key][subkey] = counts1[key][subkey] + counts2[key][subkey]
        else:
            counts[key] = counts1[key] + counts2[key]
    counts['support'] = round(counts['alt'] / float(counts['alt'] + counts['ref'] + counts['ambig'] + counts['oth']), 4)
    counts['conflict'] = round(counts['ref'] / float(counts['alt'] + counts['ref'] + counts['ambig'] + counts['oth']), 4)
    if not call1 or not call2:
        call = False
    elif call1 == call2:
        call = call1
    elif call1 == "default" or call2 == "default":
        call = "default"
    else:
        call = call1

    note = note1
    if note != "" and note2 != "":
        note += ";" + note2
    else:
        note = note2
    return counts, call, note


def get_columns(list_incompatible, long, call_all, mutations_list, single_file):
    columns = ["query", "constellations", "mrca_lineage"]
    if list_incompatible:
        columns.append("incompatible_lineages")
    if long and not call_all:
        columns.extend(["ref_count", "alt_count", "ambig_count", "other_count", "rule_count", "support", "conflict"])
    columns.append("constellation_name")
    if single_file:
        columns.extend(["call", "note"])
    if mutations_list:
        columns.extend(mutations_list)
    return columns


def call_record(record, reference_seq, constellation_names, constellation_dict, name_dict, rule_dict, mrca_lineage_dict,
                incompatible_dict, parent_lineage_dict, lineage_name_dict, output_counts, call_all, long,
                list_incompatible, mutations_list, mutation_variants, interspersion, q):
    if len(record.seq) != len(reference_seq):
        sys.stderr.write("The fasta record for %s isn't as long as the reference, is this fasta file aligned?\n"
                         % record.id)
        sys.exit(1)

    list_constellation_names = list(set([name_dict[c] for c in constellation_dict if c in constellation_names]))
    single_file = len(list_constellation_names) == 1 and output_counts

    lineages = []
    names = []
    best_constellation = None
    best_support = 0
    best_conflict = 1
    best_counts = None
    best_call = False
    scores = {}
    children = {}
    for constellation in constellation_dict:
        constellation_name = name_dict[constellation]
        if not constellation_name or constellation_name not in constellation_names:
            continue
        logging.debug("Consider constellation %s" % constellation_name)
        parents = []

        counts, call, note = count_and_classify(record.seq,
                                                constellation_dict[constellation],
                                                rule_dict[constellation])
        current_constellation = constellation
        while current_constellation in parent_lineage_dict:
            logging.debug("Current constellation %s in parent dict" % current_constellation)
            current_constellation = name_dict[lineage_name_dict[parent_lineage_dict[current_constellation]]]
            parents.append(current_constellation)
            parent_counts, parent_call, parent_note = count_and_classify(record.seq,
                                                                         constellation_dict[current_constellation],
                                                                         rule_dict[current_constellation])
            counts, call, note = combine_counts_call_notes(counts, call, note, parent_counts, parent_call, parent_note)

        for parent in parents:
            if parent not in children:
                children[parent] = []
            children[parent].append(constellation)

        if call:
            logging.debug("Have call for %s" % constellation_name)
            if call_all:
                if call != "default":
                    called_constellation_name = "%s %s" % (call, constellation_name)
                else:
                    called_constellation_name = constellation_name
                lineages.append(called_constellation_name)
                names.append(constellation)
            elif constellation in children and best_constellation in children[constellation]:
                logging.debug("Ignore as parent of best constellation")
            elif (not best_constellation) \
                    or (counts['support'] > best_support) \
                    or (counts['support'] == best_support and counts['conflict'] < best_conflict) \
                    or (counts['support'] == best_support and counts['conflict'] == best_conflict
                                                          and counts['rules'] > best_counts["rules"]) \
                    or (best_constellation in parents):
                best_constellation = constellation
                logging.debug("Set best constellation %s" % best_constellation)
                best_support = counts['support']
                best_conflict = counts['conflict']
                best_counts = counts
                best_call = call

        if interspersion:
            if counts["alt"] > 1:
                summary = constellation
                score = counts["support"]
                scores[score] = summary

        if output_counts:
            res = "%s,%i,%i,%i,%i,%i,%f,%f,%s,%s,%s\n" % (
                    record.id, counts['ref'], counts['alt'], counts['ambig'],
                    counts['oth'], counts['rules'], counts['support'],
                    counts['conflict'], call, constellation, note)
            if not single_file:
                q.put((constellation_name, res))

    if not call_all and best_constellation:
        if best_call != "default":
            best_constellation_name = "%s %s" % (best_call, name_dict[best_constellation])
        else:
            best_constellation_name = name_dict[best_constellation]
        lineages.append(best_constellation_name)
        names.append(best_constellation)

    out_entries = [record.id, "|".join(lineages), "|".join([mrca_lineage_dict[n] for n in names])]

    if list_incompatible:
        out_entries.append("|".join([incompatible_dict[n] for n in names]))
    if long and best_counts is not None:
        out_entries.append("%i,%i,%i,%i,%i,%f,%f" % (best_counts['ref'],
                                                     best_counts['alt'], best_counts['ambig'],
                                                     best_counts['oth'], best_counts['rules'],
                                                     best_counts['support'], best_counts['conflict']))
    elif long and not call_all:
        out_entries.append(",,,,,,")
    out_entries.append("|".join(names))

    if single_file:
        out_entries.append("%s,%s" % (call, note))

    if mutations_list:
        barcode_list, counts, ignore, ignore2 = generate_barcode(record.seq, mutation_variants)
        out_entries.extend(barcode_list)

    res = "%s\n" % ",".join(out_entries)
    if single_file:
        q.put((constellation_name, res))
    else:
        q.put(("summary", res))

    return res


def classify_constellations(in_fasta, list_constellation_files, constellation_names, out_csv, reference_json,
                            output_counts=False, call_all=False, long=False, label=None, list_incompatible=False,
                            mutations_list=None, dry_run=False, interspersion=False, threads=1):
    reference_seq, features_dict = load_feature_coordinates(reference_json)

    constellation_dict, name_dict, rule_dict, mrca_lineage_dict, \
    incompatible_dict, parent_lineage_dict, lineage_name_dict = load_constellations(list_constellation_files,
                                                                                    constellation_names,
                                                                                    reference_seq,
                                                                                    features_dict,
                                                                                    label,
                                                                                    include_ancestral=True,
                                                                                    rules_required=True,
                                                                                    ignore_fails=True)

    if not constellation_names:
        constellation_names = [name_dict[c] for c in constellation_dict]
    logging.debug("parent_dict: %s" % parent_lineage_dict)
    logging.debug("lineage_name_dict: %s" % lineage_name_dict)

    mutation_variants = None
    if mutations_list:
        mutations_list, mutation_variants = parse_mutations_list(mutations_list, reference_seq, features_dict)
    if dry_run:
        return

    # setup manager
    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(threads)

    # put listener to work first
    list_constellation_names = list(set([name_dict[c] for c in constellation_dict if c in constellation_names]))
    if not output_counts:
        list_constellation_names = []

    if threads > 1:
        pool.apply_async(file_writer, (out_csv, list_constellation_names, q))

    # write summary file header
    single_file = len(list_constellation_names) == 1 and output_counts
    if single_file:
        long = True

    columns = get_columns(list_incompatible, long, call_all, mutations_list, single_file)
    header = "%s\n" % ",".join(columns)
    if len(list_constellation_names) == 1:
        q.put((list_constellation_names[0], header))
    else:
        q.put(("summary", header))

        header = "query,ref_count,alt_count,ambig_count,other_count,rule_count,support,conflict,call,constellation_name,note\n"
        for constellation_name in list_constellation_names:
            q.put((constellation_name, header))

    # fire off workers
    jobs = []
    with open(in_fasta, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            job = pool.apply_async(call_record, (record, reference_seq, constellation_names, constellation_dict,
                                                 name_dict, rule_dict, mrca_lineage_dict, incompatible_dict,
                                                 parent_lineage_dict, lineage_name_dict, output_counts, call_all,
                                                 long, list_incompatible, mutations_list, mutation_variants,
                                                 interspersion, q))
            jobs.append(job)

    # collect results from the workers through the pool result queue
    for job in jobs:
        job.get()

    #now we are done, kill the listener
    if not single_file:
        q.put(("summary", "kill"))
    for constellation_name in list_constellation_names:
        q.put((constellation_name, "kill"))

    if threads == 1:
        pool.apply_async(file_writer, (out_csv, list_constellation_names, q))

    pool.close()
    pool.join()


def list_constellations(list_constellation_files, constellation_names, reference_json, label=None):
    reference_seq, features_dict = load_feature_coordinates(reference_json)

    list_of_constellations = set()
    for constellation_file in list_constellation_files:
        constellation, output_name, variants, ignore, mrca_lineage, \
        incompatible_lineage_calls, parent_lineage, lineage_name = \
            parse_variants_in(reference_seq, features_dict, constellation_file, constellation_names, label=label)
        if not constellation:
            continue
        if constellation_names and constellation not in constellation_names and output_name not in constellation_names:
            continue
        if len(variants) > 0 and mrca_lineage:
            list_of_constellations.add(output_name)
            list_of_constellations.add(mrca_lineage)
        elif len(variants) > 0:
            list_of_constellations.add(output_name)
    print("\n".join(list_of_constellations))


def parse_args():
    parser = argparse.ArgumentParser(description="""Type an alignment at specific sites and classify with a barcode.""",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--fasta-in', dest='in_fasta', help='alignment to type, in fasta format')
    parser.add_argument('--variants-config', dest='variants_in', nargs='+', help="""
    List of config files containing variants to type. This can be a one-per-line-list in a text file, 
    a json with key "sites" or a csv with columns "id" and "gene"
    Format of allowed variants include:
            snp:T6954C
            nuc:T6954C
            del:11288:9
            aa:orf1ab:T1001I
            aa:orf1ab:T1001del
            aa:orf1ab:T1001-
            orf1ab:T1001I

    snp and del positions are 1-based nucleotide position relative to the alignment
    aa position is 1-based codon position relative to the cds.
    """)
    parser.add_argument('--reference_json', help='JSON file containing keys "genome" with reference sequence '
                                                 'and "proteins", "features" or "genes" with features of interest'
                                                 ' and their coordinates')
    parser.add_argument('--variants-out', dest='variants_out', help='csv file to write')
    parser.add_argument('--append-genotypes', dest='append_genotypes', action='store_true',
                        help='if invoked, write the genotype for each variant in the config file to the output')

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    type_constellations(args.in_fasta, args.variants_in, None, args.variants_out, args.reference_json, ref_char='-')


if __name__ == '__main__':
    main()
