#!/usr/bin/env python3

import csv
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import sys
import json
import re

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
                    print("Found reference feature %s with coordinates" % name, features_dict[name])
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


def variant_to_variant_record(l, refseq, features_dict):
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
    lsplit = l.split(":")
    info = {}

    if "+" in l:
        m = re.match(r'[aa:]*(?P<cds>\w+):(?P<pos>\d+)\+(?P<alt_allele>[a-zA-Z]+)', l)
        if not m:
            sys.stderr.write("Warning: couldn't parse the following string: %s - ignoring\n" % l)
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
        print("Warning: found variant of type insertion, which will be ignored during typing")
    elif lsplit[0] in ["snp", "nuc"]:
        info = {"name": l, "type": "snp"}
        m = re.match(r'(?P<ref_allele>[ACGTUN]+)(?P<ref_start>\d+)(?P<alt_allele>[AGCTUN]*)', l[4:])
        if not m:
            sys.stderr.write("Warning: couldn't parse the following string: %s - ignoring\n" % l)
            sys.exit(1)
        info.update(m.groupdict())
        info["ref_start"] = int(info["ref_start"])
        ref_allele_check = refseq[info["ref_start"] - 1]

        if info["ref_allele"] != '?' and info["ref_allele"] != ref_allele_check:
            sys.stderr.write(
                "variants file says reference nucleotide at position %d is %s, but reference sequence has %s, "
                "context %s\n" % (info["ref_start"], info["ref_allele"], ref_allele_check, refseq[info["ref_start"] - 4:info["ref_start"] + 3]))
            sys.exit(1)

    elif lsplit[0] == "del":
        length = int(lsplit[2])
        info = {"name": l, "type": lsplit[0], "ref_start": int(lsplit[1]), "length": length,
                "ref_allele": refseq[int(lsplit[1]) - 1:int(lsplit[1]) + length - 1], "space": "nuc"}

    else:
        m = re.match(r'[aa:]*(?P<cds>\w+):(?P<ref_allele>[a-zA-Z-*]+)(?P<aa_pos>\d+)(?P<alt_allele>[a-zA-Z-*]*)', l)
        if not m:
            sys.stderr.write("Warning: couldn't parse the following string: %s - ignoring\n" % l)
            # sys.exit(1)
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


def parse_json_in(refseq, features_dict, variants_file, constellation_names=None, include_ancestral=False, label=None):
    """
    returns variant_list name and rules
    """
    variant_list = []
    name = None
    rules = None
    mrca_lineage = ""

    in_json = open(variants_file, 'r')
    json_dict = json.load(in_json, strict=False)

    if "type" in json_dict and json_dict["type"] in json_dict and "mrca_lineage" in json_dict[json_dict["type"]]:
        mrca_lineage = json_dict[json_dict["type"]]["mrca_lineage"]

    if label:
        if "type" in json_dict and json_dict["type"] in json_dict and label in json_dict[json_dict["type"]]:
            name = json_dict[json_dict["type"]][label]
    elif "label" in json_dict:
        name = json_dict["label"]
    elif "name" in json_dict:
        name = json_dict["name"]
    else:
        name = parse_name_from_file(variants_file)

    if not name:
        return variant_list, name, rules, mrca_lineage
    if constellation_names and name not in constellation_names:
        return variant_list, name, rules, mrca_lineage

    print("\nParsing constellation JSON file %s" % variants_file)

    if "sites" in json_dict:
        for site in json_dict["sites"]:
            record = variant_to_variant_record(site, refseq, features_dict)
            if record != {}:
                variant_list.append(record)
    if include_ancestral and "ancestral" in json_dict:
        for site in json_dict["ancestral"]:
            record = variant_to_variant_record(site, refseq, features_dict)
            if record != {}:
                variant_list.append(record)

    if "rules" in json_dict:
        rules = json_dict["rules"]

    in_json.close()

    return variant_list, name, rules, mrca_lineage


def parse_csv_in(refseq, features_dict, variants_file, constellation_names=None):
    """
    returns variant_list and name
    """
    variant_list = []
    compulsory = []

    name = parse_name_from_file(variants_file)
    if constellation_names and name not in constellation_names:
        return variant_list, name, compulsory

    print("\nParsing constellation CSV file %s" % variants_file)

    csv_in = open("%s" % variants_file, 'r')
    reader = csv.DictReader(csv_in, delimiter=",")

    if "id" not in reader.fieldnames:
        csv_in.close()
        csv_in = open("%s" % variants_file, 'r', encoding="utf-8-sig")
        reader = csv.DictReader(csv_in, delimiter=",")

        if "id" not in reader.fieldnames:
            csv_in.close()
            print("Warning: CSV headerline does not contain 'id': %s - ignoring" % reader.fieldnames)
            return variant_list, name, compulsory

    for row in reader:
        if ":" not in row["id"] and "gene" in reader.fieldnames:
            var = "%s:%s" % (row["gene"], row["id"])
        else:
            var = row["id"]
        record = variant_to_variant_record(var, refseq, features_dict)
        if record != {}:
            variant_list.append(record)
        if "compulsory" in reader.fieldnames and row["compulsory"] in ["True", True, "Y", "y", "Yes", "yes", "YES"]:
            compulsory.append(record["name"])

    csv_in.close()
    rules = None
    if len(compulsory) > 0:
        rules = {}
        for var in compulsory:
            rules[var] = "alt"
    return variant_list, name, rules


def parse_textfile_in(refseq, features_dict, variants_file, constellation_names=None):
    """
    returns variant_list and name
    """
    variant_list = []

    name = parse_name_from_file(variants_file)
    if constellation_names and name not in constellation_names:
        return variant_list, name

    print("\nParsing constellation text file %s" % variants_file)

    with open("%s" % variants_file, "r") as f:
        for line in f:
            l = line.split("#")[0].strip()  # remove comments from the line
            if len(l) > 0:  # skip blank lines (or comment only lines)
                record = variant_to_variant_record(l, refseq, features_dict)
                if record != {}:
                    variant_list.append(record)

    return variant_list, name


def parse_variants_in(refseq, features_dict, variants_file, constellation_names=None, include_ancestral=False, label=None):
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
    mrca_lineage = ""

    if variants_file.endswith(".json"):
        variant_list, name, rule_dict, mrca_lineage = parse_json_in(refseq, features_dict, variants_file, constellation_names, include_ancestral=include_ancestral,label=label)
    elif variants_file.endswith(".csv"):
        variant_list, name, rule_dict = parse_csv_in(refseq, features_dict, variants_file, constellation_names)

    if len(variant_list) == 0 and not variants_file.endswith(".json"):
        variant_list, name = parse_textfile_in(refseq, features_dict, variants_file, constellation_names)

    return name, variant_list, rule_dict, mrca_lineage


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
            query_allele = query.translate()
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
        query_allele = record_seq.upper()[var["ref_start"] - 4:var["ref_start"] + 3*var["length"] + 2]
        query_allele = query_allele.replace("-","")
        if len(query_allele) % 3 != 0:
            print("Warning: while typing variant %s (before,ref,after) = (%s,%s,%s) found sequence with query allele %s. Handling by adding Ns which will result in ambiguous calls" %(var["name"], var["before"], var["ref_allele"], var["after"], record_seq.upper()[var["ref_start"] - 4:var["ref_start"] + 3*var["length"] + 2]))
        while len(query_allele) % 3 != 0:
            query_allele += "N"
        query_allele = query_allele.translate()
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
        elif "N" in query_allele:
            call = 'ambig'
            if not oth_char:
                query_allele = "X"
            else:
                query_allele = "N"
        elif query_allele == "-" * var["length"]:
            call = 'alt'
            query_allele = int(var["length"] / 3)
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

def counts_follow_rules(counts, rules):
    # rules allowed include "max_ref", "min_alt"
    is_rule_follower = True
    for rule in rules:
        if ":" in rule:
            continue
        elif str(rule).startswith("min") or str(rule).startswith("max"):
            rule_parts = rule.split("_")
            if len(rule_parts) <= 1:
                continue
            if rule_parts[0] == "min" and counts[rule_parts[1]] < rules[rule]:
                is_rule_follower = False
            elif rule_parts[0] == "max" and counts[rule_parts[1]] > rules[rule]:
                is_rule_follower = False
            else:
                counts["rules"] += 1
        else:
            print("Warning: Ignoring rule %s:%s" % (rule, str(rules[rule])))
    return is_rule_follower

def count_and_classify(record_seq, variant_list, rules):
    assert rules is not None
    counts = {'ref': 0, 'alt': 0, 'ambig': 0, 'oth': 0, 'rules':0}
    is_rule_follower = True

    for var in variant_list:
        call, query_allele = call_variant_from_fasta(record_seq, var)
        #print(var, call, query_allele)
        counts[call] += 1
        if var["name"] in rules:
            if var_follows_rules(call, rules[var["name"]]):
                counts['rules'] += 1
            elif is_rule_follower:
                is_rule_follower = False

    counts['support'] = round(counts['alt']/float(counts['alt'] + counts['ref'] + counts['ambig'] + counts['oth']),4)
    counts['conflict'] = round(counts['ref'] /float(counts['alt'] + counts['ref'] + counts['ambig'] + counts['oth']),4)

    if not is_rule_follower:
        return counts, False
    else:
        return counts, counts_follow_rules(counts, rules)


def generate_barcode(record_seq, variant_list, ref_char=None, ins_char="?", oth_char="X"):
    barcode_list = []
    counts = {'ref': 0, 'alt': 0, 'ambig': 0, 'oth': 0}

    for var in variant_list:
        call, query_allele = call_variant_from_fasta(record_seq, var, ins_char, oth_char)
        # print(var, call, query_allele)
        counts[call] += 1
        if ref_char is not None and call == 'ref':
            barcode_list.append(str(ref_char))
        else:
            barcode_list.append(str(query_allele))

    counts['support'] = round(counts['alt'] / float(counts['alt'] + counts['ref'] + counts['ambig'] + counts['oth']), 4)
    counts['conflict'] = round(counts['ref'] / float(counts['alt'] + counts['ref'] + counts['ambig'] + counts['oth']), 4)

    return ''.join(barcode_list), counts


def type_constellations(in_fasta, list_constellation_files, constellation_names, out_csv, reference_json, ref_char=None,
                        output_counts=False, label=None):
    reference_seq, features_dict = load_feature_coordinates(reference_json)

    constellation_dict = {}
    for constellation_file in list_constellation_files:
        constellation, variants, ignore, mrca_lineage = parse_variants_in(reference_seq, features_dict, constellation_file, label=label)
        if not constellation:
            continue
        if len(variants) > 0:
            constellation_dict[constellation] = variants
            print("Found file %s for constellation %s containing %i variants" % (
                    constellation_file, constellation, len([v["name"] for v in variants])))
        else:
            print("Warning: %s is not a valid constellation file - ignoring" % constellation_file)

    variants_out = open(out_csv, "w")
    variants_out.write("query,%s\n" % ",".join(list(constellation_dict.keys())))

    counts_out = {}
    if output_counts:
        for constellation in constellation_dict:
            clean_name = re.sub("[^a-zA-Z0-9_\-.]", "_", constellation)
            counts_out[constellation] = open("%s.%s_counts.csv" % (out_csv.replace(".csv", ""), clean_name), "w")
            counts_out[constellation].write("query,ref_count,alt_count,ambig_count,other_count,support,conflict\n")

    with open(in_fasta, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            if len(record.seq) != len(reference_seq):
                sys.stderr.write("The fasta record for %s isn't as long as the reference, is this fasta file aligned?\n"
                                 % record.id)
                sys.exit(1)

            out_list = [record.id]
            for constellation in constellation_dict:
                barcode, counts = generate_barcode(record.seq, constellation_dict[constellation], ref_char)
                if output_counts:
                    counts_out[constellation].write("%s,%i,%i,%i,%i,%f,%f\n" % (record.id, counts['ref'], counts['alt'],
                                                 counts['ambig'], counts['oth'], counts['support'], counts['conflict']))
                out_list.append(barcode)

            variants_out.write("%s\n" % ",".join(out_list))

    variants_out.close()
    for constellation in counts_out:
        if counts_out[constellation]:
            counts_out[constellation].close()


def classify_constellations(in_fasta, list_constellation_files, constellation_names, out_csv, reference_json,
                            output_counts=False, call_all=False, long=False, label=None):

    reference_seq, features_dict = load_feature_coordinates(reference_json)

    constellation_dict = {}
    rule_dict = {}
    mrca_lineage_dict = {}
    for constellation_file in list_constellation_files:
        constellation, variants, rules, mrca_lineage = parse_variants_in(reference_seq, features_dict, constellation_file,
                                                           constellation_names, include_ancestral=True, label=label)
        if constellation_names and constellation not in constellation_names:
            continue
        if not rules:
            print("Warning: No rules provided to classify %s - ignoring" % constellation)
            continue
        else:
            rule_dict[constellation] = rules
        if len(variants) > 0:
            constellation_dict[constellation] = variants
            print("Found file %s for constellation %s containing %i variants" % (
            constellation_file, constellation, len([v["name"] for v in variants])))
            print("Rules", rule_dict[constellation])
            mrca_lineage_dict[constellation] = mrca_lineage
        else:
            print("Warning: %s is not a valid constellation file - ignoring" % constellation_file)

    variants_out = open(out_csv, "w")
    if long and not call_all:
        variants_out.write("query,constellations,mrca_lineage,ref_count,alt_count,ambig_count,other_count,rule_count,support,"
                           "conflict\n")
    else:
        variants_out.write("query,constellations,mrca_lineage\n")

    counts_out = {}
    if output_counts:
        for constellation in constellation_dict:
            clean_name = re.sub("[^a-zA-Z0-9_\-.]","_",constellation)
            counts_out[constellation] = open("%s.%s_counts.csv" % (out_csv.replace(".csv", ""), clean_name), "w")
            counts_out[constellation].write("query,ref_count,alt_count,ambig_count,other_count,rule_count,support,"
                                            "conflict,call\n")

    with open(in_fasta, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            if len(record.seq) != len(reference_seq):
                sys.stderr.write("The fasta record for %s isn't as long as the reference, is this fasta file aligned?\n"
                                 % record.id)
                sys.exit(1)

            lineages = []
            best_constellation = None
            best_support = 0
            best_conflict = 1
            best_counts = None
            for constellation in constellation_dict:
                counts, call = count_and_classify(record.seq,
                                                  constellation_dict[constellation],
                                                  rule_dict[constellation])
                if call:
                    if call_all:
                        lineages.append(constellation)
                    elif (not best_constellation) \
                        or (counts['support'] > best_support) \
                        or (counts['support'] == best_support and counts['conflict'] < best_conflict)\
                        or (counts['support'] == best_support and counts['conflict'] == best_conflict and counts['rules'] > best_counts["rules"]):
                        best_constellation = constellation
                        best_support = counts['support']
                        best_conflict = counts['conflict']
                        best_counts = counts
                elif len(constellation_dict) == 1:
                    best_counts = counts

                if output_counts:
                    counts_out[constellation].write(
                        "%s,%i,%i,%i,%i,%i,%f,%f,%s\n" % (record.id, counts['ref'], counts['alt'], counts['ambig'],
                                                       counts['oth'], counts['rules'], counts['support'],
                                                          counts['conflict'], call))
            if not call_all and best_constellation:
                lineages.append(best_constellation)
            if long and best_counts is not None:
                variants_out.write("%s,%s,%s,%i,%i,%i,%i,%i,%f,%f\n" % (record.id, "|".join(lineages),
                                                                     "|".join([mrca_lineage_dict[l] for l in lineages]),
                                                                     best_counts['ref'],
                                                                     best_counts['alt'], best_counts['ambig'],
                                                                     best_counts['oth'], best_counts['rules'],
                                                                     best_counts['support'], best_counts['conflict']))
            elif long and not call_all:
                variants_out.write("%s,%s,%s,,,,,,,\n" % (record.id, "|".join(lineages),
                                   "|".join([mrca_lineage_dict[l] for l in lineages])))
            else:
                variants_out.write("%s,%s,%s\n" % (record.id, "|".join(lineages),
                                                   "|".join([mrca_lineage_dict[l] for l in lineages])))

    variants_out.close()
    for constellation in counts_out:
        if counts_out[constellation]:
            counts_out[constellation].close()


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
