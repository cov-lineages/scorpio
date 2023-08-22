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

from .definitions import Reference, Constellation, Variant

### idea: keep constellationd in json dict form and work directly with them insead of other dicts

if sys.version_info[0] < 3:
    raise Exception("Python 3 or a more recent version is required.")



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


def parse_mutations(reference, mutations_list, ignore_fails=False):
    """
    Parse the mutations specified on command line and make a mutations constellation for them

    returns variants which is a list of dicts of snps, aas and dels,
    one dict per variant. format of subdict varies by variant type
    """
    variants = []
    problematic = []

    for mutation in mutations_list:
        record = Variant(mutation, reference, ignore_fails=ignore_fails)
        if record != {}:
            variants.append(record)
        else:
            problematic.append(mutation)

    if len(problematic) > 0:
        sys.stderr.write("The following mutations were not provided in an acceptable format: %s\n\n"
                         % ",".join(problematic))
        sys.exit(1)
    return variants


def call_variant_from_fasta(record_seq, var, ins_char="?", oth_char=None, codon=False):
    #print("Call variant for ", var)
    call = None
    query_allele = None

    if var.type == "ins":
        call = "oth"
        query_allele = ins_char
    elif var.type == "snp":
        query_allele = record_seq.upper()[var.ref_start - 1]
        if query_allele == var.ref_allele:
            call = 'ref'
        elif query_allele == "N":
            call = 'ambig'
        elif query_allele == var.alt_allele or var.alt_allele == "":
            call = 'alt'
        else:
            call = 'oth'

    elif var.type == "aa":
        try:
            query = record_seq.upper()[var.ref_start - 1:var.ref_start - 1 + 3 * len(var.ref_allele)]
            query_allele = query.translate(gap = "-")
            if query_allele == var.ref_allele:
                call = 'ref'
            elif query_allele == "X":
                call = 'ambig'
            elif query_allele == var.alt_allele or var.fuzzy:
                call = 'alt'
            else:
                call = 'oth'
        except:
            call = 'oth'

    elif var.type == "del" and var.space == "aa":
        query_allele = str(record_seq.upper()[var.ref_start - 4:var.ref_start + 3*var.length + 2])
        query = query_allele.replace("-","")
        if len(query) % 3 != 0:
            query = query.replace("N","")
        if len(query) % 3 != 0:
            query = query_allele.replace("-","N")
        if len(query) % 3 != 0:
            logging.warning("Warning: while typing variant %s (before,ref,after) = (%s,%s,%s) found sequence with query allele %s treated as %s. Handling by adding Ns which will result in ambiguous calls" %(var.name, var.before, var.ref_allele, var.after, query_allele, query))
        query_allele = query
        while len(query_allele) % 3 != 0:
            query_allele += "N"
        query_allele = Seq(query_allele).translate()
        if query_allele == var.before + var.ref_allele + var.after:
            call = 'ref'
            query_allele = 0
        elif query_allele == var.before + var.after:
            call = 'alt'
            query_allele = var.length
        elif "X" in query_allele:
            call = 'ambig'
            query_allele = "X"
        else:
            call = 'oth'
            if not oth_char:
                query_allele = "X"
    elif var.type == "del" and var.space == "nuc":
        query_allele = record_seq.upper()[var.ref_start - 1:var.ref_start + var.length - 1]

        if query_allele == var.ref_allele:
            call = 'ref'
            query_allele = 0
        elif query_allele == "-" * var.length or query_allele == "N" * var.length:
            call = 'alt'
            query_allele = max(int(var.length / 3), 1)
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

    if call == "oth" and var.type != "ins" and oth_char:
        query_allele = oth_char

    return call, query_allele


def var_follows_rule(call, rule):
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


def count_and_classify(record_seq, constellation):
    assert constellation.rules is not None
    counts = {'ref': 0, 'alt': 0, 'ambig': 0, 'oth': 0, 'rules': {},
              'substitution': {'ref': 0, 'alt': 0, 'ambig': 0, 'oth': 0},
              'indel': {'ref': 0, 'alt': 0, 'ambig': 0, 'oth': 0}}
    is_rule_follower_dict = {}
    for key in constellation.rules:
        is_rule_follower_dict[key] = True
        counts["rules"][key] = 0

    for var in constellation.variants:
        call, query_allele = call_variant_from_fasta(record_seq, var)
        #print(var, call, query_allele)
        counts[call] += 1
        if var.type in ["aa", "snp"]:
            counts["substitution"][call] += 1
        elif var.type in ["ins", "del"]:
            counts["indel"][call] += 1
        for key in constellation.rules:
            if var.name in constellation.rules[key]:
                if var_follows_rule(call, constellation.rules[key][var.name]):
                    counts['rules'][key] += 1
                elif is_rule_follower_dict[key]:
                    is_rule_follower_dict[key] = False

    counts['support'] = round(counts['alt']/float(counts['alt'] + counts['ref'] + counts['ambig'] + counts['oth']),4)
    counts['conflict'] = round(counts['ref'] /float(counts['alt'] + counts['ref'] + counts['ambig'] + counts['oth']),4)

    for key in constellation.rules:
        if not is_rule_follower_dict[key]:
            continue
        else:
            call, note = counts_follow_rules(counts, constellation.rules, key)
            if call:
                counts["rules"] = counts["rules"][key]
                call = key
                return counts, call, note
    counts["rules"] = counts["rules"]["default"]
    return counts, False, ""


def generate_barcode(record_seq, variants, ref_char=None, ins_char="?", oth_char="X",constellation_count_dict=None):
    barcode_list = []
    counts = {'ref': 0, 'alt': 0, 'ambig': 0, 'oth': 0}
    sorted_alt_sites = []

    count = 1
    for var in variants:
        call, query_allele = call_variant_from_fasta(record_seq, var, ins_char, oth_char)
        # #print(var, call, query_allele)
        counts[call] += 1
        if constellation_count_dict and "constellations" in var:
            for constellation in var.constellations:
                constellation_count_dict[constellation][call] += 1
            if call == "alt":
                sorted_alt_sites.append(var.constellations)

        if ref_char is not None and call == 'ref':
            barcode_list.append(str(ref_char))
        else:
            barcode_list.append(str(query_allele))
        if not len(barcode_list) == count:
            print(count, call, query_allele)
        count += 1

    counts['support'] = round(counts['alt'] / float(counts['alt'] + counts['ref'] + counts['ambig'] + counts['oth']), 4)
    counts['conflict'] = round(counts['ref'] / float(counts['alt'] + counts['ref'] + counts['ambig'] + counts['oth']), 4)

    return barcode_list, counts, constellation_count_dict, sorted_alt_sites


def load_constellations(list_constellation_files, constellation_names, reference, label=None,
                        include_ancestral=True, rules_required=False, ignore_fails=False):
    constellation_dict = {}
    lineage_name_dict = {}

    for constellation_file in list_constellation_files:
        constellation = Constellation(from_file=True, reference=reference, variants_file=constellation_file,
                                      include_ancestral=include_ancestral, label=label, ignore_fails=ignore_fails)
        if constellation.lineage_name:
            lineage_name_dict[constellation.lineage_name] = constellation

        if not constellation.name:
            continue

        if rules_required and not constellation.rules:
            logging.warning("Warning: No rules provided to classify %s - ignoring" % constellation)
            continue

        if len(constellation.variants) > 0:
            constellation_dict[constellation.name] = constellation
            if constellation_names and constellation.name not in constellation_names \
                    and constellation.output_name not in constellation_names:
                continue
            logging.info("\n")
            logging.info("Found file %s for constellation %s containing %i defining mutations" % (
                constellation_file, constellation.name, len([v.name for v in constellation.variants])))
            if rules_required:
                logging.info("Rules %s" % constellation.rules)
        else:
            logging.warning("Warning: %s is not a valid constellation file - ignoring" % constellation_file)

    return constellation_dict, lineage_name_dict


def parse_mutations_list(mutations_list, reference):
    new_mutations_list = []
    for entry in mutations_list:
        if '.' in entry:
            new_mutations_list.extend(parse_mutations_in(entry)) # this is a file
        else:
            new_mutations_list.append(entry)
    mutations_list = new_mutations_list
    mutation_variants = parse_mutations(reference, mutations_list)
    return mutations_list, mutation_variants


def combine_constellations(constellation_dict):
    variant_dict = {}
    constellation_count_dict = {}
    for constellation in constellation_dict.values():
        constellation_count_dict[constellation] = {"total": len(constellation.variants), 'ref': 0, 'alt': 0,
                                                   'ambig': 0, 'oth': 0}
        for variant in constellation.variants:
            if variant.name in variant_dict:
                variant_dict[variant.name]["constellations"].append(constellation)
            else:
                variant_dict[variant.name] = variant
                variant_dict[variant.name]["constellations"] = [constellation]
    union_constellation = Constellation(from_file=False, name="union", variants_list=variant_dict.values())
    return union_constellation, constellation_count_dict


def combine_constellations_by_name(constellation_dict, reference):
    new_constellation_dict = {}
    for constellation in constellation_dict.values():
        if not constellation.output_name:
            continue
        elif constellation.output_name == "None":
            constellation.output_name = constellation.name

        if constellation.output_name not in new_constellation_dict:
            new_constellation_dict[constellation.output_name] = constellation
        else:
            for variant in constellation.variants:
                found = False
                for existing_variant in new_constellation_dict[constellation.output_name].variants:
                    if variant.equals(existing_variant, reference):
                        found = True
                        break
                if not found:
                    new_constellation_dict[constellation.output_name].variants.append(variant)
    return new_constellation_dict


def file_writer(out_csv, list_of_constellations, single_file_named_summary, q):
    """listens for messages on the q, writes to files. """
    handle_dict = {}
    if len(list_of_constellations) == 1 and not single_file_named_summary:
        handle_dict[list_of_constellations[0]] = open(out_csv, "w")
    elif len(list_of_constellations) == 1 and single_file_named_summary:
        handle_dict["summary"] = open(out_csv, "w")
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


def type_record(record, reference, constellation_names, constellation_dict, constellation_count_dict, output_counts,
                append_genotypes, combination, ref_char, q):
    if len(record.seq) != len(reference.refseq):
        sys.stderr.write("The fasta record for %s isn't as long as the reference, is this fasta file aligned?\n"
                         % record.id)
        sys.exit(1)

    out_list = [record.id]

    write_extended_files = output_counts or append_genotypes
    single_file = len(constellation_names) == 1

    for constellation_name in constellation_names:
        for constellation in constellation_dict.values():
            if constellation.output_name != constellation_name:
                continue
            else:
                logging.debug("Consider constellation %s" % constellation.name)
                barcode_list, counts, sample_constellation_count_dict, sorted_alt_sites = generate_barcode(record.seq,
                                                                                                           constellation.variants,
                                                                                                           ref_char,
                                                                                                           constellation_count_dict=copy.deepcopy(constellation_count_dict))
                if write_extended_files:
                    columns = [record.id]
                    columns.append(''.join(barcode_list))
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
                    q.put((constellation.output_name, res))

                out_list.append(''.join(barcode_list))
                break

    res = "%s\n" % ",".join(out_list)
    if len(constellation_names) > 1 or not write_extended_files:
        q.put(("summary", res))
    return res


def type_constellations(in_fasta, list_constellation_files, constellation_names, out_csv, reference_json, ref_char=None,
                        output_counts=False, label=None, append_genotypes=False, mutations_list=None, dry_run=False,
                        combination=False, interspersion=False, threads=1):
    reference = Reference(reference_json)

    constellation_dict, lineage_name_dict = load_constellations(list_constellation_files,
                                                                constellation_names,
                                                                reference,
                                                                label=label,
                                                                include_ancestral=False,
                                                                rules_required=False)

    if combination:
        constellation_names = ["union"]
        constellation_dict, constellation_count_dict = combine_constellations(constellation_dict)
    else:
        constellation_count_dict = None

    if mutations_list:
        mutations_list, mutation_variants = parse_mutations_list(mutations_list, reference)
        if len(constellation_dict.keys()) == 1 and "mutations" not in constellation_dict.keys():
            constellation = list(constellation_dict.keys())[0]
            new_constellation = "%s+%s" % (constellation, '|'.join(mutations_list))
            constellation_dict[new_constellation] = constellation_dict[constellation] + mutation_variants
            del constellation_dict[constellation]
        else:
            constellation_dict["mutations"] = Constellation(from_file=False, name="mutations",
                                                            variants_list=mutation_variants)

    if dry_run:
        return

    logging.info("\n")
    logging.info("Update constellation dict")
    constellation_dict = combine_constellations_by_name(constellation_dict, reference)

    if not constellation_names:
        constellation_names = constellation_dict.keys()
    if mutations_list:
        constellation_names.append("mutations")
    logging.debug(constellation_dict)
    logging.info("Have %i constellations to type: %s" % (len(constellation_dict), list(constellation_dict.keys())))
    if constellation_count_dict:
        logging.info("Have %i candidate constellations to collect counts: %s" % (len(constellation_count_dict),
                                                                                 list(constellation_count_dict.keys())))

    # setup manager
    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(threads)

    write_extended_files = append_genotypes or output_counts

    single_file_named_summary = False
    if len(constellation_names) == 1 and not write_extended_files:
        single_file_named_summary = True

    if len(constellation_names) > 1 or single_file_named_summary:
        header = "query,%s\n" % ",".join(constellation_names)
        q.put(("summary", header))

    # put listener to work first
    list_constellation_names = constellation_names
    if not write_extended_files:
        list_constellation_names = []
    if threads > 1:
        pool.apply_async(file_writer, (out_csv, list_constellation_names, single_file_named_summary, q))

    if write_extended_files:
        for constellation in constellation_dict.values():
            if constellation.output_name not in constellation_names:
                continue
            columns = ["query"]
            if not single_file_named_summary:
                columns.append("barcode")
            if output_counts:
                columns.append("ref_count,alt_count,ambig_count,other_count,support,conflict")
            if append_genotypes:
                columns.extend([var.name for var in constellation.variants])
            if combination:
                columns.append("notes")
            header = "%s\n" % ','.join(columns)
            if len(constellation_dict) > 1:
                q.put((constellation.output_name, header))
            else:
                q.put(("summary", header))

    # fire off workers
    jobs = []
    with open(in_fasta, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            job = pool.apply_async(type_record, (record, reference, constellation_names, constellation_dict,
                                                 constellation_count_dict, output_counts, append_genotypes,
                                                 combination, ref_char, q))
            jobs.append(job)

    # collect results from the workers through the pool result queue
    for job in jobs:
        job.get()

    # now we are done, kill the listener
    if len(constellation_names) > 1 or single_file_named_summary:
        q.put(("summary", 'kill'))

    if len(constellation_names) > 1 or not single_file_named_summary:
        for constellation_name in list_constellation_names:
            q.put((constellation_name, "kill"))

    if threads == 1:
        pool.apply_async(file_writer, (out_csv, list_constellation_names, single_file_named_summary, q))

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


def classify_record(record, reference, constellation_names, constellation_dict, lineage_name_dict,
                    output_counts, call_all, long, list_incompatible, mutations_list, mutation_variants, interspersion, q):
    if len(record.seq) != len(reference.refseq):
        sys.stderr.write("The fasta record for %s isn't as long as the reference, is this fasta file aligned?\n"
                         % record.id)
        sys.exit(1)

    list_constellation_names = list(set(constellation_names))
    single_file = len(list_constellation_names) == 1 and output_counts

    constellation_calls = []  # eg Probable Omicron
    constellations_called = []  # eg the omicron constellation class
    best_constellation = None
    best_support = 0
    best_conflict = 1
    best_counts = None
    best_call = False
    scores = {}
    children = {}
    for constellation in constellation_dict.values():
        if not constellation.output_name or (constellation.output_name not in list_constellation_names and constellation.name not in constellation_names):
            continue
        logging.debug("Consider constellation %s" % constellation.name)
        parents = []

        counts, call, note = count_and_classify(record.seq,
                                                constellation)
        current_constellation = constellation
        while current_constellation.parent_lineage and current_constellation.parent_lineage in lineage_name_dict:
            logging.debug("Current %s in parent dict" % current_constellation.name)
            #print("consider parent", current_constellation.parent_lineage)
            #print(lineage_name_dict.keys())
            #print(lineage_name_dict[current_constellation.parent_lineage])
            current_constellation = lineage_name_dict[current_constellation.parent_lineage]
            parents.append(current_constellation.name)
            parent_counts, parent_call, parent_note = count_and_classify(record.seq,
                                                                         current_constellation)
            #print("parent call", parent_call)
            counts, call, note = combine_counts_call_notes(counts, call, note, parent_counts, parent_call, parent_note)
            #print("overall call", call)

        for parent in parents:
            if parent not in children:
                children[parent] = []
            children[parent].append(constellation.name)

        if call:
            #print("have call", call)
            #print("current best", best_constellation)
            #print("children", children)
            #print("parents", parents)
            logging.debug("Have call for %s" % constellation.output_name)
            if call_all:
                if call != "default" and constellation.output_name not in ["",None,"None"]:
                    called_constellation_name = "%s %s" % (call, constellation.output_name)
                elif constellation.output_name in [None,"None"]:
                    called_constellation_name = ""
                else:
                    called_constellation_name = constellation.output_name
                constellation_calls.append(called_constellation_name)
                constellations_called.append(constellation)
            elif constellation.name in children \
                    and best_constellation \
                    and best_constellation.name in children[constellation.name]:
                logging.debug("Ignore as parent of best constellation")
            elif (not best_constellation) \
                    or (counts['support'] > best_support) \
                    or (counts['support'] == best_support and counts['conflict'] < best_conflict) \
                    or (counts['support'] == best_support and counts['conflict'] == best_conflict
                                                          and counts['rules'] > best_counts["rules"]) \
                    or (best_constellation.name in parents):
                if constellation.mrca_lineage:
                    best_constellation = constellation
                    logging.debug("Set best constellation %s" % best_constellation.output_name)
                    #print("Set best constellation", best_constellation.output_name)
                    best_support = counts['support']
                    best_conflict = counts['conflict']
                    best_counts = counts
                    best_call = call

        if interspersion:
            if counts["alt"] > 1:
                summary = constellation.output_name
                score = counts["support"]
                scores[score] = summary

        if output_counts:
            out_entries = ["%s,%i,%i,%i,%i,%i,%f,%f,%s,%s,%s" % (
                    record.id, counts['ref'], counts['alt'], counts['ambig'],
                    counts['oth'], counts['rules'], counts['support'],
                    counts['conflict'], call, constellation.output_name, note)]
            if mutations_list and single_file:
                barcode_list, counts, ignore, ignore2 = generate_barcode(record.seq, mutation_variants)
                out_entries.extend(barcode_list)
            res = "%s\n" % ','.join(out_entries)
            q.put((constellation.output_name, res))

    if not call_all and best_constellation:
        if best_call != "default" and best_constellation.output_name not in ["",None,"None"]:
            best_constellation_name = "%s %s" % (best_call, best_constellation.output_name)
        elif best_constellation.output_name in [None, "None"]:
            best_constellation_name = ""
        else:
            best_constellation_name = best_constellation.output_name
        constellation_calls.append(best_constellation_name)
        constellations_called.append(best_constellation)

    if len(constellation_calls) > 1:
        zipped = zip(constellation_calls, constellations_called)
        zipped = sorted(zipped, key=lambda x: x[1].name)
        constellation_calls = [x for x, y in zipped]
        constellations_called = [y for x,y in zipped]

    out_entries = [record.id, "|".join(constellation_calls), "|".join([n.mrca_lineage for n in constellations_called])]

    if list_incompatible:
        out_entries.append("|".join([n.incompatible_lineage_calls for n in constellations_called]))
    if long and best_counts is not None:
        out_entries.append("%i,%i,%i,%i,%i,%f,%f" % (best_counts['ref'],
                                                     best_counts['alt'], best_counts['ambig'],
                                                     best_counts['oth'], best_counts['rules'],
                                                     best_counts['support'], best_counts['conflict']))
    elif long and not call_all:
        out_entries.append(",,,,,,")
    out_entries.append("|".join([n.name for n in constellations_called]))

    if single_file:
        out_entries.append("%s,%s" % (call, note))

    if mutations_list:
        barcode_list, counts, ignore, ignore2 = generate_barcode(record.seq, mutation_variants)
        out_entries.extend(barcode_list)

    res = "%s\n" % ",".join(out_entries)
    if not single_file:
        q.put(("summary", res))

    return res


def classify_constellations(in_fasta, list_constellation_files, constellation_names, out_csv, reference_json,
                            output_counts=False, call_all=False, long=False, label=None, list_incompatible=False,
                            mutations_list=None, dry_run=False, interspersion=False, threads=1):
    reference = Reference(reference_json)

    constellation_dict, lineage_name_dict = load_constellations(list_constellation_files,
                                                                                    constellation_names,
                                                                                    reference,
                                                                                    label=label,
                                                                                    include_ancestral=True,
                                                                                    rules_required=True,
                                                                                    ignore_fails=True)

    if not constellation_names:
        constellation_names = [c.output_name for c in constellation_dict.values() if c.output_name not in ["None", "", None]]
    constellation_names = list(set(constellation_names))
    logging.debug("lineage_name_dict: %s" % lineage_name_dict)

    mutation_variants = None
    if mutations_list:
        mutations_list, mutation_variants = parse_mutations_list(mutations_list, reference)
    if dry_run:
        return

    # setup manager
    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(threads)

    # put listener to work first
    list_constellation_names = list(set(constellation_names))
    single_file = len(list_constellation_names) == 1 and output_counts
    if single_file:
        long = True
    if not output_counts:
        list_constellation_names = []

    if threads > 1:
        pool.apply_async(file_writer, (out_csv, list_constellation_names, single_file, q))

    # write summary file header
    columns = get_columns(list_incompatible, long, call_all, mutations_list, single_file)
    header = "%s\n" % ",".join(columns)
    if single_file:
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
            job = pool.apply_async(classify_record, (record, reference, constellation_names, constellation_dict,
                                                     lineage_name_dict, output_counts, call_all,
                                                     long, list_incompatible, mutations_list, mutation_variants,
                                                     interspersion, q))
            jobs.append(job)

    # collect results from the workers through the pool result queue
    for job in jobs:
        job.get()

    #now we are done, kill the listener
    if not single_file:
        q.put(("summary", "kill"))
    for constellation_name in constellation_names:
        q.put((constellation_name, "kill"))

    if threads == 1:
        pool.apply_async(file_writer, (out_csv, list_constellation_names, False, q))

    pool.close()
    pool.join()


def list_constellations(list_constellation_files, constellation_names, reference_json, label=None):
    reference = Reference(reference_json)

    list_of_constellations = set()
    for constellation_file in list_constellation_files:
        constellation = Constellation(from_file=True, reference=reference, variants_file=constellation_file, label=label)
        if not constellation.output_name:
            continue
        if constellation_names \
                and constellation.output_name not in constellation_names \
                and constellation.name not in constellation_names:
            continue
        if len(constellation.variants) > 0 and constellation.mrca_lineage:
            list_of_constellations.add(constellation.output_name)
            list_of_constellations.add(constellation.mrca_lineage)
        elif len(constellation.variants) > 0:
            list_of_constellations.add(constellation.output_name)
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
