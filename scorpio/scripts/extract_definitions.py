#!/usr/bin/env python3

import sys
import csv
import operator
import argparse
import json
import logging
from Bio.Seq import Seq
from operator import itemgetter

from .type_constellations import load_feature_coordinates, resolve_ambiguous_cds

def parse_args():
    parser = argparse.ArgumentParser(description="""Pick a representative sample for each unique sequence""",
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--in-variants', dest = 'in_variants', required=True, help='CSV of containing sequence_name and nucleotide_variants or nucleotide_mutations columns, the latter being | separated list of variants')
    parser.add_argument('--in-groups', dest = 'in_groups', required=False, help='CSV of containing sequence_name and a column defining groups ')
    parser.add_argument('--group-column', dest = 'group_column', required=False, default='lineage', help='Column name defining the groups')
    parser.add_argument('--index-column', dest = 'index_column', required=False, default='sequence_name', help='Taxon column name')
    parser.add_argument('--reference_json', help='JSON file containing keys "genome" with reference sequence '
                                                 'and "proteins", "features" or "genes" with features of interest'
                                                 ' and their coordinates')
    parser.add_argument('--prefix', dest = 'prefix', default=".", help='Output directory for generated JSON')

    args = parser.parse_args()
    return args


def parse_outgroups(outgroup_file):
    """
    input is CSV with columns lineage,outgroups. Can include multiple outgroup sequence_names
    separated by a pipe.
    Returns a dictionary of outgroup sequence_name : lineage(s) for which it is an outgroup
    """
    outgroup_dict = {}
    all_outgroups = set()
    if not outgroup_file:
        return outgroup_dict, all_outgroups
    with open(outgroup_file, "r") as outgroup_handle:
        line = outgroup_handle.readline()
        while line:
            try:
                lineage, outgroups = line.strip().split(",")
                outgroups = outgroups.split("|")
                all_outgroups.update(outgroups)
                for outgroup in outgroups:
                    if outgroup in outgroup_dict:
                        outgroup_dict[outgroup].append(lineage)
                    else:
                        outgroup_dict[outgroup] = [lineage]
            except:
                continue
            line = outgroup_handle.readline()
    return outgroup_dict, all_outgroups


def get_group_dict(in_variants, group_column, index_column, subset):
    group_dict = {}
    groups = set()

    with open(in_variants,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if group_column in row and index_column in row:
                if subset and row[group_column] not in subset:
                    continue
                if row[index_column] in group_dict:
                    logging.warning("%s is a duplicate in group CSV, keeping first")
                else:
                    group_dict[row[index_column]] = row[group_column]
                    groups.add(row[group_column])

    logging.info("Found %d groups" % len(groups))

    return group_dict


def update_var_dict(var_dict, group, variants):
    if group not in var_dict:
        var_dict[group] = {"total": 1}
    else:
        var_dict[group]["total"] += 1

    for var in variants.split("|"):
        if var in var_dict[group]:
            var_dict[group][var] += 1
        else:
            var_dict[group][var] = 1
    return


def update_feature_dict(feature_dict):
    for feature in feature_dict:
        if len(feature_dict[feature]) > 2:
            cds, aa_pos = resolve_ambiguous_cds(feature_dict[feature][2], feature_dict[feature][0], feature_dict)
            if aa_pos:
                feature_dict[feature] = (aa_pos, feature_dict[feature][1] + feature_dict[feature][0] - aa_pos, cds)
    return feature_dict


def get_common_mutations(var_dict, min_occurance=3, threshold_common=0.98, threshold_intermediate=0.25):
    sorted_tuples = sorted(var_dict.items(), key=operator.itemgetter(1))
    var_dict = {k: v for k, v in sorted_tuples}

    common = []
    intermediate = []
    for var in var_dict:
        #print("var", var, var_dict[var], min_occurance, var_dict["total"])
        if var != "total" and var_dict[var] == var_dict["total"]:
            common.append(var)
        elif var != "total" and var_dict[var] >= min_occurance:
            frequency = float(var_dict[var])/var_dict["total"]
            #print("var", var, var_dict[var], frequency, threshold_common, threshold_intermediate)
            if frequency > threshold_common:
                common.append(var)
            elif frequency > threshold_intermediate:
                intermediate.append("%s:%f" % (var, frequency))
    return common, intermediate


def translate_if_possible(nuc_start, nuc_ref, nuc_alt, feature_dict, reference_seq, include_protein=False):
    nuc_end = nuc_start + len(nuc_ref)
    nuc_start = int(nuc_start)
    nuc_end = int(nuc_end)
    #print(nuc_start, nuc_end, nuc_ref, nuc_alt)
    assert nuc_start > 10
    assert nuc_end > 10
    #print(nuc_start, nuc_end, reference_seq[nuc_start-1:nuc_end-1], nuc_ref, nuc_alt)
    assert reference_seq[nuc_start-1:nuc_end-1] == nuc_ref
    query_seq = reference_seq[:nuc_start-1] + nuc_alt + reference_seq[nuc_end-1:]
    #print(len(reference_seq), len(query_seq), type(int(nuc_start-5)), type(int(nuc_end+5)))
    #print(reference_seq[nuc_start-5: nuc_end+5])
    #print(query_seq[nuc_start-5: nuc_end+5])

    for feature in feature_dict:
        if len(feature_dict[feature]) > 2:
            continue # ignore nsp definitions
        if feature_dict[feature][0] <= nuc_start <= feature_dict[feature][1]:
            start, end = nuc_start, nuc_end
            while (start - feature_dict[feature][0]) % 3 != 0:
                start -= 1
            while (end - feature_dict[feature][0]) % 3 != 0:
                end += 1
            ref_allele = Seq(reference_seq[start - 1:end - 1 ]).translate()
            query_allele = Seq(query_seq[start - 1:end - 1]).translate()
            if ref_allele == query_allele:
                return "nuc:%s%i%s" % (nuc_ref, nuc_start, nuc_alt)
            aa_pos = int((start - feature_dict[feature][0]) / 3) + 1
            if include_protein:
                feature, aa_pos = translate_to_protein_if_possible(feature, aa_pos, feature_dict)
            #print(start, end, ref_allele, query_allele, aa_pos, feature)
            return "%s:%s%i%s" % (feature, ref_allele, aa_pos, query_allele)
    return "nuc:%s%i%s" % (nuc_ref, nuc_start, nuc_alt)


def translate_to_protein_if_possible(cds, aa_start, feature_dict):
    if not cds.startswith("orf"):
        return cds, aa_start

    for feature in feature_dict:
        if len(feature_dict[feature]) < 3:
            continue  # only want nsp definitions
        if feature_dict[feature][2] == cds:
            if feature_dict[feature][0] <= aa_start <= feature_dict[feature][1]:
                return feature, aa_start-feature_dict[feature][0]+1
    return cds, aa_start

def define_mutations(list_variants, feature_dict, reference_seq, include_protein=False):
    merged_list = []
    if not list_variants:
        return merged_list

    intermediate_list = []
    for variant in list_variants:
        if ":" in variant:
            var, freq = variant.split(":")
        else:
            var = variant
            freq = None
        if var.startswith("del"):
            new_var = ':'.join(var.split("_"))
            if freq:
                merged_list.append("%s:%s" %(new_var, freq))
            else:
                merged_list.append(new_var)
        elif var.startswith("ins") or var.startswith("del"):
            i, pos, add = var.split("_")
            new_var = "%s:%s+%s" % (i, pos, add)
            if freq:
                merged_list.append("%s:%s" % (new_var, freq))
            else:
                merged_list.append(new_var)
        else:
            intermediate_list.append([var[0], int(var[1:-1]), var[-1], freq])

    intermediate_list.sort(key=itemgetter(1))
    current = ["", 1, "", None]
    for new in intermediate_list:
        #print(new, current)
        if new[1] == current[1] + len(current[0]):
            #print("merge")
            current[0] += new[0]
            current[2] += new[2]
            if new[3] and current[3]:
                current[3] = min(current[3], new[3])
            elif new[3]:
                current[3] = new[3]
        elif current[0] != "":
            var = translate_if_possible(current[1], current[0], current[2], feature_dict, reference_seq, include_protein)
            if freq:
                merged_list.append("%s:%s" % (var, freq))
            else:
                merged_list.append(var)
            current = new
        else:
            current = new
    if current[0] != "":
        var = translate_if_possible(current[1], current[0], current[2], feature_dict, reference_seq, include_protein)
        if freq:
            merged_list.append("%s:%s" % (var, freq))
        else:
            merged_list.append(var)
    return merged_list


def subtract_outgroup(common, outgroup_common):
    updated_common = []
    ancestral = []
    for var in common:
        if var in outgroup_common:
            ancestral.append(var)
        else:
            updated_common.append(var)
    return updated_common, ancestral


def write_constellation(prefix, group, list_variants, list_intermediates, list_ancestral):
    group_dict = {"name": group, "sites": list_variants, "intermediate": list_intermediates,
                  "rules": {"min_alt": max(len(list_variants) - 3, min(len(list_variants), 3)), "max_ref": 3}}
    if list_ancestral:
        group_dict["ancestral"] = list_ancestral
    with open('%s/%s.json' % (prefix, group), 'w') as outfile:
        json.dump(group_dict, outfile, indent=4)


def extract_definitions(in_variants, in_groups, group_column, index_column, reference_json, prefix, subset,
                        threshold_common, threshold_intermediate, outgroup_file, include_protein):
    if not in_groups:
        in_groups = in_variants

    outgroup_dict, outgroups = parse_outgroups(outgroup_file)

    groups_to_get = None
    if subset:
        groups_to_get = set()
        groups_to_get.update(subset)
        groups_to_get.update(outgroups)

    group_dict = get_group_dict(in_groups, group_column, index_column, groups_to_get)

    reference_seq, feature_dict = load_feature_coordinates(reference_json)
    feature_dict = update_feature_dict(feature_dict)

    var_dict = {}
    outgroup_var_dict = {}

    with open(in_variants, 'r', newline = '') as csv_in:
        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        if index_column not in reader.fieldnames:
            logging.warning("Index column %s not found in %s" % (index_column, in_variants))

        if "nucleotide_variants" in reader.fieldnames:
            var_column = "nucleotide_variants"
        elif "nucleotide_mutations" in reader.fieldnames:
            var_column = "nucleotide_mutations"
        else:
            logging.warning("No nucleotide_variants or nucleotide_mutations columns found")
            sys.exit(-1)

        for row in reader:
            if index_column in row and var_column in row:
                index = row[index_column]
                variants = row[var_column]
                if index in group_dict:
                    group = group_dict[index]
                    if subset is None or group in subset:
                        update_var_dict(var_dict, group, variants)
                if index in outgroup_dict:
                    for lineage in outgroup_dict[index]:
                        update_var_dict(outgroup_var_dict, lineage, variants)
                elif index in group_dict and group_dict[index] in outgroup_dict:
                    for lineage in outgroup_dict[group_dict[index]]:
                        update_var_dict(outgroup_var_dict, lineage, variants)
            else:
                logging.warning("Index column or variants column not in row", row)

    #print("outgroup_var_dict", outgroup_var_dict)
    #print("var_dict", var_dict)

    for group in var_dict:
        common, intermediate = get_common_mutations(var_dict[group], min_occurance=2, threshold_common=threshold_common, threshold_intermediate=threshold_intermediate)
        ancestral = None
        #print("found", common, intermediate)
        if group in outgroup_var_dict:
            outgroup_common, outgroup_intermediate = get_common_mutations(outgroup_var_dict[group], min_occurance=1, threshold_common=threshold_common, threshold_intermediate=threshold_intermediate)
            common, ancestral = subtract_outgroup(common, outgroup_common)
        nice_common = define_mutations(common, feature_dict, reference_seq, include_protein)
        nice_intermediate = define_mutations(intermediate, feature_dict, reference_seq, include_protein)
        nice_ancestral = define_mutations(ancestral, feature_dict, reference_seq, include_protein)
        write_constellation(prefix, group, nice_common, nice_intermediate, nice_ancestral)


def main():
    args = parse_args()
    extract_definitions(args.in_variants, args.in_groups, args.group_column, args.index_column, args.reference_json,
                        args.prefix, args.subset, args.threshold_common, args.threshold_intermediate)


if __name__ == '__main__':
    main()
