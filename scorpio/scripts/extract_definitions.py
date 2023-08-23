#!/usr/bin/env python3

import sys
import csv
import argparse
import json
import logging
from Bio.Seq import Seq
from operator import itemgetter

from .definitions import Reference, Constellation, Variant

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


def load_outgroup_json(outgroup_json, reference):
    constellation = Constellation(from_file=True, reference=reference, variants_file=outgroup_json, include_ancestral=True)

    return constellation.variants

def get_group_dict(in_variants, group_column, index_column, subset):
    group_dict = {}
    groups = set()
    group = "new_constellation"

    with open(in_variants,"r") as f:
        reader = csv.DictReader(f)
        if group_column and group_column not in reader.fieldnames:
            logging.warning("Group column %s not found in input file %s" %(group_column, in_variants))
            sys.exit(-1)
        if not group_column:
            if "group" in reader.fieldnames:
                group_column = "group"
                logging.warning("Using group column %s" % group_column)
            elif "lineage" in reader.fieldnames:
                group_column = "lineage"
                logging.warning("Using group column %s" % group_column)
            else:
                logging.warning("No group column specified - assume all one group.")
        if index_column not in reader.fieldnames:
            index_column = reader.fieldnames[0]
            logging.info("Using index column %s" % index_column)
        for row in reader:
            if group_column and group_column in row and index_column in row:
                if subset and row[group_column] not in subset:
                    logging.info("skip")
                    continue
                if row[index_column] in group_dict:
                    logging.warning("%s is a duplicate in group CSV, keeping first" % row[index_column])
                else:
                    group_dict[row[index_column]] = row[group_column]
                    groups.add(row[group_column])
            else:
                group_dict[row[index_column]] = group
                groups.add(group)

    if len(groups) == 0:
        logging.warning("Found no groups to define")
        sys.exit(-1)
    logging.info("Found %d groups" % len(groups))

    return group_dict

def merge_and_translate_nucs(list_variants, reference, include_protein, skip_translate):
    merged_list = []

    list_variants.sort(key=itemgetter(1))
    current = ["", 1, ""]
    for new in list_variants:
        # print(new, current)
        if new[1] == current[1] + len(current[0]):
            # print("merge")
            current[0] += new[0]
            current[2] += new[2]

        elif current[0] != "":
            var = translate_if_possible(current[1], current[0], current[2], reference, include_protein, skip_translate)
            merged_list.append(var)
            current = new
        else:
            current = new
    if current[0] != "":
        #print(current)
        var = translate_if_possible(current[1], current[0], current[2], reference, include_protein, skip_translate)
        #print(var)
        merged_list.append(var)

    return merged_list

def parse_row_variants(list_variants, reference, include_protein, skip_translate):
    merged_list = []
    if not list_variants:
        return merged_list

    intermediate_list = []
    for var in list_variants:
        if var.startswith("ins") and "_" in var:
            i, pos, add = var.split("_")
            var = "%s:%s+%s" % (i, pos, add)
            merged_list.append(var)
        elif var.startswith("ins") and ":" in var:
            i, pos, add = var.split(":")
            var = "%s:%s+%s" % (i, pos, add)
            merged_list.append(var)
        elif var.startswith("del"):
            var = ':'.join(var.split("_"))
            merged_list.append(var)
        else:
            try:
                var = var.replace("nuc:","")
                intermediate_list.append([var[0], int(var[1:-1]), var[-1]])
            except:
                print("could not add var %s to intermediate list" % var)
    merged_list.extend(merge_and_translate_nucs(intermediate_list, reference, include_protein, skip_translate))
    assert len(merged_list) <= len(list_variants)

    return [Variant(v, reference) for v in merged_list]


def update_var_dict(var_dict, group, variants):
    if group not in var_dict:
        var_dict[group] = {"total": 1}
    else:
        var_dict[group]["total"] += 1

    for var in variants:
        if var.name in var_dict[group]:
            var_dict[group][var.name].count += 1
        else:
            var_dict[group][var.name] = var
    return


def get_common_mutations(var_dict, min_occurance=3, threshold_common=0.98, threshold_intermediate=0.25):
    total = var_dict["total"]
    del var_dict["total"]
    sorted_tuples = sorted(var_dict.items(), key=lambda x: x[1].ref_start)
    var_dict = {k: v for k, v in sorted_tuples}

    common = []
    intermediate = []
    for var in var_dict:
        #print("var", var, var_dict[var], min_occurance, total)
        if var_dict[var].count == total:
            common.append(var_dict[var])
        elif var_dict[var].count >= min_occurance:
            var_dict[var].frequency = float(var_dict[var].count)/total
            #print("var", var, var_dict[var], frequency, threshold_common, threshold_intermediate)
            if var_dict[var].frequency >= threshold_common:
                common.append(var_dict[var])
            elif var_dict[var].frequency >= threshold_intermediate:
                intermediate.append(var_dict[var])
    return common, intermediate


def translate_if_possible(nuc_start, nuc_ref, nuc_alt, reference, include_protein=False, skip=False):
    if skip:
         return "nuc:%s%i%s" % (nuc_ref, nuc_start, nuc_alt)
    nuc_end = nuc_start + len(nuc_ref)
    nuc_start = int(nuc_start)
    nuc_end = int(nuc_end)
    #print(nuc_start, nuc_end, nuc_ref, nuc_alt)
    assert nuc_start > 10
    assert nuc_end > 10
    #print(nuc_start, nuc_end, refseq[nuc_start-1:nuc_end-1], nuc_ref, nuc_alt)
    assert reference.refseq[nuc_start-1:nuc_end-1] == nuc_ref
    query_seq = reference.refseq[:nuc_start-1] + nuc_alt + reference.refseq[nuc_end-1:]
    #print(len(refseq), len(query_seq), type(int(nuc_start-5)), type(int(nuc_end+5)))
    #print(refseq[nuc_start-5: nuc_end+5])
    #print(query_seq[nuc_start-5: nuc_end+5])

    for feature in reference.features_dict:
        if len(reference.features_dict[feature]) > 2:
            continue # ignore nsp definitions
        if reference.features_dict[feature][0] <= nuc_start <= reference.features_dict[feature][1]:
            start, end = nuc_start, nuc_end
            while (start - reference.features_dict[feature][0]) % 3 != 0:
                start -= 1
            while (end - reference.features_dict[feature][0]) % 3 != 0:
                end += 1
            ref_allele = Seq(reference.refseq[start - 1:end - 1 ]).translate()
            query_allele = Seq(query_seq[start - 1:end - 1]).translate()
            if ref_allele == query_allele:
                return "nuc:%s%i%s" % (nuc_ref, nuc_start, nuc_alt)
            aa_pos = int((start - reference.features_dict[feature][0]) / 3) + 1
            if include_protein:
                feature, aa_pos = translate_to_protein_if_possible(feature, aa_pos, reference)
            #print(start, end, ref_allele, query_allele, aa_pos, feature)
            return "%s:%s%i%s" % (feature, ref_allele, aa_pos, query_allele)
    return "nuc:%s%i%s" % (nuc_ref, nuc_start, nuc_alt)


def translate_to_protein_if_possible(cds, aa_start, reference):
    if not cds.startswith("orf"):
        return cds, aa_start

    for feature in reference.features_dict:
        if len(reference.features_dict[feature]) < 3:
            continue  # only want nsp definitions
        if reference.features_dict[feature][2] == cds:
            if reference.features_dict[feature][0] <= aa_start <= reference.features_dict[feature][1]:
                return feature, aa_start-reference.features_dict[feature][0]+1
    return cds, aa_start


def subtract_outgroup(common, outgroup_common, intermediate, reference):
    updated_common = []
    ancestral = []
    ancestral_remainder = outgroup_common
    updated_intermediate = intermediate
    for var in common:
        for comp in outgroup_common:
            if var.equals(comp, reference):
                ancestral.append(var)
                ancestral_remainder.remove(comp)
                break
        else:
            updated_common.append(var)
    for var in intermediate:
        for comp in outgroup_common:
            if var.equals(comp, reference):
                ancestral.append(var)
                ancestral_remainder.remove(comp)
                updated_intermediate.remove(var)
                break
    print("ancestral_remainder")
    for comp in ancestral_remainder:
        comp.print()
    return updated_common, ancestral, updated_intermediate


def write_constellation(prefix, group, list_variants, list_intermediates, list_ancestral):
    group_dict = {"name": group, "sites": [var.name for var in list_variants], "intermediate": ["%s:%f" %(var.name, var.frequency) for var in list_intermediates],
                  "rules": {"min_alt": max(len(list_variants) - 8, min(len(list_variants), 3)), "max_ref": 3}}
    if list_ancestral:
        group_dict["ancestral"] = [var.name for var in list_ancestral]
    with open('%s/%s.json' % (prefix, group), 'w') as outfile:
        json.dump(group_dict, outfile, indent=4)


def extract_definitions(in_variants, in_groups, group_column, index_column, reference_json, prefix, subset,
                        threshold_common, threshold_intermediate, outgroup_file, outgroup_json, include_protein, skip_translate):
    reference = Reference(reference_json)
    reference.update_features_dict()

    if not in_groups:
        in_groups = in_variants
        logging.debug("Using input file %s for groups" % in_variants)

    outgroup_dict, outgroups = parse_outgroups(outgroup_file)
    if outgroup_file:
        logging.debug("Parsed outgroups %s" % outgroups)

    if outgroup_json:
        outgroup_common = load_outgroup_json(outgroup_json, reference)
        logging.debug("Loaded outgroup with %i variants" % len(outgroup_common))

    groups_to_get = None
    if subset:
        groups_to_get = set()
        groups_to_get.update(subset)
        groups_to_get.update(outgroups)

    group_dict = get_group_dict(in_groups, group_column, index_column, groups_to_get)

    var_dict = {}
    outgroup_var_dict = {}

    logging.debug("Process file %s" % in_variants)
    with open(in_variants, 'r', newline = '') as csv_in:
        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        if index_column not in reader.fieldnames:
            index_column = reader.fieldnames[0]
            logging.info("Using index column %s" % index_column)

        if "nucleotide_variants" in reader.fieldnames:
            var_column = "nucleotide_variants"
        elif "nucleotide_mutations" in reader.fieldnames:
            var_column = "nucleotide_mutations"
        elif "mutations" in reader.fieldnames:
            var_column = "mutations"
        else:
            logging.warning("No nucleotide_variants or nucleotide_mutations columns found")
            sys.exit(-1)

        for row in reader:
            if index_column in row and var_column in row:
                index = row[index_column]
                variants = parse_row_variants(row[var_column].split("|"), reference, include_protein, skip_translate)
                assert len(variants) <= len(row[var_column].split("|"))
                #print(index, (index in group_dict), (index in outgroup_dict), (subset is None or group in subset), (index in group_dict and group_dict[index] in outgroup_dict))
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
            common, ancestral, intermediate = subtract_outgroup(common, outgroup_common, intermediate, reference)
        elif outgroup_json:
            common, ancestral, intermediate = subtract_outgroup(common, outgroup_common, intermediate, reference)
        logging.debug("Write constellations out to prefix %s" % prefix)
        write_constellation(prefix, group, common, intermediate, ancestral)


def main():
    args = parse_args()
    extract_definitions(args.in_variants, args.in_groups, args.group_column, args.index_column, args.reference_json,
                        args.prefix, args.subset, args.threshold_common, args.threshold_intermediate)


if __name__ == '__main__':
    main()
