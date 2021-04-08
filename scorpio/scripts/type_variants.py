#!/usr/bin/env python

import csv
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import sys
import json
import re

if sys.version_info[0] < 3:
    raise Exception("Python 3 or a more recent version is required.")


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
    for feature in ["genes"]:  # , "proteins", "features"]:
        if feature in json_dict:
            for item in json_dict[feature]:
                if "coordinates" in json_dict[feature][item]:
                    if "from" in json_dict[feature][item]["coordinates"]:
                        features_dict[item.lower()] = (json_dict[feature][item]["coordinates"]["from"],
                                                       json_dict[feature][item]["coordinates"]["to"])
                    elif "start" in json_dict[feature][item]["coordinates"]:
                        features_dict[item.lower()] = (json_dict[feature][item]["coordinates"]["start"],
                                                       json_dict[feature][item]["coordinates"]["end"])
                print("Found reference feature %s with coordinates" % item, features_dict[item.lower()])
    if len(features_dict) == 0:
        sys.stderr.write("No features (keys \"genes\", \"proteins\" or \"features\" ) provided in JSON %s " %
                         reference_json)
        sys.exit(1)

    in_json.close()
    return refseq, features_dict


def resolve_ambiguous_cds(cds, aa_pos, features_dict):
    cds = cds.lower()
    if cds in features_dict:
        return cds, aa_pos
    if cds[0] in features_dict:
        return cds, aa_pos

    if not cds.startswith("orf"):
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
    if aa_pos is None:
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
    aa:orf1ab:T1001

    to a dict
    """
    lsplit = l.split(":")
    info = {}

    if "+" in l:
        m = re.match(r'[aa:]*(?P<cds>\w+):(?P<pos>\d+)\+(?P<alt_allele>[a-zA-Z]+)', l)
        if not m:
            sys.stderr.write("couldn't parse the following string: %s\n" % l)
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
        type = "snp"
        ref_allele = lsplit[1][0]
        ref_start = int(lsplit[1][1:-1])
        alt_allele = lsplit[1][-1]
        ref_allele_check = refseq[ref_start - 1]

        if ref_allele != '?' and ref_allele != ref_allele_check:
            sys.stderr.write(
                "variants file says reference nucleotide at position %d is %s, but reference sequence has %s, "
                "context %s\n" % (ref_start, ref_allele, ref_allele_check, refseq[ref_start - 4:ref_start + 3]))
            sys.exit(1)

        info = {"name": l, "type": type, "ref_start": ref_start, "ref_allele": ref_allele, "alt_allele": alt_allele}

    elif lsplit[0] == "del":
        length = int(lsplit[2])
        info = {"name": l, "type": lsplit[0], "ref_start": int(lsplit[1]), "length": length,
                "ref_allele": refseq[int(lsplit[1]) - 1:int(lsplit[1]) + length - 1]}

    else:
        m = re.match(r'[aa:]*(?P<cds>\w+):(?P<ref_allele>[a-zA-Z-*]+)(?P<aa_pos>\d+)(?P<alt_allele>[a-zA-Z-*]*)', l)
        if not m:
            sys.stderr.write("couldn't parse the following string: %s\n" % l)
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
        if info["alt_allele"] in ['', 'del'] or '-' in info["alt_allele"]:
            info["type"] = "del"
            info["length"] = 3*len(info["ref_allele"])
            info["ref_allele"] = str(ref_allele)

    print("Found variant %s of type %s" % (info["name"], info["type"]))
    return info


def parse_variants_in(refseq, features_dict, variants_file):
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
    print("Parsing constellation file %s" % variants_file)

    if variants_file.endswith(".json"):
        in_json = open(variants_file, 'r')
        json_dict = json.load(in_json, strict=False)
        if "sites" in json_dict:
            for site in json_dict["sites"]:
                record = variant_to_variant_record(site, refseq, features_dict)
                if record != {}:
                    variant_list.append(record)
    elif variants_file.endswith(".csv"):
        with open("%s" % variants_file, 'r') as csv_in:
            reader = csv.DictReader(csv_in)
            for row in reader:
                if "id" in reader.fieldnames:
                    if ":" not in row["id"] and "gene" in reader.fieldnames:
                        var = "%s:%s" % (row["gene"], row["id"])
                    else:
                        var = row["id"]
                    record = variant_to_variant_record(var, refseq, features_dict)
                    if record != {}:
                        variant_list.append(record)
                else:
                    break
    if len(variant_list) == 0 and not variants_file.endswith(".json"):
        with open("%s" % variants_file, "r") as f:
            for line in f:
                l = line.split("#")[0].strip()  # remove comments from the line
                if len(l) > 0:  # skip blank lines (or comment only lines)
                    record = variant_to_variant_record(l, refseq, features_dict)
                    if record != {}:
                        variant_list.append(record)

    return variant_list


def call_variant(record_seq, var):
    # print("Call variant for ", var)
    call = None
    query_allele = None

    if var["type"] == "ins":
        call = "oth"
        query_allele = "?"
    elif var["type"] == "snp":
        query_allele = record_seq.upper()[var["ref_start"] - 1]
        query_allele_minus = record_seq.upper()[var["ref_start"] - 2]
        query_allele_plus = record_seq.upper()[var["ref_start"]]
        # print("Found", query_allele, query_allele_minus, query_allele_plus)
        if query_allele == var["ref_allele"]:
            call = 'ref'
        elif query_allele == var["alt_allele"]:
            call = 'alt'
        else:
            call = 'oth'

    elif var["type"] == "aa":
        try:
            query_allele = record_seq.upper()[var["ref_start"] - 1:var["ref_start"] + 2].translate()
            query_allele_minus = record_seq.upper()[var["ref_start"] - 2:var["ref_start"] + 1].translate()
            query_allele_plus = record_seq.upper()[var["ref_start"]:var["ref_start"] + 3].translate()
            # print("Found", query_allele, query_allele_minus, query_allele_plus)
            if query_allele == var["ref_allele"]:
                call = 'ref'
            elif query_allele == var["alt_allele"]:
                call = 'alt'
            else:
                call = 'oth'
        except:
            # print("Except")
            call = 'oth'
            query_allele = "X"

    elif var["type"] == "del":
        query_allele = record_seq.upper()[var["ref_start"] - 1:var["ref_start"] + var["length"] - 1]
        if query_allele == var["ref_allele"]:
            call = 'ref'
        elif query_allele == "-" * var["length"]:
            call = 'alt'
            query_allele = var["length"]
        else:
            call = 'oth'
            query_allele = "X"

    return call, query_allele


def generate_barcode(record_seq, variant_list, ref_char=None):
    barcode_list = []
    counts = {'ref': 0, 'alt': 0, 'oth': 0}

    for var in variant_list:
        call, query_allele = call_variant(record_seq, var)
        # print(var, call, query_allele)
        counts[call] += 1
        if ref_char is not None and call == 'ref':
            barcode_list.append(str(ref_char))
        else:
            barcode_list.append(str(query_allele))

    return ''.join(barcode_list), counts


def type_constellations(in_fasta, list_constellation_files, out_csv, reference_json, ref_char=None,
                        output_counts=False):
    reference_seq, features_dict = load_feature_coordinates(reference_json)
    constellation_dict = {}
    for constellation_file in list_constellation_files:
        constellation = constellation_file.split('/')[-1]
        constellation = constellation.replace(".json", "").replace(".csv", "").replace(".txt", "")
        variants = parse_variants_in(reference_seq, features_dict, constellation_file)
        print("Found file %s for constellation %s containing variants %s" % (constellation_file, constellation, ",".join([v["name"] for v in variants])))
        if len(variants) > 0:
            constellation_dict[constellation] = variants
        else:
            print("%s is not a valid constellation file" % constellation_file)
    variants_out = open(out_csv, "w")
    variants_out.write("query,%s\n" % ",".join(list(constellation_dict.keys())))

    counts_out = {}
    if output_counts:
        for constellation in constellation_dict:
            counts_out[constellation] = open("%s.counts.csv" % out_csv.replace(".csv", ""))
            counts_out[constellation].write("query,ref_count,alt_count,other_count,fraction_alt\n")

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
                    out_info = [record.id, counts['ref'], counts['alt'], counts['oth'],
                                str(round(counts['alt'] / (counts['alt'] + counts['ref'] + counts['oth']), 4))]
                    counts_out[constellation].write("%s\n" % ",".join(out_info))
                out_list.append(barcode)

            variants_out.write("%s\n" % ",".join(out_list))

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
    type_constellations(args.in_fasta, args.variants_in, args.variants_out, args.reference_json, ref_char='-')


if __name__ == '__main__':
    main()
