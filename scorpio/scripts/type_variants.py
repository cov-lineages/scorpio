#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import sys
import json
import re

if sys.version_info[0] < 3:
    raise Exception("Python 3 or a more recent version is required.")


def load_feature_coordinates(json_file):
    """
    Loads a JSON file and extracts a dictionary of coordinates for features to
    be used to translate from amino acid into nucleotide coordinates

    nuc_pos is an integer which is 1-based start pos of codon
    """
    in_json = open(json_file,'r')
    json_dict = json.load(in_json, strict=False)

    if "genome" in json_dict:
        refseq = json_dict["genome"]
    else:
        sys.stderr.write("No reference sequence (key \"genome\") provided in JSON %s " % json_file)
        sys.exit(1)

    features_dict = {}
    for feature in ["genes"]:#, "proteins", "features"]:
        if feature in json_dict:
            for item in json_dict[feature]:
                print(item,json_dict[feature][item])
                if "coordinates" in json_dict[feature][item]:
                    if "from" in json_dict[feature][item]["coordinates"]:
                        features_dict[item.lower()] = (json_dict[feature][item]["coordinates"]["from"],
                                                       json_dict[feature][item]["coordinates"]["to"])
                    elif "start" in json_dict[feature][item]["coordinates"]:
                        features_dict[item.lower()] = (json_dict[feature][item]["coordinates"]["start"],
                                                       json_dict[feature][item]["coordinates"]["end"])
                print(features_dict[item.lower()])
    if len(features_dict) == 0:
        sys.stderr.write("No features (keys \"genes\", \"proteins\" or \"features\" ) provided in JSON %s " % json_file)
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

    return nuc_pos

def variant_to_variant_record(l, refseq, features_dict):
    """
    convert a variant in one of the following formats

    snp:T6954C
    nuc:T6954C
    del:11288:9
    aa:orf1ab:T1001I
    aa:orf1ab:T1001del

    to a dict
    """
    lsplit = l.split(":")

    if lsplit[0] in ["snp","nuc"]:
        type = "snp"
        ref_allele = lsplit[1][0]
        ref_start = int(lsplit[1][1:-1])
        alt_allele = lsplit[1][-1]
        ref_allele_check = refseq[ref_start - 1]

        if ref_allele != '?' and ref_allele != ref_allele_check:
            sys.stderr.write(
                "variants file says reference nucleotide at position %d is %s, but reference sequence has %s\n" % (
                ref_start, ref_allele, ref_allele_check))
            sys.exit(1)

        newrecord = {"name": l, "type": type, "ref_start": ref_start, "ref_allele": ref_allele, "alt_allele": alt_allele}
        return newrecord

    elif lsplit[0] == "del":
        length = int(lsplit[2])
        newrecord = {"name": l, "type": lsplit[0], "ref_start": int(lsplit[1]), "length": length,
                     "ref_allele": refseq[int(lsplit[1]) - 1:int(lsplit[1]) + length - 1]}
        return newrecord

    else:
        m = re.match(r'[aa:]*(?P<cds>\w+):(?P<ref_allele>[a-zA-Z-*]+)(?P<aa_pos>\d+)(?P<alt_allele>[a-zA-Z-*]+)', l)
        if not m:
            sys.stderr.write("couldn't parse the following string: %s\n" % l)
            sys.exit(1)

        info = m.groupdict()
        info["aa_pos"] = int(info["aa_pos"])
        info["name"] = "%s:%s%d%s" % (info["cds"], info["ref_allele"], info["aa_pos"], info["alt_allele"])
        info["type"] = "aa"

        cds, aa_pos = resolve_ambiguous_cds(info["cds"], info["aa_pos"], features_dict)
        ref_start = get_nuc_position_from_aa_description(cds, aa_pos, features_dict)

        print(ref_start,info["ref_allele"], refseq[ref_start - 1:ref_start - 1 + 3*len(info["ref_allele"])])
        ref_allele = Seq(refseq[ref_start - 1:ref_start - 1 + 3*len(info["ref_allele"])])
        ref_allele_check = ref_allele.translate()
        if info["ref_allele"] != '?' and info["ref_allele"] != ref_allele_check:
            sys.stderr.write("variants file says reference amino acid in CDS %s at position %d is %s, but reference sequence has %s\n" % (
                             cds, aa_pos, ref_allele, ref_allele_check))
            sys.exit(1)

        info["cds"] = cds
        info["aa_pos"] = aa_pos
        info["ref_start"] = ref_start
        if info["alt_allele"] in ['-', 'del']:
            info["type"] = "del"
            info["length"] = 3*len(info["ref_allele"])
            info["ref_allele"] = str(ref_allele)

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

    json_file requires key "sites":[]

    csv_file requires columns "id" and "gene"

    returns variant_list which is a list of dicts of snps, aas and dels,
    one dict per variant. format of subdict varies by variant type
    """
    variant_list = []

    if variants_file.endswith(".json"):
        json_dict = json.loads(json_file)
        if "sites" in json_dict:
            for site in json_dict["sites"]:
                variant_list.append(variant_to_variant_record(site, refseq, features_dict))
    elif variants_file.endswith(".csv"):
        with open(csv_file,'w') as csv_in:
            reader = csv.DictReader(csv_in)
            for row in reader:
                if "id" in fieldnames:
                    if ":" not in row["id"] and "gene" in fieldnames:
                        var = "%s:%s" % (row["gene"], row["id"])
                    else:
                        var = row["id"]
                    variant_list.append(variant_to_variant_record(var, refseq, features_dict))
                else:
                    break
    if variant_list.empty():
        with open(list_file, "r") as f:
            for line in f:
                l = line.split("#")[0].strip()  # remove comments from the line
                if len(l) > 0:  # skip blank lines (or comment only lines)
                    variant_list.append(variant_to_variant_record(l, refseq, features_dict))

    return variant_list


def type_variants(fasta_in, reference, variant_list, variants_out_handle, write_all_variants = False):
    variants_out = open(variants_out_handle, "w")
    variants_out.write("query,ref_count,alt_count,other_count,fraction_alt")
    if write_all_variants:
        variants_out.write(',')
        variants_out.write(','.join([v['name'] for v in variant_list]))
    variants_out.write("\n")

    with open(fasta_in, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            if len(record.seq) != len(reference):
                sys.stderr.write("The fasta record for %s isn't as long as the reference, is this fasta file aligned?\n" % record.id)
                sys.exit(1)

            ref_count = 0
            alt_count = 0
            oth_count = 0

            alleles_list = []

            for var in variant_list:
                if var["type"] == "snp":
                    query_allele = record.seq.upper()[var["ref_start"] - 1]
                    if query_allele == var["ref_allele"]:
                        ref_count += 1
                    elif query_allele == var["alt_allele"]:
                        alt_count += 1
                    else:
                        oth_count += 1

                    alleles_list.append(query_allele)


                if var["type"] == "aa":
                    try:
                        query_allele = record.seq.upper()[var["ref_start"] - 1:var["ref_start"] + 2].translate()
                    except:
                        oth_count += 1
                        alleles_list.append("X")
                        continue

                    if query_allele == var["ref_allele"]:
                        ref_count += 1
                    elif query_allele == var["alt_allele"]:
                        alt_count += 1
                    else:
                        oth_count += 1

                    alleles_list.append(str(query_allele))


                if var["type"] == "del":
                    query_allele = record.seq.upper()[var["ref_start"] - 1:var["ref_start"] + var["length"] - 1]
                    if query_allele == var["ref_allele"]:
                        ref_count += 1
                        alleles_list.append("ref")
                    elif query_allele == "-" * var["length"]:
                        alt_count += 1
                        alleles_list.append("del")
                    else:
                        oth_count += 1
                        alleles_list.append("X")

            out_info = [record.id,str(ref_count),str(alt_count),str(oth_count),str(round(alt_count / (alt_count + ref_count + oth_count), 4))]
            variants_out.write(",".join(out_info))

            if write_all_variants:
                variants_out.write(",")
                variants_out.write(",".join(alleles_list))

            variants_out.write("\n")

    variants_out.close()

    pass


def parse_args():
    parser = argparse.ArgumentParser(description="""Type an alignment in Wuhan-Hu-1 coordinates for variants defined in a config file

If you have a consensus fasta file containing sequences that haven't been aligned to Wuhan-Hu-1, you can make an alignment to feed to this python script using minimap2, the latest version of gofasta and the reference fasta file:

minimap2 -a -x asm5 MN908947.fa unaligned.fasta | gofasta sam toMultiAlign --reference MN908947.fa > aligned.fasta
""",
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--fasta-in', dest = 'fasta_in', help='alignment to type, in fasta format')
    parser.add_argument('--variants-config', dest = 'variants_in', help="""
    Config file containing variants to type. This can be a one-per-line-list in a text file, a json with key "sites" 
    or a csv with columns "id" and "gene"
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
    parser.add_argument('--reference_json', help='JSON file containing keys "genome" with reference sequence and '
                                                 '"proteins", "features" or "genes" with features of interest and '
                                                 'their coordinates')
    parser.add_argument('--variants-out', dest = 'variants_out', help='csv file to write')
    parser.add_argument('--append-genotypes', dest = 'append_genotypes', action = 'store_true', help='if invoked, write the genotype for each variant in the config file to the output')

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = parse_args()
    reference_seq, features_dict = load_feature_coordinates(args.reference_json)
    variant_list = parse_variants_in(reference_seq, features_dict, args.variants_in)
    type_variants(fasta_in=args.fasta_in,
                  reference=reference_seq,
                  variants_in=variant_list,
                  variants_out_handle=args.variants_out,
                  write_all_variants=args.append_genotypes)
