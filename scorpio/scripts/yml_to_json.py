#!/usr/bin/env python3

import csv
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import sys
import json
import re
import yaml

def parse_args():
    parser = argparse.ArgumentParser(description="""Convert PHE yml to json.""",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--yml', help='in yaml')
    args = parser.parse_args()

    return args

def parse_vars(var_dict):
    sites = []
    for var in var_dict:
        if "amino-acid-change" in var:
            sites.append("%s:%s" % (var["gene"], var["amino-acid-change"]))
        elif var["type"] == "deletion":
            sites.append("del:%i:%i" % (var["one-based-reference-position"],
                                        len(var["reference-base"]) - len(var["variant-base"])))
        elif var["type"] == "insertion":
            sites.append("nuc:%i+%s" % (var["one-based-reference-position"], var["variant-base"][1:]))
        elif var["type"] == "SNP":
            sites.append("nuc:%s%i%s" % (var["reference-base"], var["one-based-reference-position"],
                                         var["variant-base"]))
        else:
            print("Could not translate variant", var)
    return sites


def parse_rules(old_rules, sites):
    rules = {}
    rules["min_alt"] = old_rules["mutations-required"]
    rules["max_ref"] = len(sites) - old_rules["allowed-wildtype"]
    if "variants" in old_rules:
        rule_sites = parse_vars(old_rules["variants"])
        for site in rule_sites:
            rules["site"] = "alt"
            if site not in yml["sites"]:
                sites.append(site)
    return rules


def convert(yml_file):
    out_dict = {}
    with open(yml_file, 'r') as stream:
        try:
            yml = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    pango_lineage = yml["belongs-to-lineage"][0]["PANGO"]
    out_dict["label"] = "c%s" %pango_lineage
    out_dict["description"] = yml["description"]
    del yml["description"]

    out_dict["sources"] = yml["information-sources"]
    if out_dict["sources"] == [""]:
        out_dict["sources"] = []
    del yml["information-sources"]

    out_dict["type"] = "variant"
    out_dict["variant"] = {

    }
    out_dict["tags"] = [
        pango_lineage
    ]

    for key in ["phe-label", "who-label"]:
        if key in yml:
            out_dict["variant"][key] = yml[key]
            out_dict["tags"].append(yml[key])
            del yml[key]
    for key in["unique-id"]:
        if key in yml:
            out_dict["tags"].append(yml[key])
            del yml[key]
    for key in ["alternate-names"]:
        if key in yml:
            out_dict["tags"].extend(yml[key])
            del yml[key]
    out_dict["variant"]["pango-lineages"] = [pango_lineage]
    out_dict["variant"]["representative-genome"] = ""
    for key in yml["belongs-to-lineage"][0]:
        if key != "PANGO":
            out_dict["tags"].append(yml[key])
    del yml["belongs-to-lineage"]

    sites = parse_vars(yml["variants"])
    del yml["variants"]
    out_dict["sites"] = sites

    if "probable" in yml["calling-definition"]:
        rules = parse_rules(yml["calling-definition"]["probable"], sites)
        del yml["calling-definition"]
        out_dict["rules"] = rules
    elif "confirmed" in yml["calling-definition"]:
        rules = parse_rules(yml["calling-definition"]["confirmed"], sites)
        del yml["calling-definition"]
        out_dict["rules"] = rules
    else:
        print("Could not parse rules from", yml["calling-definition"])

    for key in yml:
        if yml[key] not in [None, "", []]:
            print("Ignored key:", key)

    with open("%s.json" %out_dict["label"], 'w') as outfile:
        json.dump(out_dict, outfile, indent=4)

def main():
    args = parse_args()
    convert(args.yml)

if __name__ == '__main__':
    main()
