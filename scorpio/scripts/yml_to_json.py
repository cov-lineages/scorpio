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
    yml = None
    with open(yml_file, 'r') as stream:
        try:
            yml = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    sites = parse_vars(yml["variants"])
    del yml["variants"]
    yml["sites"] = sites

    if "probable" in yml["calling-definition"]:
        rules = parse_rules(yml["calling-definition"]["probable"], sites)
        del yml["calling-definition"]
        yml["rules"] = rules
    elif "confirmed" in yml["calling-definition"]:
        rules = parse_rules(yml["calling-definition"]["confirmed"], sites)
        del yml["calling-definition"]
        yml["rules"] = rules
    else:
        print("Could not parse rules from", yml["calling-definition"])

    with open(yml_file.replace("yml","json"), 'w') as outfile:
        json.dump(yml, outfile, indent=4)

def main():
    args = parse_args()
    convert(args.yml)

if __name__ == '__main__':
    main()
