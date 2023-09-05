#!/usr/bin/env python3

import sys
import csv
import json
import logging
from Bio.Seq import Seq
import re

global_aliases = {"spike": "s", "s": "spike",
                  "envelope": "e", "e": "envelope",
                  "membrane": "m", "m": "membrane",
                  "nucleocapsid": "n", "n": "nucleocapsid"}


def resolve_ambiguous_cds(cds, aa_pos, features_dict):
    cds = cds.lower()
    aa_pos = int(aa_pos)

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


class Reference:
    def __init__(self, reference_json):
        self.refseq = None
        self.features_dict = {}
        self.load_feature_coordinates(reference_json)

    def load_feature_coordinates(self, reference_json):
        """
        Loads a JSON file and extracts a dictionary of coordinates for features to
        be used to translate from amino acid into nucleotide coordinates

        nuc_pos is an integer which is 1-based start pos of codon
        """
        in_json = open(reference_json, 'r')
        json_dict = json.load(in_json, strict=False)

        if "genome" in json_dict:
            self.refseq = json_dict["genome"]
        else:
            sys.stderr.write("No reference sequence (key \"genome\") provided in JSON %s " % reference_json)
            sys.exit(1)

        for feature in ["genes", "proteins"]:  #, "features"]:
            if feature in json_dict:
                for item in json_dict[feature]:
                    name = item.lower()
                    if "name" in json_dict[feature][item]:
                        name = json_dict[feature][item]["name"].lower()
                    if name in self.features_dict or name in global_aliases and global_aliases[name] in self.features_dict:
                        continue

                    if "coordinates" in json_dict[feature][item]:
                        if "from" in json_dict[feature][item]["coordinates"]:
                            start = int(json_dict[feature][item]["coordinates"]["from"])
                            end = int(json_dict[feature][item]["coordinates"]["to"])
                        elif "start" in json_dict[feature][item]["coordinates"]:
                            start = int(json_dict[feature][item]["coordinates"]["start"])
                            end = int(json_dict[feature][item]["coordinates"]["end"])

                        if "gene" in json_dict[feature][item]:
                            self.features_dict[name] = (start, end, json_dict[feature][item]["gene"])
                        else:
                            self.features_dict[name] = (start, end)
                        logging.info("Found reference feature %s with coordinates %s" % (name, self.features_dict[name]))
        if len(self.features_dict) == 0:
            sys.stderr.write("No features (keys \"genes\", \"proteins\" or \"features\" ) provided in JSON %s " %
                             reference_json)
            sys.exit(1)

        in_json.close()

    def update_features_dict(self):
        for feature in self.features_dict:
            if len(self.features_dict[feature]) > 2:
                cds, aa_pos = resolve_ambiguous_cds(self.features_dict[feature][2], self.features_dict[feature][0], self.features_dict)
                if aa_pos:
                    self.features_dict[feature] = (aa_pos, self.features_dict[feature][1] + self.features_dict[feature][0] - aa_pos, cds)


class Variant:
    def __init__(self, variant=None, reference=None, ignore_fails=False):
        self.name = variant
        self.space = None  # aa, nuc
        self.type = None  # del, ins, snp, aa
        self.ref_start = None
        self.length = None
        self.ref_allele = None
        self.alt_allele = None

        self.cds = None
        self.aa_pos = None
        self.before = None
        self.after = None
        self.fuzzy = False

        self.count = 1
        self.frequency = 1

        self.variant_to_variant_record(variant, reference, ignore_fails)
        self.get_positions(reference)
        self.check_ref_allele(reference, ignore_fails)
        self.add_context(reference)

    def print(self):
        print("{", self.name, self.space, self.type, self.ref_start, self.length, self.ref_allele, self.alt_allele, self.cds, self.aa_pos, self.before, self.after, self.fuzzy, self.count, self.frequency, "}")

    def equals(self, other, reference):
        if self.ref_start == other.ref_start or self.aa_pos == other.aa_pos:
            if self.ref_allele == other.ref_allele and self.alt_allele == other.alt_allele:
                return True
            elif self.fuzzy or other.fuzzy:
                return False
            elif self.space != other.space:
                if self.type == "del" == other.type:
                    if self.space == "nuc" and self.length == 3*other.length:
                        return True
                    elif self.space == "aa" and self.length*3 == other.length:
                        return True
                else:
                    self_alt = self.alt_allele
                    other_alt = other.alt_allele
                    if self.space == "nuc":
                        self_alt = self.translate_alt(reference)
                    else:
                        other_alt = other.translate_alt(reference)
                    if self_alt == other_alt:
                        return True
        return False


    def get_nuc_position_from_aa_description(self, reference):
        """
        given a CDS (eg. S) and the number of an amino acid in it, get the
        1-based start position of that codon in ref coordinates

        nuc_pos is an integer which is 1-based start pos of codon
        """
        if self.aa_pos is None or self.cds not in reference.features_dict.keys() and global_aliases[
            self.cds] not in reference.features_dict.keys():
            sys.stderr.write("I don't know about cds: %s \n" % self.cds)
            sys.stderr.write("please use one of: %s" % ",".join(reference.features_dict.keys()))
            sys.exit(1)

        if self.cds in reference.features_dict:
            cds_tuple = reference.features_dict[self.cds]
        else:
            cds_tuple = reference.features_dict[global_aliases[self.cds]]
        nuc_pos = cds_tuple[0] + ((self.aa_pos - 1) * 3)

        if nuc_pos > cds_tuple[1]:
            sys.stderr.write("invalid amino acid position for cds %s : %d" % (self.cds, self.aa_pos))
            sys.exit(1)

        self.ref_start = int(nuc_pos)
        return

    def get_aa_position_from_nuc_description(self, reference):
        """
        given a nucleotide position, get the cds and amino acid position

        """
        for feature in reference.features_dict:
            if len(reference.features_dict[feature]) > 2:
                continue  # ignore nsp definitions
            if reference.features_dict[feature][0] <= self.ref_start <= self.ref_start + self.length <= \
                    reference.features_dict[feature][1]:
                start, end = self.ref_start, self.ref_start + self.length
                while (start - reference.features_dict[feature][0]) % 3 != 0:
                    start -= 1
                while (end - reference.features_dict[feature][0]) % 3 != 0:
                    end += 1

                self.cds = feature
                self.aa_pos = int((start - reference.features_dict[feature][0]) / 3) + 1

        return


    def get_positions(self, reference):
        if not self.ref_start:
            assert self.cds and self.aa_pos
            self.get_nuc_position_from_aa_description(reference)
        if not self.aa_pos:
            assert self.ref_start and self.length
            self.get_aa_position_from_nuc_description(reference)

    def translate_alt(self, reference):
        query_seq = reference.refseq[:self.ref_start - 1] + self.alt_allele + reference.refseq[self.ref_start + self.length - 1:]
        start = reference.features_dict[self.cds][0] + self.aa_pos*3
        length = self.length
        while length % 3 != 0:
            length += 1
        return Seq(query_seq[start - 1:start + length - 1]).translate()

    def check_ref_allele(self, reference, ignore_fails=True):
        assert self.space and self.type
        assert self.ref_start and self.length

        if self.type == "ins":
            return

        assert self.ref_allele is not None and self.alt_allele is not None

        if self.space == "nuc":
            ref_allele_check = reference.refseq[self.ref_start - 1:self.ref_start+self.length - 1]

            if self.ref_allele != '?' and self.ref_allele != ref_allele_check:
                sys.stderr.write(
                    "variants file says reference nucleotide at position %d is %s, but reference sequence has %s, "
                    "context %s\n" % (self.ref_start, self.ref_allele, ref_allele_check,
                                      reference.refseq[self.ref_start - 1:self.ref_start+self.length - 1]))
                self.name = None
                if not ignore_fails:
                    sys.exit(1)

        if self.space == "aa":
            ref_allele = Seq(reference.refseq[self.ref_start - 1:self.ref_start - 1 + 3 * self.length])
            ref_allele_check = ref_allele.translate()
            if self.ref_allele != '?' and self.ref_allele != ref_allele_check:
                sys.stderr.write("variants file says reference amino acid in CDS %s at position %d is %s, but "
                                 "reference sequence has %s\n" % (self.cds, self.aa_pos, self.ref_allele, ref_allele_check))
                self.name = None
                if not ignore_fails:
                    sys.exit(1)

    def add_context(self, reference):
        before = Seq(reference.refseq[self.ref_start - 4:self.ref_start - 1])
        after = Seq(reference.refseq[self.ref_start - 1 + 3 * self.length:self.ref_start - 1 + 3 * self.length + 3])
        if self.space == "aa":
            before = before.translate()
            after = after.translate()
        self.before = str(before)
        self.after = str(after)





    def variant_to_variant_record(self, l, reference, ignore_fails):
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
        # print("Parsing variant %s" %l)

        if not l:
            return

        if "#" in l:
            l = l.split("#")[0].strip()
            self.name = l
            if l == "":
                return

        if not reference:
            sys.stderr.write("No refererence\n")

        lsplit = l.split(":")

        if "+" in l:
            m = re.match(r'ins:(?P<pos>\d+)\+(?P<length>\d+)', l)
            if not m:
                m = re.match(r'[aa:]*(?P<cds>\w+):(?P<pos>\d+)\+(?P<alt_allele>[a-zA-Z]+)', l)
            if not m:
                sys.stderr.write("Warning: couldn't parse the following string: %s\n" % l)
                self.name = None
                if not ignore_fails:
                    sys.exit(1)
            info = m.groupdict()
            self.type = "ins"
            self.ref_allele = ""

            if "cds" in info and info["cds"] not in ["snp", "nuc"]:
                self.space = "aa"
                self.cds, self.aa_pos = resolve_ambiguous_cds(info["cds"], info["aa_pos"], reference.features_dict)
                self.alt_allele = info["alt_allele"]
                self.length = len(self.alt_allele)*3
                self.name = "%s:%d+%s" % (self.cds, self.aa_pos, self.alt_allele)
            else:
                self.space = "nuc"
                self.ref_start = int(info["pos"])
                if "length" in info:
                    self.length = int(info["length"])
                else:
                    self.length = len(info["alt_allele"])
                self.name = "nuc:%d+%d" % (self.ref_start, self.length)
            logging.debug("Warning: found variant of type insertion, which will be ignored during typing")
        elif lsplit[0] in ["snp", "nuc"]:
            self.space = "nuc"
            self.type = "snp"
            m = re.match(r'(?P<ref_allele>[ACGTUN]+)(?P<ref_start>\d+)(?P<alt_allele>[AGCTUN]*)', l[4:])
            if not m:
                sys.stderr.write("Warning: couldn't parse the following string: %s\n" % l)
                self.name = None
                if not ignore_fails:
                    sys.exit(1)
            info = m.groupdict()
            self.ref_allele = info["ref_allele"]
            self.ref_start = int(info["ref_start"])
            self.alt_allele = info["alt_allele"]
            self.length = len(self.ref_allele)

        elif lsplit[0] == "del":
            self.type = "del"
            self.space = "nuc"
            self.ref_start = int(lsplit[1])
            self.length = int(lsplit[2])
            self.ref_allele = reference.refseq[self.ref_start - 1:self.ref_start + self.length - 1]
            self.alt_allele = ""

        else:
            m = re.match(r'[aa:]*(?P<cds>\w+):(?P<ref_allele>[a-zA-Z-*]+)(?P<aa_pos>\d+)(?P<alt_allele>[a-zA-Z-*]*)', l)
            if not m:
                sys.stderr.write("Warning: couldn't parse the following string: %s\n" % l)
                self.name = None
                if not ignore_fails:
                    sys.exit(1)
                return

            info = m.groupdict()
            self.space = "aa"
            self.ref_allele = info["ref_allele"]
            self.alt_allele = info["alt_allele"]
            self.length = len(self.ref_allele)
            self.name = "%s:%s%s%s" % (info["cds"], info["ref_allele"], info["aa_pos"], info["alt_allele"])
            self.type = "aa"

            self.cds, self.aa_pos = resolve_ambiguous_cds(info["cds"], info["aa_pos"], reference.features_dict)

            if info["alt_allele"] in ['del', '-']:
                self.type = "del"
                self.length = len(self.ref_allele)
                self.alt_allele = ""

            elif info["alt_allele"] == '':
                self.fuzzy = True


        # print("Found variant %s of type %s" % (info["name"], info["type"]))
        return


class Constellation:
    def __init__(self, from_file=False, name=None, variants_list=None, reference=None, variants_file=None,
                 include_ancestral=False, label=None, ignore_fails=False):
        self.variants = variants_list
        if not self.variants:
            self.variants = []
        self.name = name
        self.rules = None

        self.mrca_lineage = ""
        self.incompatible_lineage_calls = ""
        self.parent_lineage = None
        self.lineage_name = None
        self.output_name = None

        if from_file:
            self.parse_variants_in(reference, variants_file, include_ancestral, label, ignore_fails)
        else:
            self.variants = sorted(self.variants, key=lambda x: int(x.ref_start))
            self.output_name = name

    def parse_name_from_file(self, constellation_file):
        name = constellation_file.split('/')[-1]
        self.name = name.replace(".json", "").replace(".csv", "").replace(".txt", "")

    def parse_json_in(self, reference, variants_file, include_ancestral=False, label=None, ignore_fails=False):
        """
        returns variants name and rules
        """
        in_json = open(variants_file, 'r')
        json_dict = json.load(in_json, strict=False)

        if "type" in json_dict and json_dict["type"] in json_dict:
            if "mrca_lineage" in json_dict[json_dict["type"]]:
                m = re.match(r'[A-Z0-9.]*', json_dict[json_dict["type"]]["mrca_lineage"])
                if not m:
                    sys.stderr.write("Warning: mrca_lineage %s not in acceptable format - ignoring\n" % json_dict[json_dict["type"]]["mrca_lineage"])
                else:
                    self.mrca_lineage = json_dict[json_dict["type"]]["mrca_lineage"]
            if "incompatible_lineage_calls" in json_dict[json_dict["type"]]:
                self.incompatible_lineage_calls = "|".join(json_dict[json_dict["type"]]["incompatible_lineage_calls"])
            if "parent_lineage" in json_dict[json_dict["type"]]:
                self.parent_lineage = json_dict[json_dict["type"]]["parent_lineage"]
            if "lineage_name" in json_dict[json_dict["type"]]:
                self.lineage_name = json_dict[json_dict["type"]]["lineage_name"]

        if "name" in json_dict:
            self.name = json_dict["name"]
        elif "label" in json_dict:
            self.name = json_dict["label"]
        else:
            self.name = self.parse_name_from_file(variants_file)

        if not self.name:
            return

        if label:
            if "type" in json_dict and json_dict["type"] in json_dict and label in json_dict[json_dict["type"]]:
                self.output_name = json_dict[json_dict["type"]][label]
        elif "label" in json_dict:
            self.output_name = json_dict["label"]
        elif self.name:
            self.output_name = self.name

        logging.debug("")
        logging.debug("Parsing constellation JSON file %s" % variants_file)

        if "sites" in json_dict:
            for site in json_dict["sites"]:
                record = Variant(site, reference, ignore_fails=ignore_fails)
                if record.name:
                    self.variants.append(record)
        if include_ancestral and "ancestral" in json_dict:
            for site in json_dict["ancestral"]:
                record = Variant(site, reference, ignore_fails=ignore_fails)
                if record.name:
                    self.variants.append(record)

        if "rules" in json_dict:
            if type(json_dict["rules"]) == dict and "default" in json_dict["rules"]:
                self.rules = json_dict["rules"]
            else:
                self.rules = {"default": json_dict["rules"]}

        in_json.close()

    def parse_csv_in(self, reference, variants_file, ignore_fails=False):
        """
        returns variants and name
        """
        compulsory = []

        self.parse_name_from_file(variants_file)

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
                return

        for row in reader:
            if ":" not in row["id"] and "gene" in reader.fieldnames:
                var = "%s:%s" % (row["gene"], row["id"])
            else:
                var = row["id"]
            record = Variant(var, reference, ignore_fails=ignore_fails)
            if record.name:
                self.variants.append(record)
            if "compulsory" in reader.fieldnames and row["compulsory"] in ["True", True, "Y", "y", "Yes", "yes", "YES"]:
                compulsory.append(record["name"])

        csv_in.close()
        if len(compulsory) > 0:
            self.rules = {"default": {}}
            for var in compulsory:
                self.rules["default"][var] = "alt"
        self.variants = sorted(self.variants, key=lambda x: int(x["ref_start"]))
        self. output_name = self.name

    def parse_textfile_in(self, reference, variants_file, ignore_fails=False):
        """
        returns variants and name
        """
        self.parse_name_from_file(variants_file)

        logging.debug("\n")
        logging.debug("Parsing constellation text file %s" % variants_file)

        with open("%s" % variants_file, "r") as f:
            for line in f:
                l = line.split("#")[0].strip()  # remove comments from the line
                if len(l) > 0:  # skip blank lines (or comment only lines)
                    record = Variant(l, reference, ignore_fails=ignore_fails)
                    if record.name:
                        self.variants.append(record)
        self.variants = sorted(self.variants, key=lambda x: int(x["ref_start"]))
        self. output_name = self.name

    def parse_variants_in(self, reference, variants_file, include_ancestral=False, label=None, ignore_fails=False):
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

        returns variants which is a list of dicts of snps, aas and dels,
        one dict per variant. format of subdict varies by variant type
        """
        if not reference:
            exit("No reference given")
        if not variants_file:
            exit("No variants file")

        if variants_file.endswith(".json"):
            self.parse_json_in(reference, variants_file, include_ancestral=include_ancestral, label=label,
                               ignore_fails=ignore_fails)
        elif variants_file.endswith(".csv"):
            self.parse_csv_in(reference, variants_file, ignore_fails=ignore_fails)

        if len(self.variants) == 0 and not variants_file.endswith(".json"):
            self.parse_textfile_in(reference, variants_file, ignore_fails=ignore_fails)

        self.variants = sorted(self.variants, key=lambda x: int(x.ref_start))