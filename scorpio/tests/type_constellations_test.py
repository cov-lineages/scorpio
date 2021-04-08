#!/usr/bin/env python

import pytest
import os

from scorpio.scripts.type_variants import *

cwd = os.getcwd()
print("Running tests")


this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data', 'type_constellations')

variants_list = [
            {'alt_allele': 'C', 'name': 'snp:T6954C', 'ref_allele': 'T', 'ref_start': 6954, 'type': "snp"},
            {'alt_allele': 'C', 'name': 'nuc:T6954C', 'ref_allele': 'T', 'ref_start': 6954, 'type': "snp"},
            {'name': 'del:11288:9', 'length': 9, 'ref_start': 11288, 'type': "del", "ref_allele": "TCTGGTTTT"},
            {'cds': 'orf1a', 'ref_allele': 'T', 'aa_pos': 1001, 'alt_allele': 'I', 'name': 'orf1ab:T1001I',
             'type': 'aa', 'ref_start': 3266},
            {'cds': 'orf1a', 'ref_allele': 'ACT', 'aa_pos': 1001, 'alt_allele': 'del', 'name': 'orf1ab:T1001del',
             'type': 'del', 'ref_start': 3266, 'length':3},
            {'cds': 'orf1a', 'ref_allele': 'ACT', 'aa_pos': 1001, 'alt_allele': '-', 'name': 'orf1ab:T1001-',
             'type': 'del', 'ref_start': 3266, 'length': 3},
            {'cds': 'orf1a', 'ref_allele': 'ACT', 'aa_pos': 1001, 'alt_allele': 'del', 'name': '1ab:T1001del',
             'type': 'del', 'ref_start': 3266, 'length': 3},
            {'cds': 'orf1a', 'ref_allele': 'ACT', 'aa_pos': 1001, 'alt_allele': 'del', 'name': '1ab:T1001del',
             'type': 'del', 'ref_start': 3266, 'length': 3},
            {'cds': 'orf1a', 'ref_allele': 'ACTATT', 'aa_pos': 1001, 'alt_allele': '-', 'name': '1ab:TI1001-',
             'type': 'del', 'ref_start': 3266, 'length': 6}
    ]

def test_load_feature_coordinates():
        json_file = "%s/../SARS-CoV-2.json" %data_dir
        print(json_file)
        refseq, result = load_feature_coordinates(json_file)
        print(result)
        expect = {
                "orf1a":   (266, 13468),
                "orf1b":   (13468, 21555),
                "s":       (21563, 25384),
                "orf3a":   (25393, 26220),
                "e":       (26245, 26472),
                "m":       (26523, 27191),
                "orf6":    (27202, 27387),
                "orf7a":   (27394, 27759),
                "orf8":    (27894, 28259),
                "n":       (28274, 29533),
                "orf10" :  (29558, 29674)}
        assert result == expect

def test_resolve_ambiguous_cds():
        in_cds_normal = ["orf1a","orf1b","s","orf3a","e","m","orf6","orf7a","orf8","n","orf10"]
        in_aa_pos_normal = [1,1,1,1,1,1,1,1,1,1,1]

        in_cds_ambiguous = ["1a","1b","S","3a","E","M","6","7a","8","N","10","1ab","1ab","1ab"]
        in_aa_pos_ambiguous = [1,1,1,1,1,1,1,1,1,1,1,1,4402,99999]

        expect_cds = ["orf1a","orf1b","s","orf3a","e","m","orf6","orf7a","orf8","n","orf10","orf1a","orf1b","orf1ab"]
        expect_aa_pos = [1,1,1,1,1,1,1,1,1,1,1,1,1,None]

        json_file = "%s/../SARS-CoV-2.json" % data_dir
        refseq, features_dict = load_feature_coordinates(json_file)

        for i in range(len(in_cds_normal)):
            result_cds, result_aa_pos = resolve_ambiguous_cds(in_cds_normal[i], in_aa_pos_normal[i], features_dict)
            assert result_cds == expect_cds[i]
            assert result_aa_pos == expect_aa_pos[i]

        for i in range(len(in_cds_ambiguous)):
            result_cds, result_aa_pos = resolve_ambiguous_cds(in_cds_ambiguous[i], in_aa_pos_ambiguous[i], features_dict)
            assert result_cds == expect_cds[i]
            assert result_aa_pos == expect_aa_pos[i]

def test_variant_to_variant_record():
    in_variant = ["snp:T6954C", "nuc:T6954C", "del:11288:9", "aa:orf1ab:T1001I", "aa:orf1ab:T1001del",
                  "aa:orf1ab:T1001-", "aa:1ab:T1001del", "1ab:T1001del", "1ab:TI1001-"]
    expect = variants_list

    json_file = "%s/../SARS-CoV-2.json" % data_dir
    refseq, features_dict = load_feature_coordinates(json_file)

    for i in range(len(in_variant)):
        result = variant_to_variant_record(in_variant[i], refseq, features_dict)
        print(result)
        assert result == expect[i]


def main():
    test_load_feature_coordinates()


if __name__ == '__main__':
    main()


