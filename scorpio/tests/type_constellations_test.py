#!/usr/bin/env python

import pytest
import os
from Bio.Seq import Seq
import filecmp

from scorpio.scripts.type_constellations import *

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
             'type': 'del', 'ref_start': 3266, 'length': 6},
            {'cds': 'orf1a', 'ref_allele': '', 'aa_pos': 1001, 'alt_allele': 'AAT', 'name': '1ab:1001+AAT',
             'type': 'ins', 'ref_start': 3266}
    ]

json_file = "%s/../SARS-CoV-2.json" % data_dir
refseq, features_dict = load_feature_coordinates(json_file)

aa1 = {"name": "aa1", "type": "aa", "ref_start": 5, "ref_allele": "L", "alt_allele": "S"}
snp1 = {"name": "snp1", "type": "snp", "ref_start": 10, "ref_allele": "T", "alt_allele": "A"}
snp2 = {"name": "snp2", "type": "snp", "ref_start": 1, "ref_allele": "A", "alt_allele": "T"}
del1 = {"name": "del1", "type": "del", "ref_start": 15, "ref_allele": "AG", "alt_allele": "", "length": 2}
ins1 = {"name": "ins1", "type": "ins", "ref_start": 17, "ref_allele": "", "alt_allele": "AAC"}

def test_load_feature_coordinates():
        result = features_dict
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

        for i in range(len(in_cds_normal)):
            result_cds, result_aa_pos = resolve_ambiguous_cds(in_cds_normal[i], in_aa_pos_normal[i], features_dict)
            assert result_cds == expect_cds[i]
            assert result_aa_pos == expect_aa_pos[i]

        for i in range(len(in_cds_ambiguous)):
            result_cds, result_aa_pos = resolve_ambiguous_cds(in_cds_ambiguous[i], in_aa_pos_ambiguous[i], features_dict)
            assert result_cds == expect_cds[i]
            assert result_aa_pos == expect_aa_pos[i]


def test_get_nuc_position_from_aa_description():
    cds = ["s", "s", "s", "s", "s"]
    aa_pos = [614, 501, 484, 681, 439]
    expect = [23402, 23063, 23012, 23603, 22877]

    for i in range(len(cds)):
        result = get_nuc_position_from_aa_description(cds[i], aa_pos[i], features_dict)
        print(result)
        assert result == expect[i]

    exit_cds = ["s", "nonsense", "s"]
    exit_aa_pos = [None, 614, 6140]
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        for i in range(len(exit_cds)):
            get_nuc_position_from_aa_description(exit_cds[i], exit_aa_pos[i], features_dict)
            assert pytest_wrapped_e.type == SystemExit
            assert pytest_wrapped_e.value.code == 1


def test_variant_to_variant_record():
    in_variant = ["snp:T6954C", "nuc:T6954C", "del:11288:9", "aa:orf1ab:T1001I", "aa:orf1ab:T1001del",
                  "aa:orf1ab:T1001-", "aa:1ab:T1001del", "1ab:T1001del", "1ab:TI1001-"]
    expect = variants_list

    for i in range(len(in_variant)):
        result = variant_to_variant_record(in_variant[i], refseq, features_dict)
        print(result)
        assert result == expect[i]


def test_set_rules():
    in_pairs = [(None, None), (None, 4), (4, None), (4, 6)]
    expect_pairs = [(8, 1), (5, 4), (4, 5), (4, 6)]
    for i in range(len(in_pairs)):
        result = set_rules(variants_list, in_pairs[i][0], in_pairs[i][1])
        print(result)
        assert result == expect_pairs[i]


def test_parse_json_in():
    variants_file = "%s/lineage_X.json" % data_dir
    variant_list, name, min_alt, max_ref, compulsory = parse_json_in(refseq, features_dict, variants_file)
    assert len(variant_list) == 24
    assert len([v for v in variant_list if v["type"] == "snp"]) == 6
    assert len([v for v in variant_list if v["type"] == "del"]) == 3
    assert len([v for v in variant_list if v["type"] == "aa"]) == 15
    assert name == "Lineage_X"
    assert min_alt == 4
    assert max_ref == 6
    assert compulsory == ["s:E484K"]


def test_parse_csv_in():
    variants_file = "%s/lineage_X.csv" % data_dir
    variant_list, name, compulsory = parse_csv_in(refseq, features_dict, variants_file)
    assert len(variant_list) == 24
    assert len([v for v in variant_list if v["type"] == "snp"]) == 6
    assert len([v for v in variant_list if v["type"] == "del"]) == 3
    assert len([v for v in variant_list if v["type"] == "aa"]) == 15
    assert name == "lineage_X"
    assert compulsory == ["s:E484K"]


def test_parse_textfile_in():
    variants_file = "%s/lineage_X.txt" % data_dir
    variant_list, name = parse_textfile_in(refseq, features_dict, variants_file)
    assert len(variant_list) == 24
    assert len([v for v in variant_list if v["type"] == "snp"]) == 6
    assert len([v for v in variant_list if v["type"] == "del"]) == 3
    assert len([v for v in variant_list if v["type"] == "aa"]) == 15
    assert name == "lineage_X"


def test_parse_variants_in():
    in_files = ["%s/lineage_X.json" % data_dir, "%s/lineage_X.csv" % data_dir, "%s/lineage_X.txt" % data_dir]
    expect_names = ["Lineage_X", "lineage_X", "lineage_X"]
    rule_dict_json = {"Lineage_X": {"min_alt": 4, "max_ref": 6, "compulsory": ["s:E484K"]}}
    rule_dict_csv = {"lineage_X": {"min_alt": 23, "max_ref": 1, "compulsory": ["s:E484K"]}}
    rule_dict_txt = {"lineage_X": {"min_alt": 23, "max_ref": 1, "compulsory": []}}
    expect_rules = [rule_dict_json, rule_dict_csv, rule_dict_txt]

    results = []
    for i in range(len(in_files)):
        name, variant_list, rule_dict = parse_variants_in(refseq, features_dict, in_files[i], {})
        assert expect_names[i] == name
        assert expect_rules[i] == rule_dict
        results.append(variant_list)
    assert results[0] == results[1] == results[2]


def test_call_variant_from_fasta():
    variants = [aa1, snp1, del1]

    ref_string = "aaaattagctcgtaagctcgcaatag"
    alt_string = "aaaatcagcacgta--ctcgcaatag"
    oth_string = "aaaattcgcccgta-gctcgcaatag"

    for var in variants:
        call, query_allele = call_variant_from_fasta(Seq(ref_string), var)
        assert call == "ref"
        assert query_allele == var["ref_allele"] or query_allele == 0

    for var in variants:
        call, query_allele = call_variant_from_fasta(Seq(alt_string), var)
        assert call == "alt"
        assert query_allele == var["alt_allele"] or query_allele == var["length"]

    for var in variants:
        call, query_allele = call_variant_from_fasta(Seq(oth_string), var, oth_char="X")
        assert call == "oth"
        assert query_allele == "X"
        call, query_allele = call_variant_from_fasta(Seq(oth_string), var)
        assert call == "oth"
        assert query_allele != var["alt_allele"] and query_allele != var["ref_allele"]

    for seq in [Seq(ref_string), Seq(alt_string), Seq(oth_string)]:
        call, query_allele = call_variant_from_fasta(seq, ins1)
        assert call == "oth"
        assert query_allele == "?"


def test_count_and_classify():
    variants = [aa1, snp1, snp2, del1, ins1]

    ref_string = "aaaattagctcgtaagctcgcaatag"
    alt_string = "aaaatcagcacgta--ctcgcaatag"
    alt_plus_string = "taaatcagcacgta--ctcgcaatag"
    oth_string = "gaaattcgcccgta-gctcgcaatag"
    seqs = [Seq(ref_string), Seq(alt_string), Seq(alt_plus_string), Seq(oth_string)]

    rules = {"min_alt": 1, "max_ref": 1, "compulsory": ["snp2"]}
    expect_classify = [False, False, True, False]
    expect_counts = [{"ref": 4, "alt": 0, "oth": 1, "compulsory": []},
                     {"ref": 1, "alt": 3, "oth": 1, "compulsory": []},
                     {"ref": 0, "alt": 4, "oth": 1, "compulsory": ["snp2"]},
                     {"ref": 0, "alt": 0, "oth": 5, "compulsory": []}]
    for i in range(len(seqs)):
        counts, classify = count_and_classify(seqs[i], variants, rules)
        print(i, counts, classify)
        assert classify == expect_classify[i]
        assert counts == expect_counts[i]

def test_generate_barcode():
    variants = [aa1, snp1, snp2, del1, ins1]

    ref_string = "aaaattagctcgtaagctcgcaatag"
    alt_string = "aaaatcagcacgta--ctcgcaatag"
    alt_plus_string = "taaatcagcacgta--ctcgcaatag"
    oth_string = "gaaattcgcccgta-gctcgcaatag"
    seqs = [Seq(ref_string), Seq(alt_string), Seq(alt_plus_string), Seq(oth_string)]

    expect_barcode_dash = ["----?", "SA-2?", "SAT2?", "XXXX?"]
    expect_barcode_ref = ["LTA0?", "SAA2?", "SAT2?", "XXXX?"]
    expect_barcode_ref_oth = ["LTA0?", "SAA2?", "SAT2?", "FCGX?"]
    expect_barcode_ins = ["----$", "SA-2$", "SAT2$", "XXXX$"]
    expect_counts = [{"ref": 4, "alt": 0, "oth": 1},
                     {"ref": 1, "alt": 3, "oth": 1},
                     {"ref": 0, "alt": 4, "oth": 1},
                     {"ref": 0, "alt": 0, "oth": 5}]
    for i in range(len(seqs)):
        barcode, counts = generate_barcode(seqs[i], variants, ref_char="-", ins_char="?", oth_char="X")
        print(i, barcode, counts)
        assert barcode == expect_barcode_dash[i]
        assert counts == expect_counts[i]

        barcode, counts = generate_barcode(seqs[i], variants, ref_char=None, ins_char="?", oth_char="X")
        print(i, barcode, counts)
        assert barcode == expect_barcode_ref[i]

        barcode, counts = generate_barcode(seqs[i], variants, ref_char=None, ins_char="?", oth_char=None)
        print(i, barcode, counts)
        assert barcode == expect_barcode_ref_oth[i]

        barcode, counts = generate_barcode(seqs[i], variants, ref_char="-", ins_char="$", oth_char="X")
        print(i, barcode, counts)
        assert barcode == expect_barcode_ins[i]


def test_classify_constellations():
    in_fasta = "%s/test.fa" % data_dir
    list_constellation_files = ["%s/lineage_X.json" % data_dir]
    constellation_names = None
    out_csv = "%s/tmp.classify_constellations.csv" % data_dir
    output_counts = False
    classify_constellations(in_fasta, list_constellation_files, constellation_names, out_csv, json_file,
                            output_counts=output_counts)

    expected = "%s/expected.classified.csv" % data_dir
    assert filecmp.cmp(out_csv, expected, shallow=False)
    os.unlink(out_csv)

def test_type_constellations():
    in_fasta = "%s/test.fa" % data_dir
    list_constellation_files = ["%s/lineage_X.json" % data_dir]
    constellation_names = None
    out_csv = "%s/tmp.type_constellations.csv" % data_dir
    output_counts = False
    type_constellations(in_fasta, list_constellation_files, constellation_names, out_csv, json_file, ref_char='-',
                            output_counts=output_counts)

    expected = "%s/expected.typed.csv" % data_dir
    assert filecmp.cmp(out_csv, expected, shallow=False)
    os.unlink(out_csv)

def main():
    test_load_feature_coordinates()


if __name__ == '__main__':
    main()


