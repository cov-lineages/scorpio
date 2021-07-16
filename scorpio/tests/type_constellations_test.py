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
            {'name': 'del:11288:9', 'length': 9, 'ref_start': 11288, 'type': "del", "ref_allele": "TCTGGTTTT",
             'space': 'nuc'},
            {'cds': 'orf1a', 'ref_allele': 'T', 'aa_pos': 1001, 'alt_allele': 'I', 'name': 'orf1ab:T1001I',
             'type': 'aa', 'ref_start': 3266, 'fuzzy': False},
            {'cds': 'orf1a', 'ref_allele': 'T', 'aa_pos': 1001, 'alt_allele': '', 'name': 'orf1ab:T1001',
             'type': 'aa', 'ref_start': 3266, 'fuzzy': True},
            {'cds': 'orf1a', 'ref_allele': 'TI', 'aa_pos': 1001, 'alt_allele': '', 'name': '1ab:TI1001',
             'type': 'aa', 'ref_start': 3266, 'fuzzy': True},
            {'cds': 'orf1a', 'ref_allele': 'T', 'aa_pos': 1001, 'alt_allele': 'del', 'name': 'orf1ab:T1001del',
             'type': 'del', 'ref_start': 3266, 'length': 1, 'space': 'aa', 'after': 'I', 'before': 'T'},
            {'cds': 'orf1a', 'ref_allele': 'T', 'aa_pos': 1001, 'alt_allele': '-', 'name': 'orf1ab:T1001-',
             'type': 'del', 'ref_start': 3266, 'length': 1, 'space': 'aa', 'after': 'I', 'before': 'T'},
            {'cds': 'orf1a', 'ref_allele': 'T', 'aa_pos': 1001, 'alt_allele': 'del', 'name': '1ab:T1001del',
             'type': 'del', 'ref_start': 3266, 'length': 1, 'space': 'aa', 'after': 'I', 'before': 'T'},
            {'cds': 'orf1a', 'ref_allele': 'T', 'aa_pos': 1001, 'alt_allele': 'del', 'name': '1ab:T1001del',
             'type': 'del', 'ref_start': 3266, 'length': 1, 'space': 'aa', 'after': 'I', 'before': 'T'},
            {'cds': 'orf1a', 'ref_allele': 'TI', 'aa_pos': 1001, 'alt_allele': '-', 'name': '1ab:TI1001-',
             'type': 'del', 'ref_start': 3266, 'length': 2, 'space': 'aa', 'after': 'Q', 'before': 'T'},
            {'cds': 'orf1a', 'ref_allele': '', 'aa_pos': 1001, 'alt_allele': 'AAT', 'name': '1ab:1001+AAT',
             'type': 'ins', 'ref_start': 3266}
    ]

json_file = "%s/../SARS-CoV-2.json" % data_dir
refseq, features_dict = load_feature_coordinates(json_file)

aa1 = {"name": "aa1", "type": "aa", "ref_start": 5, "ref_allele": "L", "alt_allele": "S", "fuzzy": False}
aa2 = {"name": "aa2", "type": "aa", "ref_start": 5, "ref_allele": "L", "alt_allele": "", "fuzzy": True}
snp1 = {"name": "snp1", "type": "snp", "ref_start": 10, "ref_allele": "T", "alt_allele": "A"}
snp2 = {"name": "snp2", "type": "snp", "ref_start": 1, "ref_allele": "A", "alt_allele": "T"}
del1 = {"name": "del1", "type": "del", "ref_start": 15, "ref_allele": "AG", "alt_allele": "", "length": 2, "space": "nuc"}
ins1 = {"name": "ins1", "type": "ins", "ref_start": 17, "ref_allele": "", "alt_allele": "AAC"}
del1a = {"name": "del1", "type": "del", "ref_start": 16, "ref_allele": "AAAAGC", "alt_allele": "", "length": 6, "space": "nuc"}
ins1a = {"name": "ins1", "type": "ins", "ref_start": 21, "ref_allele": "", "alt_allele": "AAC"}
del1b = {"name": "del1", "type": "del", "ref_start": 16, "ref_allele": "KS", "alt_allele": "-", "length": 2, "space": "aa", "after":"S", 'before': '*'}


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
                "orf10":  (29558, 29674),
                'nsp1': (1, 180),
                'nsp2': (181, 818),
                'nsp3': (819, 2763),
                'nsp4': (2764, 3263),
                'nsp5': (3264, 3569),
                'nsp6': (3570, 3859),
                'nsp7': (3860, 3942),
                'nsp8': (3943, 4140),
                'nsp9': (4141, 4253),
                'nsp10': (4254, 4392),
                'nsp11': (4392, 4392),
                'nsp12': (4393, 5324),
                'nsp13': (5325, 5925),
                'nsp14': (5926, 6452),
                'nsp15': (6453, 6798),
                'nsp16': (6799, 7096)
        }
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
    in_variant = ["snp:T6954C", "nuc:T6954C", "del:11288:9", "aa:orf1ab:T1001I", "aa:orf1ab:T1001",
                  "aa:1ab:TI1001", "aa:orf1ab:T1001del",
                  "aa:orf1ab:T1001-", "aa:1ab:T1001del", "1ab:T1001del", "1ab:TI1001-"]
    expect = variants_list

    for i in range(len(in_variant)):
        result = variant_to_variant_record(in_variant[i], refseq, features_dict)
        print(result)
        assert result == expect[i]


def test_parse_json_in():
    variants_file = "%s/lineage_X.json" % data_dir
    variant_list, name, rules, mrca_lineage, incompatible_lineages = parse_json_in(refseq, features_dict, variants_file)
    assert len(variant_list) == 24
    assert len([v for v in variant_list if v["type"] == "snp"]) == 6
    assert len([v for v in variant_list if v["type"] == "del"]) == 3
    assert len([v for v in variant_list if v["type"] == "aa"]) == 15
    assert name == "Lineage_X"
    assert rules["min_alt"] == 4
    assert rules["max_ref"] == 6
    assert rules["s:E484K"] == "alt"
    assert mrca_lineage == "B.1.1.7"
    assert incompatible_lineages == "A|B.1.351"


def test_parse_csv_in():
    variants_file = "%s/lineage_X.csv" % data_dir
    variant_list, name, rules = parse_csv_in(refseq, features_dict, variants_file)
    assert len(variant_list) == 24
    assert len([v for v in variant_list if v["type"] == "snp"]) == 6
    assert len([v for v in variant_list if v["type"] == "del"]) == 3
    assert len([v for v in variant_list if v["type"] == "aa"]) == 15
    assert name == "lineage_X"
    assert rules["s:E484K"] == "alt"


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
    rule_dict_json = {"min_alt": 4, "max_ref": 6, "s:E484K": "alt"}
    rule_dict_csv = {"s:E484K": "alt"}
    rule_dict_txt = None
    expect_rules = [rule_dict_json, rule_dict_csv, rule_dict_txt]

    results = []
    for i in range(len(in_files)):
        name, variant_list, rule_dict, mrca_lineage, incompatible_lineages = parse_variants_in(refseq, features_dict, in_files[i])
        assert expect_names[i] == name
        assert expect_rules[i] == rule_dict
        results.append(variant_list)
    assert results[0] == results[1] == results[2]


def test_call_variant_from_fasta():
    variants = [aa1, aa2, snp1, del1a, del1b]

    ref_string = "aaaattagctcgtaaaaaagctcgcaatag"
    alt_string = "aaaatcagcacgtaa------tcgcaatag"
    ambig_string = "nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn"
    oth_string = "aaaattcgcccgtaa---agctcgcaatag"

    print(Seq(ref_string).translate())
    print(Seq(alt_string.replace("-","")).translate())
    print(Seq(ambig_string).translate())
    print(Seq(oth_string.replace("-","")).translate())

    for var in variants:
        call, query_allele = call_variant_from_fasta(Seq(ref_string), var)
        assert call == "ref"
        assert query_allele == var["ref_allele"] or query_allele == 0

    for var in variants:
        call, query_allele = call_variant_from_fasta(Seq(alt_string), var)
        assert call == "alt"
        assert query_allele == var["alt_allele"] or (var["type"] == 'del' and isinstance(query_allele, int)) or ("fuzzy" in var and var["fuzzy"] and query_allele != var["ref_allele"])

    for var in variants:
        call, query_allele = call_variant_from_fasta(Seq(ambig_string), var)
        assert call == "ambig"
        assert query_allele in ["X", "N", "NN", "--"]

    for var in variants:
        call, query_allele = call_variant_from_fasta(Seq(oth_string), var, oth_char="X")
        assert call == "oth" or (call == "alt" and var["fuzzy"])
        assert query_allele == "X" or var["fuzzy"]
        call, query_allele = call_variant_from_fasta(Seq(oth_string), var)
        assert call == "oth" or (call == "alt" and var["fuzzy"])
        assert query_allele != var["alt_allele"] and query_allele != var["ref_allele"]

    for seq in [Seq(ref_string), Seq(alt_string), Seq(oth_string)]:
        call, query_allele = call_variant_from_fasta(seq, ins1)
        assert call == "oth"
        assert query_allele == "?"


def test_count_and_classify():
    variants = [aa1, aa2, snp1, snp2, del1, ins1]

    ref_string = "aaaattagctcgtaagctcgcaatag"
    alt_string = "aaaatcagcacgta--ctcgcaatag"
    alt_plus_string = "taaatcagcacgta--ctcgcaatag"
    oth_string = "gaaattcgcccgta-gctcgcaatag"
    seqs = [Seq(ref_string), Seq(alt_string), Seq(alt_plus_string), Seq(oth_string)]

    rules = {"min_alt": 1, "max_ref": 1, "snp2": "alt"}
    expect_classify = [False, False, True, False]
    expect_counts = [{"ref": 5, "alt": 0, "ambig": 0, "oth": 1, "rules": 0, "support": 0.0, "conflict": 0.8333},
                     {"ref": 1, "alt": 4, "ambig": 0, "oth": 1, "rules": 0, "support": 0.6667, "conflict": 0.1667},
                     {"ref": 0, "alt": 5, "ambig": 0, "oth": 1, "rules": 3, "support": 0.8333, "conflict": 0.0},
                     {"ref": 0, "alt": 1, "ambig": 0, "oth": 5, "rules": 0, "support": 0.1667, "conflict": 0.0}]
    for i in range(len(seqs)):
        counts, classify = count_and_classify(seqs[i], variants, rules)
        print(i, counts, classify)
        assert classify == expect_classify[i]
        assert counts == expect_counts[i]

def test_generate_barcode():
    variants = [aa1, aa2, snp1, snp2, del1a, ins1a]

    ref_string = "aaaattagctcgtaaaaaagctcgcaatag"
    alt_string = "aaaatcagcacgtaa------tcgcaatag"
    alt_plus_string = "taaatcagcacgtaa------tcgcaatag"
    oth_string = "gaaattcgcccgtaa---agctcgcaatag"
    seqs = [Seq(ref_string), Seq(alt_string), Seq(alt_plus_string), Seq(oth_string)]

    expect_barcode_dash = ["-----?", "SSA-2?", "SSAT2?", "XFXXX?"]
    expect_barcode_ref = ["LLTA0?", "SSAA2?", "SSAT2?", "XFXXX?"]
    expect_barcode_ref_oth = ["LLTA0?", "SSAA2?", "SSAT2?", "FFCGX?"]
    expect_barcode_ins = ["-----$", "SSA-2$", "SSAT2$", "XFXXX$"]
    expect_counts = [{"ref": 5, "alt": 0, "ambig": 0, "oth": 1, "support": 0.0, "conflict": 0.8333},
                     {"ref": 1, "alt": 4, "ambig": 0, "oth": 1, "support": 0.6667,
                      "conflict": 0.1667},
                     {"ref": 0, "alt": 5, "ambig": 0, "oth": 1, "support": 0.8333,
                      "conflict": 0.0},
                     {"ref": 0, "alt": 1, "ambig": 0, "oth": 5, "support": 0.1667, "conflict": 0.0}]

    for i in range(len(seqs)):
        barcode_list, counts = generate_barcode(seqs[i], variants, ref_char="-", ins_char="?", oth_char="X")
        barcode = ''.join(barcode_list)
        print(i, barcode, counts)
        assert barcode == expect_barcode_dash[i]
        assert counts == expect_counts[i]

        barcode_list, counts = generate_barcode(seqs[i], variants, ref_char=None, ins_char="?", oth_char="X")
        barcode = ''.join(barcode_list)
        print(i, barcode, counts)
        assert barcode == expect_barcode_ref[i]

        barcode_list, counts = generate_barcode(seqs[i], variants, ref_char=None, ins_char="?", oth_char=None)
        barcode = ''.join(barcode_list)
        print(i, barcode, counts)
        assert barcode == expect_barcode_ref_oth[i]

        barcode_list, counts = generate_barcode(seqs[i], variants, ref_char="-", ins_char="$", oth_char="X")
        barcode = ''.join(barcode_list)
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

def test_classify_constellations_incompatible():
    in_fasta = "%s/test.fa" % data_dir
    list_constellation_files = ["%s/lineage_X.json" % data_dir]
    constellation_names = None
    out_csv = "%s/tmp.classify_constellations.incompatible.csv" % data_dir
    output_counts = False
    classify_constellations(in_fasta, list_constellation_files, constellation_names, out_csv, json_file,
                            output_counts=output_counts, list_incompatible=True)

    expected = "%s/expected.classified.incompatible.csv" % data_dir
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


