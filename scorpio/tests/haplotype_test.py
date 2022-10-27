#!/usr/bin/env python

import pytest
import os
import filecmp
import glob

from scorpio.scripts.type_constellations import *
import constellations

cwd = os.getcwd()
this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
test_dir = os.path.join(this_dir, 'tests')
data_dir = os.path.join(this_dir, 'tests', 'data')

reference_json = "%s/SARS-CoV-2.json" % data_dir

constellations_dir = constellations.__path__[0]
list_constellation_files = []
constellation_subdirs = ["definitions"]
for dir in constellation_subdirs:
    c_dir = os.path.join(constellations_dir, dir)
    for r, d, f in os.walk(c_dir):
        for fn in f:
            if r.endswith('definitions') and fn.endswith(".json"):
                list_constellation_files.append(os.path.join(r, fn))

input = "%s/alignment.fa" % data_dir

default = {"names": [],
           "ref_char": "-",
           "output_counts": False,
           "label": None,
           "append_genotypes": False,
           "mutations":[],
           "dry_run": False,
           "combination": False,
           "interspersion": False,
           "threads": 1
           }

def test_haplotype_basic():
    run_name = "basic_haplotype"
    out_csv = "%s/haplotype/%s.csv" %(data_dir, run_name)

    type_constellations(input,
                        list_constellation_files,
                        default["names"],
                        out_csv,
                        reference_json,
                        default["ref_char"],
                        default["output_counts"],
                        default["label"],
                        default["append_genotypes"],
                        default["mutations"],
                        default["dry_run"],
                        default["combination"],
                        default["interspersion"],
                        default["threads"])

    expected = "%s/haplotype/expected.%s.csv" % (data_dir, run_name)
    assert filecmp.cmp(out_csv, expected, shallow=False)
    os.unlink(out_csv)

def test_haplotype_mrca_lineage():
    run_name = "label_haplotype"
    out_csv = "%s/haplotype/%s.csv" %(data_dir, run_name)

    type_constellations(input,
                        list_constellation_files,
                        default["names"],
                        out_csv,
                        reference_json,
                        default["ref_char"],
                        default["output_counts"],
                        "mrca_lineage",
                        default["append_genotypes"],
                        default["mutations"],
                        default["dry_run"],
                        default["combination"],
                        default["interspersion"],
                        default["threads"])

    expected = "%s/haplotype/expected.%s.csv" % (data_dir, run_name)
    assert filecmp.cmp(out_csv, expected, shallow=False)
    os.unlink(out_csv)

def test_haplotype_names():
    run_name = "names_haplotype"
    out_csv = "%s/haplotype/%s.csv" %(data_dir, run_name)
    names = ['B.1.617.1-like', 'Delta (AY.4.2-like)', 'Omicron (BA.1-like)', 'Zeta (P.2-like)']

    type_constellations(input,
                        list_constellation_files,
                        names,
                        out_csv,
                        reference_json,
                        default["ref_char"],
                        default["output_counts"],
                        default["label"],
                        default["append_genotypes"],
                        default["mutations"],
                        default["dry_run"],
                        default["combination"],
                        default["interspersion"],
                        default["threads"])

    expected = "%s/haplotype/expected.%s.csv" % (data_dir, run_name)
    assert filecmp.cmp(out_csv, expected, shallow=False)
    os.unlink(out_csv)

def test_haplotype_single_simple():
    run_name = "single_simple_haplotype"
    out_csv = "%s/haplotype/%s.csv" %(data_dir, run_name)
    names = ['Omicron (BA.1-like)']

    type_constellations(input,
                        list_constellation_files,
                        names,
                        out_csv,
                        reference_json,
                        default["ref_char"],
                        False,
                        default["label"],
                        default["append_genotypes"],
                        default["mutations"],
                        default["dry_run"],
                        default["combination"],
                        default["interspersion"],
                        default["threads"])

    expected = "%s/haplotype/expected.%s.csv" % (data_dir, run_name)
    assert filecmp.cmp(out_csv, expected, shallow=False)
    os.unlink(out_csv)

def test_haplotype_single_output_counts():
    run_name = "single_output_counts_haplotype"
    out_csv = "%s/haplotype/%s.csv" %(data_dir, run_name)
    names = ['Omicron (BA.1-like)']

    type_constellations(input,
                        list_constellation_files,
                        names,
                        out_csv,
                        reference_json,
                        default["ref_char"],
                        True,
                        default["label"],
                        default["append_genotypes"],
                        default["mutations"],
                        default["dry_run"],
                        default["combination"],
                        default["interspersion"],
                        default["threads"])

    expected = "%s/haplotype/expected.%s.csv" % (data_dir, run_name)
    assert filecmp.cmp(out_csv, expected, shallow=False)
    os.unlink(out_csv)

def test_haplotype_single_append_genotypes():
    run_name = "single_append_genotypes_haplotype"
    out_csv = "%s/haplotype/%s.csv" %(data_dir, run_name)
    names = ['Omicron (BA.1-like)']

    type_constellations(input,
                        list_constellation_files,
                        names,
                        out_csv,
                        reference_json,
                        default["ref_char"],
                        default["output_counts"],
                        default["label"],
                        True,
                        default["mutations"],
                        default["dry_run"],
                        default["combination"],
                        default["interspersion"],
                        default["threads"])

    expected = "%s/haplotype/expected.%s.csv" % (data_dir, run_name)
    assert filecmp.cmp(out_csv, expected, shallow=False)
    os.unlink(out_csv)

def test_haplotype_single_append_genotypes_output_counts():
    run_name = "single_append_genotypes_output_counts_haplotype"
    out_csv = "%s/haplotype/%s.csv" %(data_dir, run_name)
    names = ['Omicron (BA.1-like)']

    type_constellations(input,
                        list_constellation_files,
                        names,
                        out_csv,
                        reference_json,
                        default["ref_char"],
                        True,
                        default["label"],
                        True,
                        default["mutations"],
                        default["dry_run"],
                        default["combination"],
                        default["interspersion"],
                        default["threads"])

    expected = "%s/haplotype/expected.%s.csv" % (data_dir, run_name)
    assert filecmp.cmp(out_csv, expected, shallow=False)
    os.unlink(out_csv)

def test_haplotype_multi_output_counts():
    run_name = "multi_output_counts_haplotype"
    out_csv = "%s/haplotype/%s.csv" %(data_dir, run_name)
    names = ['B.1.617.1-like', 'Delta (AY.4.2-like)', 'Omicron (BA.1-like)', 'Zeta (P.2-like)']

    type_constellations(input,
                        list_constellation_files,
                        names,
                        out_csv,
                        reference_json,
                        default["ref_char"],
                        True,
                        default["label"],
                        default["append_genotypes"],
                        default["mutations"],
                        default["dry_run"],
                        default["combination"],
                        default["interspersion"],
                        default["threads"])

    for expected in glob.glob("%s/haplotype/expected.%s*.csv" % (data_dir, run_name)):
        outfile = expected.replace("expected.", "")
        assert filecmp.cmp(outfile, expected, shallow=False)
        os.unlink(outfile)

def test_haplotype_multi_append_genotypes():
    run_name = "multi_append_genotypes_haplotype"
    out_csv = "%s/haplotype/%s.csv" %(data_dir, run_name)
    names = ['B.1.617.1-like', 'Delta (AY.4.2-like)', 'Omicron (BA.1-like)', 'Zeta (P.2-like)']

    type_constellations(input,
                        list_constellation_files,
                        names,
                        out_csv,
                        reference_json,
                        default["ref_char"],
                        default["output_counts"],
                        default["label"],
                        True,
                        default["mutations"],
                        default["dry_run"],
                        default["combination"],
                        default["interspersion"],
                        default["threads"])

    for expected in glob.glob("%s/haplotype/expected.%s*.csv" % (data_dir, run_name)):
        outfile = expected.replace("expected.", "")
        assert filecmp.cmp(outfile, expected, shallow=False)
        os.unlink(outfile)

def test_haplotype_multi_append_genotypes_output_counts():
    run_name = "multi_append_genotypes_output_counts_haplotype"
    out_csv = "%s/haplotype/%s.csv" %(data_dir, run_name)
    names = ['B.1.617.1-like', 'Delta (AY.4.2-like)', 'Omicron (BA.1-like)', 'Zeta (P.2-like)']

    type_constellations(input,
                        list_constellation_files,
                        names,
                        out_csv,
                        reference_json,
                        default["ref_char"],
                        True,
                        default["label"],
                        True,
                        default["mutations"],
                        default["dry_run"],
                        default["combination"],
                        default["interspersion"],
                        default["threads"])

    for expected in glob.glob("%s/haplotype/expected.%s*.csv" % (data_dir, run_name)):
        outfile = expected.replace("expected.", "")
        assert filecmp.cmp(outfile, expected, shallow=False)
        os.unlink(outfile)

def test_haplotype_multi_mutations():
    run_name = "multi_mutations_haplotype"
    out_csv = "%s/haplotype/%s.csv" %(data_dir, run_name)
    names = ['B.1.617.1-like', 'Delta (AY.4.2-like)', 'Omicron (BA.1-like)', 'Zeta (P.2-like)']
    mutations = ["ORF3a:S26L", "nuc:T17040C"]

    type_constellations(input,
                        list_constellation_files,
                        names,
                        out_csv,
                        reference_json,
                        default["ref_char"],
                        default["output_counts"],
                        default["label"],
                        default["append_genotypes"],
                        mutations,
                        default["dry_run"],
                        default["combination"],
                        default["interspersion"],
                        default["threads"])

    for expected in glob.glob("%s/haplotype/expected.%s*.csv" % (data_dir, run_name)):
        outfile = expected.replace("expected.", "")
        assert filecmp.cmp(outfile, expected, shallow=False)
        os.unlink(outfile)

def test_haplotype_single_mutations():
    run_name = "single_mutations_haplotype"
    out_csv = "%s/haplotype/%s.csv" %(data_dir, run_name)
    names = ['Delta (AY.4.2-like)']
    mutations = ["ORF3a:S26L", "nuc:T17040C"]

    type_constellations(input,
                        list_constellation_files,
                        names,
                        out_csv,
                        reference_json,
                        default["ref_char"],
                        default["output_counts"],
                        default["label"],
                        default["append_genotypes"],
                        mutations,
                        default["dry_run"],
                        default["combination"],
                        default["interspersion"],
                        default["threads"])

    for expected in glob.glob("%s/haplotype/expected.%s*.csv" % (data_dir, run_name)):
        outfile = expected.replace("expected.", "")
        assert filecmp.cmp(outfile, expected, shallow=False)
        os.unlink(outfile)

def test_haplotype_single_mutations_output_counts():
    run_name = "single_mutations_output_counts_haplotype"
    out_csv = "%s/haplotype/%s.csv" %(data_dir, run_name)
    names = ['Delta (AY.4.2-like)']
    mutations = ["ORF3a:S26L", "nuc:T17040C"]

    type_constellations(input,
                        list_constellation_files,
                        names,
                        out_csv,
                        reference_json,
                        default["ref_char"],
                        True,
                        default["label"],
                        default["append_genotypes"],
                        mutations,
                        default["dry_run"],
                        default["combination"],
                        default["interspersion"],
                        default["threads"])

    for expected in glob.glob("%s/haplotype/expected.%s*.csv" % (data_dir, run_name)):
        outfile = expected.replace("expected.", "")
        print(outfile)
        assert filecmp.cmp(outfile, expected, shallow=False)
        os.unlink(outfile)

def test_haplotype_single_mutations_output_counts_append_genotypes():
    run_name = "single_mutations_output_counts_append_genotypes_haplotype"
    out_csv = "%s/haplotype/%s.csv" %(data_dir, run_name)
    names = ['Delta (AY.4.2-like)']
    mutations = ["ORF3a:S26L", "nuc:T17040C"]

    type_constellations(input,
                        list_constellation_files,
                        names,
                        out_csv,
                        reference_json,
                        default["ref_char"],
                        True,
                        default["label"],
                        True,
                        mutations,
                        default["dry_run"],
                        default["combination"],
                        default["interspersion"],
                        default["threads"])

    for expected in glob.glob("%s/haplotype/expected.%s*.csv" % (data_dir, run_name)):
        outfile = expected.replace("expected.", "")
        assert filecmp.cmp(outfile, expected, shallow=False)
        os.unlink(outfile)

def main():
    test_haplotype_basic()

if __name__ == '__main__':
    main()


