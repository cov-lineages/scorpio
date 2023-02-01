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
sorted(list_constellation_files)

input = "%s/alignment.fa" % data_dir

default = {"names": [],
           "output_counts": False,
           "call_all": False,
           "long": False,
           "label": None,
           "list_incompatible": False,
           "mutations":[],
           "dry_run": False,
           "interspersion": False,
           "threads": 1
           }

def test_classify_basic():
    run_name = "basic_classify"
    out_csv = "%s/classify/%s.csv" %(data_dir, run_name)

    classify_constellations(input,
                        list_constellation_files,
                        default["names"],
                        out_csv,
                        reference_json,
                        default["output_counts"],
                        default["call_all"],
                        default["long"],
                        default["label"],
                        default["list_incompatible"],
                        default["mutations"],
                        default["dry_run"],
                        default["interspersion"],
                        default["threads"])

    expected = "%s/classify/expected.%s.csv" % (data_dir, run_name)
    assert filecmp.cmp(out_csv, expected, shallow=False)
    os.unlink(out_csv)

def test_classify_dry_run():
    run_name = "basic_dry_run"
    out_csv = "%s/classify/%s.csv" %(data_dir, run_name)

    classify_constellations(input,
                        list_constellation_files,
                        default["names"],
                        out_csv,
                        reference_json,
                        default["output_counts"],
                        default["call_all"],
                        default["long"],
                        default["label"],
                        default["list_incompatible"],
                        default["mutations"],
                        True,
                        default["interspersion"],
                        default["threads"])

    outputs = glob.glob("%s/classify/expected.%s*.csv" % (data_dir, run_name))
    assert len(outputs) == 0

def test_classify_call_all():
    run_name = "call_all_classify"
    out_csv = "%s/classify/%s.csv" %(data_dir, run_name)

    classify_constellations(input,
                        list_constellation_files,
                        default["names"],
                        out_csv,
                        reference_json,
                        default["output_counts"],
                        True,
                        default["long"],
                        default["label"],
                        default["list_incompatible"],
                        default["mutations"],
                        default["dry_run"],
                        default["interspersion"],
                        default["threads"])

    expected = "%s/classify/expected.%s.csv" % (data_dir, run_name)
    assert filecmp.cmp(out_csv, expected, shallow=False)
    os.unlink(out_csv)

def test_classify_mrca_lineage():
    run_name = "label_classify"
    out_csv = "%s/classify/%s.csv" %(data_dir, run_name)

    classify_constellations(input,
                            list_constellation_files,
                            default["names"],
                            out_csv,
                            reference_json,
                            default["output_counts"],
                            default["call_all"],
                            default["long"],
                            "mrca_lineage",
                            default["list_incompatible"],
                            default["mutations"],
                            default["dry_run"],
                            default["interspersion"],
                            default["threads"])

    expected = "%s/classify/expected.%s.csv" % (data_dir, run_name)
    assert filecmp.cmp(out_csv, expected, shallow=False)
    os.unlink(out_csv)

def test_classify_mrca_lineage_call_all():
    run_name = "mrca_lineage_call_all_classify"
    out_csv = "%s/classify/%s.csv" %(data_dir, run_name)

    classify_constellations(input,
                            list_constellation_files,
                            default["names"],
                            out_csv,
                            reference_json,
                            default["output_counts"],
                            True,
                            default["long"],
                            "mrca_lineage",
                            default["list_incompatible"],
                            default["mutations"],
                            default["dry_run"],
                            default["interspersion"],
                            default["threads"])

    expected = "%s/classify/expected.%s.csv" % (data_dir, run_name)
    assert filecmp.cmp(out_csv, expected, shallow=False)
    os.unlink(out_csv)

def test_classify_who():
    run_name = "who_label_classify"
    out_csv = "%s/classify/%s.csv" %(data_dir, run_name)

    classify_constellations(input,
                            list_constellation_files,
                            default["names"],
                            out_csv,
                            reference_json,
                            default["output_counts"],
                            default["call_all"],
                            default["long"],
                            "WHO_label",
                            default["list_incompatible"],
                            default["mutations"],
                            default["dry_run"],
                            default["interspersion"],
                            default["threads"])

    expected = "%s/classify/expected.%s.csv" % (data_dir, run_name)
    assert filecmp.cmp(out_csv, expected, shallow=False)
    os.unlink(out_csv)

def test_classify_names():
    run_name = "names_classify"
    out_csv = "%s/classify/%s.csv" %(data_dir, run_name)
    names = ['B.1.617.1-like', 'Delta (AY.4.2-like)', 'Omicron (BA.1-like)', 'Zeta (P.2-like)']

    classify_constellations(input,
                            list_constellation_files,
                            names,
                            out_csv,
                            reference_json,
                            default["output_counts"],
                            default["call_all"],
                            default["long"],
                            default["label"],
                            default["list_incompatible"],
                            default["mutations"],
                            default["dry_run"],
                            default["interspersion"],
                            default["threads"])

    expected = "%s/classify/expected.%s.csv" % (data_dir, run_name)
    assert filecmp.cmp(out_csv, expected, shallow=False)
    os.unlink(out_csv)

def test_classify_names_including_unassigned():
    run_name = "names_classify_label"
    out_csv = "%s/classify/%s.csv" %(data_dir, run_name)
    names = ['B.1.617.1-like', 'Delta (AY.4.2-like)', 'Omicron (BA.1-like)', 'Zeta (P.2-like)', 'Omicron (Unassigned)']

    classify_constellations(input,
                            list_constellation_files,
                            names,
                            out_csv,
                            reference_json,
                            default["output_counts"],
                            default["call_all"],
                            default["long"],
                            "mrca_lineage",
                            default["list_incompatible"],
                            default["mutations"],
                            default["dry_run"],
                            default["interspersion"],
                            default["threads"])

    expected = "%s/classify/expected.%s.csv" % (data_dir, run_name)
    assert filecmp.cmp(out_csv, expected, shallow=False)
    os.unlink(out_csv)

def test_classify_single_output_counts():
    run_name = "single_output_counts_classify"
    out_csv = "%s/classify/%s.csv" %(data_dir, run_name)
    names = ['Omicron (BA.1-like)']

    classify_constellations(input,
                            list_constellation_files,
                            names,
                            out_csv,
                            reference_json,
                            True,
                            default["call_all"],
                            default["long"],
                            default["label"],
                            default["list_incompatible"],
                            default["mutations"],
                            default["dry_run"],
                            default["interspersion"],
                            default["threads"])

    expected = "%s/classify/expected.%s.csv" % (data_dir, run_name)
    assert filecmp.cmp(out_csv, expected, shallow=False)
    os.unlink(out_csv)

def test_classify_multi_output_counts():
    run_name = "multi_output_counts_classify"
    out_csv = "%s/classify/%s.csv" %(data_dir, run_name)
    names = ['B.1.617.1-like', 'Delta (AY.4.2-like)', 'Omicron (BA.1-like)', 'Zeta (P.2-like)']

    classify_constellations(input,
                            list_constellation_files,
                            names,
                            out_csv,
                            reference_json,
                            True,
                            default["call_all"],
                            default["long"],
                            default["label"],
                            default["list_incompatible"],
                            default["mutations"],
                            default["dry_run"],
                            default["interspersion"],
                            default["threads"])

    for expected in glob.glob("%s/classify/expected.%s*.csv" % (data_dir, run_name)):
        outfile = expected.replace("expected.", "")
        assert filecmp.cmp(outfile, expected, shallow=False)
        os.unlink(outfile)

def test_classify_multi_output_counts_call_all():
    run_name = "multi_output_counts_call_all_classify"
    comparison_run_name = "multi_output_counts_classify"
    out_csv = "%s/classify/%s.csv" %(data_dir, run_name)

    classify_constellations(input,
                            list_constellation_files,
                            default["names"],
                            out_csv,
                            reference_json,
                            True,
                            True,
                            default["long"],
                            default["label"],
                            default["list_incompatible"],
                            default["mutations"],
                            default["dry_run"],
                            default["interspersion"],
                            default["threads"])

    for expected in glob.glob("%s/classify/expected.%s.*.csv" % (data_dir, comparison_run_name)):
        outfile = expected.replace("expected.", "").replace("output_counts", "output_counts_call_all")
        print(expected, outfile)
        assert filecmp.cmp(outfile, expected, shallow=False)
        os.unlink(outfile)

    expected = "%s/classify/expected.call_all_classify.csv" % data_dir
    outfile = "%s/classify/%s.csv" % (data_dir, run_name)
    assert filecmp.cmp(outfile, expected, shallow=False)
    os.unlink(outfile)

    for outfile in glob.glob("%s/classify/%s.*.csv" % (data_dir, run_name)):
        os.unlink(outfile)

def test_classify_multi_mutations():
    run_name = "multi_mutations_classify"
    out_csv = "%s/classify/%s.csv" %(data_dir, run_name)
    names = ['B.1.617.1-like', 'Delta (AY.4.2-like)', 'Omicron (BA.1-like)', 'Zeta (P.2-like)']
    mutations = ["ORF3a:S26L", "nuc:T17040C"]

    classify_constellations(input,
                            list_constellation_files,
                            names,
                            out_csv,
                            reference_json,
                            default["output_counts"],
                            default["call_all"],
                            default["long"],
                            default["label"],
                            default["list_incompatible"],
                            mutations,
                            default["dry_run"],
                            default["interspersion"],
                            default["threads"])

    for expected in glob.glob("%s/classify/expected.%s*.csv" % (data_dir, run_name)):
        outfile = expected.replace("expected.", "")
        assert filecmp.cmp(outfile, expected, shallow=False)
        os.unlink(outfile)

def test_classify_single_mutations():
    run_name = "single_mutations_classify"
    out_csv = "%s/classify/%s.csv" %(data_dir, run_name)
    names = ['Delta (AY.4.2-like)']
    mutations = ["ORF3a:S26L", "nuc:T17040C"]

    classify_constellations(input,
                            list_constellation_files,
                            names,
                            out_csv,
                            reference_json,
                            default["output_counts"],
                            default["call_all"],
                            default["long"],
                            default["label"],
                            default["list_incompatible"],
                            mutations,
                            default["dry_run"],
                            default["interspersion"],
                            default["threads"])

    for expected in glob.glob("%s/classify/expected.%s*.csv" % (data_dir, run_name)):
        outfile = expected.replace("expected.", "")
        assert filecmp.cmp(outfile, expected, shallow=False)
        os.unlink(outfile)

def test_classify_multi_mutations():
    run_name = "multi_mutations_classify"
    out_csv = "%s/classify/%s.csv" %(data_dir, run_name)
    names = ['B.1.617.1-like', 'Delta (AY.4.2-like)', 'Omicron (BA.1-like)', 'Zeta (P.2-like)']
    mutations = ["ORF3a:S26L", "nuc:T17040C"]

    classify_constellations(input,
                            list_constellation_files,
                            names,
                            out_csv,
                            reference_json,
                            default["output_counts"],
                            default["call_all"],
                            default["long"],
                            default["label"],
                            default["list_incompatible"],
                            mutations,
                            default["dry_run"],
                            default["interspersion"],
                            default["threads"])

    for expected in glob.glob("%s/classify/expected.%s*.csv" % (data_dir, run_name)):
        outfile = expected.replace("expected.", "")
        assert filecmp.cmp(outfile, expected, shallow=False)
        os.unlink(outfile)

def test_classify_multi_mutations_output_counts():
    run_name = "multi_mutations_output_counts_classify"
    comparison_run_name = "multi_output_counts_classify"
    out_csv = "%s/classify/%s.csv" %(data_dir, run_name)
    names = ['B.1.617.1-like', 'Delta (AY.4.2-like)', 'Omicron (BA.1-like)', 'Zeta (P.2-like)']
    mutations = ["ORF3a:S26L", "nuc:T17040C"]

    classify_constellations(input,
                            list_constellation_files,
                            names,
                            out_csv,
                            reference_json,
                            True,
                            default["call_all"],
                            default["long"],
                            default["label"],
                            default["list_incompatible"],
                            mutations,
                            default["dry_run"],
                            default["interspersion"],
                            default["threads"])

    for expected in glob.glob("%s/classify/expected.%s.*.csv" % (data_dir, comparison_run_name)):
        outfile = expected.replace("expected.", "").replace("multi", "multi_mutations")
        assert filecmp.cmp(outfile, expected, shallow=False)
        os.unlink(outfile)

    expected = "%s/classify/expected.%s.csv" % (data_dir, run_name)
    outfile = expected.replace("expected.", "")
    assert filecmp.cmp(outfile, expected, shallow=False)
    os.unlink(outfile)


def test_classify_single_mutations_output_counts():
    run_name = "single_mutations_output_counts_classify"
    out_csv = "%s/classify/%s.csv" %(data_dir, run_name)
    names = ['Omicron (BA.1-like)']
    mutations = ["ORF3a:S26L", "nuc:T17040C"]

    classify_constellations(input,
                            list_constellation_files,
                            names,
                            out_csv,
                            reference_json,
                            True,
                            default["call_all"],
                            default["long"],
                            default["label"],
                            default["list_incompatible"],
                            mutations,
                            default["dry_run"],
                            default["interspersion"],
                            default["threads"])

    for expected in glob.glob("%s/classify/expected.%s*.csv" % (data_dir, run_name)):
        outfile = expected.replace("expected.", "")
        assert filecmp.cmp(outfile, expected, shallow=False)
        os.unlink(outfile)

def test_classify_multi_list_incompatible_long():
    run_name = "multi_list_incompatible_long"
    out_csv = "%s/classify/%s.csv" %(data_dir, run_name)
    names = ['B.1.617.1-like', 'Delta (B.1.617.2-like) +K417N', 'Delta (AY.4-like)', 'Omicron (BA.1-like)', 'Zeta (P.2-like)']

    classify_constellations(input,
                            list_constellation_files,
                            names,
                            out_csv,
                            reference_json,
                            default["output_counts"],
                            default["call_all"],
                            True,
                            default["label"],
                            True,
                            default["mutations"],
                            default["dry_run"],
                            default["interspersion"],
                            default["threads"])

    for expected in glob.glob("%s/classify/expected.%s*.csv" % (data_dir, run_name)):
        outfile = expected.replace("expected.", "")
        assert filecmp.cmp(outfile, expected, shallow=False)
        os.unlink(outfile)

def main():
    test_classify_basic()

if __name__ == '__main__':
    main()


