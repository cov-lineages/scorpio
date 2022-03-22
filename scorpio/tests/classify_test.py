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
    expected = "%s/classify/expected.multi_output_counts_classify.csv" % data_dir
    assert filecmp.cmp(out_csv, expected, shallow=False)
    os.unlink(out_csv)

def main():
    test_classify_basic()

if __name__ == '__main__':
    main()


