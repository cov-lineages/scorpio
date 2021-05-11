#!/usr/bin/env python3

import argparse
import sys
import os

import constellations

import scorpio.subcommands
from scorpio import __version__
from . import _program

thisdir = os.path.abspath(os.path.dirname(__file__))

def main(sysargs = sys.argv[1:]):
    parser = argparse.ArgumentParser(
        prog="scorpio",
        description="Miscellaneous fasta manipulation tools",
    )

    parser.add_argument("-v", "--version", action='version', version=f"scorpio {__version__}")
    parser.add_argument("-cv", "--constellations-version", action='version', version=f"constellations {constellations.__version__}", help="show constellation's version number and exit")

    subparsers = parser.add_subparsers(
        title="Available subcommands", help="", metavar=""
    )
    # _______________________________   common  _________________________________#
    common = argparse.ArgumentParser(prog=_program, add_help=False)

    io_group = common.add_argument_group('Input/output options')
    io_group.add_argument("-i", "--input", dest="input", required=True, help="Primary input file")
    io_group.add_argument("-m", "--metadata", dest="metadata", required=False, help="CSV of associated metadata")

    io_group.add_argument("-o", "--output", dest="output", required=False, help="Output file or path")
    io_group.add_argument("-p", "--prefix", dest="prefix", required=False, help="Output prefix. Default: scorpio")
    io_group.add_argument("--log-file", dest="log_file", metavar='<filename>', required=False,
                          help="Log file to use (otherwise uses stdout)")
    io_group.add_argument("--config", action="store", help="Input config file", dest="config")

    constellation_group = common.add_argument_group('Constellation options')
    constellation_group.add_argument("-c", "--constellations", dest="constellations", required=False, nargs='+',
                                     help="One or more JSON (or CSV) files specifying variants")
    constellation_group.add_argument("-n", "--names", dest="names", required=False, nargs='+',
                                     help="Names of constellations to include")

    misc_group = common.add_argument_group('Misc options')
    misc_group.add_argument('--tempdir', action="store",
                            help="Specify where you want the temporary stuff to go Default: $TMPDIR")
    misc_group.add_argument("--no-temp", action="store_true", help="Output all intermediate files")
    misc_group.add_argument("--verbose", action="store_true", help="Print lots of stuff to screen")
    misc_group.add_argument('-t', '--threads', action='store', dest="threads", type=int, help="Number of threads")

    # _______________________________  classify  __________________________________#

    subparser_classify = subparsers.add_parser(
        "classify",
        parents=[common],
        help="Takes a set of lineage-defining constellations with rules and classifies sequences by them"
    )
    subparser_classify.add_argument(
        '--reference-json', dest="reference_json", help='JSON file containing keys "genome" with reference sequence '
                                 'and "proteins", "features" or "genes" with features of interest'
                                 ' and their coordinates'
    )
    subparser_classify.add_argument(
        "--output-counts", dest="output_counts", action="store_true",
        help="Save a file per constellation of ref, alt and other counts"
    )
    subparser_classify.add_argument(
        "--call-all", dest="call_all", action="store_true",
        help="Allow multiple classifications"
    )
    subparser_classify.add_argument(
        "--long", dest="long", action="store_true",
        help="Write out summary file in long format"
    )

    subparser_classify.set_defaults(func=scorpio.subcommands.classify.run)

    # _______________________________  haplotype  __________________________________#

    subparser_haplotype = subparsers.add_parser(
        "haplotype",
        parents=[common],
        help="Takes a set of constellations and writes haplotypes (either as strings or individual columns)",
    )
    subparser_haplotype.add_argument(
        '--reference-json', dest="reference_json", help='JSON file containing keys "genome" with reference sequence '
                                 'and "proteins", "features" or "genes" with features of interest'
                                 ' and their coordinates'
    )
    subparser_haplotype.add_argument(
        "--ref-char", dest="ref_char", default='-', required=False,
        help="Symbol to use to represent reference allele"
    )
    subparser_haplotype.add_argument(
        "--output-counts", dest="output_counts", action="store_true",
        help="Save a file per constellation of ref, alt and other counts"
    )
    subparser_haplotype.set_defaults(func=scorpio.subcommands.haplotype.run)

    # _______________________________  report  __________________________________#

    subparser_report = subparsers.add_parser(
        "report",
        parents=[common],
        help="Merges two fasta files avoiding duplicates based on matches to "
             "metadata (takes the one in the first file)",
    )

    subparser_report.set_defaults(func=scorpio.subcommands.report.run)

    # _______________________________  define  __________________________________#

    subparser_define = subparsers.add_parser(
        "define",
        parents=[common],
        help="Takes a CSV column of per sample mutations, and a column defining groups and extracts group-defining "
             "mutations",
    )
    subparser_define.add_argument(
        '--reference-json', dest="reference_json", help='JSON file containing keys "genome" with reference sequence '
                                                        'and "proteins", "features" or "genes" with features of interest'
                                                        ' and their coordinates'
    )
    subparser_define.add_argument(
        '--in-groups', dest='in_groups', required=False,
        help='CSV of containing sequence_name and a column defining groups ')
    subparser_define.add_argument(
        '--group-column', dest='group_column', required=False, default='lineage',
        help='Column name defining the groups')
    subparser_define.add_argument(
        '--index-column', dest='index_column', required=False, default='sequence_name',
        help='Taxon column name')
    subparser_define.add_argument("--subset", dest="subset", required=False, nargs='+',
                                     help="Names of a subset of groups to define")
    subparser_define.add_argument("--threshold-common", dest="threshold_common", required=False, type=float,
                                  default=0.98, help="Frequency of a variant within group to be considered common")
    subparser_define.add_argument("--threshold-intermediate", dest="threshold_intermediate", required=False, type=float,
                                  default=0.25, help="Frequency of a variant within group to be reported as intermediate")
    subparser_define.add_argument(
        '--outgroups', dest='outgroups', required=False,
        help='Two column CSV with group, and pipe separated list of outgroup sequence_names for that list. '
             'Assumes outgroups will be in main input CSV')

    subparser_define.set_defaults(func=scorpio.subcommands.define.run)


    # _________________________________________________________________________#

    args = parser.parse_args()
    if len(sysargs) < 1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)

    """
    Resolve output prefix and outfile
    """
    if args.prefix:
        args.prefix = os.path.abspath(args.prefix)

    if args.prefix and args.output:
        prefix_head, prefix_tail = os.path.split(args.prefix)
        output_head, output_tail = os.path.split(args.output)
        if not output_head:
            args.output = os.path.join(args.prefix, args.output)
        elif output_head != prefix_head:
            sys.exit("Output path %s specified, but does not lie in outprefix %s" %(args.output, args.prefix))

    if not args.prefix:
        args.prefix = os.getcwd()

    if not args.output:
        args.output = os.path.join(args.prefix, "scorpio.csv")

    if not os.path.exists(args.prefix):
        os.mkdir(args.prefix)

    if not args.reference_json or not args.constellations:
        constellations_dir = constellations.__path__[0]
        data_dir = os.path.join(constellations_dir, "data")
        print(f"Looking in {data_dir} for data files...")
        reference_json = args.reference_json
        list_constellation_files = []

        for r, d, f in os.walk(data_dir):
            for fn in f:
                if fn == "SARS-CoV-2.json":
                    reference_json = os.path.join(r, fn)
                elif fn.endswith(".json"):
                    list_constellation_files.append(os.path.join(r, fn))
                elif fn.endswith(".csv"):
                    list_constellation_files.append(os.path.join(r, fn))
        if (not args.reference_json and reference_json == "") or (not args.constellations and list_constellation_files == []):
            print(sfunk.cyan(
                """Please either provide a reference JSON and constellation definition file, or check your environment 
                to make sure that constellations has been properly installed."""))
            exit(1)
        if not args.reference_json:
            args.reference_json = reference_json
            print("Found reference %s" %args.reference_json)
        if not args.constellations:
            args.constellations = list_constellation_files
            print("Found constellations:")
            for c in args.constellations:
                print(c)
            print("\n")

        if args.call_all and args.long:
            print("Cannot provide long format summary file with multiple calls, ignoring --long\n")

    """
    Exit with help menu if no args supplied
    """
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()
        sys.exit(-1)


if __name__ == "__main__":
    main()