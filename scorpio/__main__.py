#!/usr/bin/env python3

import argparse
import sys
import os
import logging

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
        title="Available subcommands", help="", metavar="", dest='command'
    )
    # _______________________________   common  _________________________________#
    common = argparse.ArgumentParser(prog=_program, add_help=False)

    io_group = common.add_argument_group('Input/output options')
    io_group.add_argument("-o", "--output", dest="output", required=False, help="Output file or path")
    io_group.add_argument("-p", "--prefix", dest="prefix", required=False, help="Output prefix. Default: scorpio")
    io_group.add_argument("--log-file", dest="log_file", metavar='<filename>', required=False,
                          help="Log file to use (otherwise uses stdout)")

    constellation_group = common.add_argument_group('Constellation options')
    constellation_group.add_argument("-c", "--constellations", dest="constellations", required=False, nargs='+',
                                     help="One or more JSON (or CSV) files specifying variants")
    constellation_group.add_argument("-n", "--names", dest="names", required=False, nargs='+',
                                     help="Names of constellations to include")
    constellation_group.add_argument("-l", "--label", dest="label", required=False,
                                     help="Use alternative label specified in JSON where possible")
    constellation_group.add_argument('--pangolin', dest='pangolin', action='store_true',
                                    help='Uses `mrca_lineage` as label and excludes constellations/data directory for'
                                         ' robustness')
    constellation_group.add_argument("--mutations", dest="mutations", required=False, nargs='+',
                                     help="Extra mutations to type")

    reference_group = common.add_argument_group('Reference options')
    reference_group.add_argument('--reference-json', dest="reference_json",
                                 help='JSON file containing keys "genome" with reference sequence '
                                 'and "proteins", "features" or "genes" with features of interest'
                                 ' and their coordinates')

    misc_group = common.add_argument_group('Misc options')
    misc_group.add_argument("--verbose", action="store_true", help="Print lots of stuff to screen")
    misc_group.add_argument("--dry-run", dest="dry_run", action="store_true", help="Quit after checking constellations and variants are AOK")
    misc_group.add_argument('-t', '--threads', action='store', dest="threads", type=int, default=1, help="Number of threads")

    # _______________________________  classify  __________________________________#

    subparser_classify = subparsers.add_parser(
        "classify",
        parents=[common],
        help="Takes a set of lineage-defining constellations with rules and classifies sequences by them. "
             "Output is a summary file with a row per sample recording its constellation classification"
    )
    subparser_classify.add_argument("-i", "--input", dest="input", required=True, help="Primary input file - an aligned FASTA")

    subparser_classify.add_argument(
        "--output-counts", dest="output_counts", action="store_true",
        help="Save a file per constellation with columns for the ref, alt, ambig and other counts as well as a True/False call for each sample"
    )
    subparser_classify.add_argument(
        "--call-all", dest="call_all", action="store_true",
        help="Allow multiple classifications and output these as a pipe separated list in the summary file"
    )
    subparser_classify.add_argument(
        "--long", dest="long", action="store_true",
        help="Write out summary file in long format, including the counts for the winning constellation"
    )
    subparser_classify.add_argument(
        "--list-incompatible", dest="list_incompatible", action="store_true",
        help="Adds column to the summary file listing incompatible lineages for the winning constellation as listed in the constellation JSON"
    )

    subparser_classify.set_defaults(func=scorpio.subcommands.classify.run)

    # _______________________________  haplotype  __________________________________#

    subparser_haplotype = subparsers.add_parser(
        "haplotype",
        parents=[common],
        help="Takes a set of constellations and writes haplotypes (either as strings or individual columns)."
             "Output is a summary file with a row per sample and columns recording haplotypes or haplotype barcodes",
    )
    subparser_haplotype.add_argument("-i", "--input", dest="input", required=True, help="Primary input file - an aligned FASTA")

    subparser_haplotype.add_argument(
        "--ref-char", dest="ref_char", required=False, default='-',
        help="Symbol to use to represent reference allele in barcode strings"
    )
    subparser_haplotype.add_argument(
        "--output-counts", dest="output_counts", action="store_true",
        help="Save a file per constellation with columns for the ref, alt, ambig and other counts as well as a True/False call for each sample"
    )
    subparser_haplotype.add_argument(
        "--append-genotypes", dest="append_genotypes", action="store_true",
        help="Output a column per variant site instead of a barcode string. If only a single constellation is being typed, this will be placed in the summary file, otherwise in a file per constellation"
    )
    subparser_haplotype.add_argument(
        "--combination", dest="combination", action="store_true",
        help="Combines the mutations for the specified constellations, and outputs a string across them all, with counts per found constellation"
    )
    subparser_haplotype.add_argument(
        "--interspersion", dest="interspersion", action="store_true",
        help="Evaluates an interspersion score (dev feature)"
    )
    subparser_haplotype.set_defaults(func=scorpio.subcommands.haplotype.run)


    # _______________________________  define  __________________________________#

    subparser_define = subparsers.add_parser(
        "define",
        parents=[common],
        help="Takes a CSV column of per sample mutations, and a column defining groups and extracts group-defining "
             "mutations",
    )
    subparser_define.add_argument("-i", "--input", dest="input", required=True, help="Primary input file - a CSV")

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
    subparser_define.add_argument(
        "--protein", dest="protein", action="store_true",
        help="Translates definition coordinates to proteins where possible"
    )
    subparser_define.add_argument(
        "--nucleotide", dest="nucleotide", action="store_true",
        help="Do not translate to amino acid coordinates"
    )

    subparser_define.set_defaults(func=scorpio.subcommands.define.run)

    # _______________________________  list  __________________________________#

    subparser_list = subparsers.add_parser(
        "list",
        parents=[common],
        help="Lists the constellations installed that would be typed/classified with the provided input options",
    )

    subparser_list.set_defaults(func=scorpio.subcommands.list.run)
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

    ## format logging
    format = '%(levelname)s: %(message)s'
    if args.verbose:
        level = logging.DEBUG
    elif args.command == 'list':
        level = logging.ERROR
    else:
        level = logging.INFO

    if args.log_file:
        logging.basicConfig(filename=args.log_file, level=level, format=format)
    else:
        logging.basicConfig(level=level, format=format)

    if not args.reference_json or not args.constellations:
        constellations_dir = constellations.__path__[0]
        reference_json = args.reference_json
        list_constellation_files = []

        constellation_subdirs = ["data", "definitions"]
        for dir in constellation_subdirs:
            data_dir = os.path.join(constellations_dir, dir)
            logging.info(f"Looking in {data_dir} for data files...")
            for r, d, f in os.walk(data_dir):
                for fn in f:
                    if fn == "SARS-CoV-2.json":
                        reference_json = os.path.join(r, fn)
                    elif args.pangolin and r.endswith('definitions') and fn.endswith(".json"):
                        list_constellation_files.append(os.path.join(r, fn))
                    elif not args.pangolin and fn.endswith(".json"):
                        list_constellation_files.append(os.path.join(r, fn))
                    elif not args.pangolin and fn.endswith(".csv"):
                        list_constellation_files.append(os.path.join(r, fn))
        if (not args.reference_json and reference_json == "") or (not args.constellations and list_constellation_files == []):
            logging.warning("""Please either provide a reference JSON and constellation definition file, or check your environment 
                to make sure that constellations has been properly installed.""")
            sys.exit(-1)
        if not args.reference_json:
            args.reference_json = reference_json
            logging.info("Found reference %s" %args.reference_json)
        if not args.constellations and args.command in ['haplotype', 'classify', 'list']:
            args.constellations = list_constellation_files
            logging.info("Found constellations:")
            for c in args.constellations:
                logging.info(c)
            logging.info("\n")

        if "call_all" in args and args.call_all and args.long and args.verbose:
            logging.info("Cannot provide long format summary file with multiple calls, ignoring --long\n")

        if "append_genotypes" in args and args.append_genotypes and not args.ref_char:
            args.ref_char = None
        elif "ref_char" in args and not args.ref_char:
            args.ref_char = '-'

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
