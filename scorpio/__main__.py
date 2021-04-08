#!/usr/bin/env python3

import argparse
import sys
import os

import scorpio.subcommands
from scorpio import __version__
from . import _program

thisdir = os.path.abspath(os.path.dirname(__file__))

def main(args=None):
    parser = argparse.ArgumentParser(
        prog="scorpio",
        description="Miscellaneous fasta manipulation tools",
    )

    parser.add_argument("-v", "--version", action='version', version=f"scorpio {__version__}")
    subparsers = parser.add_subparsers(
        title="Available subcommands", help="", metavar=""
    )
    # _______________________________   common  _________________________________#
    common = argparse.ArgumentParser(prog=_program, add_help=False)

    io_group = common.add_argument_group('Input/output options')
    io_group.add_argument("-i", "--input", dest="input", required=True, help="Primary input file (FASTA)")
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
    subparser_classify.set_defaults(func=scorpio.subcommands.haplotype.run)

    # _______________________________  haplotype  __________________________________#

    subparser_haplotype = subparsers.add_parser(
        "haplotype",
        parents=[common],
        help="Takes a set of constellations and writes haplotypes (either as strings or individual columns)",
    )
    subparser_haplotype.add_argument(
        '--reference-json', help='JSON file containing keys "genome" with reference sequence '
                                 'and "proteins", "features" or "genes" with features of interest'
                                 ' and their coordinates'
    )
    subparser_haplotype.add_argument(
        "--ref-char", dest="ref_char", default='-', required=False,
        help="Symbol to use to represent reference allele"
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

    # _________________________________________________________________________#

    args = parser.parse_args()

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