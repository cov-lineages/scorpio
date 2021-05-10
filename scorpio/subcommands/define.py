#!/usr/bin/env python3

from scorpio.scripts.extract_definitions import *


def run(options):
    extract_definitions(options.input,
                        options.in_groups,
                        options.group_column,
                        options.index_column,
                        options.reference_json,
                        options.prefix,
                        options.subset,
                        options.threshold_common,
                        options.threshold_intermediate)
