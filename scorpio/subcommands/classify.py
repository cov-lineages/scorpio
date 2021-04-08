#!/usr/bin/env python3

from scorpio.scripts.type_constellations import *


def run(options):
    classify_constellations(options.input,
                        options.constellations,
                        options.names,
                        options.output,
                        options.reference_json,
                        options.output_counts)
