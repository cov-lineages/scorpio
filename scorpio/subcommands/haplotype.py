#!/usr/bin/env python3

from scorpio.scripts.type_constellations import *


def run(options):
    type_constellations(options.input,
                        options.constellations,
                        options.names,
                        options.output,
                        options.reference_json,
                        options.ref_char,
                        options.output_counts,
                        options.label,
                        options.append_genotypes,
                        options.mutations,
                        options.dry_run,
                        options.combination,
                        options.interspersion)
