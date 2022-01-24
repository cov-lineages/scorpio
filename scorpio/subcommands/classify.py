#!/usr/bin/env python3

from scorpio.scripts.type_constellations import *


def run(options):
    classify_constellations(options.input,
                            options.constellations,
                            options.names,
                            options.output,
                            options.reference_json,
                            options.output_counts,
                            options.call_all,
                            options.long,
                            options.label,
                            options.list_incompatible,
                            options.mutations,
                            options.dry_run,
                            False,
                            options.threads)
