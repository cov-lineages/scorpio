#!/usr/bin/env python3

from scorpio.scripts.type_constellations import *


def run(options):
    list_constellations(options.constellations,
                        options.names,
                        options.reference_json,
                        options.label)
