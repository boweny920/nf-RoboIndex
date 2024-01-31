#!/usr/bin/env python
#
# Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#
"""A helper stage to determine method used for multiplexed data."""

import json

from cellranger.multi import config as multi_config

__MRO__ = """
stage MULTIPLEXING_METHOD(
    in  json multi_graph,
    out bool multiplexing_is_rtl,
    out bool multiplexing_is_cmo,
    src py   "../rna/stages/multi/determine_sample_assignments",
) using (
    mem_gb   = 1,
    threads  = 1,
    volatile = strict,
)
"""


def main(args, outs):
    if args.multi_graph is None:
        outs.multiplexing_is_rtl = False
        outs.multiplexing_is_cmo = False
        return

    # Read in the multi graph
    with open(args.multi_graph) as in_file:
        config = multi_config.CrMultiGraph.from_json_val(json.load(in_file))
        multiplexing_type = config.get_cell_multiplexing_type()
        if multiplexing_type is None:
            outs.multiplexing_is_rtl = False
            outs.multiplexing_is_cmo = False
        elif multiplexing_type.name == "CMO":
            outs.multiplexing_is_rtl = False
            outs.multiplexing_is_cmo = True
        else:
            outs.multiplexing_is_rtl = True
            outs.multiplexing_is_cmo = False
