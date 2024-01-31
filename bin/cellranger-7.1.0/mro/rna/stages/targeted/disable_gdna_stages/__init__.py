#!/usr/bin/env python
#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#

"""Figures out if gDNA stages should be run."""

import cellranger.csv_io as cr_csv_io

__MRO__ = """
stage DISABLE_GDNA_STAGES(
    in  csv  probe_set,
    out bool disable_targeted_gdna,
    src py   "stages/targeted/disable_gdna_stages",
)
"""


def main(args, outs):
    outs.disable_targeted_gdna = (args.probe_set is None) or (
        "region" not in cr_csv_io.load_csv_columnnames(args.probe_set)
    )
