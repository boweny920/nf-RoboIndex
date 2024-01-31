#!/usr/bin/env python
#
# Copyright (c) 2021 10X Genomics, Inc. All rights reserved
#
"""Summarize antibody analysis."""
import os
import shutil

__MRO__ = """
stage SUMMARIZE_ANTIBODY_ANALYSIS(
    in  csv  aggregate_barcodes,
    in  bool is_antibody,
    out path antibody_analysis,
    src py   "stages/feature/summarize_antibody_analysis",
) using (
    mem_gb = 4,
)
"""

ANTIBODY_ANALYSIS_FILE_NAMES = ["aggregate_barcodes.csv"]


def main(args, outs):
    list_of_files = [
        args.aggregate_barcodes,
    ]

    os.makedirs(outs.antibody_analysis, exist_ok=True)

    for (file_path, file_name) in zip(list_of_files, ANTIBODY_ANALYSIS_FILE_NAMES):
        if file_path is None or not os.path.isfile(file_path):
            continue
        shutil.copy(file_path, os.path.join(outs.antibody_analysis, file_name))
