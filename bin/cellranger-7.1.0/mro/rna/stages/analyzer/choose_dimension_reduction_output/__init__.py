#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#

import shutil

__MRO__ = """
stage CHOOSE_DIMENSION_REDUCTION_OUTPUT(
    in  h5[]   pca_h5_list,
    in  path[] pca_csv_list,
    out h5     pca_h5,
    out path   pca_csv,
    src py     "stages/analyzer/choose_dimension_reduction_output",
) using (
    volatile = strict,
)
"""


def main(args, outs):
    for h5, csv in zip(args.pca_h5_list, args.pca_csv_list):
        if h5 is not None and csv is not None:
            shutil.copy(h5, outs.pca_h5)
            shutil.copytree(csv, outs.pca_csv, dirs_exist_ok=True)
