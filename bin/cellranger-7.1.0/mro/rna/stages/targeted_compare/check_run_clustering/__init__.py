#!/usr/bin/env python
#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#

"""Checks the matrix is at least 2x2 and disables clustering stages if that is the case."""


from cellranger.matrix import CountMatrix

__MRO__ = """
stage CHECK_RUN_CLUSTERING(
    in  h5   matrix_h5,
    out bool skip_clustering,
    src py   "stages/targeted_compare/check_run_clustering",
) using (
    volatile = strict,
)
"""


def main(args, outs):

    rows, cols, _ = CountMatrix.load_dims_from_h5(args.matrix_h5)

    outs.skip_clustering = bool(rows < 2 or cols < 2)
