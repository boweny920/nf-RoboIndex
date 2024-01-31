#!/usr/bin/env python
#
# Copyright (c) 2021 10x Genomics, Inc. All rights reserved.
#

"""Run the core FBPCA algorithm from the filtered matrix and return the dimension reduced matrix."""


import numpy as np

import cellranger.analysis.constants as analysis_constants
import cellranger.analysis.run_fbpca as cr_fbpca
import cellranger.matrix as cr_matrix

__MRO__ = """
stage RUN_FBPCA_CORE(
    in  h5     filt_matrix_h5,
    in  int    num_pcs,
    in float   mem_gb,
    out npy    dimred_matrix,
    src py     "stages/analyzer/run_fbpca_core",
) split (
) using (
    volatile = strict,
)
"""


def split(args):
    return {"chunks": [], "join": {"__mem_gb": args.mem_gb}}


def join(args, outs, chunk_defs, chunk_outs):
    matrix = cr_matrix.CountMatrix.load_h5_file(args.filt_matrix_h5)

    # l2 norm
    matrix.m = matrix.m.astype("float64")
    cr_matrix.inplace_csc_column_normalize_l2(matrix.m)

    num_pcs = args.num_pcs
    if num_pcs is None:
        num_pcs = analysis_constants.CBC_N_COMPONENTS_DEFAULT
    dimred_matrix = cr_fbpca.fbpca_reduce_dimension(matrix, num_pcs)

    np.save(outs.dimred_matrix, dimred_matrix)
