#!/usr/bin/env python
#
# Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#
"""Run FBPCA from a cellranger CountMatrix object, which should typically contains.

    one library type (e.g. Gene Expression, Peaks)
"""
from __future__ import annotations

from typing import Any, Optional

import numpy as np
from fbpca import pca

import cellranger.analysis.constants as analysis_constants
import cellranger.matrix as cr_matrix


def fbpca_reduce_dimension(matrix, dimred, seed=0):
    """Fast Randomized PCA."""
    np.random.seed(seed)
    X = matrix.m.T
    k = min((dimred, X.shape[0], X.shape[1]))
    U, s, _ = pca(X, k=k)  # Automatically centers. pylint: disable=invalid-name
    indicies = list(range(k))
    return U[:, indicies] * s[indicies]


def run_fbpca(
    matrix: cr_matrix.CountMatrix, num_pcs: Optional[int] = None
) -> tuple[np.ndarray, dict[str, Any]]:
    """Run FBPCA and returned a dimension-reduced matrix and barcode, feature info."""

    # Take features with some minimum number of total UMIs in pseudobulk across all batches
    feature_indices = np.flatnonzero(matrix.get_counts_per_feature())
    all_bcs_matrix = matrix.select_features(feature_indices)

    # filter barcodes with zero count
    bc_indices = np.flatnonzero(all_bcs_matrix.get_counts_per_bc())
    matrix = matrix.select_barcodes(bc_indices)

    # l2 norm
    matrix.m = matrix.m.astype("float64")
    cr_matrix.inplace_csc_column_normalize_l2(matrix.m)

    if num_pcs is None:
        num_pcs = analysis_constants.CBC_N_COMPONENTS_DEFAULT
    dimred_matrix = fbpca_reduce_dimension(matrix, num_pcs)

    # restore the zero count entries to the dimred matrix
    assert matrix.get_shape()[1] == dimred_matrix.shape[0]  # number of barcodes
    restored_matrix = np.zeros((all_bcs_matrix.get_shape()[1], dimred_matrix.shape[1]), dtype=float)

    for (i, bc) in enumerate(matrix.bcs):
        restored_matrix[all_bcs_matrix.bc_to_int(bc), :] = dimred_matrix[i, :]

    bc_feature_info = {
        "barcodes": all_bcs_matrix.bcs,
        "features": all_bcs_matrix.feature_ref.feature_defs,
    }

    return restored_matrix, bc_feature_info
