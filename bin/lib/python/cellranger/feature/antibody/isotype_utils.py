#!/usr/bin/env python3
#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#

"""Utils for isotype manipulation and stats."""
from __future__ import annotations

import re
from itertools import compress
from typing import TYPE_CHECKING, Optional

import numpy as np
from martian import log_info
from scipy.stats import spearmanr

import cellranger.rna.library as rna_library

if TYPE_CHECKING:
    from cellranger.matrix import CountMatrix

# TODO: Link json_schemas from OligoDB for header validation.
ISOTYPE_HEADER = "isotype_control"
HEADER_CORRELATIONS = "id,name,pearson_r,spearman_r,spearman_pval"
# Arbitrary number of UMI required to run a feature
MIN_FILTERED_UMIS = 100


def get_isotype_ids(feature_matrix: CountMatrix) -> list[int]:
    """Check if a reference has isotypes and return their ids.

    Args:
        feature_matrix (CountMatrix): Filtered FB count matrix

    Returns:
        list[int]: List of ids for isotypes
    """
    isotype_indexes = []
    if ISOTYPE_HEADER in feature_matrix.feature_ref.all_tag_keys:
        for feature in feature_matrix.feature_ref.feature_defs:
            if re.match("true", feature.tags[ISOTYPE_HEADER], re.IGNORECASE):
                isotype_indexes.append(feature.index)

    if isotype_indexes:
        # Check that they have data to run on
        np_dense = feature_matrix.m.toarray()[isotype_indexes, :]
        keep = np.sum(np_dense, axis=1) > MIN_FILTERED_UMIS
        isotype_indexes = list(compress(isotype_indexes, keep))
    return isotype_indexes


def calculate_isotype_correlations(fb_filtered_matrix: CountMatrix) -> Optional[np.ndarray]:
    """Calculate correlation values of features with the isotypes.

    Args:
        fb_filtered_matrix (CountMatrix): Filtered matrix and ANTIBODY library specific

    Returns:
        isotype_csv_path(Optional(np.ndarray))
    """
    isotype_ids = get_isotype_ids(fb_filtered_matrix)
    if len(isotype_ids) > 1:
        # Convert to array
        np_dense = fb_filtered_matrix.m.toarray()
        np_isotypes = np_dense[isotype_ids, :]
        keep = np.where(np.sum(np_dense, axis=1) > MIN_FILTERED_UMIS)
        if len(keep[0]) > np_isotypes.shape[0]:
            # Need to filter the data because we can't run rows with only zeros in spearman
            np_dense_log = np.log10(np_dense + 1.0)[keep]
            mean_per_barcode_isotype_log = np.sum(np_isotypes, axis=0) / len(isotype_ids)

            # Calculate both pearson and spearman
            pearson_log = np.corrcoef(np_dense_log, mean_per_barcode_isotype_log, rowvar=True)
            rho, pval = spearmanr(a=np_dense_log, b=mean_per_barcode_isotype_log, axis=1)
            # We want the last row without the last column
            subset_pearson_log = pearson_log[-1, :-1]
            subset_rho = rho[-1, :-1]
            subset_pval = pval[-1, :-1]
            names = np.array(fb_filtered_matrix.feature_ref.get_feature_names())[keep[0]]
            feature_ids = np.array(
                fb_filtered_matrix.feature_ref.get_feature_ids_by_type(
                    rna_library.ANTIBODY_LIBRARY_TYPE
                )
            )[keep[0]]

            # Merge all columns
            final = np.column_stack(
                [
                    feature_ids,
                    names,
                    subset_pearson_log,
                    subset_rho,
                    subset_pval,
                ]
            )
            return final
        else:
            log_info(
                "Not enough features in the filtered matrix. Skipping antibody isotype correlations."
            )

    else:
        log_info(
            f"Isotype correlations not generated.\n"
            f"This feature requires a {ISOTYPE_HEADER} column in the feature reference to be enabled."
            f"This message can be discarded if you don't run isotypes in your experiment"
        )
    return None


def write_isotypes_to_csv(isotype_correlations_csv_path: str, isotype_correlations: np.ndarray):
    """Write out a numpy isotype array to csv path.

    Args:
        isotype_correlations_csv_path (str): Path to the outs of the correlations
        isotype_correlations (np.ndarray): isotype_correlations results array with isotype correlations.
    """
    np.savetxt(
        isotype_correlations_csv_path,
        isotype_correlations,
        delimiter=",",
        fmt="%s",
        header=HEADER_CORRELATIONS,
        comments="",
    )
