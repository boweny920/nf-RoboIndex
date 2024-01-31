#!/usr/bin/env python
#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#

"""Generates targeted and parent barcode summary csv."""


import martian
import numpy as np
from six import ensure_str

import cellranger.pandas_utils as pdu
from cellranger.molecule_counter import MoleculeCounter
from cellranger.pandas_utils import FEATURE_DF_COUNT_COL, FEATURE_DF_DUP_COL, FEATURE_DF_UMI_COL
from cellranger.rna.library import GENE_EXPRESSION_LIBRARY_TYPE
from cellranger.targeted.targeted_compare_constants import (
    BOTH,
    CONTROL_LABEL,
    COUNT_MODE_ALL,
    COUNT_MODE_TARGETED,
    COUNT_MODES_TARGETED_COMPARISON,
    NEITHER,
    PARENT_ONLY,
    TARGETED_COMPARE_CELL_CALL_GRP_COL,
    TARGETED_COMPARE_LABELS,
    TARGETED_LABEL,
    TARGETED_ONLY,
)

__MRO__ = """
stage COMPARE_BARCODE_RANK_PLOTS(
    in  h5  targeted_molecule_info,
    in  h5  parent_molecule_info,
    out csv barcode_summary_csv,
    src py  "stages/targeted_compare/compare_barcode_rank_plots",
) split (
) using (
    volatile = strict,
)
"""


def split(args):
    return {
        "chunks": [],
        "join": {"__mem_gb": 4},
    }


# pylint: disable=too-many-locals
def make_merged_barcode_df(molecule_info_fns, targeted_feature_indices, count_mode):
    """Merge barcode info of targeted and parent experiments."""

    # make a dict of dataframes (targeted and parent keys)
    barcodes_per_sample = {}
    for mol_info_idx, mol_info_fn in enumerate(molecule_info_fns):

        sample_label = TARGETED_COMPARE_LABELS[mol_info_idx]

        with MoleculeCounter.open(mol_info_fn, "r") as mc:
            gex_library_indices = mc.get_library_indices_by_type()[GENE_EXPRESSION_LIBRARY_TYPE]

            if count_mode == COUNT_MODE_TARGETED:
                filter_feature_idx = targeted_feature_indices
            elif count_mode == COUNT_MODE_ALL:
                filter_feature_idx = None

            barcode_summary = pdu.collapse_barcode_counts(
                mc,
                filter_library_idx=gex_library_indices,
                filter_feature_idx=filter_feature_idx,
            )
            barcode_summary.drop([FEATURE_DF_DUP_COL], axis=1, inplace=True)

        barcodes_per_sample[sample_label] = barcode_summary

    # FIXME (maria.alexis): not sure this will work correctly with multiplexing
    merged_barcode_info = barcodes_per_sample[TARGETED_LABEL].merge(
        barcodes_per_sample[CONTROL_LABEL],
        on=pdu.FEATURE_DF_BARCODE_SEQ_COL,
        suffixes=[f"_{TARGETED_LABEL}", f"_{CONTROL_LABEL}"],
        how="outer",
    )
    merged_barcode_info.fillna(0, inplace=True)

    del barcodes_per_sample

    return merged_barcode_info


def join(args, outs, chunk_defs, chunk_outs):

    molecule_info_fns = [args.targeted_molecule_info, args.parent_molecule_info]

    merged_barcode_df = None
    merged_barcode_df_index_cols = [
        pdu.FEATURE_DF_BARCODE_SEQ_COL,
        f"{pdu.MOL_INFO_CELL_COL}_{TARGETED_LABEL}",
        f"{pdu.MOL_INFO_CELL_COL}_{CONTROL_LABEL}",
    ]

    with MoleculeCounter.open(args.targeted_molecule_info, "r") as mc:
        targeted_feature_indices = mc.feature_reference.get_target_feature_indices()

    # merge counts in all genes and targeted-genes
    count_modes = COUNT_MODES_TARGETED_COMPARISON
    for count_mode in count_modes:

        barcode_info_df = make_merged_barcode_df(
            molecule_info_fns, targeted_feature_indices, count_mode
        )

        # rename columns with suffixes for what we're tallying
        barcode_info_df.columns = [
            col if col in merged_barcode_df_index_cols else ".".join([col, count_mode])
            for col in barcode_info_df.columns
        ]

        if merged_barcode_df is None:
            merged_barcode_df = barcode_info_df
        else:
            merged_barcode_df = merged_barcode_df.merge(
                barcode_info_df,
                on=merged_barcode_df_index_cols,
                how="outer",
            )
            merged_barcode_df.fillna(0, inplace=True)
        for col in [
            f"{pdu.MOL_INFO_CELL_COL}_{TARGETED_LABEL}",
            f"{pdu.MOL_INFO_CELL_COL}_{CONTROL_LABEL}",
        ]:
            merged_barcode_df[col] = merged_barcode_df[col].replace({np.nan: False}).astype(bool)

    assert merged_barcode_df is not None
    for col in [FEATURE_DF_COUNT_COL, FEATURE_DF_UMI_COL]:
        for label in [TARGETED_LABEL, CONTROL_LABEL]:
            for count_mode in count_modes:
                merged_barcode_df[f"{col}_{label}.{count_mode}"] = merged_barcode_df[
                    f"{col}_{label}.{count_mode}"
                ].astype(int)

    def _get_cell_group_helper(row):
        targeted_call, parent_call = row
        if targeted_call and parent_call:
            return BOTH
        elif targeted_call:
            return TARGETED_ONLY
        elif parent_call:
            return PARENT_ONLY
        else:
            return NEITHER

    # classify by whether it's a cell in one, both, or neither experiment
    merged_barcode_df[TARGETED_COMPARE_CELL_CALL_GRP_COL] = merged_barcode_df[
        [
            f"{pdu.MOL_INFO_CELL_COL}_{TARGETED_LABEL}",
            f"{pdu.MOL_INFO_CELL_COL}_{CONTROL_LABEL}",
        ]
    ].apply(_get_cell_group_helper, axis=1)

    # write to CSV
    outs.barcode_summary_csv = ensure_str(martian.make_path("barcode_summary.csv"))
    pdu.sanitize_dataframe(merged_barcode_df).to_csv(outs.barcode_summary_csv, index=False)
