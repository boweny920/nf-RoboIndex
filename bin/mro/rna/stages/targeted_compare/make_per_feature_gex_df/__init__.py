#!/usr/bin/env python
#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#

"""Make per-gene GEX table for targeted and parent experiments."""


import martian
import numpy as np
import pandas as pd
from six import ensure_str

import cellranger.pandas_utils as pdu
import cellranger.targeted.utils as cr_tgt_utils
from cellranger.molecule_counter import MoleculeCounter
from cellranger.pandas_utils import FEATURE_DF_COUNT_COL, FEATURE_DF_DUP_COL, FEATURE_DF_UMI_COL
from cellranger.rna.library import GENE_EXPRESSION_LIBRARY_TYPE
from cellranger.targeted.targeted_compare_constants import (
    CONTROL_LABEL,
    NORM_MODE_NONE,
    NORM_MODE_RAW,
    NORM_MODE_TARGETED,
    NORM_MODES_TARGETED_COMPARISON,
    TARGETED_COMPARE_ENRICHMENT_COL,
    TARGETED_COMPARE_IS_ENRICHED_COL,
    TARGETED_COMPARE_IS_TGT_GENE_COL,
    TARGETED_COMPARE_LABELS,
    TARGETED_COMPARE_RECOVERY_COL,
    TARGETED_LABEL,
)
from tenkit.stats import robust_divide

__MRO__ = """
stage MAKE_PER_FEATURE_GEX_DF(
    in  h5  targeted_molecule_info,
    in  h5  parent_molecule_info,
    out csv feature_summary_csv,
    src py  "stages/targeted_compare/make_per_feature_gex_df",
) split (
) using (
    volatile = strict,
)
"""

FEATURE_ID_COL = "feature_id"
FEATURE_NAME_COL = "feature_name"
FEATURE_DF_COUNT_CELLS_COL = FEATURE_DF_COUNT_COL + "_cells"
FEATURE_DF_UMI_CELLS_COL = FEATURE_DF_UMI_COL + "_cells"
FEATURE_DF_COLUMNS = [
    FEATURE_ID_COL,
    FEATURE_NAME_COL,
    FEATURE_DF_COUNT_COL,
    FEATURE_DF_UMI_COL,
    FEATURE_DF_DUP_COL,
    FEATURE_DF_COUNT_CELLS_COL,
    FEATURE_DF_UMI_CELLS_COL,
]

# pylint: disable=singleton-comparison


def split(args):
    return {"chunks": [], "join": {"__mem_gb": 4}}


def make_merged_feature_df(molecule_info_fns, targeted_feature_ids, norm_mode):
    """Merge targeted and parent gene expression counts per targeted gene and rescale counts if need be,."""

    features_per_sample = {}
    for mol_info_idx, mol_info_fn in enumerate(molecule_info_fns):

        sample_label = TARGETED_COMPARE_LABELS[mol_info_idx]

        with MoleculeCounter.open(mol_info_fn, "r") as mc:
            gex_library_indices = mc.get_library_indices_by_type()[GENE_EXPRESSION_LIBRARY_TYPE]

            feature_summary = pdu.collapse_feature_counts(
                mc, filter_library_idx=gex_library_indices
            )
            feature_summary = feature_summary[
                feature_summary["feature_type"] == GENE_EXPRESSION_LIBRARY_TYPE
            ]
        feature_summary = feature_summary[FEATURE_DF_COLUMNS]
        features_per_sample[sample_label] = feature_summary

    merged_feature_info = features_per_sample[TARGETED_LABEL].merge(
        features_per_sample[CONTROL_LABEL],
        on=[FEATURE_ID_COL, FEATURE_NAME_COL],
        suffixes=[f"_{TARGETED_LABEL}", f"_{CONTROL_LABEL}"],
        how="outer",
    )

    merged_feature_info[TARGETED_COMPARE_IS_TGT_GENE_COL] = merged_feature_info[
        FEATURE_ID_COL
    ].isin(targeted_feature_ids)

    # rescale reads and UMIs by total count
    cols_to_renorm = [
        FEATURE_DF_UMI_COL,
        FEATURE_DF_COUNT_COL,
        FEATURE_DF_UMI_COL + "_cells",
        FEATURE_DF_COUNT_COL + "_cells",
    ]

    if norm_mode in [NORM_MODE_RAW, NORM_MODE_TARGETED]:

        if norm_mode == NORM_MODE_TARGETED:
            # we don't care about off-target genes if just rescaling to targeted genes, set to NA to simplify rescaling
            for col in cols_to_renorm:
                for sample_label in [TARGETED_LABEL, CONTROL_LABEL]:
                    merged_feature_info.loc[
                        merged_feature_info[TARGETED_COMPARE_IS_TGT_GENE_COL] == False,
                        f"{col}_{sample_label}",
                    ] = np.nan

        for col in cols_to_renorm:
            min_count = min(
                merged_feature_info[f"{col}_{TARGETED_LABEL}"].sum(),
                merged_feature_info[f"{col}_{CONTROL_LABEL}"].sum(),
            )
            for label in [TARGETED_LABEL, CONTROL_LABEL]:
                total_count = merged_feature_info[f"{col}_{label}"].sum()
                merged_feature_info[f"{col}_{label}"] = merged_feature_info[
                    f"{col}_{label}"
                ].multiply(robust_divide(min_count, total_count))

    # drop these columns as they don't make sense post-rescaling
    merged_feature_info.drop(
        [
            f"{FEATURE_DF_DUP_COL}_{TARGETED_LABEL}",
            f"{FEATURE_DF_DUP_COL}_{CONTROL_LABEL}",
        ],
        axis=1,
        inplace=True,
    )

    return merged_feature_info


def join(args, outs, chunk_defs, chunk_outs):

    with MoleculeCounter.open(args.targeted_molecule_info, "r") as mc:
        targeted_feature_ids = mc.feature_reference.get_target_feature_ids()

    molecule_info_fns = [args.targeted_molecule_info, args.parent_molecule_info]

    merged_per_feature_df = None

    # merge gene counts at different rescalings
    for norm_mode in NORM_MODES_TARGETED_COMPARISON:

        feature_info_df = make_merged_feature_df(molecule_info_fns, targeted_feature_ids, norm_mode)
        # only add enrichments once, when normalized to raw reads
        if norm_mode == NORM_MODE_RAW:
            feature_info_df[TARGETED_COMPARE_ENRICHMENT_COL] = np.log2(
                (feature_info_df[f"{FEATURE_DF_COUNT_COL}_cells_{TARGETED_LABEL}"]).divide(
                    feature_info_df[f"{FEATURE_DF_COUNT_COL}_cells_{CONTROL_LABEL}"]
                )
            )
            feature_info_df_subset = feature_info_df[
                np.isfinite(feature_info_df[TARGETED_COMPARE_ENRICHMENT_COL])
            ]
            enrichment_params, _ = cr_tgt_utils.fit_enrichments(
                feature_info_df_subset[TARGETED_COMPARE_IS_TGT_GENE_COL].to_numpy(copy=True),
                feature_info_df_subset[TARGETED_COMPARE_ENRICHMENT_COL].to_numpy(copy=True),
                method=cr_tgt_utils.BOTH_TIED,
            )
            read_enrichment_threshold = enrichment_params.log_rpu_threshold
            print("Read enrichment threshold: %s" % read_enrichment_threshold)
            feature_info_df[TARGETED_COMPARE_IS_ENRICHED_COL] = pd.notna(
                feature_info_df[TARGETED_COMPARE_ENRICHMENT_COL]
            ) & (feature_info_df[TARGETED_COMPARE_ENRICHMENT_COL] > read_enrichment_threshold)
        # only add UMI ratios once, at full depth
        elif norm_mode == NORM_MODE_NONE:
            feature_info_df[TARGETED_COMPARE_RECOVERY_COL] = (
                feature_info_df[f"{FEATURE_DF_UMI_COL}_cells_{TARGETED_LABEL}"]
            ).divide(feature_info_df[f"{FEATURE_DF_UMI_COL}_cells_{CONTROL_LABEL}"])

        # rename columns with suffixes for that norm mode, exclude constant index columns
        feature_info_df.columns = [
            col
            if col
            in [
                FEATURE_ID_COL,
                FEATURE_NAME_COL,
                TARGETED_COMPARE_ENRICHMENT_COL,
                TARGETED_COMPARE_IS_ENRICHED_COL,
                TARGETED_COMPARE_RECOVERY_COL,
                TARGETED_COMPARE_IS_TGT_GENE_COL,
            ]
            else ".".join([col, norm_mode])
            for col in feature_info_df.columns
        ]

        if merged_per_feature_df is None:
            merged_per_feature_df = feature_info_df
        else:
            merged_per_feature_df = merged_per_feature_df.merge(
                feature_info_df,
                on=[FEATURE_ID_COL, FEATURE_NAME_COL, TARGETED_COMPARE_IS_TGT_GENE_COL],
                how="outer",
            )

    # write to CSV
    outs.feature_summary_csv = ensure_str(martian.make_path("feature_summary.csv"))
    pdu.sanitize_dataframe(merged_per_feature_df).to_csv(outs.feature_summary_csv, index=False)
