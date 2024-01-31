#!/usr/bin/env python3
#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#

"""Takes the output from the targeted comparison tool and produced a.

websummary of the results and a json of relevant targeting metrics.
"""

from __future__ import annotations

import json
import math
import os
from collections import Counter, namedtuple

import h5py
import martian
import numpy as np
import pandas as pd
from six import ensure_binary, ensure_str
from sklearn.metrics import adjusted_rand_score

import cellranger.matrix as cr_matrix
import cellranger.targeted.targeted_compare_utils as cr_tgt_cmp_utils
import tenkit.safe_json as tk_safe_json
from cellranger.molecule_counter import MoleculeCounter
from cellranger.pandas_utils import (
    FEATURE_DF_BARCODE_SEQ_COL,
    FEATURE_DF_COUNT_COL,
    FEATURE_DF_UMI_COL,
    MOL_INFO_CELL_COL,
    sanitize_dataframe,
)
from cellranger.rna.library import GENE_EXPRESSION_LIBRARY_TYPE
from cellranger.targeted import simple_utils
from cellranger.targeted.targeted_compare_constants import (
    BOTH,
    CONTROL_LABEL,
    CONTROL_WTA_LABEL,
    COUNT_MODE_TARGETED,
    NORM_MODE_NONE,
    NORM_MODE_RAW,
    NORM_MODE_TARGETED,
    PARENT_ONLY,
    TARGETED_COMPARE_CELL_CALL_GRP_COL,
    TARGETED_COMPARE_ENRICHMENT_COL,
    TARGETED_COMPARE_IS_ENRICHED_COL,
    TARGETED_COMPARE_IS_TGT_GENE_COL,
    TARGETED_COMPARE_LABELS,
    TARGETED_COMPARE_RECOVERY_COL,
    TARGETED_LABEL,
    TARGETED_ONLY,
    TargetCompareConstants,
)
from cellranger.targeted.targeted_compare_utils import TARGETED_COMPARE_WS, TargetedControlCompareSC
from cellranger.version import get_version
from cellranger.websummary.metrics import (
    SpatialTargetedCompareMetricAnnotations,
    TargetedCompareMetricAnnotations,
)
from tenkit.stats import robust_divide
from websummary import summarize

__MRO__ = """
stage SUMMARIZE_COMPARISON(
    in  string sample_id,
    in  string sample_desc,
    in  h5     targeted_molecule_info,
    in  h5     parent_molecule_info,
    in  csv    target_set,
    in  csv    feature_summary_csv,
    in  csv    barcode_summary_csv,
    in  csv    sc_rna_corr_df_by_barcode,
    in  csv    sc_rna_corr_df_by_gene,
    in  h5     combined_matrix_h5,
    in  h5     sc_rna_tsne_h5,
    out csv    metrics_summary_csv,
    out json   summary_metrics_json,
    out html   web_summary,
    out json   alarms_json,
    src py     "stages/targeted_compare/summarize_comparison",
) using (
    mem_gb = 4,
    volatile = strict,
)
"""


ClusteringConfig = namedtuple(
    "ClusteringConfig", ["matrix_h5", "tsne_h5", "clusters_h5", "skip_clustering"]
)


def split(args):
    with MoleculeCounter.open(args.targeted_molecule_info, "r") as in_mc:
        mol_info_mem_gb_tgt = int(math.ceil(in_mc.estimate_mem_gb(in_mc.nrows(), scale=4.0)))
    with MoleculeCounter.open(args.parent_molecule_info, "r") as in_mc:
        mol_info_mem_gb_ctrl = int(math.ceil(in_mc.estimate_mem_gb(in_mc.nrows(), scale=4.0)))
    return {"chunks": [], "join": {"__mem_gb": max(mol_info_mem_gb_tgt, mol_info_mem_gb_ctrl)}}


def _basic_metrics(mol_info_fns, feature_summary_df, barcode_summary_df):
    """Fetches and/or calculates basic metrics at full depth on both samples.

    Metrics computed here
    are:
        Estimated number of cells, Total Read Pairs, Fraction of Reads Confidently Mapped to the
        Targeted Transcriptome, Mean Reads per Cell, Median UMIs per Cell, Mean Targeted Reads per Cell,
        Median Targeted UMIs per Cell, Median Genes per Cell, Median Targeted Genes per Cell, Number of
        Targeted Genes Detected Exclusively in either sample.
    """
    metrics = {}

    for (mol_info_fn, label) in zip(mol_info_fns, TARGETED_COMPARE_LABELS):
        # number of reads, number of cells
        with MoleculeCounter.open(mol_info_fn, "r") as mc:
            gex_library_indices = mc.get_library_indices_by_type()[GENE_EXPRESSION_LIBRARY_TYPE]
            num_cells = len(
                set(
                    mc.get_filtered_barcodes(
                        mc.get_barcode_info(),
                        mc.get_library_info(),
                        mc.get_barcodes(),
                        library_type=GENE_EXPRESSION_LIBRARY_TYPE,
                    )
                )
            )
            raw_reads_per_lib = mc.get_raw_read_pairs_per_library()
            total_reads = np.sum([raw_reads_per_lib[idx] for idx in gex_library_indices])
        metrics[f"num_cells_{label}"] = num_cells
        metrics[f"total_read_pairs_{label}"] = total_reads

        total_targeted_reads = feature_summary_df[
            feature_summary_df[TARGETED_COMPARE_IS_TGT_GENE_COL]
        ][f"{FEATURE_DF_COUNT_COL}_{label}.{NORM_MODE_NONE}"].sum()
        metrics[f"total_targeted_read_pairs_{label}"] = total_targeted_reads

        # frac reads on target
        frac_reads_on_target = robust_divide(
            metrics[f"total_targeted_read_pairs_{label}"],
            metrics[f"total_read_pairs_{label}"],
        )
        metrics[f"fraction_reads_on_target_{label}"] = frac_reads_on_target

        # reads and UMIs per cell, on target reads and UMIs per cell
        if label == TARGETED_LABEL:
            groups_to_include = [BOTH, TARGETED_ONLY]
        elif label == CONTROL_LABEL:
            groups_to_include = [BOTH, PARENT_ONLY]
        mean_reads_per_cell = np.divide(total_reads, num_cells)
        metrics[f"mean_reads_per_cell_{label}"] = mean_reads_per_cell
        mean_reads_per_cell = np.divide(total_targeted_reads, num_cells)
        metrics[f"mean_targeted_reads_per_cell_{label}"] = mean_reads_per_cell
        median_umis_per_cell = np.nanmedian(
            barcode_summary_df[
                barcode_summary_df[TARGETED_COMPARE_CELL_CALL_GRP_COL].isin(groups_to_include)
            ][f"{FEATURE_DF_UMI_COL}_{label}.{COUNT_MODE_TARGETED}"]
        )
        metrics[f"median_targeted_umis_per_cell_{label}"] = median_umis_per_cell

    subset_df = barcode_summary_df[barcode_summary_df[TARGETED_COMPARE_CELL_CALL_GRP_COL] == BOTH]
    cell_targeted_depth_factor = (
        (subset_df[f"{FEATURE_DF_COUNT_COL}_{TARGETED_LABEL}.{COUNT_MODE_TARGETED}"] + 1).divide(
            subset_df[f"{FEATURE_DF_COUNT_COL}_{CONTROL_LABEL}.{COUNT_MODE_TARGETED}"] + 1
        )
    ).mean()
    metrics["cell_targeted_depth_factor"] = cell_targeted_depth_factor
    metrics["num_panel_genes"] = feature_summary_df[TARGETED_COMPARE_IS_TGT_GENE_COL].sum()

    # cell-call groups
    metrics[f"num_cell_calls_{BOTH}"] = Counter(
        barcode_summary_df[TARGETED_COMPARE_CELL_CALL_GRP_COL]
    )[BOTH]
    metrics[f"num_cell_calls_{TARGETED_ONLY}"] = Counter(
        barcode_summary_df[TARGETED_COMPARE_CELL_CALL_GRP_COL]
    )["targeted-only"]
    metrics[f"num_cell_calls_{PARENT_ONLY}"] = Counter(
        barcode_summary_df[TARGETED_COMPARE_CELL_CALL_GRP_COL]
    )["parent-only"]
    metrics["frac_cell_barcode_overlap"] = robust_divide(
        metrics["num_cell_calls_both"],
        min(
            metrics[f"num_cells_{TARGETED_LABEL}"],
            metrics[f"num_cells_{CONTROL_LABEL}"],
        ),
    )

    return metrics


# pylint: disable=too-many-locals, invalid-name
def _calculate_metrics(
    mol_info_fns,
    feature_summary_df,
    barcode_summary_df,
):
    """Calculate a few summary metrics useful to display in WS table."""

    metrics = _basic_metrics(mol_info_fns, feature_summary_df, barcode_summary_df)

    # just subset to targeted genes
    tgt_summary_df = feature_summary_df[feature_summary_df[TARGETED_COMPARE_IS_TGT_GENE_COL]]

    # prpc, uses raw
    reads_in_cells_targeted = "{}_cells_{}.{}".format(
        FEATURE_DF_COUNT_COL, TARGETED_LABEL, NORM_MODE_RAW
    )
    reads_in_cells_parent = "{}_cells_{}.{}".format(
        FEATURE_DF_COUNT_COL, CONTROL_LABEL, NORM_MODE_RAW
    )
    # num targeted genes detected, uses raw
    num_panel_genes_detected = tgt_summary_df[
        (tgt_summary_df[reads_in_cells_targeted] > 0) | (tgt_summary_df[reads_in_cells_parent] > 0)
    ].shape[0]
    metrics["num_targeted_genes_detected"] = num_panel_genes_detected
    num_panel_genes_detected_parent = tgt_summary_df[
        (tgt_summary_df[reads_in_cells_parent] > 0)
    ].shape[0]
    num_panel_genes_detected_targeted = tgt_summary_df[
        (tgt_summary_df[reads_in_cells_targeted] > 0)
    ].shape[0]
    metrics[f"num_targeted_genes_detected_{CONTROL_LABEL}"] = num_panel_genes_detected_parent
    metrics[f"total_targeted_genes_detected_{TARGETED_LABEL}"] = num_panel_genes_detected_targeted
    metrics[f"total_targeted_genes_detected_{CONTROL_LABEL}"] = num_panel_genes_detected_parent

    num_panel_genes_exclusive_targeted = tgt_summary_df[
        (tgt_summary_df[reads_in_cells_targeted] > 0) & (tgt_summary_df[reads_in_cells_parent] == 0)
    ].shape[0]
    num_panel_genes_exclusive_parent = tgt_summary_df[
        (tgt_summary_df[reads_in_cells_parent] > 0) & (tgt_summary_df[reads_in_cells_targeted] == 0)
    ].shape[0]
    metrics[
        f"num_targeted_genes_detected_exclusive_{TARGETED_LABEL}"
    ] = num_panel_genes_exclusive_targeted
    metrics[
        f"num_targeted_genes_detected_exclusive_{CONTROL_LABEL}"
    ] = num_panel_genes_exclusive_parent

    # r-squared and mean ratio of relative read and UMI counts.
    # For reads, uses NORM_MORE_RAW; for UMIs, uses NORM_MODE_NONE
    metrics["targeted_gene_read_rsquared"] = (
        np.log10(tgt_summary_df[reads_in_cells_targeted] + 1).corr(
            np.log10(tgt_summary_df[reads_in_cells_parent] + 1)
        )
        ** 2
    )
    umis_in_cells_targeted = "{}_cells_{}.{}".format(
        FEATURE_DF_UMI_COL, TARGETED_LABEL, NORM_MODE_NONE
    )
    umis_in_cells_parent = "{}_cells_{}.{}".format(
        FEATURE_DF_UMI_COL, CONTROL_LABEL, NORM_MODE_NONE
    )
    metrics["targeted_gene_umi_rsquared"] = (
        np.log10(tgt_summary_df[umis_in_cells_targeted] + 1).corr(
            np.log10(tgt_summary_df[umis_in_cells_parent] + 1)
        )
        ** 2
    )

    targeted_enrichments = tgt_summary_df[tgt_summary_df[reads_in_cells_parent] > 0][
        [
            TARGETED_COMPARE_ENRICHMENT_COL,
            TARGETED_COMPARE_IS_ENRICHED_COL,
            TARGETED_COMPARE_RECOVERY_COL,
        ]
    ]

    is_enriched = targeted_enrichments[TARGETED_COMPARE_IS_ENRICHED_COL]
    metrics["num_targeted_genes_enriched"] = is_enriched.sum()
    metrics["frac_targeted_genes_enriched"] = robust_divide(
        metrics["num_targeted_genes_enriched"],
        metrics[f"num_targeted_genes_detected_{CONTROL_LABEL}"],
    )
    metrics["frac_non_targeted_genes_enriched"] = robust_divide(
        feature_summary_df[
            (~feature_summary_df[TARGETED_COMPARE_IS_TGT_GENE_COL])
            & (feature_summary_df[reads_in_cells_parent] > 0)
        ][TARGETED_COMPARE_IS_ENRICHED_COL].sum(),
        feature_summary_df[
            (~feature_summary_df[TARGETED_COMPARE_IS_TGT_GENE_COL])
            & (feature_summary_df[reads_in_cells_parent] > 0)
        ].shape[0],
    )

    metrics["mean_frac_UMIs_recovered_per_gene"] = np.mean(
        targeted_enrichments[TARGETED_COMPARE_RECOVERY_COL]
    )

    finite_enrichments = targeted_enrichments[
        np.isfinite(targeted_enrichments[TARGETED_COMPARE_ENRICHMENT_COL])
    ][TARGETED_COMPARE_ENRICHMENT_COL]
    metrics["mean_read_enrichment"] = np.power(2, np.mean(finite_enrichments))

    return metrics


# pylint: disable=too-many-arguments
def compare_cell_clustering(
    targeted_clustering_config,
    parent_subset_clustering_config,
    parent_wta_clustering_config,
    barcode_summary_df,
    plot_same_cluster_k=False,
):
    """Creates a dataframe of all barcodes in all 3 samples.

    Columns contain the cluster
    each cell belongs to as well as its coordinates for t-SNE plotting in 2D.
    """
    if (
        targeted_clustering_config.skip_clustering
        or parent_subset_clustering_config.skip_clustering
        or parent_wta_clustering_config.skip_clustering
    ):
        return {}, None

    clustering_stats = {}
    use_clustering_key = None

    parent_wta_df = _load_tsne_projection_coords(
        parent_wta_clustering_config.matrix_h5, parent_wta_clustering_config.tsne_h5
    )
    parent_wta_df["clusters"] = _load_clustering(
        parent_wta_clustering_config.clusters_h5, n_clusters=use_clustering_key
    )
    clustering_stats[f"{CONTROL_WTA_LABEL}_clusters_best"] = parent_wta_df["clusters"].max()
    parent_wta_df.columns = parent_wta_df.columns.map(
        lambda x: x if x == FEATURE_DF_BARCODE_SEQ_COL else x + "." + CONTROL_WTA_LABEL
    )

    # if plot same K for all clusters, use the parent WTA clustering as a reference for the other two
    if plot_same_cluster_k:
        use_clustering_key = "_{}".format(clustering_stats[f"{CONTROL_WTA_LABEL}_clusters_best"])
    else:
        use_clustering_key = None

    targeted_df = _load_tsne_projection_coords(
        targeted_clustering_config.matrix_h5, targeted_clustering_config.tsne_h5
    )
    targeted_df["clusters"] = _load_clustering(
        targeted_clustering_config.clusters_h5, n_clusters=use_clustering_key
    )
    clustering_stats[f"{TARGETED_LABEL}_clusters_best"] = targeted_df["clusters"].max()
    targeted_df.columns = targeted_df.columns.map(
        lambda x: x if x == FEATURE_DF_BARCODE_SEQ_COL else x + "." + TARGETED_LABEL
    )

    parent_subset_df = _load_tsne_projection_coords(
        parent_subset_clustering_config.matrix_h5, parent_subset_clustering_config.tsne_h5
    )
    parent_subset_df["clusters"] = _load_clustering(
        parent_subset_clustering_config.clusters_h5, n_clusters=use_clustering_key
    )
    clustering_stats[f"{CONTROL_LABEL}_clusters_best"] = parent_subset_df["clusters"].max()
    parent_subset_df.columns = parent_subset_df.columns.map(
        lambda x: x if x == FEATURE_DF_BARCODE_SEQ_COL else x + "." + CONTROL_LABEL
    )

    # merge the 3 dataframes on barcode
    merged_df = targeted_df.merge(parent_wta_df, on=FEATURE_DF_BARCODE_SEQ_COL, how="outer")
    merged_df = merged_df.merge(parent_subset_df, on=FEATURE_DF_BARCODE_SEQ_COL, how="left")

    # hacky but NaN -> -1 for cell barcodes only in sample, and make these columns ints again
    for col in [
        f"clusters.{TARGETED_LABEL}",
        f"clusters.{CONTROL_LABEL}",
        f"clusters.{CONTROL_WTA_LABEL}",
    ]:
        merged_df[col] = (merged_df[col].replace(np.nan, -1)).astype(int)

    # reassign cluster indices to maximize agreement in per-barcode assignment between numbered clusters
    # across experiments. e.g. if clusters 1 and 2 are swapped in the targeted experiment, swap the
    # cluster numbers before plotting tSNEs so their colors are consitent.
    merged_df = relabel_clusters(merged_df, CONTROL_LABEL, TARGETED_LABEL)
    merged_df = relabel_clusters(merged_df, CONTROL_LABEL, CONTROL_WTA_LABEL)

    # get clustering similarity stats just for bookkeeping
    subset_df = merged_df[
        (merged_df[f"clusters.{TARGETED_LABEL}"] != -1)
        & (merged_df[f"clusters.{CONTROL_LABEL}"] != -1)
    ]
    clustering_stats["adj_rand_score_ctrl_tgt"] = adjusted_rand_score(
        subset_df[f"clusters.{TARGETED_LABEL}"],
        subset_df[f"clusters.{CONTROL_LABEL}"],
    )

    subset_df = merged_df[
        (merged_df[f"clusters.{CONTROL_LABEL}"] != -1)
        & (merged_df[f"clusters.{CONTROL_WTA_LABEL}"] != -1)
    ]
    clustering_stats["adj_rand_score_ctrl_ctrl_wta"] = adjusted_rand_score(
        subset_df[f"clusters.{CONTROL_LABEL}"],
        subset_df[f"clusters.{CONTROL_WTA_LABEL}"],
    )

    # set non-cells to -1
    for label in [TARGETED_LABEL, CONTROL_LABEL, CONTROL_WTA_LABEL]:
        if label == TARGETED_LABEL:
            cells = barcode_summary_df[barcode_summary_df[f"{MOL_INFO_CELL_COL}_{TARGETED_LABEL}"]][
                FEATURE_DF_BARCODE_SEQ_COL
            ].tolist()
        else:
            cells = barcode_summary_df[barcode_summary_df[f"{MOL_INFO_CELL_COL}_{CONTROL_LABEL}"]][
                FEATURE_DF_BARCODE_SEQ_COL
            ].tolist()
        merged_df.loc[~merged_df[FEATURE_DF_BARCODE_SEQ_COL].isin(cells), f"clusters.{label}"] = -1

    return clustering_stats, merged_df


def relabel_clusters(merged_df, label1, label2):
    """Greedily reassigns/renumber clusters according to maximum fraction of overlapping.

    cell barcodes.
    """
    max_clusters_1 = merged_df[f"clusters.{label1}"].max()
    max_clusters_2 = merged_df[f"clusters.{label2}"].max()

    def _get_union(row, subset_df):
        cluster1, cluster2 = row[:2]
        return subset_df[
            (subset_df.iloc[:, 0] == cluster1) | (subset_df.iloc[:, 1] == cluster2)
        ].shape[0]

    # get intersection and union of cell barcodes assigned to all pairs of clusters
    subset_df = merged_df[
        (merged_df[f"clusters.{label1}"] != -1) & (merged_df[f"clusters.{label2}"] != -1)
    ]

    df = subset_df.groupby([f"clusters.{label1}", f"clusters.{label2}"]).size().reset_index()
    df.rename(columns={0: "intersection"}, inplace=True)

    n_in_clusters_1 = subset_df.groupby(f"clusters.{label1}").size().to_dict()
    n_in_clusters_2 = subset_df.groupby(f"clusters.{label2}").size().to_dict()
    df["n_in_cluster_1"] = df[f"clusters.{label1}"].map(n_in_clusters_1)
    df["n_in_cluster_2"] = df[f"clusters.{label2}"].map(n_in_clusters_2)
    df["min_n_in_cluster"] = df[[f"clusters.{label1}", f"clusters.{label2}"]].min(axis=1)

    df["frac_overlap"] = df["intersection"].divide(df["min_n_in_cluster"])
    df = df.loc[df["frac_overlap"] >= 0.5]
    df.sort_values("frac_overlap", ascending=False, inplace=True)

    # do a greedy reassignment of pairings of clusters across samples, by assigning clusters with the highest overlap to one another
    cluster_assignments_1 = {}
    cluster_assignments_2 = {}
    matching_clusters = set()
    while len(cluster_assignments_1) < min(max_clusters_1, max_clusters_2) and df.shape[0] > 0:
        cluster1, cluster2 = (
            df[f"clusters.{label1}"].iloc[0],
            df[f"clusters.{label2}"].iloc[0],
        )
        # keep cluster numbering from first group
        cluster_assignments_1[cluster1] = cluster1
        cluster_assignments_2[cluster2] = cluster1
        matching_clusters.add(cluster1)
        # remove both those clusters from the running
        df = df[(df[f"clusters.{label1}"] != cluster1) & (df[f"clusters.{label2}"] != cluster2)]

    # if some clusters were unmatchable, just reassign those to the smallest possible non-colliding value
    if len(cluster_assignments_1) < max_clusters_1:
        cluster_assignments_values = set(cluster_assignments_1.values())
        for cluster1 in range(1, max_clusters_1 + 1):
            if cluster1 not in cluster_assignments_1:
                new_cluster_idx = [
                    i for i in range(1, max_clusters_1 + 1) if i not in cluster_assignments_values
                ][0]
                cluster_assignments_1[cluster1] = new_cluster_idx
                cluster_assignments_values.add(new_cluster_idx)
    if len(cluster_assignments_2) < max_clusters_2:
        cluster_assignments_values = set(cluster_assignments_2.values())
        for cluster2 in range(1, max_clusters_2 + 1):
            if cluster2 not in cluster_assignments_2:
                new_cluster_idx = [
                    i for i in range(1, max_clusters_2 + 1) if i not in cluster_assignments_values
                ][0]
                cluster_assignments_2[cluster2] = new_cluster_idx
                cluster_assignments_values.add(new_cluster_idx)

    cluster_assignments_1[-1] = -1
    cluster_assignments_2[-1] = -1
    merged_df[f"clusters.{label1}"] = np.vectorize(cluster_assignments_1.get)(
        merged_df[f"clusters.{label1}"]
    )
    merged_df[f"clusters.{label2}"] = np.vectorize(cluster_assignments_2.get)(
        merged_df[f"clusters.{label2}"]
    )

    # keep track of sample-specific clusters and shared clusters for coloring purposes
    # label1 is considered the reference here and will get default colors so no clusters considered sample-specific
    merged_df[f"specific_cluster.{label1}"] = False
    merged_df[f"specific_cluster.{label2}"] = ~(
        merged_df[f"clusters.{label2}"].isin(matching_clusters)
    )

    return merged_df


def _load_clustering(clustering_h5, n_clusters=None):
    with h5py.File(clustering_h5, "r") as in_f:
        if n_clusters is None:
            n_clusters_keys = list(in_f["kmeans"].keys())
            n_clusters = n_clusters_keys[
                np.argmin([np.squeeze(in_f["kmeans"][k]["cluster_score"]) for k in n_clusters_keys])
            ]
        clusters = np.array(in_f["kmeans"][n_clusters]["clusters"])
    return clusters


def _load_tsne_projection_coords(matrix_h5, tsne_h5):
    bcs = cr_matrix.CountMatrix.load_bcs_from_h5(matrix_h5)
    with h5py.File(tsne_h5, "r") as in_f:
        tsne_df = pd.DataFrame(
            np.array(in_f["tsne"]["_2"]["transformed_tsne_matrix"]), columns=["x", "y"]
        )
    tsne_df["barcode"] = bcs
    return tsne_df


def _get_median_gene_metrics(mol_info_fns, targeted_gene_ids):

    genes_per_cell_metrics = {}

    for (mol_info_fn, label) in zip(mol_info_fns, TARGETED_COMPARE_LABELS):
        mat = TargetedControlCompareSC.mol_info_to_coo(mol_info_fn)
        genes_per_cell = TargetedControlCompareSC.get_genes_per_cell_from_coo_matrix(mat)
        del mat
        genes_per_cell_metrics[f"median_genes_per_cell_{label}"] = np.nanmedian(genes_per_cell)
        del genes_per_cell
        mat = TargetedControlCompareSC.mol_info_to_coo(mol_info_fn, targeted_gene_ids)
        genes_per_cell = TargetedControlCompareSC.get_genes_per_cell_from_coo_matrix(mat)
        del mat
        genes_per_cell_metrics[f"median_targeted_genes_per_cell_{label}"] = np.nanmedian(
            genes_per_cell
        )
        del genes_per_cell

    return genes_per_cell_metrics


def prune_output_dfs(barcode_summary_df, feature_summary_df, is_spatial):
    """Prunes the dataframes to just the raw depth counts which will be shown to the customer."""

    def _reorder_columns(df, cols_to_move):
        # Move barcode/feature id columns to front
        assert all(col in df.columns for col in cols_to_move)
        reordered_cols = cols_to_move + [
            col for col in df.columns.tolist() if col not in cols_to_move
        ]
        df = df.reindex(columns=reordered_cols)
        return df

    def _trans_to_readable_names(colnames):
        new_colnames = []
        for col in colnames:
            new_col = " ".join(
                [word.capitalize() if word not in ["in"] else word for word in col.split("_")]
            )
            new_col = new_col.replace("Num", "Number of")
            new_colnames.append(new_col)
        return new_colnames

    def _trans_to_spatial_names(colnames):
        new_colnames = []
        for col in colnames:
            new_col = " ".join(
                [word.capitalize() if word not in ["in"] else word for word in col.split("_")]
            )
            new_col = new_col.replace("Cell", "Tissue-Covered Spot")
            new_colnames.append(new_col)
        return new_colnames

    barcode_summary_df = _reorder_columns(
        barcode_summary_df, cols_to_move=[FEATURE_DF_BARCODE_SEQ_COL]
    )
    feature_summary_df = _reorder_columns(
        feature_summary_df, cols_to_move=["feature_id", "feature_name"]
    )

    # prune barcode_summary_df -> replace "." separator with "_"
    new_colnames = [col.replace(".", "_") for col in barcode_summary_df.columns]
    barcode_summary_df.columns = _trans_to_readable_names(new_colnames)

    if is_spatial:
        new_colnames_spatial = [
            col.replace("Cell", "Tissue-Covered Spot") for col in barcode_summary_df.columns
        ]
        barcode_summary_df.columns = _trans_to_spatial_names(new_colnames_spatial)

    # prune feature_summary_df -> remove all downsamplings and rescalings, remove ".None" suffix for raw depth
    cols_to_keep = [
        col
        for col in feature_summary_df.columns
        if col.split(".")[-1] not in [NORM_MODE_RAW, NORM_MODE_TARGETED]
    ]
    feature_summary_df = feature_summary_df[cols_to_keep]
    new_colnames = [col.replace("." + NORM_MODE_NONE, "") for col in feature_summary_df.columns]
    feature_summary_df.columns = _trans_to_readable_names(new_colnames)

    if is_spatial:
        new_colnames_spatial = [
            col.replace("Cell", "Tissue-Covered Spot") for col in feature_summary_df.columns
        ]
        feature_summary_df.columns = _trans_to_spatial_names(new_colnames_spatial)

    return barcode_summary_df, feature_summary_df


def join(args, outs, chunk_defs, chunk_outs):
    # Is the sample a spatial one?
    IS_SPATIAL = MoleculeCounter.open(args.targeted_molecule_info, "r").is_spatial_data()

    mol_info_fns = [args.targeted_molecule_info, args.parent_molecule_info]
    with MoleculeCounter.open(args.targeted_molecule_info, "r") as mc:
        species_list = mc.feature_reference.get_genomes()
    _, targeted_gene_ids, _ = simple_utils.parse_target_csv(args.target_set)

    barcode_summary_df = pd.read_csv(
        ensure_str(args.barcode_summary_csv), converters={FEATURE_DF_BARCODE_SEQ_COL: ensure_binary}
    )
    feature_summary_df = pd.read_csv(ensure_str(args.feature_summary_csv))

    # calculate metrics needed for websummary
    metrics = _calculate_metrics(mol_info_fns, feature_summary_df, barcode_summary_df)
    metrics.update(_get_median_gene_metrics(mol_info_fns, targeted_gene_ids))

    targeted_clustering_config = ClusteringConfig(
        args.targeted_matrix_h5,
        args.tsne_targeted_h5,
        args.targeted_clusters_h5,
        args.skip_clustering_targeted,
    )
    parent_subset_clustering_config = ClusteringConfig(
        args.parent_matrix_subset_h5,
        args.tsne_parent_subset_h5,
        args.parent_clusters_subset_h5,
        args.skip_clustering_parent_subset,
    )
    parent_wta_clustering_config = ClusteringConfig(
        args.parent_matrix_wta_h5,
        args.tsne_parent_wta_h5,
        args.parent_clusters_wta_h5,
        args.skip_clustering_parent_wta,
    )
    clustering_metrics, clustering_df = compare_cell_clustering(
        targeted_clustering_config,
        parent_subset_clustering_config,
        parent_wta_clustering_config,
        barcode_summary_df,
        args.plot_same_cluster_k,
    )
    metrics.update(clustering_metrics)

    # Determine if spatial or not
    if IS_SPATIAL:
        target_compare_constants = TargetCompareConstants(
            cell_or_spot="tissue covered spot",
            caps_cell_or_spot="Tissue Covered Spot",
            alt_caps_cell_or_spot="Tissue covered spot",
        )
        metadata = SpatialTargetedCompareMetricAnnotations()
        name = "Space Ranger"
    else:
        target_compare_constants = TargetCompareConstants(
            cell_or_spot="cell", caps_cell_or_spot="Cell", alt_caps_cell_or_spot="Cell"
        )
        metadata = TargetedCompareMetricAnnotations()
        name = "Cell Ranger"
    websummary_data = {
        TARGETED_COMPARE_WS: {
            "sample": {
                "id": args.sample_id,
                "description": args.sample_desc,
                "command": name,
                "subcommand": "Targeted Compare",
                "pipeline_version": get_version(),
            }
        },
        # FIXME: can probably just rename TARGETED_COMPARE_WS to summary
        # for now, hack by duplicating arguments so title-tag shows up correctly
        "summary": {
            "sample": {
                "id": args.sample_id,
                "description": args.sample_desc,
            }
        },
    }

    # add tables
    cr_tgt_cmp_utils.add_tables(
        websummary_data[TARGETED_COMPARE_WS],
        args.target_set,
        mol_info_fns,
        metrics,
        metadata,
        species_list,
    )

    # add plots
    cr_tgt_cmp_utils.add_plots(
        websummary_data[TARGETED_COMPARE_WS],
        metrics,
        barcode_summary_df,
        feature_summary_df,
        clustering_df,
        IS_SPATIAL,
        target_compare_constants,
    )

    # Pull up the correct template information
    template_path = os.path.dirname(__file__)
    template_file = os.path.join(template_path, "targeted-compare.html")
    with open(template_file) as infile:
        template = infile.read()

    # sanitize websummary
    websummary_data = tk_safe_json.json_sanitize(websummary_data)
    with open(outs.web_summary, "w") as outfile:
        summarize.generate_html_summary(websummary_data, template, template_path, outfile)

    assert "alarms" in websummary_data[TARGETED_COMPARE_WS]
    outs.alarms_json = martian.make_path("alarms.json")
    with open(outs.alarms_json, "w") as out_f:
        json.dump(
            websummary_data[TARGETED_COMPARE_WS]["alarms"],
            out_f,
            indent=4,
            sort_keys=True,
        )

    outs.summary_metrics_json = martian.make_path("summary_metrics.json")
    with open(outs.summary_metrics_json, "w") as out_f:
        tk_safe_json.dump_numpy(metrics, out_f, indent=4, sort_keys=True)

    outs.metrics_summary_csv = martian.make_path("metrics_summary.csv")

    cr_tgt_cmp_utils.build_tgtcmp_metrics_csv(
        metrics, outs.metrics_summary_csv, is_spatial=IS_SPATIAL
    )

    barcode_summary_df, feature_summary_df = prune_output_dfs(
        barcode_summary_df, feature_summary_df, is_spatial=IS_SPATIAL
    )
    # merge with clustering df
    if clustering_df is not None:
        barcode_summary_df = barcode_summary_df.merge(
            clustering_df[
                [
                    FEATURE_DF_BARCODE_SEQ_COL,
                    f"clusters.{TARGETED_LABEL}",
                    f"clusters.{CONTROL_LABEL}",
                    f"clusters.{CONTROL_WTA_LABEL}",
                ]
            ],
            left_on=FEATURE_DF_BARCODE_SEQ_COL.capitalize(),
            right_on=FEATURE_DF_BARCODE_SEQ_COL,
            how="left",
        )
        for label in [TARGETED_LABEL, CONTROL_LABEL, CONTROL_WTA_LABEL]:
            barcode_summary_df[f"clusters.{label}"] = (
                barcode_summary_df[f"clusters.{label}"].replace(np.nan, -1).astype(int)
            )
        barcode_summary_df.rename(
            columns={
                f"clusters.{TARGETED_LABEL}": "Cluster in Targeted Sample",
                "clusters.{}".format(
                    CONTROL_LABEL
                ): "Cluster in Parent Sample Subset to Targeted Genes",
                f"clusters.{CONTROL_WTA_LABEL}": "Cluster in Parent Sample",
            },
            inplace=True,
        )
        # pandas is annoying and duplicates this column
        barcode_summary_df.drop(FEATURE_DF_BARCODE_SEQ_COL, inplace=True, axis=1)

    outs.barcode_summary_csv = martian.make_path("barcode_summary.csv")
    outs.feature_summary_csv = martian.make_path("feature_summary.csv")
    sanitize_dataframe(barcode_summary_df).to_csv(ensure_str(outs.barcode_summary_csv), index=False)
    sanitize_dataframe(feature_summary_df).to_csv(ensure_str(outs.feature_summary_csv), index=False)
