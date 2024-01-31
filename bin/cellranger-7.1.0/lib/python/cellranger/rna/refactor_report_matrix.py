#!/usr/bin/env python
#
# Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
#
"""Functions for refatoring report_matrix.py.

right now, computes only median confidently mapped reads per singlet per genome
"""
from __future__ import annotations

# pylint: disable=too-many-format-args
import numpy as np
import pandas as pd

import cellranger.feature.utils as feature_utils
import cellranger.utils as cr_utils


def set_default_no_cells(genomes, antibody_present, crispr_present, custom_present):
    """Set default values for metrics."""
    metrics = {}
    metrics["median_total_reads_per_singlet"] = 0

    # rna/report_matrix is a convoluted mess, so handling edge case here
    # for some metrics computed there

    for genome in genomes:
        metrics[f"{genome}_filtered_bcs_median_unique_genes_detected"] = 0
        metrics[f"{genome}_filtered_bcs_mean_unique_genes_detected"] = 0
        metrics[f"{genome}_filtered_bcs_total_unique_genes_detected"] = 0
        metrics[f"{genome}_filtered_bcs_median_counts"] = 0
        metrics[f"{genome}_filtered_bcs_mean_counts"] = 0

    if antibody_present:
        metrics["ANTIBODY_multi_filtered_bcs"] = 0
        metrics["ANTIBODY_filtered_bcs_transcriptome_union"] = 0
        metrics["ANTIBODY_multi_filtered_bcs_median_counts"] = 0
        metrics["ANTIBODY_multi_filtered_bcs_mean_counts"] = 0
        metrics["ANTIBODY_multi_usable_reads_per_filtered_bc"] = 0

    if crispr_present:
        metrics["CRISPR_multi_filtered_bcs"] = 0
        metrics["CRISPR_multi_filtered_bcs_median_counts"] = 0
        metrics["CRISPR_multi_filtered_bcs_mean_counts"] = 0
        metrics["CRISPR_multi_usable_reads_per_filtered_bc"] = 0

    if custom_present:
        metrics["Custom_multi_filtered_bcs"] = 0
        metrics["Custom_multi_filtered_bcs_median_counts"] = 0
        metrics["Custom_multi_filtered_bcs_mean_counts"] = 0
        metrics["Custom_multi_usable_reads_per_filtered_bc"] = 0
    return metrics


def set_default_no_barcode_metrics():
    """Set default values for metrics."""
    metrics = {}
    # metrics["median_total_reads_per_singlet"] = 0
    return metrics


def compute_per_cell_metrics(
    filtered_barcodes_path,
    per_barcode_metrics_path,
    genomes,
    antibody_present=False,
    crispr_present=False,
    custom_present=False,
    is_targeted=False,
):
    """Right now, computes only median confidently mapped reads per singlet per genome."""
    # input validation
    input_files = [per_barcode_metrics_path, filtered_barcodes_path]
    input_files_present = feature_utils.all_files_present(input_files)
    if not input_files_present:
        raise ValueError("Per barcode metrics or filtered_barcodes CSV is not present")

    filtered_barcodes = [
        x.decode("utf8") for x in feature_utils.get_gex_cell_list(filtered_barcodes_path)
    ]
    num_cells = len(filtered_barcodes)

    try:
        per_barcode_metrics = pd.read_csv(per_barcode_metrics_path)
    except pd.errors.EmptyDataError:
        return set_default_no_barcode_metrics()

    if num_cells == 0:
        return set_default_no_cells(
            genomes=genomes,
            antibody_present=antibody_present,
            crispr_present=crispr_present,
            custom_present=custom_present,
        )

    per_filtered_bc_metrics = per_barcode_metrics[
        per_barcode_metrics["barcode"].isin(filtered_barcodes)
    ]

    metrics = {}

    metrics["median_total_reads_per_singlet"] = np.median(
        per_filtered_bc_metrics["raw_reads"].values
    )

    genome_to_barcodes = cr_utils.load_barcode_csv(filtered_barcodes_path)

    for genome in genomes:

        # must get the per-barcode metrics for only barcodes from this genome
        genome_barcodes = [x.decode("utf8") for x in genome_to_barcodes[genome.encode()]]
        per_filtered_bc_metrics_genome = per_barcode_metrics[
            per_barcode_metrics["barcode"].isin(genome_barcodes)
        ]

        # pylint: disable=cell-var-from-loop
        def add_median_metric(metric_name, barcode_metric_key):
            """Calculates the median over genome-assigned-barcodes of a given per-barcode-metric and adds it under metric_name (adds genome prefix)."""
            key = f"{genome}_{metric_name}"
            value = np.median(
                per_filtered_bc_metrics_genome[f"{barcode_metric_key}_{genome}"].values
            )
            metrics[key] = value

        add_median_metric("median_reads_per_singlet", "mapped_reads")
        add_median_metric("median_conf_reads_per_singlet", "conf_reads")

        if is_targeted:
            add_median_metric("median_reads_per_singlet_ontarget", "ontarget_reads")
            add_median_metric("median_reads_per_singlet_offtarget", "offtarget_reads")
            add_median_metric("median_conf_reads_per_singlet_ontarget", "conf_ontarget_reads")
            add_median_metric("median_conf_reads_per_singlet_offtarget", "conf_offtarget_reads")

    return metrics
