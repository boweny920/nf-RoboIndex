#!/usr/bin/env python
#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#

"""Computes per-cell and per-gene GEX correlations for targeted-comparison tool."""


import math

import martian
import numpy as np
import pandas as pd
from six import ensure_binary, ensure_str

import cellranger.targeted.targeted_compare_utils as cr_tgt_cmp_utils
from cellranger.molecule_counter import MoleculeCounter
from cellranger.targeted.targeted_compare_constants import (
    CONTROL_LABEL,
    NORM_MODE_RAW,
    TARGETED_COMPARE_CELL_CALL_GRP_COL,
    TARGETED_LABEL,
)
from cellranger.targeted.targeted_compare_utils import TargetedControlCompareSC

__MRO__ = """
stage COMPUTE_SC_CORRELATIONS(
    in  h5     targeted_molecule_info,
    in  h5     parent_molecule_info,
    in  csv    barcode_summary_csv,
    out csv    sc_rna_corr_df_by_barcode,
    out csv    sc_rna_corr_df_by_gene,
    out h5     targeted_matrix_h5,
    out h5     parent_matrix_h5,
    out h5     combined_matrix_h5,
    src py     "stages/targeted_compare/compute_sc_correlations",
) split (
) using (
    volatile = strict,
)
"""


def split(args):
    with MoleculeCounter.open(args.targeted_molecule_info, "r") as mc:
        # 4.0 empirically determined using memory_profiler in a jupyter notebook
        tgt_mol_info_mem_gb = int(math.ceil(mc.estimate_mem_gb(mc.nrows(), scale=4.0)))
        n_genes = len(mc.feature_reference.get_target_feature_ids())
    with MoleculeCounter.open(args.parent_molecule_info, "r") as mc:
        # 4.0 empirically determined using memory_profiler in a jupyter notebook
        ctrl_mol_info_mem_gb = int(math.ceil(mc.estimate_mem_gb(mc.nrows(), scale=4.0)))

    barcodes = pd.read_csv(
        ensure_str(args.barcode_summary_csv), usecols=[TARGETED_COMPARE_CELL_CALL_GRP_COL]
    )
    n_bc = barcodes[barcodes[TARGETED_COMPARE_CELL_CALL_GRP_COL] != "neither"].shape[0]

    # memory consumed by matrix and some buffer for table and data structures
    csc_mem_gb = int(math.ceil(8 * n_bc * n_genes * 3 / 1e9)) + 1

    return {
        "chunks": [],
        "join": {
            # max of a couple of combinations of objects that are stored simulatenously, plus buffer
            "__mem_gb": max(csc_mem_gb + tgt_mol_info_mem_gb, csc_mem_gb + ctrl_mol_info_mem_gb)
        },
    }


def calculate_scrna_correlations(
    mol_info_fns,
    downsample_rates,
    target_genes_ids_csv,
    barcode_summary_csv,
    targeted_matrix_h5,
    parent_matrix_subset_h5,
    parent_matrix_wta_h5,
):
    """Calculate correlations: 1) between targeted and parent cells using targeted genes;

    2) between targeted genes using counts across cells.
    """

    np.random.seed(0)

    sc_rna_comparison = TargetedControlCompareSC(
        tgt_mol_info_fn=mol_info_fns[0],
        ctrl_mol_info_fn=mol_info_fns[1],
        target_gene_ids=target_genes_ids_csv,
        barcode_summary_fn=barcode_summary_csv,
        tgt_downsample=downsample_rates[0],
        ctrl_downsample=downsample_rates[1],
    )
    sc_rna_comparison.scrna_mat_to_h5(targeted_matrix_h5, sample_label=TARGETED_LABEL)
    sc_rna_comparison.scrna_mat_to_h5(parent_matrix_subset_h5, sample_label=CONTROL_LABEL)
    del sc_rna_comparison

    # write original parent matrix not subset to targeted genes
    sc_rna_comparison = TargetedControlCompareSC(
        tgt_mol_info_fn=None,
        ctrl_mol_info_fn=mol_info_fns[1],
        target_gene_ids=None,
        barcode_summary_fn=barcode_summary_csv,
        tgt_downsample=downsample_rates[0],
        ctrl_downsample=downsample_rates[1],
    )
    sc_rna_comparison.scrna_mat_to_h5(parent_matrix_wta_h5, sample_label=CONTROL_LABEL)


def join(args, outs, chunk_defs, chunk_outs):

    molecule_info_fns = [args.targeted_molecule_info, args.parent_molecule_info]
    with MoleculeCounter.open(args.targeted_molecule_info, "r") as mc:
        target_gene_ids = mc.feature_reference.get_target_feature_ids()
        target_gene_ids = [ensure_binary(gene_id) for gene_id in target_gene_ids]

    # subsample to raw reads, use full depth if downsample_sc_data is False
    # FIXME (maria.alexis): downsampling not currently in use, but should eventually subsample to rrpc
    if args.downsample_sc_data:
        downsample_rates = cr_tgt_cmp_utils.get_downsampling_rate_per_sample(
            molecule_info_fns, NORM_MODE_RAW
        )
    else:
        downsample_rates = [None, None]

    # nomenclature: parent_subset is parent subset to tgt genes, parent_wta is parent all genes
    targeted_matrix_h5 = martian.make_path(outs.targeted_matrix_h5)
    parent_matrix_subset_h5 = martian.make_path(outs.parent_matrix_subset_h5)
    parent_matrix_wta_h5 = martian.make_path(outs.parent_matrix_wta_h5)

    calculate_scrna_correlations(
        molecule_info_fns,
        downsample_rates,
        target_gene_ids,
        args.barcode_summary_csv,
        targeted_matrix_h5,
        parent_matrix_subset_h5,
        parent_matrix_wta_h5,
    )
