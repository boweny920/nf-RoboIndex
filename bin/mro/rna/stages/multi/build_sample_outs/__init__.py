#!/usr/bin/env python
#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#
"""A helper stage to build CS outs tailored to the multiplexing strategy."""
import os
from cellranger.cr_io import hard_link

__MRO__ = """
stage BUILD_SAMPLE_OUTS(
    in  SampleSlfeOuts  sample_slfe_outs,
    in  path            rna_analysis,
    in  path            crispr_analysis,
    in  path            antibody_analysis,
    in  path            antigen_analysis,
    in  cloupe          cloupe,
    in  html            web_summary,
    in  csv             metrics_summary_csv,
    in  VdjOutputsCS    vdj_b_outs,
    in  VdjOutputsCS    vdj_t_outs,
    in  VdjOutputsCS    vdj_t_gd_outs,
    in  bool            is_rtl_multiplexed,
    in  BEAM_ANALYZER   beam_analyzer,
    out SampleOutputsCS sample_outs,
    src py              "stages/multi/build_sample_outs",
)
"""


def main(args, outs):
    count = {}
    count["analysis"] = hard_link(args.rna_analysis)
    # Both of these aggregate barcodes files are equivalent
    if args.antibody_analysis is not None and os.path.exists(
        os.path.join(args.antibody_analysis, "aggregate_barcodes.csv")
    ):
        count["aggregate_barcodes"] = hard_link(
            os.path.join(args.antibody_analysis, "aggregate_barcodes.csv")
        )
    elif args.antigen_analysis is not None and os.path.exists(
        os.path.join(args.antigen_analysis, "aggregate_barcodes.csv")
    ):
        count["aggregate_barcodes"] = hard_link(
            os.path.join(args.antigen_analysis, "aggregate_barcodes.csv")
        )
    else:
        count["aggregate_barcodes"] = None
    count["sample_cloupe"] = hard_link(args.cloupe)
    count["crispr_analysis"] = hard_link(args.crispr_analysis)

    sample_slfe_outs = args.sample_slfe_outs or {}
    count["feature_reference_csv"] = hard_link(sample_slfe_outs.get("feature_reference"))
    count["probe_set"] = hard_link(sample_slfe_outs.get("probe_set"))
    count["sample_alignments"] = hard_link(sample_slfe_outs.get("bam_file"))
    count["sample_alignments_index_bai"] = hard_link(sample_slfe_outs.get("bai_index_file"))
    count["sample_alignments_index_csi"] = hard_link(sample_slfe_outs.get("csi_index_file"))
    count["sample_filtered_barcodes_csv"] = hard_link(sample_slfe_outs.get("filtered_barcodes"))
    count["sample_filtered_feature_bc_matrix"] = hard_link(
        sample_slfe_outs.get("filtered_matrix_h5")
    )
    count["sample_filtered_feature_bc_matrix_mex"] = hard_link(
        sample_slfe_outs.get("filtered_matrix_mex")
    )
    count["sample_raw_feature_bc_matrix"] = (
        hard_link(sample_slfe_outs.get("raw_matrix_h5")) if args.is_rtl_multiplexed else None
    )
    count["sample_raw_feature_bc_matrix_mex"] = (
        hard_link(sample_slfe_outs.get("raw_matrix_mex")) if args.is_rtl_multiplexed else None
    )
    count["sample_molecule_info"] = hard_link(sample_slfe_outs.get("molecule_info"))
    count["target_panel"] = hard_link(sample_slfe_outs.get("target_panel"))

    sample_outs = {}
    sample_outs["count"] = count
    sample_outs["metrics_summary"] = hard_link(args.metrics_summary_csv)
    sample_outs["vdj_b"] = args.vdj_b_outs
    sample_outs["vdj_t"] = args.vdj_t_outs
    sample_outs["vdj_t_gd"] = args.vdj_t_gd_outs
    sample_outs["web_summary"] = hard_link(args.web_summary)

    if args.beam_analyzer is not None:
        sample_outs["antigen_analysis"] = {}
        for key in ["antigen_specificity_scores", "per_barcode"]:
            sample_outs["antigen_analysis"][key] = hard_link(args.beam_analyzer.get(key))
    else:
        sample_outs["antigen_analysis"] = None

    outs.sample_outs = sample_outs
