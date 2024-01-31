#!/usr/bin/env python
#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#
"""Calculate antigen specificity for Beam runs."""

import csv
import json
import os
from collections import defaultdict

import cellranger.matrix as cr_matrix
import cellranger.rna.library as rna_library
import tenkit.safe_json as tk_safe_json
from cellranger.feature.antigen.specificity import (
    AntigenAssigner,
    BarcodeAS,
    CellsCategory,
    ClonalGroupLevel,
)

__MRO__ = """
stage CALCULATE_ANTIGEN_SPECIFICITY(
    in  h5     filtered_feature_counts_matrix,
    in  string beam_mode,
    in  csv    filtered_contig_annotations,
    in  csv    clonotypes_csv,
    in  json   vdj_cell_barcodes,
    out csv    antigen_specificity_scores,
    out csv    antigen_assignment,
    out csv    clonotype_concordance,
    out csv    exact_subclonotype_concordance,
    out json   summary,
    src py     "stages/feature/antigen_specificity",
)
"""

TARGETING_ANTIGEN = "targeting_antigen"
ALLELE = "mhc_allele"
BEAM_AB = "beam_ab"
BEAM_T = "beam_t"


def get_control_antigens(fr_obj, beam_mode):
    """Build a map of the control feature id for each allele and a dict of feature id to allele."""
    control_by_allele = {}
    antigen_to_allele = {}
    for feature_def in fr_obj.feature_defs:
        if ALLELE not in fr_obj.all_tag_keys or beam_mode == BEAM_AB:
            allele = ""
        else:
            allele = feature_def.tags.get(ALLELE, "")
        if feature_def.tags.get(TARGETING_ANTIGEN, "").lower() == "false":
            control_by_allele[allele] = feature_def.id
        antigen_to_allele[feature_def.id] = allele
    return control_by_allele, antigen_to_allele


def build_clonotype_maps(filtered_contig_annotations):
    """Build a map of vdj cell barcode to clonotype_id / exact_subclonotype_id."""
    bc_clonotype_map = {}
    bc_exact_subclonotype_map = {}
    with open(filtered_contig_annotations) as csvfile:
        csvreader = csv.reader(csvfile)
        seen_header = False
        for row in csvreader:
            if not seen_header:
                assert row[0] == "barcode"
                assert row[28] == "raw_clonotype_id"
                assert row[30] == "exact_subclonotype_id"
                seen_header = True
                continue

            if row[0] in bc_clonotype_map and (row[28] is None or row[28] == ""):
                continue
            elif row[0] in bc_clonotype_map:
                assert bc_clonotype_map[row[0]] == row[28]
                assert bc_exact_subclonotype_map[row[0]] == row[30]
            else:
                bc_clonotype_map[row[0]] = row[28]
                bc_exact_subclonotype_map[row[0]] = row[30]
    return bc_clonotype_map, bc_exact_subclonotype_map


# pylint: disable=too-many-locals
def main(args, outs):
    filt_ag_bc_mat = cr_matrix.CountMatrix.load_h5_file(
        args.filtered_feature_counts_matrix
    ).select_features_by_type(rna_library.ANTIGEN_LIBRARY_TYPE)

    control_by_allele, antigen_to_allele = get_control_antigens(
        filt_ag_bc_mat.feature_ref, args.beam_mode
    )

    # exclude any on-target antigens that have zero UMI counts
    nonzero_features = list(filt_ag_bc_mat.get_counts_per_feature() > 0)
    feature_ids = filt_ag_bc_mat.get_feature_ids_by_type(rna_library.ANTIGEN_LIBRARY_TYPE)
    filtered_ids = [
        i for (i, v) in zip(feature_ids, nonzero_features) if v or i in control_by_allele.values()
    ]
    antigen_to_allele = {k: v for k, v in antigen_to_allele.items() if k in filtered_ids}
    filt_ag_bc_mat = filt_ag_bc_mat.select_features_by_ids(filtered_ids)

    # check that filtered_contig_annotations.csv is generated
    annotation_file = args.filtered_contig_annotations is not None or os.path.exists(
        args.filtered_contig_annotations
    )

    # Don't compute specificity scores if there are no off target antigens
    # or allele is not defined for beam-t or there are no non-zero on-target antigens
    # pylint: disable=too-many-boolean-expressions
    if (
        args.beam_mode is None
        or not control_by_allele
        or (args.beam_mode == BEAM_T and all(allele == "" for allele in control_by_allele))
        or not annotation_file
        or not any(i not in control_by_allele.values() for i in filtered_ids)
    ):
        outs.antigen_specificity_scores = None
        outs.antigen_assignment = None
        outs.clonotype_concordance = None
        outs.exact_subclonotype_concordance = None
        outs.summary = None
        return

    # load clonotype information for each barcode
    bc_clonotype_map, bc_exact_subclonotype_map = build_clonotype_maps(
        args.filtered_contig_annotations
    )

    barcode_map = defaultdict(BarcodeAS)
    for allele, control in control_by_allele.items():
        antigens = [k for k, v in antigen_to_allele.items() if v == allele]
        antigen_id_map = filt_ag_bc_mat.select_features_by_ids(antigens).feature_ids_map
        for bc in filt_ag_bc_mat.bcs:
            counts = (
                filt_ag_bc_mat.select_barcodes_by_seq([bc])
                .select_features_by_ids(antigens)
                .get_counts_per_feature()
            )
            antigen_umi_map = {
                f_id.decode(): counts[f_idx] for f_id, f_idx in antigen_id_map.items()
            }
            control_umi = antigen_umi_map.pop(control.decode())
            if bc.decode() in barcode_map:
                barcode_map[bc.decode()].update_barcode(
                    control={control.decode(): control_umi}, antigens=antigen_umi_map, allele=allele
                )
            else:
                clonotype_id = bc_clonotype_map.get(bc.decode(), "None")
                exact_subclonotype_id = (
                    "None" if clonotype_id == "None" else bc_exact_subclonotype_map[bc.decode()]
                )
                barcode_map[bc.decode()] = BarcodeAS(
                    barcode=bc.decode(),
                    clonotype_id=clonotype_id,
                    exact_subclonotype_id=exact_subclonotype_id,
                    control={control.decode(): control_umi},
                    antigens=antigen_umi_map,
                    allele=allele,
                )

    metrics = {}
    antigen_assignments = AntigenAssigner(barcode_map)
    antigen_assignments.write_antigen_specificity_csv(outs.antigen_specificity_scores)
    antigen_assignments.write_antigen_assignment_csv(outs.antigen_assignment)
    metrics.update(antigen_assignments.get_antigen_assignment_metrics(CellsCategory.ALL_CELLS))
    metrics.update(antigen_assignments.get_antigen_assignment_metrics(CellsCategory.GEX_ONLY_CELLS))

    clonotype_concordance = antigen_assignments.generate_cells_by_clonotype(
        ClonalGroupLevel.CLONOTYPE, args.clonotypes_csv
    )
    clonotype_concordance.write_clonotype_concordance_csv(outs.clonotype_concordance)
    metrics.update(clonotype_concordance.get_clonotype_concordance_metrics())

    exact_subclonotype_concordance = antigen_assignments.generate_cells_by_clonotype(
        ClonalGroupLevel.EXACT_SUBCLONOTYPE, args.clonotypes_csv
    )
    exact_subclonotype_concordance.write_clonotype_concordance_csv(
        outs.exact_subclonotype_concordance
    )
    metrics.update(exact_subclonotype_concordance.get_clonotype_concordance_metrics())

    with open(outs.summary, "w") as outfile:
        json.dump(tk_safe_json.json_sanitize(metrics), outfile, sort_keys=True, indent=4)
