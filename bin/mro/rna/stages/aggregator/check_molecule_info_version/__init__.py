#!/usr/bin/env python
#
# Copyright (c) 2018 10x Genomics, Inc. All rights reserved.
#
"""Aggr preflight check + convert legacy molecule info h5 to current version."""

import os
import shutil

import martian

import cellranger.constants as cr_constants
import cellranger.molecule_counter as cr_mol_counter
from cellranger.aggr.preprocessing import check_molecule_info_version_split
from cellranger.analysis.constants import CBC_MAX_NCELLS
from cellranger.molecule_counter_converter import convert_v2_to_v4, upgrade_file

__MRO__ = """
stage CHECK_MOLECULE_INFO_VERSION(
    in  map[] sample_defs,
    out map[] updated_sample_defs,
    src py    "stages/aggregator/check_molecule_info_version",
) split using (
    in  int   mol_h5_version,
    in  map   sample_def,
    out map   updated_sample_def,
)
"""
SP_PRODUCT = "sp"


def split(args):
    return check_molecule_info_version_split(args)


def main(args, outs):
    outs.updated_sample_def = args.sample_def.copy()
    outs.updated_sample_def[cr_constants.AGG_ORIGINAL_H5_FIELD] = args.sample_def[
        cr_constants.AGG_H5_FIELD
    ]

    if args.mol_h5_version == cr_mol_counter.CURR_FILE_VERSION:
        # avoid copying, pass it along
        return

    v2_mole_info_h5 = args.sample_def[cr_constants.AGG_H5_FIELD]
    v2_file_basename = os.path.basename(v2_mole_info_h5)
    out_mole_info_h5 = martian.make_path(v2_file_basename)

    if args.mol_h5_version == 2:
        convert_v2_to_v4(v2_mole_info_h5, out_mole_info_h5)
    else:
        shutil.copy(v2_mole_info_h5, out_mole_info_h5)

    try:
        with cr_mol_counter.MoleculeCounter.open(out_mole_info_h5, "r+") as mc:
            upgrade_file(mc)
    except ValueError as err:
        martian.exit(str(err))

    outs.updated_sample_def[cr_constants.AGG_H5_FIELD] = out_mole_info_h5


def join(args, outs, chunk_defs, chunk_outs):
    outs.updated_sample_defs = [chunk_out.updated_sample_def for chunk_out in chunk_outs]

    if any("batch" in sample_def for sample_def in outs.updated_sample_defs):
        ncells = 0
        for sample_def in outs.updated_sample_defs:
            with cr_mol_counter.MoleculeCounter.open(
                sample_def[cr_constants.AGG_H5_FIELD], "r"
            ) as mc:
                ncells += mc.get_num_filtered_barcodes_for_library(0)
        if ncells > CBC_MAX_NCELLS:
            martian.exit(
                "You provided {:,} cells in total, but chemistry batch correction only supports up to {:,} cells.".format(
                    ncells, CBC_MAX_NCELLS
                )
            )
