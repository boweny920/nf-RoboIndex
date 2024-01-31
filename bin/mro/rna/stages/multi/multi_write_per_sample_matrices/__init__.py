#!/usr/bin/env python
#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#
"""Slice the filtered matrix into per sample matrices."""

import csv
import json

import martian

import cellranger.matrix as cr_matrix
import cellranger.rna.matrix as rna_matrix
import cellranger.utils as cr_utils

__MRO__ = """
struct SampleMatrices(
    string sample,
    h5     filtered_matrix_h5,
    path   filtered_matrix_mex,
    h5     raw_matrix_h5,
    path   raw_matrix_mex,
    csv    filtered_barcodes,
)

stage MULTI_WRITE_PER_SAMPLE_MATRICES(
    in  h5               matrix_h5,
    in  h5               raw_matrix_h5,
    in  csv              filtered_barcodes,
    in  json             sample_barcodes,
    in  json             sample_cell_barcodes,
    out SampleMatrices[] sample_matrices,
    src py               "stages/multi/multi_write_per_sample_matrices",
) split (
    in  string           sample,
    in  string[]         barcodes,
    in  string[]         cell_barcodes,
) using (
    mem_gb   = 4,
    volatile = strict,
)
"""


def split(args):
    mem_gb = 3 + cr_matrix.CountMatrix.get_mem_gb_from_matrix_h5(args.raw_matrix_h5, scale=1.5)

    with open(args.sample_cell_barcodes) as inf:
        sample_cell_barcodes = json.loads(inf.read())
    with open(args.sample_barcodes) as inf:
        sample_barcodes = json.loads(inf.read())
    samples = list(sample_barcodes.keys())

    chunks = []
    for sample in samples:
        barcodes = sample_barcodes[sample]
        cell_barcodes = sample_cell_barcodes[sample]

        chunks.append(
            {
                "sample": sample,
                "barcodes": barcodes,
                "cell_barcodes": cell_barcodes,
                "__mem_gb": mem_gb,
                "__threads": 1,
            }
        )

    return {
        "chunks": chunks,
        "join": {
            "__mem_gb": 3,
            "__threads": 1,
        },
    }


def main(args, outs):

    # make paths for new per-sample matrices (these are already filtered for target panel in case of targeted)
    sample_filtered_matrix_h5 = martian.make_path(
        f"{args.sample}_filtered_feature_barcode_matrix.h5"
    ).decode("utf8")
    sample_filtered_matrix_mex = martian.make_path(
        f"{args.sample}_filtered_feature_barcode_matrix"
    ).decode("utf8")
    sample_filtered_barcodes_csv = martian.make_path(f"{args.sample}_filtered_barcodes.csv").decode(
        "utf8"
    )

    # make paths for new per-sample (all genes) matrices with all barcodes,
    # for RTL-type multiplexing this includes non-cell barcodes
    sample_raw_matrix_h5 = martian.make_path(f"{args.sample}_raw_feature_barcode_matrix.h5").decode(
        "utf8"
    )
    sample_raw_matrix_mex = martian.make_path(f"{args.sample}_raw_feature_barcode_matrix").decode(
        "utf8"
    )

    # create the per-sample matrices
    filtered_matrix = cr_matrix.CountMatrix.load_h5_file(args.matrix_h5)
    raw_matrix = cr_matrix.CountMatrix.load_h5_file(args.raw_matrix_h5)

    # filter matrix barcodes for cell-associated sample barcodes, maintain original order
    cell_barcodes = {b.encode() for b in args.cell_barcodes}
    filtered_bcs = [b for b in filtered_matrix.bcs if b in cell_barcodes]
    filtered_matrix = filtered_matrix.select_barcodes_by_seq(filtered_bcs)

    # filter matrix barcodes for sample barcodes, maintain original order
    barcodes = {b.encode() for b in args.barcodes}
    sample_bcs = [b for b in raw_matrix.bcs if b in barcodes]
    raw_matrix = raw_matrix.select_barcodes_by_seq(sample_bcs)

    chemistry = cr_matrix.CountMatrix.load_chemistry_from_h5(args.matrix_h5)
    gem_groups = list({cr_utils.split_barcode_seq(bc)[1] for bc in barcodes})
    matrix_attrs = cr_matrix.make_matrix_attrs_count(args.sample, gem_groups, chemistry)

    # write matrices
    filtered_matrix.save_h5_file(
        sample_filtered_matrix_h5,
        extra_attrs=matrix_attrs,
        sw_version=martian.get_pipelines_version(),
    )
    raw_matrix.save_h5_file(
        sample_raw_matrix_h5,
        extra_attrs=matrix_attrs,
        sw_version=martian.get_pipelines_version(),
    )

    rna_matrix.save_mex(
        filtered_matrix, sample_filtered_matrix_mex, martian.get_pipelines_version()
    )
    rna_matrix.save_mex(raw_matrix, sample_raw_matrix_mex, martian.get_pipelines_version())

    # write the filtered barcodes file for sample
    sample_barcode_set = set(args.cell_barcodes)
    with open(args.filtered_barcodes) as csvfile, open(
        sample_filtered_barcodes_csv, "w"
    ) as outfile:
        for row in csv.reader(csvfile, delimiter=","):
            barcode = row[1]
            if barcode in sample_barcode_set:
                print(",".join(row), file=outfile)

    outs.sample_matrices = [
        {
            "sample": args.sample,
            "filtered_matrix_h5": sample_filtered_matrix_h5,
            "filtered_matrix_mex": sample_filtered_matrix_mex,
            "raw_matrix_h5": sample_raw_matrix_h5,
            "raw_matrix_mex": sample_raw_matrix_mex,
            "filtered_barcodes": sample_filtered_barcodes_csv,
        }
    ]


def join(args, outs, chunk_defs, chunk_outs):
    outs.sample_matrices = [c.sample_matrices[0] for c in chunk_outs]
