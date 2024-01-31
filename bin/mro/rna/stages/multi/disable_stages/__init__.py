#!/usr/bin/env python
#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#
"""A helper stage to determine whether or not to disable the legacy bam file (holding all reads)."""

from cellranger.matrix import CountMatrix
from cellranger.rna.library import MULTIPLEXING_LIBRARY_TYPE

__MRO__ = """
stage DISABLE_STAGES(
    in  bool  no_bam,
    in  bool  disable_multi,
    in  bool  is_pd,
    in  h5    raw_feature_bc_matrix,
    out bool  disable_legacy_bam,
    out bool  disable_sample_bams,
    out bool  disable_assign_tags,
    src py    "stages/multi/disable_stages",
) using (
    volatile = strict,
)
"""


def main(args, outs):
    # Don't make an uber bam in multi for CS code
    outs.disable_legacy_bam = args.no_bam or ((not args.is_pd) and (not args.disable_multi))

    # disable the new sample bam if no_bam was specified, or if it's not a multi run
    outs.disable_sample_bams = args.no_bam or args.disable_multi

    # Do we have multiplexing data present? Used to enable multiplexing runs in COUNT_PD pipelines for
    # R&D purposes
    feature_ref = CountMatrix.load_feature_ref_from_h5_file(args.raw_feature_bc_matrix)
    no_cmos = feature_ref.get_count_of_feature_type(MULTIPLEXING_LIBRARY_TYPE) == 0
    outs.disable_assign_tags = args.disable_multi and no_cmos
