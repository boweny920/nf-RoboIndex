#!/usr/bin/env python
#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#
"""Constants used for targeted-compare."""


from __future__ import annotations

from collections import namedtuple

import cellranger.websummary.plotly_tools as pltly

# Column names, labels, and colors used by targeted-compare
TARGETED_LABEL = "in_targeted_sample"
CONTROL_LABEL = "in_parent_sample"
CONTROL_WTA_LABEL = "parent_all_genes"

TARGETED_COMPARE_LABELS = [TARGETED_LABEL, CONTROL_LABEL]
TARGETED_COMPARE_COLORS = {TARGETED_LABEL: "red", CONTROL_LABEL: "blue"}

TARGETED_COMPARE_ENRICHMENT_COL = "read_enrichment"
TARGETED_COMPARE_IS_ENRICHED_COL = "is_enriched"
TARGETED_COMPARE_RECOVERY_COL = "umi_recovery"
TARGETED_COMPARE_IS_TGT_GENE_COL = "is_targeted"
TARGETED_COMPARE_CELL_CALL_GRP_COL = "cell_call_category"

# targeted-compare normalization modes used in pipeline
NORM_MODE_RAW = "raw_counts"
NORM_MODE_NONE = "none"
NORM_MODE_TARGETED = "targeted_counts"
NORM_MODES_TARGETED_COMPARISON = [NORM_MODE_TARGETED, NORM_MODE_RAW, NORM_MODE_NONE]
COUNT_MODE_ALL = "in_all_genes"
COUNT_MODE_TARGETED = "in_targeted_genes"
COUNT_MODES_TARGETED_COMPARISON = [COUNT_MODE_TARGETED, COUNT_MODE_ALL]

# cell-calling groups used for targeted-compare
PARENT_ONLY = "parent-only"
TARGETED_ONLY = "targeted-only"
NEITHER = "neither"
BOTH = "both"
CellCallGroup = namedtuple("cell_call_group", ["name", "label", "color", "title"])
CELL_CALL_CATEGORIES = [
    (CellCallGroup(BOTH, "TGT + PAR", "purple", "Barcodes in Both Samples")),
    (CellCallGroup(TARGETED_ONLY, "TGT only", "red", "Barcodes in Targeted Only")),
    (CellCallGroup(PARENT_ONLY, "PAR only", "blue", "Barcodes in Parent Only")),
    (CellCallGroup(NEITHER, "neither", "purple", "Barcodes Unobserved")),
]

TARGETED_COMPARE_CLUSTER_COLORS = [
    # 10 default plotly colors
    "#1F77B4",
    "#FF7F0E",
    "#2CA02C",
    "#D62728",
    "#9467BD",
    "#8C564B",
    "#E377C2",
    "#7F7F7F",
    "#BCBD22",
    "#17BECF",
    # rainbow colors for novel clusters
    "#FF0000",
    "#CCFF00",
    "#00FF66",
    "#0066FF",
    "#CC00FF",
    "#FF4D00",
    "#80FF00",
    "#00FFB2",
    "#001AFF",
    "#FF00E6",
    "#FF9900",
    "#33FF00",
    "#00FFFF",
    "#3300FF",
    "#FF0099",
    "#FFE500",
    "#00FF19",
    "#00B3FF",
    "#7F00FF",
    "#FF004D",
]

# Help text for plots and tables
MIN_COUNTS = 1

_READ_ENRICHMENT_HELP_TEXT = {
    "title": "Per Gene Read Enrichment",
    "data": [
        [
            "Left: Correlation of per Gene Read Counts in Targeted and Parent Samples",
            [
                "The plot shows the number of {cell_or_spot}-associated reads per targeted gene in the targeted and parent samples. "
                "Pearson correlation of log10 read counts is shown in the upper left. "
                "The targeted and parent samples are normalized by the total number of read pairs across samples to account for sequencing depth."
            ],
        ],
        [
            "Right: Per Gene Read Enrichment Histogram",
            [
                "Histogram of the read enrichment (log2) per targeted gene. "
                "Enrichment is defined as the number of reads observed in the targeted sample, divided by that in the parent sample, after normalizing samples "
                "to the total number of read pairs in either sample. "
                "Only genes with at least "
                + str(MIN_COUNTS)
                + " {cell_or_spot}-associated UMI in the parent and targeted samples are shown."
            ],
        ],
    ],
}

_UMI_RECOVERY_HELP_TEXT = {
    "title": "Per Gene Sensitivity",
    "data": [
        [
            "Left: Correlation of per Gene UMI Counts in Targeted and Parent Samples",
            [
                "The plot shows the number of {cell_or_spot}-associated UMIs per targeted gene in the targeted and parent samples. "
                "Pearson correlation of log10 UMI counts is shown in the upper left."
            ],
        ],
        [
            "Right: Histogram of the Ratio of Targeted UMIs Recovered per Gene Relative to Parent",
            [
                "Histogram of the number of UMIs in the targeted sample divided by the number of UMIs in the parent sample, per targeted gene. "
                "Only genes with at least "
                + str(MIN_COUNTS)
                + " {cell_or_spot}-associated UMI in the parent sample are shown."
            ],
        ],
    ],
}

_TSNE_SELECTOR_HELP = {
    "title": "t-SNE Projection of {caps_cell_or_spot}s, Colored by Cluster",
    "data": [
        [
            "",
            [
                "The axes correspond to the 2-dimensional embedding produced by the t-SNE algorithm. In this space, pairs of {cell_or_spot}s that are close to each other have "
                "more similar gene expression profiles than {cell_or_spot}s that are distant from each other. t-SNE projections are drawn according to the selection in the "
                "dropdown on the right. If selecting 'Parent - All Genes', 'Parent - Targeted Genes', or 'Targeted', the same t-SNE projection will be shown in all three plots, while if selecting"
                " 'Per Sample', the three different t-SNEs will be shown side-by-side. {alt_caps_cell_or_spot}s are always colored by their clustering in that sample, irrespective of xy coordinates."
                " Barcodes identified as {cell_or_spot}s in either sample are shown. {alt_caps_cell_or_spot}s in the parent sample with zero targeted UMI counts are excluded from the 'Parent - Targeted genes' results. The display is limited to a random subset of {cell_or_spot}s."
            ],
        ],
        ["", ["Left: Parent sample using all genes"]],
        ["", ["Middle: Parent sample using only targeted genes"]],
        ["", ["Right: Targeted sample"]],
    ],
}


def _recursive_str_format(value, kwargs):
    """Recursively format a thing.

    Args:
        value: What values to change
        kwargs: The change to be made

    Returns:
        Strings "Cell" or "Tissue covered spot" depending on the pipeline being run
    """

    if isinstance(value, str):
        return value.format(**kwargs)
    elif isinstance(value, dict):
        keys = value.keys()
        for k in keys:
            value[k] = _recursive_str_format(value[k], kwargs)
        return value
    elif isinstance(value, list):
        for i, j in enumerate(value):
            value[i] = _recursive_str_format(j, kwargs)
        return value
    else:
        return value


class TargetCompareConstants:  # pylint: disable=too-few-public-methods
    """A collection of constants."""

    def __init__(self, cell_or_spot="cell", caps_cell_or_spot="Cell", alt_caps_cell_or_spot="Cell"):
        kwargs = {
            "cell_or_spot": cell_or_spot,
            "caps_cell_or_spot": caps_cell_or_spot,
            "alt_caps_cell_or_spot": alt_caps_cell_or_spot,
        }
        self.read_enrichment_help_text = _recursive_str_format(_READ_ENRICHMENT_HELP_TEXT, kwargs)
        self.umi_recovery_help_text = _recursive_str_format(_UMI_RECOVERY_HELP_TEXT, kwargs)
        self.tsne_selector_help = _recursive_str_format(_TSNE_SELECTOR_HELP, kwargs)


CELL_CALLING_HELP_TEXT = {
    "title": "Cell Calling Comparison",
    "data": [
        [
            "Left: Barcode Rank Plot",
            [
                "The plot shows the count of UMIs in targeted genes mapped to each barcode. "
                "The targeted and parent samples are normalized to the total number of targeted gene expression read pairs. "
                "Tooltip shows the fraction of cell barcodes in a plot region that are called as cells in 1) the targeted and parent samples (TGT + PAR), "
                "2) the targeted sample only (TGT only), 3) the parent sample only (PAR only), or 4) neither sample (neither)."
            ],
        ],
        [
            "Right: UMIs per Barcode in the Targeted and Parent Samples",
            [
                "The plot shows the number of UMIs mapped to each barcode in the targeted and parent samples. "
                "Barcodes are colored by whether they are called as cells in one sample, both samples, or neither."
            ],
        ],
    ],
}

SPATIAL_UMI_COMPARISON_HELP_TEXT = {
    "title": "UMI Comparison",
    "data": [
        [
            "UMIs per Barcode in the Targeted and Parent Samples",
            [
                "The plot shows the number of UMIs mapped to each barcode in the targeted and parent samples. "
                "Barcodes are colored by whether they are called as tissue-covered spot in one sample, both samples, or neither."
            ],
        ],
    ],
}

TSNE_PLOT_COMMON = {
    "config": pltly.PLOT_CONFIG,
    "layout": {
        "hovermode": False,
        "autosize": True,
        "margin": {"l": 5, "r": 5},
        "width": 300,
        "height": 430,
        "legend": {"font": {"size": 11}, "orientation": "h"},
        "xaxis": {"showticklabels": False, "fixedrange": False},
        "yaxis": {"showticklabels": False, "fixedrange": False},
    },
}
