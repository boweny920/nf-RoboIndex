#!/usr/bin/env python
#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#
"""Utilities for dealing with targeted and features."""
from __future__ import annotations

import csv
import math
from collections import Counter, OrderedDict
from copy import deepcopy
from typing import Optional  # pylint: disable=import-error,unused-import

import numpy as np
import pandas as pd
import scipy.sparse as sp_sparse
from six import ensure_binary, ensure_str

import cellranger.matrix as cr_matrix
import cellranger.pandas_utils as pdu
import cellranger.simple_csv as simple_csv
import cellranger.targeted.simple_utils as simple_utils
import cellranger.utils as cr_utils
import cellranger.webshim.common as cr_webshim
import cellranger.webshim.constants.shared as shared_constants
import cellranger.websummary.plotly_tools as pltly
import tenkit.stats as tk_stats
from cellranger.molecule_counter import (
    BARCODE_IDX_COL_NAME,
    FEATURE_IDX_COL_NAME,
    GEM_GROUP_COL_NAME,
    MoleculeCounter,
)
from cellranger.pandas_utils import (
    FEATURE_DF_BARCODE_SEQ_COL,
    FEATURE_DF_COUNT_COL,
    FEATURE_DF_UMI_COL,
)
from cellranger.rna.library import GENE_EXPRESSION_LIBRARY_TYPE
from cellranger.targeted.targeted_compare_constants import (
    BOTH,
    CELL_CALL_CATEGORIES,
    CELL_CALLING_HELP_TEXT,
    CONTROL_LABEL,
    CONTROL_WTA_LABEL,
    COUNT_MODE_ALL,
    COUNT_MODE_TARGETED,
    NEITHER,
    NORM_MODE_NONE,
    NORM_MODE_RAW,
    SPATIAL_UMI_COMPARISON_HELP_TEXT,
    TARGETED_COMPARE_CELL_CALL_GRP_COL,
    TARGETED_COMPARE_CLUSTER_COLORS,
    TARGETED_COMPARE_COLORS,
    TARGETED_COMPARE_ENRICHMENT_COL,
    TARGETED_COMPARE_IS_TGT_GENE_COL,
    TARGETED_COMPARE_RECOVERY_COL,
    TARGETED_LABEL,
    TSNE_PLOT_COMMON,
    TargetCompareConstants,
)
from cellranger.webshim.common import convert_numpy_array_to_line_chart
from cellranger.webshim.constants.targetedcompare import tc_metrics
from cellranger.websummary.react_components import React10X

# pylint: disable=invalid-name
######################################################################
# Functions related to TargetedCompare Websummary
######################################################################

# used for making certain plots along both dimensions
BY_BARCODE = "barcode"
BY_GENE = "gene"
# make of WS key
TARGETED_COMPARE_WS = "targeted_compare_ws"
# bins for histograms
NUM_BINS = 40
# minimum counts to ensure meaningful enrichments
MIN_COUNTS = 1
TSNE_LABEL_TO_NAME_DICT = {
    CONTROL_WTA_LABEL: "Parent - All Genes",
    CONTROL_LABEL: "Parent - Targeted Genes",
    TARGETED_LABEL: "Targeted",
    "each": "Per Sample",
}


class TsneData(React10X):  # pylint: disable=too-few-public-methods
    """POD class for data on a tSNE."""

    def __init__(self, key, name, plot_json):
        super().__init__()
        self.key = key
        self.name = name
        self.plot_json = plot_json


class Tsnes(React10X):  # pylint: disable=too-few-public-methods
    """Stores similar data for a collection of tSNEs plot jsons."""

    def __init__(self):
        super().__init__()
        self.tsnes = []

    def add_tsne(self, tsne_data):
        assert isinstance(tsne_data, TsneData)
        self.tsnes.append(tsne_data)

    def to_dict_for_json(self):
        return {"tsnes": [tsne.to_dict_for_json() for tsne in self.tsnes]}


def make_new_tsne_plot(
    clusters_df,
    which_coords_label,
    which_clustering_label,
    cr_tgt_cmp_constants,
    is_spatial,
    novel_cluster_color_idx=0,
):  # pylint: disable=unused-argument
    """Builds the plot representing the dropdown of the tSNEs for the parent, parent subset.

    to targeted genes, and targeted samples.

    Args:
        clusters_df (pd.DataFrame): dataframe with xy coordinates and clustering of cells
            according done for each individual sample.
        which_coords_label (str): which tSNE (i.e. xy coordinates) to plot
        which_clustering_label (str): which clustering (i.e. colors overlaid) to plot
        cr_tgt_cmp_constants: ignored
        is_spatial: Whether to plot as spatial data.
        novel_cluster_color_idx: What color to start at.

    Returns:
        plot (dict): plotly plot that can be ingested by webshim.
    """

    if which_coords_label == "each":
        which_coords_label = which_clustering_label
    max_clusters = clusters_df[f"clusters.{which_clustering_label}"].max()
    n_cells = clusters_df.shape[0]
    data = []

    for cluster in range(-1, max_clusters + 1):
        if cluster == 0:
            continue
        subset_df = clusters_df[clusters_df[f"clusters.{which_clustering_label}"] == cluster]
        frac_cells_in_cluster = tk_stats.robust_divide(subset_df.shape[0], n_cells)
        x_colname = f"x.{which_coords_label}"
        y_colname = f"y.{which_coords_label}"

        subset_df = subset_df.drop_duplicates(subset=[x_colname, y_colname])
        np.random.seed(0)
        if subset_df.shape[0] == 0:
            continue
        subset_df = subset_df.sample(n=min(500, subset_df.shape[0]))

        if cluster == -1:
            # non-cell
            cluster_color = "lightgray"
        elif (
            which_clustering_label != CONTROL_LABEL
            and subset_df.iloc[0][f"specific_cluster.{which_clustering_label}"]
        ):
            # not a shared cluster, use a new color
            cluster_color = TARGETED_COMPARE_CLUSTER_COLORS[novel_cluster_color_idx]
            novel_cluster_color_idx += 1
        else:
            cluster_color = TARGETED_COMPARE_CLUSTER_COLORS[cluster - 1]

        if cluster != -1:
            name = f"{cluster} ({frac_cells_in_cluster:.1%})"
        else:
            if is_spatial:
                name = f"Discordant<br>Spots ({frac_cells_in_cluster:.1%})"
            else:
                name = f"Non-cell ({frac_cells_in_cluster:.1%})"

        data.append(
            {
                "name": name,
                "x": np.round(subset_df[x_colname], 2).tolist(),
                "y": np.round(subset_df[y_colname], 2).tolist(),
                "type": "scattergl",
                "mode": "markers",
                "marker": {
                    "opacity": 0.7,
                    "size": 4,
                    "color": cluster_color if cluster != -1 else "lightgray",
                    "line": {"width": 0, "color": cluster},
                },
                "showlegend": True,
            }
        )
    plot = deepcopy(TSNE_PLOT_COMMON)
    plot["layout"]["title"] = TSNE_LABEL_TO_NAME_DICT[which_clustering_label]
    plot["data"] = data
    return plot


class TsneSelector(React10X):
    """Mimics the data requirements and model of the React Component in TsneSelector.js."""

    def __init__(self, plot_help_txt):
        """Initialize the data for a TsneSelector Component, which displays 3 tSNEs in a row: parent.

        sample (left), parent sample subset to targeted genes (middle), targeted sample (right).

        :param plot_help_txt: dict of data needed for the React component `DynamicHelptext` that
        will be displayed above the plots
        """
        super().__init__()
        self._left_plots = None
        self._middle_plots = None
        self._right_plots = None
        self.plot_help_txt = plot_help_txt

    def _validate_same_tsnes(self):
        """The plots/tables should all be arrays with the same tSNE represented at each position in them."""
        to_check = [
            x
            for x in [self._left_plots, self._middle_plots, self._right_plots]
            if x and isinstance(x, Tsnes)
        ]
        if len(to_check) > 1:
            baseline = to_check[0].tsnes
            for alt in to_check[1:]:
                if any(x.key != y.key or x.name != y.name for x, y in zip(baseline, alt.tsnes)):
                    raise ValueError(
                        "The tsnes for the plots of a tsne selector must be the same and "
                        "in order."
                    )

    @property
    def left_plots(self):
        return self._left_plots

    @left_plots.setter
    def left_plots(self, value):
        """Parent (all genes) tSNEs."""
        assert isinstance(value, Tsnes)
        self._left_plots = value
        self._validate_same_tsnes()

    @property
    def middle_plots(self):
        """Parent subset to targeted genes tSNEs."""
        return self._right_plots

    @middle_plots.setter
    def middle_plots(self, value):
        assert isinstance(value, Tsnes)
        self._middle_plots = value
        self._validate_same_tsnes()

    @property
    def right_plots(self):
        return self._right_plots

    @right_plots.setter
    def right_plots(self, value):
        """Targeted tSNEs."""
        assert isinstance(value, Tsnes)
        self._right_plots = value
        self._validate_same_tsnes()

    def to_dict_for_json(self):
        l_plots = self._left_plots.to_dict_for_json()
        m_plots = self._middle_plots.to_dict_for_json()
        r_plots = self._right_plots.to_dict_for_json()
        return {
            "left_plots": l_plots,
            "middle_plots": m_plots,
            "right_plots": r_plots,
            "plot_help": self.plot_help_txt,
        }


def build_tsne_dropdown(clusters_df, cr_tgt_cmp_constants, help_text, is_spatial):
    """Builds the plot representing the dropdown of the tSNEs.

    Plots for the parent, parent subset to targeted genes, and targeted samples.

    Args:
        clusters_df (pd.DataFrame): dataframe with xy coordinates and clustering of cells
            according done for each individual sample.
        cr_tgt_cmp_constants: ignored
        help_text: dict of data needed for the React component `DynamicHelptext` that
            will be displayed above the plots
        is_spatial: is this a spatial sample or not

    Returns:
        plot (dict): plotly plot that can be ingested by webshim.
    """

    # init 3 different plot collections for the left (parent), middle (parent subset), and right
    # (targeted) tSNEs. Each Tsne will contain the plot info for all plots using the different
    # xy-coordinates that can be selected in the dropdown menu.
    left_tsnes = Tsnes()
    middle_tsnes = Tsnes()
    right_tsnes = Tsnes()

    # get cluster color index for parent_wta and targeted samples (parent subset sample will use
    # default colors regardless of whether it's a novel cluster). parent_wta novel cluster colors
    # will start at 10 (i.e. no offset from end of default colors), while targeted cluster colors
    # will start at 10 + num_novel_clusters_control_wta
    novel_cluster_color_index_control_wta = 10
    num_novel_clusters_control_wta = len(
        set(
            clusters_df.loc[clusters_df[f"specific_cluster.{CONTROL_WTA_LABEL}"]][
                f"clusters.{CONTROL_WTA_LABEL}"
            ]
        ).difference({-1})
    )
    novel_cluster_color_index_targeted = 10 + num_novel_clusters_control_wta

    for coords in [CONTROL_WTA_LABEL, CONTROL_LABEL, TARGETED_LABEL, "each"]:
        new_tsne = TsneData(
            coords,  # key
            TSNE_LABEL_TO_NAME_DICT[coords],  # name
            make_new_tsne_plot(
                clusters_df,
                coords,
                CONTROL_WTA_LABEL,
                cr_tgt_cmp_constants,
                is_spatial,
                novel_cluster_color_idx=novel_cluster_color_index_control_wta,
            ),  # plot_json
        )
        left_tsnes.add_tsne(new_tsne)
        new_tsne = TsneData(
            coords,  # key
            TSNE_LABEL_TO_NAME_DICT[coords],  # name
            make_new_tsne_plot(
                clusters_df,
                coords,
                CONTROL_LABEL,
                cr_tgt_cmp_constants,
                is_spatial,
                novel_cluster_color_idx=0,
            ),  # plot_json
        )
        middle_tsnes.add_tsne(new_tsne)
        new_tsne = TsneData(
            coords,  # key
            TSNE_LABEL_TO_NAME_DICT[coords],  # name
            make_new_tsne_plot(
                clusters_df,
                coords,
                TARGETED_LABEL,
                cr_tgt_cmp_constants,
                is_spatial,
                novel_cluster_color_idx=novel_cluster_color_index_targeted,
            ),  # plot_json
        )
        right_tsnes.add_tsne(new_tsne)

    tsne_select = TsneSelector(help_text)

    tsne_select.left_plots = left_tsnes
    tsne_select.middle_plots = middle_tsnes
    tsne_select.right_plots = right_tsnes

    return tsne_select.to_dict_for_json()


def build_tgtcmp_metrics_csv(
    websummary_data, fn, is_spatial, style_func=lambda *args: "", target_func=lambda *args: ""
):
    """Makes a CSV of targeted compare summary metrics."""

    tables = []
    metrics = tc_metrics(is_spatial)
    for table_dict in metrics:
        table_name = table_dict["name"]
        table_metrics = table_dict["metrics"]

        rows = []
        for metric_dict in table_metrics:
            name = metric_dict["name"]
            cr_webshim.add_table_rows(
                websummary_data, name, metric_dict, rows, style_func, target_func
            )

        if len(rows) > 0:
            tables.append({"name": table_name, "rows": rows, "headers": []})

    csv_metrics = OrderedDict()
    for table in tables:
        if not table:
            continue
        for metric, _, value in table["rows"]:
            if isinstance(metric, dict):
                metric = metric["v"]
            if isinstance(value, dict):
                value = value["v"]
            if metric not in csv_metrics:
                csv_metrics[metric] = value

    with open(fn, "w") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(list(csv_metrics.keys()))
        writer.writerow(list(csv_metrics.values()))


def make_targeted_metrics_table(summary_data, metadata, species_list):
    """Make a table of paired metrics for targeted and control sample for TargetedCompare WS."""

    metric_keys = [
        "num_targeted_genes_detected_in_parent_sample",
        "num_targeted_genes_enriched",
        "cell_targeted_depth_factor",
    ]

    metrics = metadata.gen_metric_list(summary_data, metric_keys, species_list)
    helptext = metadata.gen_metric_helptext(metric_keys)

    enrichment_table_rows = [[metric.name, metric.value_string] for metric in metrics]

    return {
        "help": {"title": "Targeted Enrichment", "data": helptext},
        "table": {"rows": enrichment_table_rows},
    }


def make_paired_metrics_table(summary_data, metadata, species_list):
    """Make a table of paired metrics for targeted and control sample for TargetedCompare WS."""

    metric_keys = [
        "total_read_pairs",
        "num_cells",
        "fraction_reads_on_target",
        "mean_reads_per_cell",
        "mean_targeted_reads_per_cell",
        "median_targeted_genes_per_cell",
        "total_targeted_genes_detected",
        "num_targeted_genes_detected_exclusive",
        "median_targeted_umis_per_cell",
    ]

    metrics_table_rows = []
    for metric in metric_keys:
        targeted_metric = metadata.gen_metric_list(
            summary_data, ["_".join([metric, TARGETED_LABEL])], species_list
        )
        control_metric = metadata.gen_metric_list(
            summary_data, ["_".join([metric, CONTROL_LABEL])], species_list
        )
        metric_name = (
            targeted_metric[0].name if len(targeted_metric) > 0 else control_metric[0].name
        )
        value_on_target = targeted_metric[0].value_string if len(targeted_metric) > 0 else ""
        value_off_target = control_metric[0].value_string if len(control_metric) > 0 else ""
        metrics_table_rows.append([metric_name, value_on_target, value_off_target])

    helptext = metadata.gen_metric_helptext(
        ["_".join([metric, TARGETED_LABEL]) for metric in metric_keys]
    )

    metrics_table = {
        "header": [
            "",
            "Targeted&nbsp;Sample",
            "Parent&nbsp;Sample",
        ],
        "rows": metrics_table_rows,
    }
    return {"help": {"title": "Paired Metrics", "data": helptext}, "table": metrics_table}


def add_barcode_scatter(merged_barcode_df, count_mode=COUNT_MODE_ALL):
    """Makes scatterplot of targeted vs control cell-barcodes.

    Points are colored by whether they are
    called as cells in one, both, or neither sample.

    Args:
        merged_barcode_df (pd.DataFrame): dataframe with counts per barcode in targeted
            and control samples.
        count_mode (str): one of targeted_umis or all_umis, specifies which set of UMIs to
            plot.

    Returns:
        plot (dict): plotly plot that can be ingested by webshim.
    """

    data = []

    # reverse trace plotting order
    for cell_call_category in CELL_CALL_CATEGORIES[::-1]:
        category = cell_call_category.name
        color = cell_call_category.color
        if category == NEITHER:
            continue
        subset_df = merged_barcode_df.loc[
            merged_barcode_df[TARGETED_COMPARE_CELL_CALL_GRP_COL] == category
        ]
        n = subset_df.shape[0]

        targeted_counts_colname = f"{FEATURE_DF_UMI_COL}_{TARGETED_LABEL}.{count_mode}"
        control_counts_colname = f"{FEATURE_DF_UMI_COL}_{CONTROL_LABEL}.{count_mode}"
        # avoid plotting duplicates -- useless, and triggers a plotly bug with toggling traces
        subset_df.drop_duplicates(
            subset=[targeted_counts_colname, control_counts_colname], inplace=True
        )
        data.append(
            {
                "name": category + " (n = " + str(n) + ")",
                "x": subset_df[control_counts_colname].tolist(),
                "y": subset_df[targeted_counts_colname].tolist(),
                "text": subset_df[FEATURE_DF_BARCODE_SEQ_COL].tolist(),
                "type": "scattergl",
                "mode": "markers",
                "marker": {"size": 5, "color": [color] * subset_df.shape[0], "opacity": 0.25},
            }
        )

    if count_mode == COUNT_MODE_TARGETED:
        x_axis_title = "UMIs in Targeted Genes per Barcode in Parent Sample"
        y_axis_title = "UMIs in Targeted Genes per Barcode in Targeted Sample."
    elif count_mode == COUNT_MODE_ALL:
        x_axis_title = "UMIs per Barcode in Parent Sample"
        y_axis_title = "UMIs per Barcode in Targeted Sample"
    else:
        raise ValueError(f"Unsupported count mode {count_mode}")

    plot = {
        "config": pltly.PLOT_CONFIG,
        "layout": {
            "showlegend": True,
            "legend": {"orientation": "h", "y": -0.25},
            "hovermode": "closest",
            "xaxis": {"type": "log", "title": x_axis_title, "fixedrange": False},
            "yaxis": {"type": "log", "title": y_axis_title, "fixedrange": False},
        },
        "data": data,
    }
    return plot


def add_barcode_violin(merged_barcode_df, count_mode=COUNT_MODE_ALL):
    """Makes scatterplot of targeted vs control cell-barcodes.

    Points are colored by whether they are
    called as cells in one, both, or neither sample.

    Args:
        merged_barcode_df (pd.DataFrame): dataframe with counts per barcode in targeted
            and control samples.
        count_mode (str): one of targeted_umis or all_umis, specifies which set of UMIs to
            plot.

    Returns:
        plot (dict): plotly plot that can be ingested by webshim.
    """

    data = []

    if count_mode == COUNT_MODE_TARGETED:
        y_axis_title = "UMIs in Targeted Genes per Barcode"
    elif count_mode == COUNT_MODE_ALL:
        y_axis_title = "UMIs per Barcode"
    else:
        raise ValueError(f"Unsupported count mode {count_mode}")

    # reverse trace plotting order
    for cell_call_category in CELL_CALL_CATEGORIES[::-1]:
        category = cell_call_category.name
        if category == NEITHER:
            continue
        subset_df = merged_barcode_df[
            merged_barcode_df[TARGETED_COMPARE_CELL_CALL_GRP_COL] == category
        ]

        targeted_counts_colname = f"{FEATURE_DF_UMI_COL}_{TARGETED_LABEL}.{count_mode}"
        control_counts_colname = f"{FEATURE_DF_UMI_COL}_{CONTROL_LABEL}.{count_mode}"
        # avoid plotting duplicates -- useless, and triggers a plotly bug with toggling traces
        subset_df.drop_duplicates(
            subset=[targeted_counts_colname, control_counts_colname], inplace=True
        )

        data.append(
            {
                "name": "Targeted",
                "y": subset_df[targeted_counts_colname].tolist(),
                "text": subset_df[FEATURE_DF_BARCODE_SEQ_COL].tolist(),
                "type": "violin",
                "box": {"visible": True},
                "line": {"color": "black"},
                "meanline_visible": True,
                "fillcolor": "#e74c3c",
                "opacity": 0.6,
            }
        )

        data.append(
            {
                "name": "Parent",
                "y": subset_df[control_counts_colname].tolist(),
                "text": subset_df[FEATURE_DF_BARCODE_SEQ_COL].tolist(),
                "type": "violin",
                "box": {"visible": True},
                "line": {"color": "black"},
                "meanline_visible": True,
                "fillcolor": "#3498db",
                "opacity": 0.6,
            }
        )

    plot = {
        "config": pltly.PLOT_CONFIG,
        "layout": {
            "title": {"text": "Distribution of UMIs Detected"},
            "showlegend": True,
            "xaxis": {"visible": False},
            "hovermode": "closest",
            "yaxis": {"title": y_axis_title},
        },
        "data": data,
    }
    return plot


def add_barcode_difference_table(barcode_summary_df):
    """Make a table of how many barcodes are in one sample.

    the other sample, both or none.
    """

    ref_table_rows = []

    category_counts = Counter(barcode_summary_df[TARGETED_COMPARE_CELL_CALL_GRP_COL])

    for cell_call_category in CELL_CALL_CATEGORIES:
        category = cell_call_category.name
        if category == NEITHER:
            continue
        counts = category_counts[category]
        title = cell_call_category.title
        ref_table_rows.append([title, counts])

    ref_table = {"header": ["Barcode Concordance"], "rows": ref_table_rows}
    return ref_table


def add_plots(
    ws_data,
    metrics,
    barcode_summary_df,
    feature_summary_df,
    clustering_df,
    is_spatial,
    cr_tgt_cmp_constants,
):
    """Add plots that will constitute the bulk of the TargetedCompare WS."""
    # initalize the class TargetCompareConstants
    if is_spatial:
        help_text_constants = TargetCompareConstants(
            cell_or_spot="tissue covered spot",
            caps_cell_or_spot="Tissue Covered Spot",
            alt_caps_cell_or_spot="Tissue covered spot",
        )
    else:
        help_text_constants = TargetCompareConstants(
            cell_or_spot="cell", caps_cell_or_spot="Cell", alt_caps_cell_or_spot="Cell"
        )

    ws_data["read_enrichment"] = {
        "left_plot": add_gex_scatterplot(
            feature_summary_df,
            rsquared=metrics["targeted_gene_read_rsquared"],
            normalization_mode=NORM_MODE_RAW,
            to_plot="Reads",
        ),
        "right_plot": add_enrichment_histogram(
            feature_summary_df,
            mean_ratio=metrics["mean_read_enrichment"],
            xlab="Ratio of read counts in targeted sample relative to parent sample, log2",
            to_plot="Reads",
        ),
        "help_text": help_text_constants.read_enrichment_help_text,
    }
    ws_data["umi_enrichment"] = {
        "left_plot": add_gex_scatterplot(
            feature_summary_df,
            rsquared=metrics["targeted_gene_umi_rsquared"],
            normalization_mode=NORM_MODE_NONE,
            to_plot="UMIs",
        ),
        "right_plot": add_enrichment_histogram(
            feature_summary_df,
            mean_ratio=metrics["mean_frac_UMIs_recovered_per_gene"],
            xlab="Ratio of UMI counts in targeted sample relative to parent sample",
            to_plot="UMIs",
        ),
        "help_text": help_text_constants.umi_recovery_help_text,
    }

    if is_spatial:
        ws_data["spatial_umi_comparison"] = {
            "plot": add_barcode_scatter(barcode_summary_df),
            "help_text": SPATIAL_UMI_COMPARISON_HELP_TEXT,
        }

    else:
        ws_data["cell_calling_comparison"] = {
            "left_plot": add_barcode_rank_plot(barcode_summary_df),
            "right_plot": add_barcode_scatter(barcode_summary_df),
            "help_text": CELL_CALLING_HELP_TEXT,
        }

    # if not enough cells and single-cell clustering was disabled
    if clustering_df is not None:
        ws_data["tsne_selector"] = build_tsne_dropdown(
            clustering_df,
            cr_tgt_cmp_constants,
            help_text_constants.tsne_selector_help,
            is_spatial,
        )
    return ws_data


def generate_alarms(metrics, metadata, species_list):
    """Add alarms for TargetedCompare WS."""
    alarm_keys = [
        "frac_cell_barcode_overlap",
        "frac_targeted_genes_enriched",
        "cell_targeted_depth_factor",
    ]

    alarms = metadata.gen_metric_list(metrics, alarm_keys, species_list)
    return [metric.alarm_dict for metric in alarms if metric.alarm_dict]


def add_tables(ws_data, target_set, mol_info_fns, metrics, metadata, species_list):
    """Wrapper adding tables to TargetedCompare WS."""

    sample_id = ws_data["sample"]["id"]
    sample_desc = ws_data["sample"]["description"]
    pipeline_version = ws_data["sample"]["pipeline_version"]

    ws_data["run_info_table"] = make_run_info_table(
        sample_id, sample_desc, pipeline_version, mol_info_fns, target_set
    )
    ws_data["paired_metrics_table"] = make_paired_metrics_table(metrics, metadata, species_list)
    ws_data["targeted_metrics"] = make_targeted_metrics_table(metrics, metadata, species_list)

    alarms = generate_alarms(metrics, metadata, species_list)
    ws_data["alarms"] = {"alarms": alarms}
    return ws_data


# pylint: disable=too-many-locals
def add_barcode_rank_plot(merged_barcode_df, count_mode=COUNT_MODE_ALL):
    """Plots cell-calling comparison plot.

    Plots targeted knee plot as usual, with tooltip specifying the fraction of cells that are called
    in both samples, one sample only, or neither along varying regions of the barcode rank plot.
    """

    # sort by UMI count and cell status in targeted sample
    merged_barcode_df.sort_values(
        by=[
            f"{pdu.FEATURE_DF_UMI_COL}_{TARGETED_LABEL}.{count_mode}",
            f"{pdu.MOL_INFO_CELL_COL}_{TARGETED_LABEL}",
        ],
        ascending=[False, False],
        inplace=True,
    )
    merged_barcode_df.dropna(
        subset=[f"{pdu.FEATURE_DF_UMI_COL}_{TARGETED_LABEL}.{count_mode}"], inplace=True
    )
    merged_barcode_df.reset_index(inplace=True, drop=True)

    targeted_counts = merged_barcode_df[
        f"{pdu.FEATURE_DF_UMI_COL}_{TARGETED_LABEL}.{count_mode}"
    ].values
    targeted_cell_calls = merged_barcode_df[f"{pdu.MOL_INFO_CELL_COL}_{TARGETED_LABEL}"]
    cell_call_groups = merged_barcode_df[TARGETED_COMPARE_CELL_CALL_GRP_COL].values

    data = []

    last_tgt_cell_idx = targeted_cell_calls[targeted_cell_calls].last_valid_index()
    last_tgt_barcode_idx = targeted_counts.shape[0]

    def _format_cell_call_group_props(cell_call_groups):
        n_per_category = Counter(cell_call_groups)
        n_total = float(np.sum(list(n_per_category.values())))
        hovertext = ""
        for cell_call_category in CELL_CALL_CATEGORIES:
            category = cell_call_category.name
            label = cell_call_category.label
            hovertext += f"{100 * n_per_category[category] / n_total:.0f}% {label}<br>"
        return hovertext

    segments = [
        (
            "cells",
            np.unique(
                np.round(
                    np.subtract(np.logspace(np.log10(1), np.log10(last_tgt_cell_idx + 1), 10), 1)
                ).astype(int)
            ),
        ),
        (
            "background",
            np.unique(
                np.round(
                    np.logspace(np.log10(last_tgt_cell_idx), np.log10(last_tgt_barcode_idx), 10)
                ).astype(int)
            ),
        ),
    ]
    for (group, subsegments) in segments:
        for i in range(len(subsegments) - 1):
            plot_rows = convert_numpy_array_to_line_chart(
                targeted_counts[subsegments[i] : subsegments[i + 1] + 1], int
            )
            # follows build_plot_data_dict logic, +1 because logplot
            data.append(
                {
                    "y": [row[1] for row in plot_rows],
                    "x": [row[0] + 1 + subsegments[i] for row in plot_rows],
                    "name": "Targeted Cells" if group == "cells" else "Background",
                    "legendgroup": "Targeted Cells" if group == "cells" else "Background",
                    "hoverinfo": "text",
                    "text": _format_cell_call_group_props(
                        cell_call_groups[subsegments[i] : subsegments[i + 1]]
                    ),
                    "type": "scattergl",
                    "mode": "lines",
                    "line": {
                        "color": TARGETED_COMPARE_COLORS[TARGETED_LABEL]
                        if group == "cells"
                        else "lightgray",
                        "width": shared_constants.BC_RANK_PLOT_LINE_WIDTH,
                    },
                    "showlegend": i == 0,
                }
            )

    plot = {
        "config": pltly.PLOT_CONFIG,
        "data": data,
        "layout": {
            "showlegend": True,
            "hovermode": "closest",
            "xaxis": {"type": "log", "title": "Barcodes", "fixedrange": False},
            "yaxis": {"type": "log", "title": "UMI counts", "fixedrange": False},
        },
    }

    return plot


def add_enrichment_histogram(merged_counts, mean_ratio: float, xlab, to_plot="Reads"):
    """Makes histogram of per-gene enrichments.

    Args:
        merged_counts (pd.DataFrame): dataframe with counts per gene in targeted
            and control samples, downsampled to various rates.
        mean_ratio: Ratio being plotted.
        xlab: x axis title.
        to_plot (str): specifies the normalization scheme used to plot
            GEX values and along which to subset the GEX table, default is raw_reads.

    Returns:
        plot (dict): plotly plot that can be ingested by webshim.
    """
    if to_plot == "Reads":
        enrichments = merged_counts[merged_counts[TARGETED_COMPARE_IS_TGT_GENE_COL]][
            TARGETED_COMPARE_ENRICHMENT_COL
        ]
        unit = "read enrichment"
    elif to_plot == "UMIs":
        enrichments = merged_counts[merged_counts[TARGETED_COMPARE_IS_TGT_GENE_COL]][
            TARGETED_COMPARE_RECOVERY_COL
        ]
        unit = "UMI ratio"
    else:
        raise ValueError(f"Unsupported plot key {to_plot}")

    # this will filter out genes only observed in tgt and control samples with +inf/-inf enrichments
    enrichments = enrichments[np.isfinite(enrichments)]

    # make histogram -- there might be some weird plotly issues happening but using the histogram
    # utility seems to work fine in a notebook but not in our websummaries
    bins = np.linspace(enrichments.min(), enrichments.max(), NUM_BINS + 1)
    counts = [0.0] * NUM_BINS
    for i in range(len(bins) - 1):
        counts[i] = enrichments[(enrichments >= bins[i]) & (enrichments < bins[i + 1])].shape[0]
    bin_midpoints = [np.mean(bins[i : i + 2]) for i in range(len(bins) - 1)]

    enrichment_hist = {
        "config": pltly.PLOT_CONFIG,
        "data": [
            {
                "type": "bar",
                "x": bin_midpoints,
                "y": counts,
                "text": [
                    "%.2f - %.2f (%d)" % (bins[i], bins[i + 1], counts[i])
                    for i in range(len(counts) - 1)
                ],
            }
        ],
        "layout": {
            "showlegend": False,
            "bargap": 0.1,
            "hovermode": "closest",
            "xaxis": {"title": xlab},
            "yaxis": {"title": "Number of targeted genes"},
            "title": f"Mean {unit} = {mean_ratio:.3f}",
        },
    }

    return enrichment_hist


def add_gex_scatterplot(
    merged_counts, rsquared: float, normalization_mode: str = NORM_MODE_RAW, to_plot: str = "Reads"
):
    """Makes scatterplot of gene normalized gene expression values in control and.

    targeted samples.

    Args:
        merged_counts (pd.DataFrame): dataframe with counts per gene in targeted and
                                      control samples, downsampled to various rates.
        rsquared (float): [TODO]
        normalization_mode (str): specifies the normalization scheme used to plot GEX values and
                                  along which to subset the GEX table, default is targeted_reads.
        to_plot (str): What to plot ("Reads" or "UMIs").

    Returns:
        plot (dict): plotly plot that can be ingested by webshim.
    """
    if to_plot == "Reads":
        COL = f"{FEATURE_DF_COUNT_COL}_cells"
    elif to_plot == "UMIs":
        COL = f"{FEATURE_DF_UMI_COL}_cells"
    else:
        raise ValueError(f"Unsupported plot key {to_plot}")

    merged_counts = merged_counts[merged_counts[TARGETED_COMPARE_IS_TGT_GENE_COL]]

    targeted_counts_col = f"{COL}_{TARGETED_LABEL}.{normalization_mode}"
    control_counts_col = f"{COL}_{CONTROL_LABEL}.{normalization_mode}"
    feature_name_col = "feature_name"

    # drop duplicates due to plotly issues
    subset_df = merged_counts[
        [targeted_counts_col, control_counts_col, feature_name_col]
    ].drop_duplicates(subset=[targeted_counts_col, control_counts_col])
    y = np.round(np.log10(subset_df[targeted_counts_col] + 1), 2).tolist()
    x = np.round(np.log10(subset_df[control_counts_col] + 1), 2).tolist()
    gene_names = subset_df[feature_name_col].tolist()

    xmin = np.min(x) if len(x) > 0 else 0
    ymin = np.min(y) if len(y) > 0 else 0
    xmax = np.max(x) if len(x) > 0 else 1
    ymax = np.max(y) if len(y) > 0 else 1

    gex_scatterplot = {
        "config": pltly.PLOT_CONFIG,
        "data": [
            {"type": "scattergl", "mode": "markers", "name": "", "x": x, "y": y, "text": gene_names}
        ],
        "layout": {
            "hovermode": "closest",
            "showlegend": False,
            "xaxis": {
                "title": f"{to_plot} per gene in parent sample (log10+1)",
                "fixedrange": False,
                "showline": False,
            },
            "yaxis": {
                "title": f"{to_plot} per gene in targeted sample (log10+1)",
                "fixedrange": False,
                "showline": False,
            },
            "annotations": [
                {
                    "x": (xmax - xmin) * 0.2 + xmin,
                    "y": (ymax - ymin) * 0.8 + ymin,
                    "text": "R<sup>2</sup> = %.3f" % rsquared,
                    "showarrow": False,
                    "font": {"size": 14},
                }
            ],
            "shapes": [
                {  # add x=y line
                    "type": "line",
                    "x0": min(xmin, ymin),
                    "y0": min(xmin, ymin),
                    "x1": max(xmax, ymax),
                    "y1": max(xmax, ymax),
                    "line": {"color": "gray", "dash": "dash", "width": 2},
                }
            ],
        },
    }
    return gex_scatterplot


def make_run_info_table(sample_id, sample_desc, pipeline_version, mol_info_fns, target_set_fn):
    """Make a table of basic run information such as target set and.

    reference transcriptome. Table is directly ingestedin websummary.
    """

    ref_table_rows = [
        ["Sample ID", sample_id],
        ["Sample Description", sample_desc],
        ["Targeted h5 Path", mol_info_fns[0]],
        ["Parent h5 Path", mol_info_fns[1]],
        ["Pipeline Version", pipeline_version],
    ]

    with MoleculeCounter.open(mol_info_fns[0], "r") as mc:
        reference_name = ",".join([g for g in mc.feature_reference.get_genomes() if g != ""])
        ref_table_rows.append(["Reference Name", reference_name])

    target_set_header, target_gene_ids, _ = simple_utils.parse_target_csv(target_set_fn)
    ref_table_rows.append(["Target Set", target_set_header["panel_name"]])
    ref_table_rows.append(["Number of Targeted Genes", len(target_gene_ids)])

    ref_table = {"header": ["Sample"], "rows": ref_table_rows}
    return ref_table


######################################################################
# Other utilities, mostly used for TargetedCompare pipeline
######################################################################

UNION = "union"
INTERSECTION = "intersection"


class TargetedControlCompareSC:
    """Useful class for targeted-comparison that does some comparisons using the.

    filtered-feature barcode matrix in both samples. Does basic operations to
    downsample and subset these matrices approrpiately and calculate necessary
    correlations.
    """

    def __init__(
        self,
        tgt_mol_info_fn,
        ctrl_mol_info_fn,
        target_gene_ids,
        barcode_summary_fn,
        tgt_downsample=None,
        ctrl_downsample=None,
    ):
        self.tgt_mol_info_fn = tgt_mol_info_fn
        self.ctrl_mol_info_fn = ctrl_mol_info_fn
        self.tgt_downsample = tgt_downsample
        self.ctrl_downsample = ctrl_downsample

        self.barcode_seqs = TargetedControlCompareSC.load_barcode_seqs(
            barcode_summary_fn, group=UNION
        )
        self.target_gene_ids = target_gene_ids

        self.tgt_sc_matrix = TargetedControlCompareSC.mol_info_to_coo(
            self.tgt_mol_info_fn,
            self.target_gene_ids,
            self.barcode_seqs,
            downsample=self.tgt_downsample,
        )
        self.ctrl_sc_matrix = TargetedControlCompareSC.mol_info_to_coo(
            self.ctrl_mol_info_fn,
            self.target_gene_ids,
            self.barcode_seqs,
            downsample=self.ctrl_downsample,
        )
        if tgt_mol_info_fn is not None and ctrl_mol_info_fn is not None:
            assert self.tgt_sc_matrix.shape == self.ctrl_sc_matrix.shape

    @staticmethod
    def load_barcode_seqs(barcode_summary_fn, group=INTERSECTION):
        """Read the merged barcode summary file. Select and store the cell-barcodes on.

        which to slice the matrices (either the union or the intersection across both
        samples).
        """
        if barcode_summary_fn is None:
            return None

        try:
            barcode_summary_df = simple_csv.read_simple_csv(
                ensure_str(barcode_summary_fn),
                dtype=bytes,
                usecols=[TARGETED_COMPARE_CELL_CALL_GRP_COL, FEATURE_DF_BARCODE_SEQ_COL],
            )
        except Exception as ex:
            raise RuntimeError("Couldn't read a properly formatted barcode csv") from ex
        try:
            if group == UNION:
                bcs_in_common = barcode_summary_df[
                    barcode_summary_df[TARGETED_COMPARE_CELL_CALL_GRP_COL] != b"neither"
                ]
            elif group == INTERSECTION:
                bcs_in_common = barcode_summary_df[
                    barcode_summary_df[TARGETED_COMPARE_CELL_CALL_GRP_COL] == b"both"
                ]
            else:
                raise ValueError(f"invalid group mode {group}")
        except Exception as ex:
            raise RuntimeError("Couldn't get barcode sets from barcode csv") from ex

        result = {ensure_binary(bc) for bc in bcs_in_common[FEATURE_DF_BARCODE_SEQ_COL]}
        return sorted(result)

    @staticmethod
    def mol_info_to_coo(
        mol_info_fn,
        target_gene_ids=None,
        barcode_seqs=None,
        downsample=None,
        return_mat_type="umis",
    ):
        """Downsample a molecule_info.h5 and return a sparse matrix subsetted.

        to the selected cell barcodes and targeted gene ids (in that order).
        """
        if mol_info_fn is None:
            return None

        with MoleculeCounter.open(mol_info_fn, "r") as mc:
            filter_library_idx = mc.get_library_indices_by_type()[GENE_EXPRESSION_LIBRARY_TYPE]

            feature_id_to_index_map = mc.feature_reference.id_map

            gem_wells = sorted(set(mc.get_gem_groups()))
            num_barcodes = len(mc.get_barcodes())
            all_barcodes = cr_utils.format_barcode_seqs(mc.get_barcodes(), gem_wells)
            barcode_seq_to_index_map = {bc: i for i, bc in enumerate(all_barcodes)}
            del all_barcodes

            # default to all genes and cell barcodes if none provided
            if barcode_seqs is None:
                barcode_seqs = mc.get_filtered_barcodes(
                    mc.get_barcode_info(),
                    mc.get_library_info(),
                    mc.get_barcodes(),
                    library_type=GENE_EXPRESSION_LIBRARY_TYPE,
                )
            if target_gene_ids is None:
                target_gene_ids = [
                    ensure_binary(fd.id)
                    for fd in mc.feature_reference.feature_defs
                    if fd.feature_type == GENE_EXPRESSION_LIBRARY_TYPE
                ]

            # get mapping of original gene and barcode indices to new ones
            orig_new_feature_index_map = {}
            for (idx, target_gene_id) in enumerate(target_gene_ids):
                if not isinstance(target_gene_id, bytes):
                    raise ValueError(
                        f"target_gene_id must be bytes, but was {type(target_gene_id)}"
                    )
                feat = feature_id_to_index_map[target_gene_id].index
                orig_new_feature_index_map[feat] = idx

            orig_new_barcode_index_map = {}
            for bc_idx, bc_seq in enumerate(barcode_seqs):
                if not isinstance(bc_seq, bytes):
                    raise ValueError("Barcode sequence was not of bytes type.")
                if bc_seq in barcode_seq_to_index_map:
                    orig_new_barcode_index_map[barcode_seq_to_index_map[bc_seq]] = bc_idx

        mol_info_df = pdu.mol_info_from_h5(
            mol_info_fn,
            with_umi=False,
            with_gem_group=True,
            with_library_idx=False,
            filter_library_idx=filter_library_idx,
            filter_feature_idx=list(orig_new_feature_index_map),
            filter_barcode_idx=list(orig_new_barcode_index_map),
            with_cell_call=False,
            downsample=downsample,
        )

        # offset barcode indices by gem wells
        mol_info_df[BARCODE_IDX_COL_NAME] = np.add(
            mol_info_df[BARCODE_IDX_COL_NAME],
            np.multiply(np.subtract(mol_info_df[GEM_GROUP_COL_NAME], 1), num_barcodes),
        )

        # remap barcode and gene indices in molecule info to their new ones
        if mol_info_df.shape[0] > 0:
            mol_info_df[FEATURE_IDX_COL_NAME] = remap_indices(
                mol_info_df[FEATURE_IDX_COL_NAME], orig_new_feature_index_map
            )
            mol_info_df[BARCODE_IDX_COL_NAME] = remap_indices(
                mol_info_df[BARCODE_IDX_COL_NAME], orig_new_barcode_index_map
            )

        if return_mat_type == "umis":
            ones = np.ones(mol_info_df.shape[0], dtype=cr_matrix.DEFAULT_DATA_DTYPE)
            matrix = sp_sparse.coo_matrix(
                (ones, (mol_info_df[FEATURE_IDX_COL_NAME], mol_info_df[BARCODE_IDX_COL_NAME])),
                shape=(len(orig_new_feature_index_map), len(barcode_seqs)),
            )
            del ones
        elif return_mat_type == "reads":
            matrix = sp_sparse.coo_matrix(
                (
                    mol_info_df[FEATURE_DF_COUNT_COL],
                    (mol_info_df[FEATURE_IDX_COL_NAME], mol_info_df[BARCODE_IDX_COL_NAME]),
                ),
                shape=(len(orig_new_feature_index_map), len(barcode_seqs)),
            )
        else:
            raise ValueError(f"Invalid count type for return matrix {return_mat_type}")
        del mol_info_df

        return matrix

    def compare_sc_expression(self, group_by=BY_BARCODE):
        """Computes the gene expression correlation.

        Args:
            group_by (str): Either `barcode` for per cell-barcode or `gene` for per gene.

        Returns:
            pd.DataFrame: The correlations.
        """
        assert group_by in [BY_BARCODE, BY_GENE]

        # init these variables in case no mol_info entries
        dense_mat_tgt, dense_mat_ctrl = None, None

        # axis over which to do matrix computations
        if group_by == BY_BARCODE:
            axis = 0
        elif group_by == BY_GENE:
            axis = 1
        else:
            raise ValueError(f"Unsupported grouping mode {group_by}")

        def vectorized_corr(X, Y, axis=axis):
            """Compute correlation between all pairs of rows (X_i, Y_i) or cols (X_j, Y_j),.

            dependent on axis.
            """
            X_mu = X.mean(axis=axis)
            Y_mu = Y.mean(axis=axis)
            numerator = np.mean(np.multiply(np.subtract(X, X_mu), np.subtract(Y, Y_mu)), axis=axis)
            denominator = np.multiply(X.std(axis=axis), Y.std(axis=axis))
            return np.squeeze(np.asarray(np.divide(numerator, denominator))).transpose()

        # corr has to go through a dense matrix, control the total num of entries in the m x n matrix
        # we "densify" in case there are lots of cells and this would get pretty big
        max_entries_to_densify = 5000000
        chunk_size = max(
            1, int(math.floor(float(max_entries_to_densify) / self.tgt_sc_matrix.shape[axis]))
        )
        chunk_starts = np.arange(0, self.tgt_sc_matrix.shape[abs(axis - 1)], chunk_size).tolist()
        correlations_df = []

        for chunk_start in chunk_starts:

            chunk_stop = chunk_start + chunk_size

            if group_by == BY_BARCODE:
                dense_mat_tgt = self.tgt_sc_matrix.tocsr()[:, chunk_start:chunk_stop].todense()
                dense_mat_ctrl = self.ctrl_sc_matrix.tocsr()[:, chunk_start:chunk_stop].todense()
            elif group_by == BY_GENE:
                dense_mat_tgt = self.tgt_sc_matrix.tocsr()[chunk_start:chunk_stop, :].todense()
                dense_mat_ctrl = self.ctrl_sc_matrix.tocsr()[chunk_start:chunk_stop, :].todense()

            correlations = pd.DataFrame(
                vectorized_corr(dense_mat_tgt, dense_mat_ctrl, axis=axis), columns=["correlations"]
            )
            if group_by == BY_BARCODE:
                correlations["barcode_seq"] = self.barcode_seqs[chunk_start:chunk_stop]
                correlations["UMIs_targeted"] = np.squeeze(np.asarray(dense_mat_tgt.sum(axis=axis)))
                correlations["UMIs_control"] = np.squeeze(np.asarray(dense_mat_ctrl.sum(axis=axis)))
            elif group_by == BY_GENE:
                correlations["gene_id"] = self.target_gene_ids[chunk_start:chunk_stop]
                correlations["UMIs_targeted"] = dense_mat_tgt.sum(axis=axis)
                correlations["UMIs_control"] = dense_mat_ctrl.sum(axis=axis)

            correlations_df.append(correlations)

        del dense_mat_tgt, dense_mat_ctrl
        correlation_df = pd.concat(correlations_df)

        return correlation_df

    @staticmethod
    def matrix_to_h5(csc_matrix, barcodes, feature_ref, outfile):
        """Saves a csc matrix to h5, removing zeroed out axes."""
        combined_matrix = cr_matrix.CountMatrix(feature_ref, barcodes, csc_matrix)
        combined_matrix, _, _ = combined_matrix.select_nonzero_axes()
        combined_matrix.save_h5_file(outfile)

    def scrna_mat_to_h5(self, outfile, sample_label):
        """Saves a CountMatrix object to h5.

        sample_label argument can be one of TARGETED_LABEL,
        CONTROL_LABEL, or "both" and specifies which matrix should be saved to file.
        """
        if sample_label == TARGETED_LABEL:
            concat_matrix = self.tgt_sc_matrix.tocsc()
            concat_barcodes = self.barcode_seqs
            with MoleculeCounter.open(self.tgt_mol_info_fn, "r") as mc:
                orig_feature_ref = mc.get_feature_ref()
        elif sample_label == CONTROL_LABEL:
            concat_matrix = self.ctrl_sc_matrix.tocsc()
            concat_barcodes = self.barcode_seqs
            with MoleculeCounter.open(self.ctrl_mol_info_fn, "r") as mc:
                orig_feature_ref = mc.get_feature_ref()
        elif sample_label == BOTH:
            concat_matrix = sp_sparse.hstack((self.tgt_sc_matrix, self.ctrl_sc_matrix))
            concat_matrix = concat_matrix.tocsc()
            # append additional suffix to barcodes to know which sample they are from
            concat_barcodes = np.hstack(
                [
                    ["-".join([bc, TARGETED_LABEL]) for bc in self.barcode_seqs],
                    ["-".join([bc, CONTROL_LABEL]) for bc in self.barcode_seqs],
                ]
            )
            with MoleculeCounter.open(self.tgt_mol_info_fn, "r") as mc:
                orig_feature_ref = mc.get_feature_ref()
        else:
            raise RuntimeError("can't save matrix %s to h5" % sample_label)

        if self.target_gene_ids is not None:
            target_gene_indices = [
                fd.index for fd in orig_feature_ref.feature_defs if fd.id in self.target_gene_ids
            ]
            new_feature_ref = orig_feature_ref.select_features(target_gene_indices)
        else:
            new_feature_ref = orig_feature_ref

        TargetedControlCompareSC.matrix_to_h5(
            concat_matrix, concat_barcodes, new_feature_ref, outfile
        )

    @staticmethod
    def get_genes_per_cell_from_coo_matrix(mat):
        """Gets detected genes per cell from sparse matrix."""
        # rows are now barcodes
        mat = mat.transpose().tocsr()
        genes_per_cell = []
        for row in mat:
            genes_per_cell.append(row.nonzero()[1].shape[0])
        return genes_per_cell


def remap_indices(orig_array, remap_dict):
    """Remaps a numpy array according to the key-value pairs in remap_dict.

    This is a more memory-efficient solution relative to np.vectorize, taken from
    https://stackoverflow.com/questions/47171356
    """
    remap_keys = np.array(list(remap_dict.keys()))
    remap_values = np.array(list(remap_dict.values()))
    sorted_idx = remap_keys.argsort()
    sorted_keys = remap_keys[sorted_idx]
    sorted_values = remap_values[sorted_idx]
    new_array = sorted_values[np.searchsorted(sorted_keys, orig_array)]
    return new_array


def get_downsampling_rate_per_sample(
    mol_info_fns: list[bytes], normalization_mode: str
) -> list[Optional[float]]:
    """Get downsampling rate for each sample from the mol_info, dependent on the normalization.

    mode.
    """
    if normalization_mode == NORM_MODE_NONE:
        return [1.0] * len(mol_info_fns)
    elif normalization_mode != NORM_MODE_RAW:
        raise ValueError(f"unsupported normalization mode {normalization_mode}")
    reads_per_sample: list[int] = []
    for mol_info_fn in mol_info_fns:
        with MoleculeCounter.open(mol_info_fn, "r") as mc:
            gex_library_indices = mc.get_library_indices_by_type()[GENE_EXPRESSION_LIBRARY_TYPE]
            raw_reads_per_lib = mc.get_raw_read_pairs_per_library()
            total_reads = np.sum(
                [reads for (li, reads) in enumerate(raw_reads_per_lib) if li in gex_library_indices]
            )
        reads_per_sample.append(total_reads)

    min_reads_per_sample = np.nanmin(reads_per_sample)
    if normalization_mode == NORM_MODE_NONE:
        downsampling_rate_by_sample = [1.0] * len(mol_info_fns)
    else:
        downsampling_rate_by_sample = [
            tk_stats.robust_divide(min_reads_per_sample, reads) for reads in reads_per_sample
        ]
    #  convert 1.0 to None to make sure pdu skips downsampling
    downsampling_rate_by_sample = [
        None if ds_rate == 1.0 else ds_rate for ds_rate in downsampling_rate_by_sample
    ]
    return downsampling_rate_by_sample
