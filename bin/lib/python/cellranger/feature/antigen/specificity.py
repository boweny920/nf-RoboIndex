#!/usr/bin/env python
#
# Copyright (c) 2022 10X Genomics, Inc. All rights reserved.
#
#
"""Generate antigen specificity scores and antigen assignments."""
from __future__ import annotations

import csv
from ast import literal_eval
from collections import Counter, OrderedDict, defaultdict, namedtuple
from enum import Enum
from typing import AnyStr, Optional

import numpy as np
from scipy.stats import beta

import cellranger.rna.library as rna_library
import tenkit.stats as tk_stats
from cellranger.report import PercentMetric  # pylint: disable=no-name-in-module

SIGNAL_PRIOR = 1
NOISE_PRIOR = 3

FEATURE_SEPARATOR = "|"
UNASSIGNED = "Unassigned"
BLANK = "Blank"

CANONICAL_VDJ_GENE_PAIRS = ["TRA_TRB", "IGH_IGK", "IGH_IGL"]

ANTIGEN_SPECIFICITY_CSV_HEADER = [
    "barcode",
    "antigen",
    "antigen_umi",
    "control",
    "control_umi",
    "antigen_specificity_score",
    "mhc_allele",
    "raw_clonotype_id",
    "exact_subclonotype_id",
]

AssignmentPerCell = namedtuple(
    "AssignmentPerCell",
    [
        "barcode",
        "num_antigens",
        "assignment",
        "specificity_score",
        "antigen_umi",
        "control_umi",
        "clonotype_id",
        "exact_subclonotype_id",
    ],
)

Concordance = namedtuple(
    "Concordance",
    [
        "clonotype_key",
        "size",
        "canonical_pair",
        "assigned_antigen",
        "num_bcs_with_assigned_antigen",
        "concordance",
    ],
)


class ClonalGroupLevel(Enum):
    """Specifies grouping level for clonotypes."""

    # pylint:disable=invalid-name
    CLONOTYPE = "clonotypes"
    EXACT_SUBCLONOTYPE = "exact_subclonotypes"


class CellsCategory(Enum):
    """Specifies category of cells."""

    # pylint:disable=invalid-name
    ALL_CELLS = "cells"
    GEX_ONLY_CELLS = "gex_only_cells"


# pylint: disable=too-few-public-methods
class AssignmentMetadata:
    """Metadata for feature assignments."""

    def __init__(
        self,
        cells_name: CellsCategory,
        num_cells,
        num_cells_blanks,
        num_cells_unassigned,
        num_cells_singlets,
        num_cells_multiplets,
    ):
        self.cells_name = cells_name
        self.num_cells = num_cells
        self.num_cells_blanks = num_cells_blanks
        self.frac_cells_without_antigen_umis = tk_stats.robust_divide(
            self.num_cells_blanks, self.num_cells
        )
        self.num_cells_unassigned = num_cells_unassigned
        self.frac_cells_unassigned = tk_stats.robust_divide(
            self.num_cells_unassigned, self.num_cells
        )
        self.num_cells_singlets = num_cells_singlets
        self.frac_singlets = tk_stats.robust_divide(self.num_cells_singlets, self.num_cells)
        self.num_cells_multiplets = num_cells_multiplets
        self.frac_multiplets = tk_stats.robust_divide(self.num_cells_multiplets, self.num_cells)


class AntigenAssigner:
    """Assigns antigens to barcodes."""

    def __init__(
        self,
        data: dict[str, BarcodeAS],
    ):
        self.data = OrderedDict(sorted(data.items()))
        self._assignment_per_cell = None

    @property
    def assignment_per_cell(self):
        """Returns a dict indexed by bc containing metadata on features assigned."""

        if self._assignment_per_cell is None:
            result: dict[str, AssignmentPerCell] = OrderedDict()
            for bc_as in self.data.values():
                assignment = bc_as.get_assigned_antigen()
                if assignment in [UNASSIGNED, BLANK]:
                    result[bc_as.barcode] = AssignmentPerCell(
                        bc_as.barcode,
                        0,
                        assignment,
                        np.nan,
                        np.nan,
                        np.nan,
                        bc_as.clonotype_id,
                        bc_as.exact_subclonotype_id,
                    )
                else:
                    antigens = assignment.split(FEATURE_SEPARATOR)
                    num_antigens = len(antigens)
                    specificity_scores = FEATURE_SEPARATOR.join(
                        [
                            str(round(s, 3))
                            for a, s in bc_as.specificity_scores.items()
                            if a in antigens
                        ]
                    )
                    antigen_umi = FEATURE_SEPARATOR.join(
                        [str(u) for a, u in bc_as.antigens.items() if a in antigens]
                    )
                    ctrl_umi = FEATURE_SEPARATOR.join(
                        [
                            str(bc_as.controls[c])
                            for a, c in bc_as.antigen_to_control.items()
                            if a in antigens
                        ]
                    )
                    result[bc_as.barcode] = AssignmentPerCell(
                        bc_as.barcode,
                        num_antigens,
                        assignment,
                        specificity_scores,
                        antigen_umi,
                        ctrl_umi,
                        bc_as.clonotype_id,
                        bc_as.exact_subclonotype_id,
                    )
            self._assignment_per_cell = result
        return self._assignment_per_cell

    def get_assignment_metadata(self, cells_name: CellsCategory):
        """Returns metadata for feature assignments for cells or gex_only cells."""
        if cells_name == CellsCategory.ALL_CELLS:
            data = self.data
        elif cells_name == CellsCategory.GEX_ONLY_CELLS:
            data = {k: v for k, v in self.data.items() if v.clonotype_id == "None"}

        num_cells = len(data)
        num_cells_blanks = len(
            [v for k, v in self.assignment_per_cell.items() if v.assignment == BLANK and k in data]
        )
        num_cells_unassigned = len(
            [
                v
                for k, v in self.assignment_per_cell.items()
                if v.assignment == UNASSIGNED and k in data
            ]
        )
        num_cells_singlets = len(
            [v for k, v in self.assignment_per_cell.items() if v.num_antigens == 1 and k in data]
        )
        num_cells_multiplets = (
            num_cells - num_cells_blanks - num_cells_unassigned - num_cells_singlets
        )
        return AssignmentMetadata(
            cells_name,
            num_cells,
            num_cells_blanks,
            num_cells_unassigned,
            num_cells_singlets,
            num_cells_multiplets,
        )

    def get_antigen_assignment_metrics(
        self,
        cells_name: CellsCategory,
    ) -> dict[str, float]:
        """Compute summary metrics on antigen assignments."""

        report_prefix = rna_library.get_library_type_metric_prefix(rna_library.ANTIGEN_LIBRARY_TYPE)
        assignment_metadata = self.get_assignment_metadata(cells_name)

        name_fr_ce_w_no_features = report_prefix + f"frac_{cells_name.value}_blank_antigen"

        frac_cells_with_no_features = assignment_metadata.frac_cells_without_antigen_umis

        name_fr_ce_no_features = report_prefix + "frac_{}_unassigned_antigen".format(
            cells_name.value
        )

        frac_cells_unassigned_features = assignment_metadata.frac_cells_unassigned

        name_fr_ce_w_sg_features = report_prefix + "frac_{}_with_single_antigen".format(
            cells_name.value
        )
        frac_cells_with_single_features = assignment_metadata.frac_singlets

        name_fr_ce_w_mult_features = report_prefix + "frac_{}_with_multiple_antigen".format(
            cells_name.value
        )
        frac_cells_multiple_features = assignment_metadata.frac_multiplets

        res = {
            name_fr_ce_w_no_features: frac_cells_with_no_features,
            name_fr_ce_no_features: frac_cells_unassigned_features,
            name_fr_ce_w_sg_features: frac_cells_with_single_features,
            name_fr_ce_w_mult_features: frac_cells_multiple_features,
            f"{report_prefix}{cells_name.value}_num_blanks": assignment_metadata.num_cells_blanks,
            f"{report_prefix}{cells_name.value}_num_unassigned": assignment_metadata.num_cells_unassigned,
            f"{report_prefix}{cells_name.value}_num_singlets": assignment_metadata.num_cells_singlets,
            f"{report_prefix}{cells_name.value}_num_multiplets": assignment_metadata.num_cells_multiplets,
        }
        return res

    @staticmethod
    def load_from_file(filename: AnyStr):
        """Load the data from an antigen_specificity_scores.csv file.

        Args:
            filename: the file to load from

        Returns:
            a new AntigenAssigner object
        """
        barcode_map = defaultdict(BarcodeAS)
        CsvRow = namedtuple("CsvRow", ANTIGEN_SPECIFICITY_CSV_HEADER)
        if filename is None:
            return AntigenAssigner({})
        with open(filename) as csvfile:
            csvreader = csv.reader(csvfile)
            next(csvreader)
            for row in map(CsvRow._make, csvreader):
                if row.barcode in barcode_map:
                    barcode_map[row.barcode].update_barcode(
                        control={row.control: int(row.control_umi)},
                        antigens={row.antigen: int(row.antigen_umi)},
                        scores={row.antigen: float(row.antigen_specificity_score)},
                        allele=row.mhc_allele,
                    )
                else:
                    exact_subclonotype_id = (
                        "None" if row.raw_clonotype_id == "None" else row.exact_subclonotype_id
                    )
                    barcode_map[row.barcode] = BarcodeAS(
                        barcode=row.barcode,
                        clonotype_id=row.raw_clonotype_id,
                        exact_subclonotype_id=exact_subclonotype_id,
                        control={row.control: int(row.control_umi)},
                        antigens={row.antigen: int(row.antigen_umi)},
                        specificity_scores={row.antigen: float(row.antigen_specificity_score)},
                        allele=row.mhc_allele,
                    )
            return AntigenAssigner(barcode_map)

    def write_antigen_specificity_csv(self, path: AnyStr):
        """Create antigen specificity CSV."""
        with open(path, "w") as csv_handle:
            csv_writer = csv.writer(csv_handle, delimiter=",")
            csv_writer.writerow(ANTIGEN_SPECIFICITY_CSV_HEADER)
            for bc_as in self.data.values():
                bc = [bc_as.barcode] * len(bc_as.antigens)
                antigen = bc_as.antigens.keys()
                ag_umi = bc_as.antigens.values()
                ctrl = bc_as.antigen_to_control.values()
                ctrl_umi = [bc_as.controls[c] for c in ctrl]
                allele = [bc_as.seen_antigens[a] for a in bc_as.antigens]
                score = [round(score, 3) for score in bc_as.specificity_scores.values()]
                c_id = [bc_as.clonotype_id] * len(bc_as.antigens)
                esc_id = [bc_as.exact_subclonotype_id] * len(bc_as.antigens)
                csv_writer.writerows(
                    zip(bc, antigen, ag_umi, ctrl, ctrl_umi, score, allele, c_id, esc_id)
                )

    def write_antigen_assignment_csv(self, path: AnyStr):
        """Create antigen assignemnt CSV."""
        with open(path, "w") as csv_handle:
            csv_writer = csv.writer(csv_handle, delimiter=",")
            csv_writer.writerow(AssignmentPerCell._fields)
            csv_writer.writerows(self.assignment_per_cell.values())

    def generate_cells_by_clonotype(self, grouped_by: ClonalGroupLevel, clonotypes_csv: AnyStr):
        """Returns CellsPerClonotype grouped by clonotype_id."""
        per_clonotype_dict = defaultdict(list[BarcodeAS])
        for barcode_as in self.data.values():
            if grouped_by == ClonalGroupLevel.CLONOTYPE:
                key = barcode_as.clonotype_id
            elif grouped_by == ClonalGroupLevel.EXACT_SUBCLONOTYPE:
                key = (
                    "None"
                    if barcode_as.clonotype_id == "None"
                    else "_".join([barcode_as.clonotype_id, barcode_as.exact_subclonotype_id])
                )
            per_clonotype_dict[key].append(barcode_as)
        return CellsPerClonotype.from_dictionary(per_clonotype_dict, grouped_by, clonotypes_csv)


class BarcodeAS:
    """Data structure representing antigen counts and specificity scores for a single barcode."""

    control_by_allele: dict[str, str] = defaultdict(str)  # allele:control
    seen_antigens: dict[str, str] = defaultdict(str)  # antigen:allele

    __slots__ = (
        "barcode",
        "clonotype_id",
        "exact_subclonotype_id",
        "controls",
        "_antigens",
        "_specificity_scores",
        "_assignments",
    )

    def __init__(
        self,
        barcode: bytes,
        clonotype_id: str,
        exact_subclonotype_id: str,
        control: dict[str, int],  # key:value = control_id:umi
        antigens: dict[str, int],  # key:value = feature_id:umi
        allele: str,
        specificity_scores: Optional[dict[str, float]] = None,  # key:value = feature_id:score
    ):
        self.barcode = barcode
        self.clonotype_id = clonotype_id
        self.exact_subclonotype_id = exact_subclonotype_id

        # sanity checks
        assert len(control) == 1
        assert set(control.keys()).isdisjoint(set(antigens.keys()))
        if specificity_scores:
            assert antigens.keys() == specificity_scores.keys()
        self.controls = control
        self._antigens = antigens
        self._specificity_scores = specificity_scores
        self._assignments = None

        self.update_controls(allele, control)
        self.update_antigens(allele, antigens)

    @classmethod
    def update_controls(cls, allele, control):
        """Update class variable: control_by_allele."""

        if allele in cls.control_by_allele:
            assert all(c == cls.control_by_allele[allele] for c in control)
        else:
            cls.control_by_allele[allele] = list(control.keys())[0]

    @classmethod
    def update_antigens(cls, allele, antigens):
        """Update class variable: seen_antigens."""
        for antigen in antigens:
            if antigen in cls.seen_antigens:
                assert cls.seen_antigens[antigen] == allele
            else:
                cls.seen_antigens[antigen] = allele

    @property
    def antigen_to_control(self):
        ags = OrderedDict(sorted(self.seen_antigens.items()))
        return {antigen: self.control_by_allele[allele] for antigen, allele in ags.items()}

    @property
    def antigens(self):
        """Generates an orderd dict of antigens:umi_counts for all observed antigens."""
        ags = OrderedDict((a, 0) for a in self.antigen_to_control)
        if list(self._antigens.keys()) != list(ags.keys()):
            for antigen, umi_count in self._antigens.items():
                ags[antigen] = umi_count
            self._antigens = ags
        return self._antigens

    @property
    def specificity_scores(self):
        """Generates an orderd dict of antigens:scores for all observed antigens."""
        if self._specificity_scores is None:
            self._specificity_scores = self.calculate_antigen_specificity()
        elif list(self._specificity_scores.keys()) != list(self.antigens.keys()):
            scores = OrderedDict((a, 0.0) for a in self.antigens)
            for antigen, score in self._specificity_scores.items():
                scores[antigen] = score
            self._specificity_scores = scores
        return self._specificity_scores

    def assignments(self, threshold=75):
        """Assigns antigens with specificity score above a threshold to this barcode."""
        if self._assignments is None:
            self._assignments = {k: v >= threshold for k, v in self.specificity_scores.items()}
        return self._assignments

    def calculate_antigen_specificity(self):
        """Calculate specificity scores for each antigen from the beta-distribution."""
        noise = [self.controls[c] for c in self.antigen_to_control.values()]
        signal = self.antigens.values()
        scores = [
            (1 - beta.cdf(0.925, S + SIGNAL_PRIOR, N + NOISE_PRIOR)) * 100
            for S, N in zip(signal, noise)
        ]
        return dict(zip(list(self.antigens.keys()), scores))

    def get_assigned_antigen(self):
        """Returns antigen(s) assigned to this barcode as a string."""
        if not any(self.assignments().values()):
            if sum(self.antigens.values()) == 0:
                return BLANK
            else:
                return UNASSIGNED
        else:
            return FEATURE_SEPARATOR.join([k for k, v in self.assignments().items() if v])

    def update_barcode(self, control, antigens, allele, scores=None):
        """Updates a BarcodeAS instance with additional antigen data."""
        assert len(control) == 1
        assert all(a not in self._antigens.keys() for a in antigens)
        assert set(control.keys()).isdisjoint(set(antigens.keys()))
        if scores:
            assert antigens.keys() == scores.keys()

        self._antigens.update(antigens)
        self.controls.update(control)
        if scores:
            self._specificity_scores.update(scores)
        self.update_controls(allele, control)
        self.update_antigens(allele, antigens)


class CellsPerClonotype:
    """Stores a list of BarcodeAS objects per Clonotype or Exact subclonotype.

    Basically a dictionary like:

            "Clonotype1": [BarcodeAS(object),BarcodeAS(object),...],
            "Clonotype2": [BarcodeAS(object),BarcodeAS(object),...],
    """

    def __init__(self, grouped_by: ClonalGroupLevel, clonotypes_csv: AnyStr):
        """Initialize like a dictionary."""
        self.grouped_by = grouped_by
        self.clonotypes_csv = clonotypes_csv
        self._data: dict[str, list[BarcodeAS]] = {}
        self._num_chains_per_clonotype = None
        self._concordance_per_clonotype = None

    def __getitem__(self, item: str):
        """Returns a list for a given assignment."""
        return self._data[item]

    def __setitem__(self, key: str, value: list[BarcodeAS]):
        assert isinstance(value, list)
        self._data[key] = value

    def __delitem__(self, key: str):
        del self._data[key]

    def __iter__(self):
        return iter(self._data)

    def __len__(self):
        return len(self._data)

    def clonotype_size(self, key):
        return len(self._data[key])

    def clonotype_concordance(self, key):
        """Returns clonotype concordance for a given clonotype."""
        assignments = []
        for bc in self._data[key]:
            assignments.append(bc.get_assigned_antigen())
        assignments = [UNASSIGNED if a == BLANK else a for a in assignments]
        max_antigen = Counter(sorted(assignments)).most_common(1)[0]
        size = self.clonotype_size(key)
        concordance = tk_stats.robust_divide(max_antigen[1], size)

        clonotype_id = (
            key.split("_")[0] if self.grouped_by == ClonalGroupLevel.EXACT_SUBCLONOTYPE else key
        )
        is_canonical_pair = (
            self.num_chains_per_clonotype[clonotype_id] if clonotype_id != "None" else False
        )
        return Concordance(
            key, size, is_canonical_pair, max_antigen[0], max_antigen[1], concordance
        )

    @property
    def num_chains_per_clonotype(self):
        """Returns a dict indexed by clonotype containing metadata on num chains."""
        if self._num_chains_per_clonotype is None:
            result: dict[str, bool] = {}
            with open(self.clonotypes_csv) as infile:
                csvreader = csv.reader(infile)
                seen_header = False
                for row in csvreader:
                    if not seen_header:
                        assert row[0] == "clonotype_id"
                        assert row[3] == "cdr3s_aa"
                        seen_header = True
                        continue
                    chains = "_".join(sorted(genes.split(":")[0] for genes in row[3].split(";")))
                    result[row[0]] = chains in CANONICAL_VDJ_GENE_PAIRS
            self._num_chains_per_clonotype = result
        return self._num_chains_per_clonotype

    @property
    def concordance_per_clonotype(self):
        """Returns a dict indexed by clonotype containing metadata on antigen assigned and concordance."""

        if self._concordance_per_clonotype is None:
            result: dict[str, Concordance] = {}
            for clonotype in self._data:
                result.update({clonotype: self.clonotype_concordance(clonotype)})
            order_by = {
                k: [
                    literal_eval(i) if i != "None" else 0
                    for i in k.removeprefix("clonotype").split("_")
                ]
                for k in result
            }
            result = OrderedDict(sorted(result.items(), key=lambda kv: order_by[kv[0]]))
            self._concordance_per_clonotype = result
        return self._concordance_per_clonotype

    def get_clonotype_concordance_metrics(self):
        """Compute summary metrics on clonotype (or exact_subclonotype) concordance."""
        report_prefix = rna_library.get_library_type_metric_prefix(rna_library.ANTIGEN_LIBRARY_TYPE)

        name_median_concordance_gt9 = report_prefix + "median_concordance_of_{}_size_gt9".format(
            self.grouped_by.value
        )
        name_min_concordance_gt9 = report_prefix + "lowest_concordance_of_{}_size_gt9".format(
            self.grouped_by.value
        )
        concordance_gt9 = [
            v.concordance
            for k, v in self.concordance_per_clonotype.items()
            if v.size > 9 and k != "None"
        ]

        name_median_concordance_gt9_canonical_pair = (
            report_prefix + f"median_concordance_of_{self.grouped_by.value}_size_gt9_canonical_pair"
        )
        name_min_concordance_gt9_canonical_pair = (
            report_prefix + f"lowest_concordance_of_{self.grouped_by.value}_size_gt9_canonical_pair"
        )
        concordance_gt9_canonical_pair = [
            v.concordance
            for k, v in self.concordance_per_clonotype.items()
            if v.size > 9 and k != "None" and v.canonical_pair
        ]

        name_aggregate_concordance = report_prefix + "aggregate_concordance_of_{}".format(
            self.grouped_by.value
        )

        name_aggregate_concordance_canonical_pair = (
            report_prefix + f"aggregate_concordance_of_{self.grouped_by.value}_canonical_pair"
        )

        aggregate_concordance = PercentMetric()
        aggregate_concordance_canonical_pair = PercentMetric()
        for clonotype, value in self.concordance_per_clonotype.items():
            if clonotype != "None":
                aggregate_concordance.add_value(value.num_bcs_with_assigned_antigen, value.size)
                if value.canonical_pair:
                    aggregate_concordance_canonical_pair.add_value(
                        value.num_bcs_with_assigned_antigen, value.size
                    )

        return {
            name_median_concordance_gt9: np.median(concordance_gt9) if concordance_gt9 else np.nan,
            name_min_concordance_gt9: np.min(concordance_gt9) if concordance_gt9 else np.nan,
            name_median_concordance_gt9_canonical_pair: np.median(concordance_gt9_canonical_pair)
            if concordance_gt9_canonical_pair
            else np.nan,
            name_min_concordance_gt9_canonical_pair: np.min(concordance_gt9_canonical_pair)
            if concordance_gt9_canonical_pair
            else np.nan,
            name_aggregate_concordance: aggregate_concordance.report(),
            name_aggregate_concordance_canonical_pair: aggregate_concordance_canonical_pair.report(),
        }

    @staticmethod
    def from_dictionary(
        data: dict[str, list[BarcodeAS]], grouped_by: ClonalGroupLevel, clonotypes_csv: AnyStr
    ):
        """Loads from a dictionary."""
        return_value = CellsPerClonotype(grouped_by, clonotypes_csv)
        for clonotype, bc_list in data.items():
            return_value[clonotype] = bc_list
        return return_value

    def write_clonotype_concordance_csv(self, path: AnyStr):
        """Create clonotype concordance CSV."""
        with open(path, "w") as csv_handle:
            csv_writer = csv.writer(csv_handle, delimiter=",")
            csv_writer.writerow(Concordance._fields)
            for clonotype in self.concordance_per_clonotype.values():
                csv_writer.writerow(clonotype)
