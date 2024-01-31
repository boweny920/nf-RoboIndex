#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#


from __future__ import annotations

from typing import TYPE_CHECKING, Union

import numpy as np
import pandas as pd
from six import ensure_binary, ensure_str

if TYPE_CHECKING:
    Path = Union[bytes, str]


def load_csv_columnnames(path: Path) -> list[str]:
    """Gets the headers of csv files with # as comments.

    Args:
        path (Union[str, bytes]): path of the csv

    Returns:
        List[str]: List of column names of probe_set.csv
    """
    return list(pd.read_csv(path, comment="#").columns)


def load_csv_rownames(csv_file):
    rownames = np.atleast_1d(
        pd.read_csv(
            ensure_str(csv_file), usecols=[0], converters={0: ensure_binary}, na_filter=False
        ).values.squeeze()
    )
    return rownames


def write_target_features_csv(path: Path, feature_ints, header=None):
    """Write a list of ints to a file."""
    with open(path, "w") as outf:
        if header is not None:
            outf.write(str(header) + "\n")
        for i in feature_ints:
            outf.write(f"{i}\n")


def read_target_features_csv(path, header=None):
    """Read feature IDs from a file."""
    if header is not None:
        feature_ids = pd.read_csv(ensure_str(path))
        feature_ids = feature_ids[header].tolist()
    else:
        feature_ids = pd.read_csv(ensure_str(path), usecols=[0], header=None)
        feature_ids = feature_ids[0].tolist()
    return feature_ids
