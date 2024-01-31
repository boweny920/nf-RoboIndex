#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#
"""Data loading utilities which have limited dependencies on anything else."""

from __future__ import annotations

import json
import os

import pandas as pd

from cellranger import constants as cr_constants
from cellranger.spatial.pipeline_mode import PipelineMode, Product, SlideType

# String and other constants used in our spatial assay.
IMAGEX_LOWRES = "pxl_col_in_lowres"
IMAGEY_LOWRES = "pxl_row_in_lowres"

TISSUE_POSITIONS_HEADER = [
    "barcode",
    "in_tissue",
    "array_row",
    "array_col",
    "pxl_row_in_fullres",
    "pxl_col_in_fullres",
]

TISSUE_POSITIONS_HEADER_TYPES = {
    "barcode": "str",
    "in_tissue": "int32",
    "array_row": "int32",
    "array_col": "int32",
    "pxl_row_in_fullres": "int32",
    "pxl_col_in_fullres": "int32",
}

# Maximum high and low image dimension for display and, sometimes, processing
HIRES_MAX_DIM_DEFAULT = 2000
HIRES_MAX_DIM_DICT = {
    PipelineMode(Product.VISIUM, SlideType.XL): 4000,
    PipelineMode(Product.VISIUM_HD, SlideType.VISIUM_HD): 4000,
}
LORES_MAX_DIM = 600

# constants for dark_images mro parameter
DARK_IMAGES_NONE = 0
DARK_IMAGES_CHANNELS = 1
DARK_IMAGES_COLORIZED = 2
DISALLOWED_DARK_IMAGES_EXTENSION = [".png"]


def parse_slide_sample_area_id(slide_sample_area_id):
    """Given an input to the pipeline like V19L01-006-B1,.

    parse out slide sample id and area id
    """

    slide_sample_id, area_id = slide_sample_area_id[:-3], slide_sample_area_id[-2:]
    return slide_sample_id, area_id


def get_galfile_path(barcode_whitelist: str) -> str:
    """Given a barcode whitelist, return the path to the corresponding GAL file."""

    path_to_galfile = os.path.join(cr_constants.BARCODE_WHITELIST_PATH, barcode_whitelist + ".gal")
    return path_to_galfile


def read_from_json(filename):
    """Read from a given json file."""

    with open(filename) as json_file:
        data = json.load(json_file)

    return data


def get_scalefactors(scalefactors_fn: str) -> dict[str, float]:
    """Load the scale factors.

    Args:
        scalefactors_fn:

    Returns:
    """
    with open(scalefactors_fn) as scalefactors:
        return json.load(scalefactors)


def get_lowres_coordinates(tissue_positions_csv: str, scalefactors_json: str) -> pd.DataFrame:
    """Return a pandas data frame that is just like the tissue_positions_csv but has the lowres scaled image coordinates.

    Args:
        tissue_positions_csv (str): Path to the tissue_positions.csv
        scalefactors_json (str): Path to the scalefactors_json.json

    Returns:
        pd.DataFrame:
    """
    coords = read_tissue_positions_csv(tissue_positions_csv)

    # read in scalefactors json and adjust coords for downsampled image
    scalef = get_scalefactors(scalefactors_json)["tissue_lowres_scalef"]
    coords[IMAGEY_LOWRES] = coords["pxl_row_in_fullres"] * scalef
    coords[IMAGEX_LOWRES] = coords["pxl_col_in_fullres"] * scalef
    return coords


def read_tissue_positions_csv(tissue_positions_fn) -> pd.DataFrame:
    # output dir to search for a file name
    # raw data
    # file name
    """Read the tissue positions csv as a pandas dataframe.

    Args:
        tissue_positions_fn (str): Filename

    Returns:
        pd.DataFrame: Csv as a dataframe
    """
    # For backwards compatibility
    ## First check if the file has a header. If there are digits there is no header
    with open(tissue_positions_fn) as f:
        first_line = f.readline()

    no_header = any(map(str.isdigit, first_line))

    # Set the kwargs according to the header state
    kwargs = {"names": TISSUE_POSITIONS_HEADER} if no_header else {"header": 0}

    coords = pd.read_csv(
        tissue_positions_fn,
        **kwargs,
        dtype=TISSUE_POSITIONS_HEADER_TYPES,
        sep=",",
    )
    coords = coords.set_index("barcode")
    coords.index = coords.index.astype(bytes, copy=True)
    return coords


def get_remove_image_pages(loupe_alignment_file: str) -> set:
    """Set of image pages which should be skipped during tiling."""
    if loupe_alignment_file is None:
        return set()

    alignment_data = read_from_json(loupe_alignment_file)
    pages = alignment_data.get("removeImagePages", [])
    return set(pages)


def get_image_page_names(loupe_alignment_file: str) -> list:
    """Return channel names from alignment json."""
    if loupe_alignment_file is None:
        return []
    alignment_data = read_from_json(loupe_alignment_file)
    return alignment_data.get("imagePageNames", [])
