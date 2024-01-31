#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#


from __future__ import annotations

import json
from collections import OrderedDict
from typing import Any

import cellranger.rna.library as rna_library


def get_library_type(sample_def):
    """Note: function modified to substitute library type rna_library.FEATURETEST_LIBRARY_TYPE with library type rna_library.MULTIPLEXING_LIBRARY_TYPE."""
    if sample_def.get(rna_library.LIBRARY_TYPE) == rna_library.FEATURETEST_LIBRARY_TYPE:
        return rna_library.MULTIPLEXING_LIBRARY_TYPE
    ltype = sample_def.get(rna_library.LIBRARY_TYPE)
    return ltype or None


def get_subsample_rate(sample_def):
    return sample_def.get("subsample_rate", 1.0)


def get_target_set_name(sample_def):
    # for untargeted case, target_set_name is set to null
    name = sample_def.get(rna_library.TARGET_SET_KEY)
    if name is None:
        name = rna_library.DEFAULT_TARGET_SETS[0]
    return name


def get_target_set(sample_def):
    return sample_def.get("target_set")


# Return pipeline_parameters_json from this sample bag.
# Migrate pipeline parameters whose name has changed
# for backward compatibility with older sample bags.
def get_params_json_overrides(sample_def: dict[str, Any]) -> dict[str, Any]:
    if "metadata" in sample_def:
        params_json = sample_def["metadata"].get("pipeline_parameters_json")
    else:
        params_json = sample_def.get("pipeline_parameters_json")
    if not params_json:
        return {}

    params = json.loads(params_json) if isinstance(params_json, str) else params_json

    # Migrate pipeline parameters whose name has changed.
    # The key is the previous name, and the value is the current name.
    MIGRATE_PIPELINE_PARAMETERS = {"enforce_library_concordance": "check_library_compatibility"}
    return {MIGRATE_PIPELINE_PARAMETERS.get(key, key): value for key, value in params.items()}


def get_gem_well(sample_def):
    """Get the GEM well from a sample def.

    Defaults to 1 if the gem_group is specified as None.

    Args:
      sample_def (dict): Sample def
    Returns:
      int: GEM well
    """
    gg = sample_def["gem_group"] or 1
    return int(gg)


def assign_library_ids(sample_defs, default_library_type=rna_library.DEFAULT_LIBRARY_TYPE):
    """Assign library ids to a list of sample defs if not given already.

    Each (gem group, library_type) must have a distinct library id. If missing,
    the library id becomes an increasing integer.

    Args:
      sample_defs (list of dict): Sample defs, some possibly already containing a library_id key
      default_library_type (str): Use if a library_type is None or empty
    Returns:
      list of str: Library IDs, one per sample def.
    """
    # Resulting library ID for each sample def
    sd_libs = []

    # Map tuple of (gem_group (int), library_type (str)) to library ID (int)
    library_to_id = OrderedDict()

    # Assign a library ID to each sample def (most commonly by generating a new, unique one)
    for sd in sample_defs:
        gem_group = get_gem_well(sd)
        lib_type = get_library_type(sd) or default_library_type
        lib_tuple = (gem_group, lib_type)  # library
        sd_libs.append(lib_tuple)  # library for each sample def

        # If a library ID wasn't given by the user, assign it an integer that maps
        # uniquely to (gem_group, library_type)
        default_id = library_to_id.get(lib_tuple, str(len(library_to_id)))

        lib_id = sd.get("library_id", default_id)

        # Check if this library ID is already assigned to a different library.
        if lib_id in set(library_to_id.values()):
            # Find the first colliding library
            other_lib = next(
                other_lib
                for other_lib, other_lib_id in library_to_id.items()
                if other_lib_id == lib_id
            )
            if other_lib != lib_tuple:
                raise ValueError(
                    'Library ID "%s" is already associated with library "%s." A library ID must identify exactly one library.'
                    % (lib_id, other_lib)
                )

        if lib_tuple not in library_to_id:
            library_to_id[lib_tuple] = lib_id

        elif library_to_id[lib_tuple] != lib_id:
            # Library already has an ID
            raise ValueError(
                'Library "%s" already has ID "%s." Cannot assign it the different ID "%s."'
                % (str(lib_tuple), library_to_id[lib_tuple], lib_id)
            )

    # remap library_ids to ints
    library_to_id = OrderedDict((tup, i) for i, tup in enumerate(library_to_id))

    # Map sample defs to library ids
    return [library_to_id[lib] for lib in sd_libs]
