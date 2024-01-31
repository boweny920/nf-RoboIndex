#!/usr/bin/env python
#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#

"""Does targeted-compare preflight checks; checks for identical references and that both.

molecule_info h5 files come from a targeted and a non-targeted parent run respectively.
"""


import martian

import cellranger.env as cr_env
import cellranger.preflight as cr_preflight
import cellranger.reference as cr_reference
from cellranger.molecule_counter import (
    ANALYSIS_PARAMETERS_METRIC,
    INTRON_MODE_HISTORIC_DEFAULT,
    INTRON_MODE_PARAM,
    TOTAL_READS_METRIC,
    USABLE_READS_METRIC,
    MoleculeCounter,
)
from cellranger.targeted import simple_utils
from cellranger.targeted.targeted_constants import TARGETING_METHOD_HC

__MRO__ = """
stage TARGETED_COMPARE_PREFLIGHT(
    in  h5   targeted_molecule_info,
    in  h5   parent_molecule_info,
    in  json target_set,
    src py   "stages/targeted_compare/targeted_compare_preflight",
) using (
    volatile = strict,
)
"""


def run_preflight_checks(args):
    """Checks both molecule info files exist and that one is from a targeted.

    run and the other is not. Checks the integrity of the target panel file
    and requires that the target panel file and the metadata of the targeted
    molecule_info specify the same target feature ids.
    """

    cr_preflight.validate_targeted_compare_mol_info(
        args.targeted_molecule_info,
        expecting_targeted=True,
        required_metrics=["reference_fasta_hash", "reference_gtf_hash", "chemistry_endedness"],
        required_library_metrics=[TOTAL_READS_METRIC, USABLE_READS_METRIC],
    )
    cr_preflight.validate_targeted_compare_mol_info(
        args.parent_molecule_info,
        expecting_targeted=False,
        required_metrics=["reference_fasta_hash", "reference_gtf_hash", "chemistry_endedness"],
        required_library_metrics=[TOTAL_READS_METRIC, USABLE_READS_METRIC],
    )

    if args.targeted_molecule_info == args.parent_molecule_info:
        raise cr_preflight.PreflightException(
            "The same file was specified for both the targeted and parent input files."
        )

    # check targeted and parent are consistent in terms of chemistry and reference
    targeted_mol_info = MoleculeCounter.open(args.targeted_molecule_info, "r")
    parent_mol_info = MoleculeCounter.open(args.parent_molecule_info, "r")

    for (key, readable_name) in zip(
        ["reference_fasta_hash", "reference_gtf_hash", "chemistry_endedness"],
        ["genome reference", "annotation GTF", "endedness"],
    ):
        if targeted_mol_info.get_metric(key) != parent_mol_info.get_metric(key):
            if key == "chemistry_endedness":
                msg = (
                    "The targeted and parent molecule files have different "
                    "{readable_name}, but are required to be the same."
                ).format(readable_name=readable_name)
            else:
                msg = (
                    "The targeted and parent molecule files have different {readable_name}s."
                    "\nPlease re-run '{product} count'/'{product} count --target-panel' with uniform {readable_name} in"
                    " order to run targeted-compare."
                ).format(readable_name=readable_name, product=cr_env.product())
            raise cr_preflight.PreflightException(msg)

    targeted_mol_info_intron_mode = (
        targeted_mol_info.get_all_metrics()
        .get(ANALYSIS_PARAMETERS_METRIC, {})
        .get(INTRON_MODE_PARAM, INTRON_MODE_HISTORIC_DEFAULT)
    )
    parent_mol_info_intron_mode = (
        parent_mol_info.get_all_metrics()
        .get(ANALYSIS_PARAMETERS_METRIC, {})
        .get(INTRON_MODE_PARAM, INTRON_MODE_HISTORIC_DEFAULT)
    )
    if targeted_mol_info_intron_mode != parent_mol_info_intron_mode:
        msg = (
            "The targeted and parent molecule files were produced with different include-introns settings, respectively:"
            " {} and {}. Please re-run '{product} count'/'{product} count --target-panel' with uniform include-introns "
            "parameters in order to run targeted-compare.".format(
                targeted_mol_info_intron_mode, parent_mol_info_intron_mode, product=cr_env.product()
            )
        )
        raise cr_preflight.PreflightException(msg)
    # Load feature IDs from target panel for validation (also validates target panel file)
    try:
        _, _, _ = simple_utils.parse_target_csv(
            args.target_set, expected_targeting_method=TARGETING_METHOD_HC
        )
    except Exception as err:
        raise cr_preflight.PreflightException(str(err))

    mol_info_target_panel_hash = targeted_mol_info.get_metric("target_panel_hash")

    targeted_mol_info.close()
    parent_mol_info.close()

    # check target panel csv and targeted molecule_info panel hashes are the same
    target_panel_csv_hash = cr_reference.compute_hash_of_file(args.target_set)
    if mol_info_target_panel_hash != target_panel_csv_hash:
        msg = (
            "Target panel csv hash %s is not the same as the molecule info target panel hash %s. "
            % (
                target_panel_csv_hash,
                mol_info_target_panel_hash,
            )
        )
        msg += "Please make sure to use the same target-panel that was when running {} count --target-panel.".format(
            cr_env.product()
        )
        raise cr_preflight.PreflightException(msg)


# pylint: disable=unused-argument
def main(args, outs):
    try:
        run_preflight_checks(args)
    except cr_preflight.PreflightException as err:
        martian.exit(err.msg)

    cr_preflight.record_package_versions()
