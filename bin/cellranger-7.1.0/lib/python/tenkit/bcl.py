#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#

from __future__ import annotations

import collections
import glob
import os
import re
import subprocess
import xml.etree.ElementTree as etree
from collections.abc import Iterable
from typing import Optional, TypedDict, Union

import martian

BCL2FASTQ_V1 = "OLD"
BCL2FASTQ_V2 = "NEW"


def get_bcl2fastq_v1(hostname: str) -> tuple[Optional[str], Optional[str]]:
    try:
        subprocess.check_call(["which", "configureBclToFastq.pl"])

        try:
            subprocess.check_call(["which", "perl"])
        except subprocess.CalledProcessError:
            msg = (
                "On machine: %s, perl not found on PATH. (Required for configureBclToFastq.pl)"
                % hostname
            )
            return (None, msg)

        return ("1.8.4", None)
    except subprocess.CalledProcessError:
        msg = "On machine: %s, configureBclToFastq.pl not found on PATH." % hostname
        return (None, msg)


def get_bcl2fastq_v2(hostname: str) -> tuple[Optional[bytes], Optional[str]]:
    try:
        subprocess.check_call(["which", "bcl2fastq"])
        # Restore the LD_LIBRARY_PATH set aside by sourceme.bash/shell10x.
        # Required for some installations of bcl2fastq.
        new_environ = dict(os.environ)
        new_environ["LD_LIBRARY_PATH"] = os.environ.get("_TENX_LD_LIBRARY_PATH", "")
        output = subprocess.check_output(
            ["bcl2fastq", "--version"], env=new_environ, stderr=subprocess.STDOUT
        )
        match = None
        for l in output.split(b"\n"):
            match = re.match(b"bcl2fastq v([0-9.]+)", l)
            if match is not None:
                return (match.groups()[0], None)

        return (
            None,
            "bcl2fastq version not recognized -- please check the output of bcl2fastq --version",
        )
    except subprocess.CalledProcessError:
        msg = "On machine: %s, bcl2fastq not found on PATH." % hostname
        return (None, msg)


def check_bcl2fastq(hostname: str, rta_str: str) -> tuple[str, Union[bytes, str]]:
    x, y, z = [int(xi) for xi in rta_str.split(".")][0:3]

    # RTA <1.18.54
    # must run 1.8.4
    if x == 1 and ((y < 18) or ((y == 18) and z < 54)):
        v1, _ = get_bcl2fastq_v1(hostname)
        if v1 is not None:
            return (BCL2FASTQ_V1, v1)
        else:
            msg = "mkfastq requires bcl2fastq 1.8.4 for RTA version: %s" % rta_str
            martian.exit(msg)
            raise SystemExit()

    # RTA >= 1.18.54
    # run 2.17 or higher
    else:
        v2, msg = get_bcl2fastq_v2(hostname)

        if v2 is not None:
            v2x, v2y = [int(v2_part) for v2_part in v2.split(b".")][0:2]
        else:
            msg = (
                "No valid bcl2fastq found on path. Recommended version of bcl2fastq is v2.20.\n\n%s"
                % msg
            )
            martian.exit(msg)
            raise SystemExit()

        if v2x == 2 and v2y >= 17:
            return (BCL2FASTQ_V2, v2)
        else:
            msg = (
                "Incorrect bcl2fastq version found: %s. Recommended version of bcl2fastq is v2.20."
                % v2
            )
            martian.exit(msg)
            raise SystemExit()


def get_rta_version(input_path: str) -> tuple[str, Optional[bool], dict[str, str]]:
    """Query the BCL folder for the RTA version of the run.

    Also finds out whether the I2 read needs to be reverse
    complemented.
    """

    rp_nextseq = os.path.join(input_path, "RunParameters.xml")
    rp_other = os.path.join(input_path, "runParameters.xml")

    if os.path.exists(rp_nextseq):
        run_parameters_xml = rp_nextseq
    else:
        run_parameters_xml = rp_other

    tree = etree.parse(run_parameters_xml)

    # Do we need to RC the I2 read?
    # Our current understanding is that NextSeq and HiSeq X / 4000 require it
    rc_i2 = False

    # TENKIT-60 NovaSeq runParameters.xml doesn't have the "Setup" node
    setup_node = tree.getroot().find("Setup")
    if setup_node is not None:
        application = tree.getroot().find("Setup").find("ApplicationName").text
        application_version = tree.getroot().find("Setup").find("ApplicationVersion").text
    else:
        # CSR-477 iSeq has a new runParameters variant!
        application_name = tree.getroot().find("ApplicationName")
        if application_name is not None:
            application = application_name.text
        else:
            application = tree.getroot().find("Application").text
        application_version = tree.getroot().find("ApplicationVersion").text

    assert application is not None
    assert application_version
    if application.find("NextSeq") >= 0:
        rc_i2 = True
    elif application.find("MiSeq") >= 0:
        rc_i2 = False
    # according to https://support.illumina.com/content/dam/illumina-support/documents/documentation/system_documentation/miseq/indexed-sequencing-overview-guide-15057455-03.pdf
    # NovaSeq follows MiSeq/HiSeq 2500 workflow for I5 index; this is effectively
    # a noop thanks to rc_i2 being False by default but let's make it explicit
    elif application.find("NovaSeq") >= 0:
        rc_i2 = detect_rc_i2_via_recipe_novaseq(input_path, application_version)
    elif application.find("HiSeq") >= 0:
        # Hiseq 4000 has version 3 -- hopefully this is stable??
        app_str = re.search(r"[\d\.]+", application_version).group()
        main_app_ver = int(app_str.split(".")[0])
        if main_app_ver > 2:
            rc_i2 = True
        else:
            rc_i2 = False
    # according to https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/experiment-design/illumina-adapter-sequences-1000000002694-07.pdf
    # iSeq 100 = rc.  index has to be exactly 0, to avoid MiSeq conflict.
    elif application.find("iSeq") == 0:
        rc_i2 = True

    # Can we run bcl2fastq 2, or do we need to use the old version?
    if setup_node is not None:
        rta_tag = tree.getroot().find("Setup").find("RTAVersion")
        if rta_tag is None:
            rta_tag = tree.getroot().find("RTAVersion")
    # TENKIT-60 NovaSeq lowercases the tag
    else:
        rta_tag = tree.getroot().find("RtaVersion")

    rta_version = rta_tag.text
    assert rta_version is not None
    if rta_version.startswith("v"):
        rta_version = rta_version[1:]

    params = {"ApplicationName": application, "ApplicationVersion": application_version}
    return (rta_version, rc_i2, params)


def get_sequencer_type(input_path: str) -> Optional[str]:
    """Returns the sequencer type from runParameters.xml in the input path.

    TODO: Perhaps look at the flowcell ID instead (see Preyas'
    Illumina identification code to see if that's more robust)
    """
    _, _, params = get_rta_version(input_path)
    if "ApplicationName" in params:
        return params["ApplicationName"].split()[0]
    else:
        return None


def detect_rc_i2_via_recipe_novaseq(flowcell_path: str, application_version) -> Optional[bool]:
    """Determine if I2 is RC by inspecting the Recipe xml file.

    Finds the recipe xml file given the flowcell_path. Return true if I2 is RC,
    false if RC is fwd-complement, and None if the recipe xml couldn't be
    found / read.
    """

    app_str = re.search(r"[\d\.]+", application_version).group()
    app_ver = [int(x) for x in app_str.split(".")]

    # We use a different detection scheme in ICS 1.8.0 and above.
    new_novaseq_detect_mode = False
    if app_ver >= [1, 8, 0]:
        new_novaseq_detect_mode = True

    run_info_xml = os.path.join(flowcell_path, "RunInfo.xml")
    if not os.path.exists(run_info_xml):
        return None
    (_, flowcell) = load_run_info(run_info_xml)

    recipe_file = None

    novaseq_recipe = os.path.join(flowcell_path, "Recipe", flowcell + ".xml")
    if os.path.exists(novaseq_recipe):
        recipe_file = novaseq_recipe

    if recipe_file is None:
        # Other sequencers (eg NextSeq) have a different convention
        # try and find any recipe file.
        xmls = glob.glob(os.path.join(flowcell_path, "Recipe", "*.xml"))
        if len(xmls) > 0:
            # We have no idea what >1 xml file here means
            # so we just pick one
            recipe_file = xmls[0]

    if recipe_file is not None:
        if new_novaseq_detect_mode:
            return new_novaseq_detect_rc_i2_from_recipe_xml(recipe_file)
        else:
            return old_novaseq_detect_rc_i2_from_recipe_xml(recipe_file)
    else:
        return None


def old_novaseq_detect_rc_i2_from_recipe_xml(recipe_xml: str) -> bool:
    """Determine if I2 is RC by inspecting the Recipe xml file.

    Based on a scheme from Illumina, if the "IndexPreparation-i5" or "Index2Preparation" ChemistryStep
    exists and uses the "BP14" reagent, the the I2 read is RC.
    """

    tree = etree.parse(recipe_xml)
    chem_steps = [x for x in tree.iter() if x.tag == "ChemistryStep"]
    i5_steps = [x for x in chem_steps if x.attrib.get("Description").startswith("Index")]

    reagents = set()
    for step in i5_steps:
        for el in step.iter():
            r = el.attrib.get("ReagentName")
            if r:
                reagents.add(r)

    return "BP14" in reagents


def new_novaseq_detect_rc_i2_from_recipe_xml(recipe_xml):
    """New scheme for detecting workflow mode in NovaSeq SW version 1.8.0 and above.

    Parse if  `<ChemistryRef Description="PETurnaround" ChemistryName="PETurnaround" />`
    comes before or after `<ReadRef Description="IndexRead i5" ReadName="IndexRead2" />`.
    If PETuraround comes after, then you're in workflow A, if it comes before IndexRead2workflow B.
    """
    tree = etree.parse(recipe_xml)
    protocol = [x for x in tree.iter() if x.tag == "Protocol"][0]

    # find the 2 key steps in the protocol list
    index_read2 = [
        x
        for x in enumerate(protocol.iter())
        if x[1].tag == "ReadRef" and x[1].attrib.get("ReadName") == "IndexRead2"
    ]
    turnaround = [
        x
        for x in enumerate(protocol.iter())
        if x[1].tag == "ChemistryRef" and x[1].attrib.get("ChemistryName") == "PETurnaround"
    ]

    # Check which step comes first
    if len(index_read2) == 1 and len(turnaround) == 1:
        if turnaround > index_read2:
            # workflow A
            return False
        else:
            # workflow B
            return True

    else:
        # We don't recognize the steps in the protocol, so bail
        return None


class ReadInfo(TypedDict):
    read_length: int
    read_name: str
    index_read: bool
    original_read_name: str


def load_run_info(run_info_xml: str) -> tuple[list[ReadInfo], Optional[str]]:
    """Get the read names and read lengths from the Illumina RunInfo.xml file."""
    tree = etree.parse(run_info_xml)
    reads_node = tree.getroot().find("Run").find("Reads")
    reads = reads_node.findall("Read")

    # Now we follow some hard-coded conventions on order in which reads appear.
    # R1 is always first
    # R2 is always last
    # Index reads (if they exist) are in the middle, in numerical order

    # NOTE -- if there is only one index read it is assumed to be I1.
    # BclToFastq should give us this
    # NOTE -- we assume paired end reads!
    read_info = [
        ReadInfo(
            read_length=x,
            index_read=idx != 0 and idx != len(reads) - 1,
            read_name="R1" if idx == 0 else "R2" if idx == len(reads) - 1 else "I" + str(idx),
            original_read_name="R" + str(idx + 1),
        )
        for idx, x in enumerate(int(read.attrib["NumCycles"]) for read in reads)
    ]
    flowcell = tree.getroot().find("Run").find("Flowcell").text

    # NB: currently you have to comment out the next two lines to get
    # nosetests to run correctly outside of a stage.
    martian.log_info("Read Info: %s" % read_info)
    martian.log_info("Flowcell ID: %s" % flowcell)
    return (read_info, flowcell)


def make_bases_mask_val(
    read_info: Iterable[ReadInfo],
    barcode_read: Optional[str] = None,
    sample_index_read: Optional[str] = None,
    dual_indexed: bool = False,
    ignore_dual_index: bool = False,
) -> str:
    """Undocumented.

    Args:
        read_info: The ReadInfo block from RunInfo.xml
        barcode_read: The read to use as the barcode. Can be an index.
        sample_index_read: The ReadInfo read (I1, I2) to use as the sample index
        dual_indexed: If the input BCLs were dual-indexed intentionally, then
                      preserve the index status of the 2nd, non-sample index index
                      in the mask.  Too often, people inadvertently ran dual-indexed
                      values for the barcode read (Cell Ranger v1, GemCode), so the
                      default behavior is to treat the off-SI indexed read as a
                      normal, non-index read.
        ignore_dual_index: Stub out any dual index with Ns.
    """
    # We will emit all reads in RunInfo.xml
    def base_mask(read):
        if read["read_name"][0] == "R":
            return "Y" + str(read["read_length"])
        elif read["read_name"][0] == "I":
            if read["read_name"] == barcode_read:
                return "Y" + str(read["read_length"])
            elif read["read_name"] == sample_index_read:
                return "I" + str(read["read_length"])
            elif dual_indexed:
                if ignore_dual_index:
                    return "N" + str(read["read_length"])
                else:
                    return "I" + str(read["read_length"])
            else:
                return "Y" + str(read["read_length"])
        else:
            martian.throw("read name was not recognized: %s" % read["read_name"])
            raise SystemExit()

    # Special hack to convert the bases_mask
    # to only give 8bp in I1, when in dual indexing
    # mode but using `ignore_dual_index` (aka --filter-single-index)
    # This has the following effect on the a bases mask:
    # Y28,I10,I10,Y90 -> Y28,I8,N12,Y90.
    # See details in CELLRANGER-3909
    if dual_indexed and ignore_dual_index:
        print("original read_info: ", read_info)
        rr = collections.OrderedDict((r["read_name"], r) for r in read_info)
        old_i1_length = rr["I1"]["read_length"]
        rr["I1"]["read_length"] = 8
        rr["I2"]["read_length"] = rr["I2"]["read_length"] + (old_i1_length - 8)
        read_info = list(rr.values())

        print("edited read info for i1 only: ", read_info)

    masks = [base_mask(r) for r in read_info]
    return ",".join(masks)


def get_bcl2fastq_read_type_map(
    read_info: Iterable[ReadInfo],
    barcode_read: Optional[str] = None,
    sample_index_read: Optional[str] = None,
    dual_indexed: bool = False,
    ignore_dual_index: bool = False,
) -> dict[str, str]:
    """Get a mapping between read name and output file.

    Read name from ReadInfo (R1,I1,I2,R2) and bcl2fastq
    output file naming (R1/R2/R3/R4/I1)

    The guarantee here is that the 10X sample index will always be on I1,
    if generated.  If dual-indexing is specified, the secondary index will
    be on I2.  Upstream pipestances can expect read1 to be in R1, sample
    indexes to be on I1, and read2 to be on R2.

    Args:
        read_info: The ReadInfo block from RunInfo.xml
        barcode_read: The ReadInfo read to use as the barcode (can be I2)
        sample_index_read: The ReadInfo read (I1, I2) to use as the sample index
    """
    # read_names = [r["read_name"] for r in read_info]
    read_map = {}
    reads_counted = 0
    for r in read_info:
        read_name = r["read_name"]
        assert isinstance(read_name, str)
        if read_name == sample_index_read:
            read_map[read_name] = "I1"
        elif dual_indexed and r["index_read"]:
            # handle I2 == barcode_read case (ATAC)
            # ignore_dual_index will be True since ATAC sample
            # index is I7-only, so this goes first
            if read_name == barcode_read:
                reads_counted += 1
                read_map[read_name] = "R%d" % reads_counted
            # ignore I2 if ignore_dual_index specified
            elif ignore_dual_index:
                continue
            else:
                read_map[read_name] = "I2"
        else:
            reads_counted += 1
            read_map[read_name] = "R%d" % reads_counted
    return read_map


if __name__ == "__main__":
    v = get_bcl2fastq_v2("host")
    print(v)
    check_bcl2fastq("host", "2.3.4")
