#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#

"""Things that are terrible and should never be used.

.. deprecated::
   Use of the methods in this module suggests the author of the code
   misunderstood something about how cellranger works.
"""


from __future__ import annotations

# NOTE: MEM_GB_PER_THREAD is only used for a few stages where we've encountered memory oversubscription issues
# on clusters without memory reservations. As of 3/15/2017, it's used by:
# - RUN_PCA
# - RUN_DIFFERENTIAL_EXPRESSION
# - RUN_GRAPH_CLUSTERING
MEM_GB_PER_THREAD = 8


def get_thread_request_from_mem_gb(mem_gb):
    """A hack to deal with users who don't understand the --mempercore option.

    For systems without memory reservations, reserve multiple threads if necessary to
    avoid running out of memory.

    .. deprecated::
       Don't use this.  This is what the --mempercore option to mrp is for.
    """
    est_threads = round(float(mem_gb) / MEM_GB_PER_THREAD)
    # make sure it's 1, 2, or 4
    for threads in [1, 2, 4]:
        if est_threads <= threads:
            return threads
    return 4
