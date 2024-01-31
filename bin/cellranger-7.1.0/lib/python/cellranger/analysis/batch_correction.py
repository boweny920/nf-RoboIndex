#!/usr/bin/env python
#
# Copyright (c) 2018 10x Genomics, Inc. All rights reserved.
#
"""Library functions for performing batch correction."""


from __future__ import annotations

import struct
from collections import Counter
from typing import Union

import numpy as np
import sklearn.neighbors as sk_neighbors
from sklearn.metrics.pairwise import rbf_kernel

DEFAULT_BALLTREE_LEAFSIZE = 40


def batch_effect_score(
    dimred_matrix: np.ndarray[tuple[int, int], np.dtype[np.float64]],
    batch_ids: np.ndarray[int, np.dtype[np.string_]],
    knn_method: Union[str, int] = "sqrt",
    subsample_frac: float = 1.0,
):
    """For each cell, search KNN and calculate the proportion of cells from.

    the same batch. Then compute the ratio of the proportion to the
    percentage (number of cells) of this batch. The batch_effect_score is
    defined as the average of the ratio. Closer to 1 means no batch effect.
    """
    num_bcs = dimred_matrix.shape[0]
    assert num_bcs == len(batch_ids)

    # Do subsampling before computing batch percentages to avoid introducing noise
    np.random.seed(0)
    if subsample_frac < 1.0:
        select_bc_idx = np.array(
            [i for i in range(num_bcs) if np.random.uniform() < subsample_frac]
        )
        num_bcs_final = len(select_bc_idx)
        dimred_matrix = dimred_matrix[select_bc_idx]
        batch_ids = batch_ids[select_bc_idx]
    else:
        num_bcs_final = num_bcs

    # Return batch score of 1 if there is only 1 cell
    if num_bcs_final < 2:
        return 1

    if knn_method == "sqrt":
        num_nn = int(np.sqrt(num_bcs_final))
    elif isinstance(knn_method, int):
        num_nn = min(num_bcs_final - 1, knn_method)
    else:
        raise ValueError('knn_method must be "sqrt" or an integer')

    # For a cell in a given batch, what fraction of other cells share the same batch?
    # This is the null expectation of same_batch_frac below, given perfect batch mixing
    counter = Counter(batch_ids)
    batch_to_frac = {batch: (count - 1) / (num_bcs_final - 1) for batch, count in counter.items()}
    null_same_batch_frac = np.fromiter((batch_to_frac[i] for i in batch_ids), dtype=np.float64)

    balltree = sk_neighbors.BallTree(dimred_matrix, leaf_size=DEFAULT_BALLTREE_LEAFSIZE)
    knn_idx = balltree.query(dimred_matrix, k=num_nn + 1, return_distance=False)

    same_batch_frac = np.mean(batch_ids[:, None] == batch_ids[knn_idx[:, 1:]], axis=1)
    same_batch_ratio = same_batch_frac / null_same_batch_frac

    return np.mean(same_batch_ratio)


def find_knn(curr_matrix, ref_matrix, knn):
    """For each row in curr_matrix, find k nearest neighbors in ref_matrix,.

    return an array of shape=[curr_matrix.shape[0] * knn, ], which stores
    the index of nearest neighbors in ref_matrix
    """
    balltree = sk_neighbors.BallTree(ref_matrix, leaf_size=DEFAULT_BALLTREE_LEAFSIZE)
    num_nn = min(ref_matrix.shape[0], knn)
    nn_idx = balltree.query(curr_matrix, k=num_nn, return_distance=False)
    return nn_idx.ravel().astype(int)


def serialize_batch_nearest_neighbor(fp, batch_nearest_neighbor):
    for (a, b), s in batch_nearest_neighbor.items():
        fp.write(struct.pack("qqQ", a, b, len(s)))
        for i, j in s:
            fp.write(struct.pack("qq", i, j))


def deserialize_batch_nearest_neighbor(fp):
    """>>> from cStringIO import StringIO.

    >>> batch1 = dict()
    >>> batch1[(0, 1)] = set([(1, 2), (3, 4), (5, 6)])
    >>> batch1[(1, 2)] = set([(7, 8), (9, 10)])
    >>> batch1[(3, 4)] = set([(11, 12)])
    >>> fp = StringIO()
    >>> serialize_batch_nearest_neighbor(fp, batch1)
    >>> fp.seek(0)
    >>> batch2 = deserialize_batch_nearest_neighbor(fp)
    >>> batch1 == batch2
    True
    """
    batch_nearest_neighbor = {}
    while True:
        fmt = "qqQ"
        sz = struct.calcsize("qqQ")
        buf = fp.read(sz)
        if len(buf) == 0:
            break
        elif len(buf) != sz:
            raise RuntimeError("corrupted batch_nearest_neighbor stream (key)")
        a, b, slen = struct.unpack(fmt, buf)
        fmt = "qq"
        sz = struct.calcsize("qq")
        s = set()
        for _ in range(slen):
            buf = fp.read(sz)
            if len(buf) != sz:
                raise RuntimeError("corrupted batch_nearest_neighbor stream (set)")
            i, j = struct.unpack(fmt, buf)
            s.add((i, j))
        batch_nearest_neighbor[(a, b)] = s
    return batch_nearest_neighbor


def correction_vector(dimred_matrix, cur_submatrix_idx, mnn_cur_idx, mnn_ref_idx, sigma):
    """Compute the batch-correction vector.

    1. For each MNN pair in current dataset and the reference, a pair-specific
    batch-correction vector is computed as the vector difference between the
    paired cells.
    2. For each barcode in cur dataset, a batch-correction vector is calculated
    as a weighted average of these pair-specific vectors, as computed with a
    Gaussian kernel.
    """
    num_pcs = dimred_matrix.shape[1]
    corr_vector = np.zeros((0, num_pcs))

    # the number of mnn and submatrix dim might be very large, process by chunk to save memory
    cur_submatrix_size = len(cur_submatrix_idx)
    mnn_size = len(mnn_cur_idx)
    # based on empirical testing
    cur_submatrix_chunk_size = int(1e6 / num_pcs)
    mnn_chunk_size = int(2e7 / num_pcs)

    for i in range(0, cur_submatrix_size, cur_submatrix_chunk_size):
        cur_submatrix_chunk = cur_submatrix_idx[i : i + cur_submatrix_chunk_size]
        cur_submatrix = dimred_matrix[cur_submatrix_chunk]

        weighted_sum, weights_sum = np.zeros(cur_submatrix.shape), np.zeros(cur_submatrix.shape)

        for j in range(0, mnn_size, mnn_chunk_size):
            mnn_cur_chunk = mnn_cur_idx[j : j + mnn_chunk_size]
            mnn_ref_chunk = mnn_ref_idx[j : j + mnn_chunk_size]

            mnn_cur = dimred_matrix[mnn_cur_chunk]
            weights = rbf_kernel(cur_submatrix, mnn_cur, gamma=0.5 * sigma)
            bias = dimred_matrix[mnn_ref_chunk] - mnn_cur
            weighted_sum += np.dot(weights, bias)
            weights_sum += np.tile(np.sum(weights, axis=1), (num_pcs, 1)).T

        corr_vector = np.vstack((corr_vector, weighted_sum / weights_sum))

    return corr_vector
