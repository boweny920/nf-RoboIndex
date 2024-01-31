#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#

from __future__ import annotations

import errno
import gzip
import io
import os
import shutil
import sys
from typing import IO, TYPE_CHECKING, Any, Literal, TextIO, Union, overload

import lz4.frame as lz4
import martian

import cellranger.h5_constants as h5_constants

if TYPE_CHECKING:
    import subprocess

    Path = Union[bytes, str]


def fixpath(path: Path) -> Path:
    return os.path.abspath(os.path.expandvars(os.path.expanduser(path)))


def get_input_path(oldpath: Path, is_dir: bool = False) -> Path:
    if not isinstance(oldpath, str):
        sys.exit(f"'{oldpath}' is not a valid string type and is an invalid path")
    path = fixpath(oldpath)
    if not os.path.exists(path):
        sys.exit("Input file does not exist: %s" % path)
    if not os.access(path, os.R_OK):
        sys.exit(f"Input file path {path} does not have read permissions")
    if is_dir:
        if not os.path.isdir(path):
            sys.exit("Please provide a directory, not a file: %s" % path)
    else:
        if not os.path.isfile(path):
            sys.exit("Please provide a file, not a directory: %s" % path)
    return path


def get_input_paths(paths: list[Path]) -> list[Path]:
    return [get_input_path(path) for path in paths]


def get_output_path(oldpath: Path) -> Path:
    path = fixpath(oldpath)
    dirname = os.path.dirname(path)
    if not os.path.exists(dirname):
        sys.exit("Output directory does not exist: %s" % dirname)
    if not os.path.isdir(dirname):
        sys.exit("Please provide a directory, not a file: %s" % dirname)
    return path


@overload
def open_maybe_gzip(filename: Path, mode: Literal["r"] = "r") -> TextIO:
    ...


@overload
def open_maybe_gzip(filename: Path, mode: Literal["w"] = ...) -> TextIO:
    ...


@overload
def open_maybe_gzip(filename: Path, mode: Literal["rb"] = ...) -> io.BufferedReader:
    ...


@overload
def open_maybe_gzip(filename: Path, mode: Literal["wb"] = ...) -> io.BufferedReader:
    ...


def open_maybe_gzip(filename: Path, mode: str = "r") -> IO[Any]:
    # this _must_ be a bytes
    if not isinstance(filename, bytes):
        filename = str(filename).encode()
    if filename.endswith(h5_constants.GZIP_SUFFIX):
        raw = gzip.open(filename, mode, 2)
    elif filename.endswith(h5_constants.LZ4_SUFFIX):
        raw = lz4.open(filename, mode)
    else:
        return open(filename, mode)

    bufsize = 1024 * 1024  # 1MB of buffering
    if mode == "r":
        return io.TextIOWrapper(io.BufferedReader(raw, buffer_size=bufsize))
    elif mode == "w":
        return io.TextIOWrapper(io.BufferedWriter(raw, buffer_size=bufsize))
    elif mode == "rb":
        return io.BufferedReader(raw, buffer_size=bufsize)
    elif mode == "wb":
        return io.BufferedWriter(raw, buffer_size=bufsize)

    else:
        raise ValueError("Unsupported mode for compression: %s" % mode)


class CRCalledProcessError(Exception):
    def __init__(self, msg: str):
        super().__init__(msg)
        self.msg = msg

    def __str__(self) -> str:
        return self.msg


def check_completed_process(p: subprocess.CompletedProcess, cmd: str):
    """Raises an exception if the completed process failed.

    Args:
        p:   Subprocess
        cmd: Command that was run

    Raises:
        CRCalledProcessError
    """
    if p.returncode is None:
        raise CRCalledProcessError("Process did not finish: %s ." % cmd)
    elif p.returncode != 0:
        raise CRCalledProcessError("Process returned error code %d: %s ." % (p.returncode, cmd))


def mkdir(path: Path, exist_ok: bool = True):
    """Create a directory.

    By default succeed if it already exists.

    Useful because transient NFS server issues may induce double creation attempts.
    """
    if exist_ok:
        os.makedirs(path, exist_ok=True)
    else:
        os.mkdir(path)


def remove(f: Path, nonexistent_ok: bool = True):
    """Delete a file. By default succeed if it doesn't exist.

    Useful because transient NFS server issues may induce double deletion attempts.
    """
    if nonexistent_ok:
        try:
            os.remove(f)
        except OSError as e:
            if e.errno == errno.ENOENT:
                pass
            else:
                raise
    else:
        os.remove(f)


def hardlink_with_fallback(src: Path, dst: Path):
    """Hard-links src to dst, falling back to copy if it fails.

    If `src` is a directory, it will attempt to recursively hardlink.
    """

    def _hardlink_file_with_fallback(src: Path, dst: Path):
        """Hardlink a file, fallback to copy if fail."""
        try:
            os.link(src, dst)
        except OSError as ex:
            if ex.errno in [errno.EPERM, errno.EOPNOTSUPP, errno.EXDEV, errno.EMLINK, errno.EEXIST]:
                try:
                    shutil.copy2(src, dst)
                except shutil.SameFileError:
                    pass
            else:
                raise

    if os.path.isdir(src):
        # recursively hardlink a path, fallback to copy if fail
        shutil.copytree(src, dst, copy_function=_hardlink_file_with_fallback, dirs_exist_ok=True)
    else:
        _hardlink_file_with_fallback(src, dst)


def hard_link(f, relative_path=None):
    """Make a new hard link in a stage directory to the file f, defaulting to the basename of f."""
    if not f:
        return None

    if relative_path is None:
        relative_path = os.path.basename(f)

    new_path = martian.make_path(relative_path).decode("utf8")
    hardlink_with_fallback(f, new_path)

    return new_path


def concatenate_files(out_path: Path, in_paths: list[Path], mode: str = ""):
    with open(out_path, "w" + mode) as out_file:
        for in_path in in_paths:
            with open(in_path, "r" + mode) as in_file:
                shutil.copyfileobj(in_file, out_file)


def concatenate_headered_files(out_path: Path, in_paths: list[Path], mode: str = ""):
    """Concatenate files, taking the first line of the first file.

    and skipping the first line for subsequent files.
    Asserts that all header lines are equal.
    """
    with open(out_path, "w" + mode) as out_file:
        if len(in_paths) > 0:
            # Write first file
            with open(in_paths[0], "r" + mode) as in_file:
                header = in_file.readline()
                out_file.write(header)
                shutil.copyfileobj(in_file, out_file)

        # Write remaining files
        for in_path in in_paths[1:]:
            with open(in_path, "r" + mode) as in_file:
                this_header = in_file.readline()
                assert this_header == header
                shutil.copyfileobj(in_file, out_file)


def write_empty_json(filename: Path):
    with open(filename, "wb") as f:
        f.write(b"{}")


def touch(path: Path):
    """Create an empty file."""
    fh = open(path, "w")
    fh.close()
