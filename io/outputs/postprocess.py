"""
postprocess.py
==============
Error computation for GPR receiver-line comparisons.

All error functions follow the convention used throughout the notebooks:
  - amplitude / real / imaginary  →  normalised:  |test − ref| / |ref|
  - phase                         →  absolute:    test − ref

If the test and reference datasets have different receiver grids, the
reference is interpolated onto the test positions before computing the error.

Functions
---------
field_error(ref, test, qty_idx)  – single-quantity error array
all_errors(ref, test)            – all four quantities + r in one call
error_stats(err)                 – (mean, std, max_abs) of a finite error array
"""

from __future__ import annotations
import numpy as np
from fieldreader import GPRDataset


def field_error(ref: GPRDataset, test: GPRDataset, qty_idx: int) -> np.ndarray:
    """
    Compute the error between *test* and *ref* for one quantity.

    Parameters
    ----------
    ref      : reference dataset (e.g. analytical solution)
    test     : dataset to evaluate (e.g. elfe3D run)
    qty_idx  : 0 = amplitude, 1 = phase, 2 = real, 3 = imaginary

    Returns
    -------
    err : np.ndarray, same length as test.r
    """
    ref_data  = ref.field(qty_idx)
    test_data = test.field(qty_idx)

    # Interpolate reference onto test grid if receiver positions differ
    if len(ref.r) != len(test.r) or not np.allclose(ref.r, test.r):
        ref_data = np.interp(test.r, ref.r, ref_data)

    if qty_idx in (0, 2, 3):   # normalised
        with np.errstate(divide="ignore", invalid="ignore"):
            err = np.abs(test_data - ref_data) / np.abs(ref_data)
            err[~np.isfinite(err)] = np.nan
    else:                       # absolute (phase)
        err = test_data - ref_data

    return err


def all_errors(
    ref: GPRDataset,
    test: GPRDataset,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute errors for all four quantities.

    Returns
    -------
    (amp_err, phase_err, real_err, imag_err, r)  where r = test.r
    """
    return (
        field_error(ref, test, 0),
        field_error(ref, test, 1),
        field_error(ref, test, 2),
        field_error(ref, test, 3),
        test.r,
    )


def error_stats(err: np.ndarray) -> tuple[float, float, float]:
    """
    Return (mean, std, max_abs) computed over the finite values of *err*.
    """
    finite = err[np.isfinite(err)]
    return (
        float(np.nanmean(finite)),
        float(np.nanstd(finite)),
        float(np.nanmax(np.abs(finite))),
    )
