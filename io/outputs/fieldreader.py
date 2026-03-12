"""
fieldreader.py
==============
Data containers and loaders for GPR post-processing.

Classes
-------
GPRDataset          – holds one receiver-line dataset (amp, phase, real, imag, r, label)
AnalyticalLoader    – reads a semi-analytical CSV (Evert / empymod format)
ElfeLoader          – reads one electric_fields_receiver_line.txt output file

Functions
---------
load_elfe_batch     – bulk-load a list of run directories into ElfeLoader objects
"""

from __future__ import annotations
import os
import numpy as np
from dataclasses import dataclass
from typing import Optional


# ---------------------------------------------------------------------------
# Dataset container
# ---------------------------------------------------------------------------

@dataclass
class GPRDataset:
    """
    Single-orientation GPR receiver-line dataset.

    Fields
    ------
    amp         : |Ex|  absolute amplitude
    phase       : ∠Ex   phase in radians
    real        : Re(Ex)
    imag        : Im(Ex)
    orientation : 'Endfire', 'Broadside', or 'Oblique'
    r           : receiver positions (metres)
    label       : source / run name used in legends
    """
    amp:         np.ndarray
    phase:       np.ndarray
    real:        np.ndarray
    imag:        np.ndarray
    orientation: str
    r:           np.ndarray
    label:       str

    def field(self, idx: int) -> np.ndarray:
        """Return field array by quantity index: 0=amp, 1=phase, 2=real, 3=imag."""
        return (self.amp, self.phase, self.real, self.imag)[idx]

    # Keep notebook tuple-index access working (ds[0]…ds[6])
    def __getitem__(self, idx: int):
        return (self.amp, self.phase, self.real, self.imag,
                self.orientation, self.r, self.label)[idx]


# Default quantities for 4-panel plots
QUANTITIES = [
    ("Amplitude (Ex)",   "Amplitude (V/m)"),
    ("Phase (Ex)",       "Phase (rad)"),
    ("Real (Ex)",        "Real (V/m)"),
    ("Imaginary (Ex)",   "Imaginary (V/m)"),
]


# ---------------------------------------------------------------------------
# Analytical CSV loader
# ---------------------------------------------------------------------------

class AnalyticalLoader:
    """
    Load a semi-analytical CSV (Evert quadrature or empymod DLF output).

    Expected column layout:
        0   : r  (source–receiver separation, m)
        1   : |Ex| Endfire
        3,4 : Re, Im  Endfire
        5   : |Ex| Broadside
        7,8 : Re, Im  Broadside   (sign-flipped by default)
        9   : |Ex| Oblique
        11,12: Re, Im  Oblique    (sign-flipped by default)

    Parameters
    ----------
    filepath         : path to the CSV file
    label            : legend label for all datasets from this file
    sign_flip_bs_ob  : apply −1 to Broadside/Oblique Re and Im (default True)
    """

    def __init__(self, filepath: str, label: str, sign_flip_bs_ob: bool = True):
        self.filepath        = filepath
        self.label           = label
        self.sign_flip_bs_ob = sign_flip_bs_ob
        self._data: Optional[np.ndarray] = None

    def _load(self) -> np.ndarray:
        if self._data is None:
            self._data = np.loadtxt(self.filepath, delimiter=",", skiprows=1)
        return self._data

    def load_all(self) -> list[GPRDataset]:
        """Return [Endfire, Broadside, Oblique]."""
        d   = self._load()
        r   = d[:, 0]
        sgn = -1.0 if self.sign_flip_bs_ob else 1.0

        ef_re, ef_im = d[:, 3], d[:, 4]
        bs_re, bs_im = sgn * d[:, 7], sgn * d[:, 8]
        ob_re, ob_im = sgn * d[:, 11], sgn * d[:, 12]

        def _ds(re, im, amp_col, orientation):
            amp   = d[:, amp_col]
            phase = np.arctan2(im, re)
            return GPRDataset(amp, phase, re, im, orientation, r, self.label)

        return [
            _ds(ef_re, ef_im, 1,  "Endfire"),
            _ds(bs_re, bs_im, 5,  "Broadside"),
            _ds(ob_re, ob_im, 9,  "Oblique"),
        ]

    def endfire(self)   -> GPRDataset: return self.load_all()[0]
    def broadside(self) -> GPRDataset: return self.load_all()[1]
    def oblique(self)   -> GPRDataset: return self.load_all()[2]


# ---------------------------------------------------------------------------
# elfe3D output loader
# ---------------------------------------------------------------------------

class ElfeLoader:
    """
    Load one electric_fields_receiver_line.txt and extract receiver slices.

    Column layout:
        0,1,2 : x, y, z positions
        4,5   : Re(Ex), Im(Ex)

    Parameters
    ----------
    filepath      : path to the text file
    label         : legend label
    num_endfire   : number of endfire receivers (starting at row 0)
    num_broadside : broadside receiver count after the endfire block (0 = absent)
    num_oblique   : oblique receiver count after the broadside block  (0 = absent)
    """

    def __init__(
        self,
        filepath:      str,
        label:         str,
        num_endfire:   int,
        num_broadside: int = 0,
        num_oblique:   int = 0,
    ):
        self.filepath      = filepath
        self.label         = label
        self.num_endfire   = num_endfire
        self.num_broadside = num_broadside
        self.num_oblique   = num_oblique
        self._data: Optional[np.ndarray] = None

    def _load(self) -> np.ndarray:
        if self._data is None:
            self._data = np.loadtxt(self.filepath)
        return self._data

    def _slice(self, start: int, count: int, orientation: str) -> GPRDataset:
        d     = self._load()
        rows  = d[start:start + count]
        rx    = rows[:, 0]
        re_ex = rows[:, 4]
        im_ex = rows[:, 5]
        amp   = np.abs(re_ex + 1j * im_ex)
        phase = np.angle(re_ex + 1j * im_ex)
        return GPRDataset(amp, phase, re_ex, im_ex, orientation, rx, self.label)

    def endfire(self)   -> GPRDataset:
        return self._slice(0, self.num_endfire, "Endfire")

    def broadside(self) -> GPRDataset:
        return self._slice(self.num_endfire, self.num_broadside, "Broadside")

    def oblique(self)   -> GPRDataset:
        return self._slice(self.num_endfire + self.num_broadside,
                           self.num_oblique, "Oblique")


# ---------------------------------------------------------------------------
# Batch helper
# ---------------------------------------------------------------------------

def load_elfe_batch(
    base_folder:   str,
    run_names:     list[str],
    labels:        list[str],
    num_endfire:   int,
    num_broadside: int = 0,
    num_oblique:   int = 0,
    filename:      str = "electric_fields_receiver_line.txt",
) -> list[ElfeLoader]:
    """
    Create one ElfeLoader per run directory.

    Parameters
    ----------
    base_folder   : parent directory containing the run sub-folders
    run_names     : list of sub-folder names
    labels        : corresponding legend labels
    num_endfire   : endfire receiver count (same for all runs)
    """
    return [
        ElfeLoader(
            os.path.join(base_folder, run, filename),
            label, num_endfire, num_broadside, num_oblique,
        )
        for run, label in zip(run_names, labels)
    ]
