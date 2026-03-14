"""
workflow.py
===========
Single class that orchestrates the full GPR simulation pipeline:

    Step 1  –  write_inputs()   create all FEM input files
    Step 2  –  run_tetgen()     call TetGen to mesh the geometry
    Step 3  –  run_solver()     call elfe3d_gpr
    Step 4  –  load_results()   read electric_fields_receiver_line.txt
    Step 5  –  plot_*()         visualise and compare with reference

Intended use: import into a Jupyter notebook (see workflow.ipynb).

Path configuration
------------------
  base_dir      : root of your project, e.g. r"C:\Projects\GPR"
  run_name      : name of the sub-folder that holds input + mesh files,
                    e.g. "in_air_only"  →  base_dir/in_air_only/
  exec_rel_path : path to elfe3d_gpr *relative to base_dir*,
                    e.g. "../elfe3d_gpr"  or  "../../solver/elfe3d_gpr"

WSL note
--------
  If you run this notebook from Windows (Jupyter / VS Code) set
  use_wsl=True.  Every subprocess call is then prefixed with "wsl"
  so that Linux binaries (tetgen15, elfe3d_gpr) are reached through
  the Windows Subsystem for Linux.
  If you run the notebook from inside WSL directly, set use_wsl=False.
"""

from __future__ import annotations

import os
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath
from typing import Optional

import numpy as np

# ── output-side imports (same folder) ──────────────────────────────────────
from outputs.fieldreader import AnalyticalLoader, ElfeLoader, GPRDataset
from outputs.postprocess import field_error, error_stats
from outputs.visualize   import (
    ReceiverLinePlot,
    ReceiverLineErrorPlot,
    ReceiverLineCombined,
    ReceiverLineCombinedMulti,
    ErrorHistogramPlot,
    ErrorStatPlot,
)

# ── input-side imports (adjust the import path if your io/inputs/ lives
#    elsewhere relative to this file) ───────────────────────────────────────
# workflow.py lives in master\ — io\inputs is one level down
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "io", "inputs"))
try:
    from inputs.writeinputfiles  import write_all_input_files   # noqa: F401
    from inputs.writetetgenpoly  import write_tetgen_poly        # noqa: F401
    _INPUTS_AVAILABLE = True
except ImportError:
    _INPUTS_AVAILABLE = False
    # The workflow still runs steps 2-5 if you skip write_inputs()


# ---------------------------------------------------------------------------
# Configuration dataclass
# ---------------------------------------------------------------------------

@dataclass
class WorkflowConfig:
    """
    All paths and simulation parameters in one place.
    Pass this to GPRWorkflow().
    """

    # ── Paths ───────────────────────────────────────────────────────────────
    base_dir:        str = r"C:\Projects\GPR\master"  # project root (Windows path)
    inputs_rel_path: str = r"io\inputs"               # inputs parent, relative to base_dir
                                                       # → base_dir\io\inputs\<run_name>\
    run_name:        str = "in_air_only"               # input sub-folder name
    exec_rel_path:   str = "elfe3d_gpr"               # executable, relative to base_dir
                                                       # → base_dir\elfe3d_gpr
    poly_filename:   str = "GPR_model.poly"            # .poly file written by write_inputs()
    output_file:   str = "electric_fields_receiver_line.txt"

    # ── TetGen ──────────────────────────────────────────────────────────────
    tetgen_bin:    str = "tetgen15"
    tetgen_flags:  str = "-pq1.2kAaen"

    # ── Environment ─────────────────────────────────────────────────────────
    use_wsl:       bool = True   # True  → prefix all subprocess calls with "wsl"
                                 # False → run natively (use when inside WSL)

    # ── Receiver counts (for ElfeLoader) ────────────────────────────────────
    num_endfire:   int  = 64
    num_broadside: int  = 0
    num_oblique:   int  = 0

    # ── FEM model parameters (forwarded to write_inputs) ────────────────────
    # Add whatever your writeinputfiles functions expect.
    # These are passed as **kwargs so you can extend freely.
    model_params:  dict = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Main workflow class
# ---------------------------------------------------------------------------

class GPRWorkflow:
    """
    Orchestrates the complete GPR simulation pipeline.

    Parameters
    ----------
    config : WorkflowConfig instance with all paths and parameters set.
    """

    def __init__(self, config: WorkflowConfig):
        self.cfg = config
        self._results: Optional[list[GPRDataset]] = None  # populated by load_results()

    # ── Derived paths ───────────────────────────────────────────────────────

    @property
    @property
    def input_dir(self) -> Path:
        """base_dir / inputs_rel_path / run_name  →  master\io\inputs\in_air_only"""
        return Path(self.cfg.base_dir) / self.cfg.inputs_rel_path / self.cfg.run_name

    @property
    def poly_path(self) -> Path:
        return self.input_dir / self.cfg.poly_filename

    @property
    def output_path(self) -> Path:
        return self.input_dir / self.cfg.output_file

    @property
    def exec_path(self) -> Path:
        """Absolute path to the elfe3d_gpr executable."""
        return (Path(self.cfg.base_dir) / self.cfg.exec_rel_path).resolve()

    def _wsl_path(self, windows_path: Path) -> str:
        """
        Convert a Windows path to a WSL path string.
        e.g. C:\\Projects\\GPR\\in_air_only → /mnt/c/Projects/GPR/in_air_only
        """
        p = str(windows_path).replace("\\", "/")
        if len(p) >= 2 and p[1] == ":":
            drive = p[0].lower()
            p = f"/mnt/{drive}/{p[3:]}"
        return p

    def _cmd(self, *args) -> list[str]:
        """Prepend 'wsl' if use_wsl=True."""
        cmd = list(args)
        if self.cfg.use_wsl:
            cmd = ["wsl"] + cmd
        return cmd

    # ── Step 1: write inputs ─────────────────────────────────────────────────

    def write_inputs(self, **kwargs):
        """
        Create all FEM input files and the TetGen .poly geometry file.

        Merges config.model_params with any extra kwargs you pass here.
        Calls write_tetgen_poly() and write_all_input_files() from io/inputs/.

        If the input modules are not importable, prints a warning and skips.
        """
        if not _INPUTS_AVAILABLE:
            print("⚠  io/inputs not found — skipping write_inputs(). "
                  "Check that io/inputs/ is on the path and imports succeed.")
            return

        params = {**self.cfg.model_params, **kwargs}
        os.makedirs(self.input_dir, exist_ok=True)

        print(f"Writing input files to:  {self.input_dir}")
        write_tetgen_poly(output_dir=str(self.input_dir), **params)
        write_all_input_files(output_dir=str(self.input_dir), **params)
        print("✓  Input files written.")

    # ── Step 2: run TetGen ───────────────────────────────────────────────────

    def run_tetgen(self):
        """
        Call TetGen to mesh the geometry.

        Runs:   tetgen15 <flags> <path/to/poly>
        The mesh files are written into the same directory as the .poly file.
        """
        poly  = self._wsl_path(self.poly_path) if self.cfg.use_wsl else str(self.poly_path)
        cmd   = self._cmd(self.cfg.tetgen_bin, self.cfg.tetgen_flags, poly)

        print(f"Running TetGen:\n  {' '.join(cmd)}")
        t0 = time.time()
        result = subprocess.run(cmd, capture_output=True, text=True)
        elapsed = time.time() - t0

        if result.returncode != 0:
            print("✗  TetGen failed:")
            print(result.stderr or result.stdout)
            raise RuntimeError(f"TetGen exited with code {result.returncode}")

        print(f"✓  TetGen finished in {elapsed:.1f}s.")
        if result.stdout:
            print(result.stdout[-2000:])   # last 2 kB of output

    # ── Step 3: run solver ───────────────────────────────────────────────────

    def run_solver(self):
        """
        Call the elfe3d_gpr FEM solver.

        The solver is invoked from the directory that contains it
        (exec_path.parent) with no additional arguments — it reads all
        input files by convention.

        Runs:   ./elfe3d_gpr
        """
        exec_file = self._wsl_path(self.exec_path)   if self.cfg.use_wsl else str(self.exec_path)
        exec_dir  = self._wsl_path(self.exec_path.parent) if self.cfg.use_wsl else str(self.exec_path.parent)

        # cd into the solver directory first so it finds its config files
        if self.cfg.use_wsl:
            cmd = ["wsl", "bash", "-c", f"cd '{exec_dir}' && ./elfe3d_gpr"]
        else:
            cmd = ["bash", "-c", f"cd '{exec_dir}' && ./elfe3d_gpr"]

        print(f"Running solver:\n  {' '.join(cmd)}")
        t0 = time.time()
        result = subprocess.run(cmd, capture_output=True, text=True)
        elapsed = time.time() - t0

        if result.returncode != 0:
            print("✗  Solver failed:")
            print(result.stderr or result.stdout)
            raise RuntimeError(f"elfe3d_gpr exited with code {result.returncode}")

        print(f"✓  Solver finished in {elapsed:.1f}s.")
        if result.stdout:
            print(result.stdout[-2000:])

    # ── Step 4: load results ──────────────────────────────────────────────────

    def load_results(self, label: str = "elfe3D") -> list[GPRDataset]:
        """
        Read electric_fields_receiver_line.txt and return a list of GPRDatasets.

        Returns a list with up to 3 elements depending on which orientations
        are configured: [endfire] or [endfire, broadside] or [endfire, broadside, oblique].

        The result is also stored in self._results for use by the plot methods.
        """
        loader = ElfeLoader(
            filepath=str(self.output_path),
            label=label,
            num_endfire=self.cfg.num_endfire,
            num_broadside=self.cfg.num_broadside,
            num_oblique=self.cfg.num_oblique,
        )

        results = [loader.endfire()]
        if self.cfg.num_broadside > 0:
            results.append(loader.broadside())
        if self.cfg.num_oblique > 0:
            results.append(loader.oblique())

        self._results = results
        print(f"✓  Loaded {len(results)} orientation(s) from {self.output_path}")
        return results

    # ── Step 5: plotting convenience methods ─────────────────────────────────

    def _get_datasets(self, datasets):
        """Use supplied datasets or fall back to self._results."""
        if datasets is not None:
            return datasets
        if self._results is None:
            raise RuntimeError("Call load_results() first, or pass datasets explicitly.")
        return self._results

    def plot_comparison(
        self,
        analytical: GPRDataset,
        datasets:   list[GPRDataset] = None,
        suptitle:   str = "",
        output_dir: str = None,
        fname:      str = None,
    ) -> None:
        """2×2 comparison plot (amplitude log, phase, real, imag)."""
        ds = self._get_datasets(datasets)
        ReceiverLinePlot([analytical] + ds).plot(
            suptitle=suptitle, output_dir=output_dir, fname=fname,
        )

    def plot_errors(
        self,
        analytical: GPRDataset,
        datasets:   list[GPRDataset] = None,
        suptitle:   str = "",
        output_dir: str = None,
        fname:      str = None,
    ) -> None:
        """2×2 error plot vs analytical reference."""
        ds = self._get_datasets(datasets)
        ReceiverLineErrorPlot(ds, reference=analytical).plot(
            suptitle=suptitle, output_dir=output_dir, fname=fname,
        )

    def plot_combined(
        self,
        analytical: GPRDataset,
        datasets:   list[GPRDataset] = None,
        suptitle:   str = "",
        output_dir: str = None,
        fname:      str = None,
    ) -> None:
        """2×4 combined figure (fields top, errors bottom). Multi-dataset version."""
        ds = self._get_datasets(datasets)
        if len(ds) == 1:
            ReceiverLineCombined(ds[0], analytical).plot(
                suptitle=suptitle, output_dir=output_dir, fname=fname,
            )
        else:
            ReceiverLineCombinedMulti(ds, analytical).plot(
                suptitle=suptitle, output_dir=output_dir, fname=fname,
            )

    def plot_histograms(
        self,
        analytical: GPRDataset,
        datasets:   list[GPRDataset] = None,
        suptitle:   str = "",
        output_dir: str = None,
        fname:      str = None,
    ) -> None:
        """1×2 amplitude + phase error histograms."""
        ds = self._get_datasets(datasets)
        ErrorHistogramPlot(ds, reference=analytical).plot(
            suptitle=suptitle, output_dir=output_dir, fname=fname,
        )

    def print_error_summary(
        self,
        analytical: GPRDataset,
        datasets:   list[GPRDataset] = None,
    ) -> None:
        """Print mean / std / max errors for each dataset and quantity."""
        ds = self._get_datasets(datasets)
        qty_names = ["Amplitude", "Phase", "Real", "Imaginary"]
        for d in ds:
            print(f"\n── {d.label} ({d.orientation}) ──")
            for qi, name in enumerate(qty_names):
                err = field_error(analytical, d, qi)
                m, s, mx = error_stats(err)
                unit = "%" if qi != 1 else "rad"
                scale = 100 if qi != 1 else 1
                print(f"  {name:12s}:  mean={m*scale:.3f}{unit}  "
                      f"std={s*scale:.3f}{unit}  max={mx*scale:.3f}{unit}")

    # ── Full pipeline shortcut ────────────────────────────────────────────────

    def run_all(self, write_inputs_kwargs: dict = None):
        """
        Run the full pipeline: write_inputs → run_tetgen → run_solver.
        Does NOT load results or plot — call load_results() and plot_* after.
        """
        self.write_inputs(**(write_inputs_kwargs or {}))
        self.run_tetgen()
        self.run_solver()
