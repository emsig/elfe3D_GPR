"""
visualize.py
============
Plotting classes for GPR receiver-line data.

All styling logic lives here — no separate style module needed.

Classes
-------
ReceiverLinePlot          – 2×2 multi-dataset comparison (amp log, phase, re, im)
ReceiverLineErrorPlot     – 2×2 errors vs a reference dataset
ReceiverLineCombined      – 2×4 single dataset: top row = fields, bottom = errors
ReceiverLineCombinedMulti – 2×4 multi-dataset: top row = fields, bottom = errors
ErrorHistogramPlot        – 1×2 amplitude + phase error histograms
ErrorStatPlot             – 2×2 mean / std / max vs a sweep parameter
"""

from __future__ import annotations

import os
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.ticker import FixedLocator, FixedFormatter

from fieldreader import GPRDataset, QUANTITIES
from postprocess import field_error, error_stats


# ---------------------------------------------------------------------------
# Styling
# ---------------------------------------------------------------------------

_PALETTE = [
    "#1f77b4", "#d62728", "#2ca02c", "#e66100",
    "#17becf", "#9467bd", "#bcbd22", "#1fa2ff", "#ff6f61", "#40e0d0",
]
_LINESTYLES = [
    "--", "-.", (0, (1, 1)), (0, (3, 1, 1, 1)), (0, (5, 1)), ":",
]
_FIXED = {
    "Evert": {
        "color": "#000000", "linestyle": "-", "lw_factor": 1.5,
        "label": "Analytical Solution",
    },
    "empymod - 2001 DLF": {
        "color": "#4B2E05", "linestyle": "--", "lw_factor": 1.5,
        "label": "empymod - 2001 DLF",
    },
}


def _build_styles(datasets: list[GPRDataset], base_lw: float = 2.5) -> dict[str, dict]:
    """
    Return a label → style-dict mapping for a list of datasets.
    Fixed sources (Evert, empymod) always get the same heavy lines.
    All other sources are assigned palette colours in the order they appear.
    """
    styles = {}
    var_idx = 0
    for ds in datasets:
        lbl = ds.label
        if lbl in styles:
            continue
        if lbl in _FIXED:
            spec = _FIXED[lbl]
            styles[lbl] = {
                "color":     spec["color"],
                "linestyle": spec["linestyle"],
                "linewidth": base_lw * spec["lw_factor"],
                "label":     spec["label"],
                "zorder":    1,
                "is_fixed":  True,
            }
        else:
            styles[lbl] = {
                "color":     _PALETTE[var_idx % len(_PALETTE)],
                "linestyle": _LINESTYLES[var_idx % len(_LINESTYLES)],
                "linewidth": base_lw,
                "label":     lbl,
                "zorder":    3 + var_idx,
                "is_fixed":  False,
            }
            var_idx += 1
    return styles


def _auto_xticks(ax, r: np.ndarray) -> None:
    start = np.floor(r.min() * 2) / 2
    stop  = np.ceil(r.max()  * 2) / 2 + 0.01
    ax.set_xticks(np.arange(start, stop, 1.0))


def _save(fig, output_dir: Optional[str], fname: Optional[str]) -> None:
    if output_dir and fname:
        fig.savefig(os.path.join(output_dir, fname), dpi=300)
    plt.show()


# ---------------------------------------------------------------------------
# 1. ReceiverLinePlot  –  2×2 comparison
# ---------------------------------------------------------------------------

class ReceiverLinePlot:
    """
    2×2 comparison of amplitude (log), phase, real, and imaginary parts
    for multiple datasets on the same axes.

    Parameters
    ----------
    datasets   : list of GPRDataset (mix of analytical and numerical runs)
    quantities : list of 4 (title, ylabel); defaults to QUANTITIES
    base_lw    : base line width for variable sources
    font       : dict of font sizes: suptitle, label, tick, legend
    figsize    : figure size
    """

    def __init__(
        self,
        datasets:   list[GPRDataset],
        quantities: list[tuple] = None,
        base_lw:    float = 2.5,
        font:       dict = None,
        figsize:    tuple = (12, 8),
    ):
        self.datasets   = datasets
        self.quantities = quantities or QUANTITIES
        self.styles     = _build_styles(datasets, base_lw)
        self.font       = {"suptitle": 24, "label": 18, "tick": 18, "legend": 15, **(font or {})}
        self.figsize    = figsize

    def plot(self, suptitle: str = "", output_dir: str = None, fname: str = None) -> plt.Figure:
        fig, axes = plt.subplots(2, 2, figsize=self.figsize, sharex=True)

        for j, (title, ylabel) in enumerate(self.quantities):
            ax = axes[j // 2, j % 2]
            fn = ax.semilogy if j == 0 else ax.plot

            # Fixed sources drawn first (underneath)
            for ds in self.datasets:
                st = self.styles[ds.label]
                if not st["is_fixed"]:
                    continue
                fn(ds.r, ds.field(j), label=st["label"], color=st["color"],
                   linestyle=st["linestyle"], linewidth=st["linewidth"], zorder=st["zorder"])

            # Variable sources
            for ds in self.datasets:
                st = self.styles[ds.label]
                if st["is_fixed"]:
                    continue
                fn(ds.r, ds.field(j), label=st["label"], color=st["color"],
                   linestyle=st["linestyle"], linewidth=st["linewidth"],
                   marker="o", markersize=2, zorder=st["zorder"])

            ax.set_title(title, fontsize=self.font["label"], fontweight="bold")
            ax.set_ylabel(ylabel, fontsize=self.font["label"])
            ax.tick_params(labelsize=self.font["tick"])
            ax.grid(True, which="both" if j == 0 else "major", linestyle="--", linewidth=0.5)
            if j // 2 == 1:
                ax.set_xlabel("Distance (m)", fontsize=self.font["label"])
                _auto_xticks(ax, self.datasets[0].r)
            if j == 3:
                ax.legend(fontsize=self.font["legend"])

        fig.suptitle(suptitle, fontsize=self.font["suptitle"], fontweight="bold")
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        _save(fig, output_dir, fname)
        return fig


# ---------------------------------------------------------------------------
# 2. ReceiverLineErrorPlot  –  2×2 errors
# ---------------------------------------------------------------------------

class ReceiverLineErrorPlot:
    """
    2×2 grid showing errors of each dataset relative to a reference.
    Fixed sources are skipped automatically.
    """

    def __init__(
        self,
        datasets:   list[GPRDataset],
        reference:  GPRDataset,
        quantities: list[tuple] = None,
        base_lw:    float = 2.5,
        font:       dict = None,
        figsize:    tuple = (12, 8),
    ):
        self.datasets   = datasets
        self.reference  = reference
        self.quantities = quantities or QUANTITIES
        self.styles     = _build_styles(datasets, base_lw)
        self.font       = {"suptitle": 24, "label": 18, "tick": 18, "legend": 15, **(font or {})}
        self.figsize    = figsize

    def plot(self, suptitle: str = "", output_dir: str = None, fname: str = None) -> plt.Figure:
        var_ds = [ds for ds in self.datasets if not self.styles[ds.label]["is_fixed"]]

        fig, axes = plt.subplots(2, 2, figsize=self.figsize, sharex=True)

        for j, (title, _) in enumerate(self.quantities):
            ax = axes[j // 2, j % 2]
            fn = ax.semilogy if j != 1 else ax.plot

            for ds in var_ds:
                err = field_error(self.reference, ds, j)
                st  = self.styles[ds.label]
                fn(ds.r, err, label=st["label"], color=st["color"],
                   linestyle=st["linestyle"], linewidth=st["linewidth"], zorder=st["zorder"])

            ax.axhline(0, linestyle="--", linewidth=0.8, color="gray")
            ax.set_title(
                f"Normalized Error in {title}" if j != 1 else f"Error in {title}",
                fontsize=self.font["label"], fontweight="bold",
            )
            ax.set_ylabel("Error", fontsize=self.font["label"])
            ax.tick_params(labelsize=self.font["tick"])
            ax.grid(True, which="both" if j != 1 else "major", linestyle="--", linewidth=0.5)
            if j // 2 == 1:
                ax.set_xlabel("Distance (m)", fontsize=self.font["label"])
                _auto_xticks(ax, self.reference.r)
            if j == 3:
                ax.legend(fontsize=self.font["legend"])

        fig.suptitle(suptitle, fontsize=self.font["suptitle"], fontweight="bold")
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        _save(fig, output_dir, fname)
        return fig


# ---------------------------------------------------------------------------
# 3. ReceiverLineCombined  –  2×4 single dataset
# ---------------------------------------------------------------------------

class ReceiverLineCombined:
    """
    2×4 figure for one dataset vs an analytical reference.
      Row 1 : field curves  (analytical + dataset)
      Row 2 : error curves
    Phase column is placed last.
    """

    def __init__(
        self,
        dataset:          GPRDataset,
        analytical:       GPRDataset,
        quantities:       list[tuple] = None,
        ds_color:         str = "#1f77b4",
        err_color:        str = "#ff7f0e",
        analytical_label: str = "Semi-Analytical Quadrature",
        font:             dict = None,
        figsize:          tuple = None,
    ):
        self.dataset          = dataset
        self.analytical       = analytical
        self.quantities       = quantities or QUANTITIES
        self.ds_color         = ds_color
        self.err_color        = err_color
        self.analytical_label = analytical_label
        self.font             = {"suptitle": 24, "label": 18, "tick": 14, "legend": 14, **(font or {})}
        ncols = len(self.quantities)
        self.figsize = figsize or (6 * ncols, 10)

    def _col_order(self):
        phase = next((i for i, (t, _) in enumerate(self.quantities) if "Phase" in t), None)
        n = len(self.quantities)
        return ([i for i in range(n) if i != phase] + [phase]) if phase is not None else list(range(n))

    def plot(self, suptitle: str = "", output_dir: str = None, fname: str = None) -> plt.Figure:
        ncols   = len(self.quantities)
        ordered = self._col_order()
        lw      = 2.5
        fig, axes = plt.subplots(2, ncols, figsize=self.figsize, sharex="col")
        panel_labels = [f"({c})" for c in "abcdefgh"[:2 * ncols]]

        for col, j in enumerate(ordered):
            title, ylabel = self.quantities[j]

            # top row – fields
            ax = axes[0, col]
            ax.text(-0.12, 1.08, panel_labels[col], transform=ax.transAxes,
                    fontsize=self.font["label"], fontweight="bold", va="top")
            fn = ax.semilogy if j == 0 else ax.plot
            ln = fn(self.analytical.r, self.analytical.field(j),
                    color="black", linestyle="-", linewidth=lw * 1.5,
                    label=self.analytical_label, zorder=2)[0]
            ln.set_path_effects([pe.Stroke(linewidth=lw * 1.5 + 2, foreground="white"), pe.Normal()])
            fn(self.dataset.r, self.dataset.field(j),
               color=self.ds_color, linestyle="--", linewidth=lw,
               label=self.dataset.label, zorder=3)
            ax.set_title(title, fontsize=self.font["label"], fontweight="bold")
            ax.set_ylabel(ylabel, fontsize=self.font["label"])
            ax.tick_params(labelsize=self.font["tick"])
            ax.grid(True, which="both" if j == 0 else "major", linestyle="--", linewidth=0.5)
            if col == 0:
                ax.legend(fontsize=self.font["legend"])

            # bottom row – errors
            ax = axes[1, col]
            ax.text(-0.12, 1.08, panel_labels[ncols + col], transform=ax.transAxes,
                    fontsize=self.font["label"], fontweight="bold", va="top")
            err = field_error(self.analytical, self.dataset, j)
            fn2 = ax.semilogy if j in (0, 2, 3) else ax.plot
            fn2(self.dataset.r, err, color=self.err_color, linewidth=lw, label="Error")
            ax.axhline(0, linestyle="--", linewidth=1, color="gray")
            ax.set_title(
                f"Normalized Error in {title}" if j != 1 else f"Error in {title}",
                fontsize=self.font["label"], fontweight="bold",
            )
            ax.set_ylabel("Error", fontsize=self.font["label"])
            ax.set_xlabel("Distance (m)", fontsize=self.font["label"])
            ax.tick_params(labelsize=self.font["tick"])
            ax.grid(True, which="both" if j == 0 else "major", linestyle="--", linewidth=0.5)
            if col == 0:
                ax.legend(fontsize=self.font["legend"])

        fig.suptitle(suptitle, fontsize=self.font["suptitle"], fontweight="bold")
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        _save(fig, output_dir, fname)
        return fig


# ---------------------------------------------------------------------------
# 4. ReceiverLineCombinedMulti  –  2×4 multiple datasets
# ---------------------------------------------------------------------------

class ReceiverLineCombinedMulti:
    """
    2×4 figure for multiple datasets vs an analytical reference.
      Row 1 : field curves  (analytical + all datasets)
      Row 2 : error curves
    """

    def __init__(
        self,
        datasets:         list[GPRDataset],
        analytical:       GPRDataset,
        quantities:       list[tuple] = None,
        analytical_label: str = "Analytical Solution",
        base_lw:          float = 2.5,
        font:             dict = None,
        figsize:          tuple = None,
    ):
        self.datasets         = datasets
        self.analytical       = analytical
        self.quantities       = quantities or QUANTITIES
        self.analytical_label = analytical_label
        self.styles           = _build_styles(datasets, base_lw)
        self.font             = {"suptitle": 24, "label": 18, "tick": 14, "legend": 14, **(font or {})}
        ncols = len(self.quantities)
        self.figsize = figsize or (6 * ncols, 10)

    def _col_order(self):
        phase = next((i for i, (t, _) in enumerate(self.quantities) if "Phase" in t), None)
        n = len(self.quantities)
        return ([i for i in range(n) if i != phase] + [phase]) if phase is not None else list(range(n))

    def plot(self, suptitle: str = "", output_dir: str = None, fname: str = None) -> plt.Figure:
        ncols   = len(self.quantities)
        ordered = self._col_order()
        lw      = 2.5
        fig, axes = plt.subplots(2, ncols, figsize=self.figsize, sharex="col")
        panel_labels = [f"({c})" for c in "abcdefgh"[:2 * ncols]]

        for col, j in enumerate(ordered):
            title, ylabel = self.quantities[j]

            # top row
            ax = axes[0, col]
            ax.text(-0.12, 1.08, panel_labels[col], transform=ax.transAxes,
                    fontsize=self.font["label"], fontweight="bold", va="top")
            fn = ax.semilogy if j == 0 else ax.plot
            ln = fn(self.analytical.r, self.analytical.field(j),
                    color="black", linestyle="-", linewidth=lw * 1.5,
                    label=self.analytical_label, zorder=1)[0]
            ln.set_path_effects([pe.Stroke(linewidth=lw * 1.5 + 2, foreground="white"), pe.Normal()])
            for ds in self.datasets:
                st = self.styles[ds.label]
                fn(ds.r, ds.field(j), color=st["color"], linestyle=st["linestyle"],
                   linewidth=st["linewidth"], label=st["label"], zorder=st["zorder"])
            ax.set_title(title, fontsize=self.font["label"], fontweight="bold")
            ax.set_ylabel(ylabel, fontsize=self.font["label"])
            ax.tick_params(labelsize=self.font["tick"])
            ax.grid(True, which="both" if j == 0 else "major", linestyle="--", linewidth=0.5)
            if col == 0:
                ax.legend(fontsize=self.font["legend"])

            # bottom row
            ax = axes[1, col]
            ax.text(-0.12, 1.08, panel_labels[ncols + col], transform=ax.transAxes,
                    fontsize=self.font["label"], fontweight="bold", va="top")
            for ds in self.datasets:
                err = field_error(self.analytical, ds, j)
                st  = self.styles[ds.label]
                fn2 = ax.semilogy if j in (0, 2, 3) else ax.plot
                fn2(ds.r, err, color=st["color"], linestyle=st["linestyle"],
                    linewidth=st["linewidth"], label=st["label"], zorder=st["zorder"])
            ax.axhline(0, linestyle="--", linewidth=1, color="gray")
            ax.set_title(
                f"Normalized Error in {title}" if j != 1 else f"Error in {title}",
                fontsize=self.font["label"], fontweight="bold",
            )
            ax.set_ylabel("Error", fontsize=self.font["label"])
            ax.set_xlabel("Distance (m)", fontsize=self.font["label"])
            ax.tick_params(labelsize=self.font["tick"])
            ax.grid(True, which="both" if j == 0 else "major", linestyle="--", linewidth=0.5)
            if col == 0:
                ax.legend(fontsize=self.font["legend"])

        fig.suptitle(suptitle, fontsize=self.font["suptitle"], fontweight="bold")
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        _save(fig, output_dir, fname)
        return fig


# ---------------------------------------------------------------------------
# 5. ErrorHistogramPlot  –  1×2 amplitude + phase histograms
# ---------------------------------------------------------------------------

class ErrorHistogramPlot:
    """
    1×2 histogram of amplitude (log-spaced bins, %) and phase errors.
    Legend shows mean, std, and max for each dataset.
    """

    def __init__(
        self,
        datasets:  list[GPRDataset],
        reference: GPRDataset,
        base_lw:   float = 2.5,
        font:      dict = None,
        figsize:   tuple = (12, 4),
    ):
        self.datasets  = datasets
        self.reference = reference
        self.styles    = _build_styles(datasets, base_lw)
        self.font      = {"suptitle": 27, "label": 18, "tick": 18, "legend": 10, **(font or {})}
        self.figsize   = figsize

    def plot(self, suptitle: str = "", output_dir: str = None, fname: str = None) -> plt.Figure:
        fig, axes = plt.subplots(1, 2, figsize=self.figsize)

        for ds in self.datasets:
            st        = self.styles[ds.label]
            amp_err   = field_error(self.reference, ds, 0)
            phase_err = field_error(self.reference, ds, 1)

            # amplitude
            ax  = axes[0]
            fin = amp_err[np.isfinite(amp_err)]
            pos = fin[fin > 0]
            bins = np.logspace(np.log10(pos.min()), np.log10(pos.max()), 40)
            m, s, mx = error_stats(amp_err)
            lbl = f"{m*100:.2f}%,  {s*100:.2f}%,  {mx*100:.2f}%"
            ax.hist(fin, bins=bins, color=st["color"], alpha=0.7, label=lbl)
            ax.set_xscale("log")
            ax.set_title("Histogram of Abs. Norm. Error\nAmplitude (Ex)",
                         fontsize=self.font["label"], fontweight="bold")
            ax.set_xlabel("Absolute Normalised Error", fontsize=self.font["label"])
            ax.set_ylabel("Count", fontsize=self.font["label"])
            ax.tick_params(labelsize=self.font["tick"])
            ax.grid(True, linestyle="--", linewidth=0.5)
            ax.legend(fontsize=self.font["legend"], title="mean,  std,  max")

            # phase
            ax  = axes[1]
            fin = phase_err[np.isfinite(phase_err)]
            m, s, mx = error_stats(phase_err)
            lbl = f"{m:.4f},  {s:.4f},  {mx:.4f}"
            ax.hist(fin, bins=40, color=st["color"], alpha=0.7, label=lbl)
            ax.set_title("Histogram of Error\nPhase (Ex)",
                         fontsize=self.font["label"], fontweight="bold")
            ax.set_xlabel("Error (rad)", fontsize=self.font["label"])
            ax.set_ylabel("Count", fontsize=self.font["label"])
            ax.tick_params(labelsize=self.font["tick"])
            ax.grid(True, linestyle="--", linewidth=0.5)
            ax.legend(fontsize=self.font["legend"], title="mean,  std,  max")

        fig.suptitle(suptitle, fontsize=self.font["suptitle"], fontweight="bold")
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        _save(fig, output_dir, fname)
        return fig


# ---------------------------------------------------------------------------
# 6. ErrorStatPlot  –  mean / std / max vs sweep parameter
# ---------------------------------------------------------------------------

class ErrorStatPlot:
    """
    2×2 figure showing mean, std, and max absolute error as a function of a
    scalar parameter (e.g. PML thickness).

    Parameters
    ----------
    datasets     : one GPRDataset per parameter value, in order
    reference    : reference dataset
    param_values : scalar values corresponding to each dataset
    xtick_labels : optional custom labels (e.g. LaTeX λ/n strings)
    """

    def __init__(
        self,
        datasets:     list[GPRDataset],
        reference:    GPRDataset,
        param_values: list[float],
        quantities:   list[tuple] = None,
        xtick_labels: list[str] = None,
        font:         dict = None,
        figsize:      tuple = (14, 10),
    ):
        self.datasets     = datasets
        self.reference    = reference
        self.param_values = param_values
        self.quantities   = quantities or QUANTITIES
        self.xtick_labels = xtick_labels
        self.font         = {"suptitle": 27, "label": 18, "tick": 14, "legend": 12, **(font or {})}
        self.figsize      = figsize

    def plot(self, suptitle: str = "", output_dir: str = None, fname: str = None) -> plt.Figure:
        # Collect per-dataset stats for each quantity
        means = [[] for _ in range(4)]
        stds  = [[] for _ in range(4)]
        maxs  = [[] for _ in range(4)]

        for ds in self.datasets:
            for qi in range(4):
                err = field_error(self.reference, ds, qi)
                m, s, mx = error_stats(err)
                means[qi].append(m)
                stds[qi].append(s)
                maxs[qi].append(mx)

        x = np.arange(len(self.datasets))
        fig, axes = plt.subplots(2, 2, figsize=self.figsize)

        for qi, (title, _) in enumerate(self.quantities):
            ax = axes[qi // 2, qi % 2]
            ax.plot(x, means[qi], "o-", label="Mean",    linewidth=2)
            ax.plot(x, stds[qi],  "s-", label="Std",     linewidth=2)
            ax.plot(x, maxs[qi],  "^-", label="Max abs", linewidth=2)
            ax.set_title(
                f"Normalized Error in {title}" if qi != 1 else f"Error in {title}",
                fontsize=self.font["label"], fontweight="bold",
            )
            ax.set_ylabel("Error", fontsize=self.font["label"])
            ax.set_xlabel("Parameter", fontsize=self.font["label"])
            ax.tick_params(labelsize=self.font["tick"])
            ax.grid(True, linestyle="--", linewidth=0.5)
            ax.legend(fontsize=self.font["legend"])

            if self.xtick_labels:
                ax.xaxis.set_major_locator(FixedLocator(x))
                ax.xaxis.set_major_formatter(FixedFormatter(self.xtick_labels))
            else:
                ax.set_xticks(x)
                ax.set_xticklabels([str(v) for v in self.param_values],
                                    fontsize=self.font["tick"])

        fig.suptitle(suptitle, fontsize=self.font["suptitle"], fontweight="bold")
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        _save(fig, output_dir, fname)
        return fig
