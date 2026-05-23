"""
run.py
======
Example entry point showing how fieldreader / postprocess / visualize
replace the repeated boilerplate in the thesis notebooks.

Each section corresponds to one experiment from the original notebooks.
Uncomment the sections you need and set the paths at the top.
"""

import os
import numpy as np

from .fieldreader import AnalyticalLoader, ElfeLoader, load_elfe_batch
from .postprocess import field_error, all_errors, error_stats
from .visualize   import (
    ReceiverLinePlot,
    ReceiverLineErrorPlot,
    ReceiverLineCombined,
    ReceiverLineCombinedMulti,
    ErrorHistogramPlot,
    ErrorStatPlot,
)

# ---------------------------------------------------------------------------
# Paths  –  edit these
# ---------------------------------------------------------------------------

BASE        = r"F:\Projects\EMGeoInversion\Tests_Thesis"
ANALYTICAL  = os.path.join(BASE, "semi-analytic_100MHz")


# ===========================================================================
# Section 1 – Half-space CMP  (plotting_th2)
# ===========================================================================

POST = os.path.join(BASE, "6HS", "postprocess")

evert   = AnalyticalLoader(os.path.join(ANALYTICAL, "Exx_single_freq_4_100MHz.csv"),   label="Evert").endfire()
empymod = AnalyticalLoader(os.path.join(ANALYTICAL, "GPR-2001-4-dlf.csv"), label="empymod - 2001 DLF").endfire()

loaders = load_elfe_batch(
    base_folder=os.path.join(BASE, "6HS"),
    run_names=["out_HF_l1d_l2PML_CMP_BA", "out_HF_l1d_l2PML_CMP_BK", "out_HF_l1d_l2PML_CMP_F"],
    labels=[r"Uniform $k_\text{min}$", r"Uniform $k_\text{max}$", "Varying Stretch"],
    num_endfire=256,
)
runs = [l.endfire() for l in loaders]

ReceiverLinePlot([evert, empymod] + runs).plot(
    suptitle="Endfire – Half-Space CMP",
    output_dir=POST, fname="hs_cmp_comparison.png",
)
ReceiverLineErrorPlot([evert, empymod] + runs, reference=evert).plot(
    suptitle="Errors vs Evert – Half-Space CMP",
    output_dir=POST, fname="hs_cmp_error.png",
)
ReceiverLineCombinedMulti(runs, analytical=evert).plot(
    suptitle="Combined – Half-Space CMP",
    output_dir=POST, fname="hs_cmp_combined.png",
)


# ===========================================================================
# Section 2 – Two-layer CO  (plotting_th3)
# ===========================================================================

POST2 = os.path.join(BASE, "6TL", "postprocess")

evert_tl = AnalyticalLoader(
    os.path.join(ANALYTICAL, "Exx_single_freq_4_9_100MHz_NR.csv"), label="Evert"
).endfire()

loaders2 = load_elfe_batch(
    base_folder=os.path.join(BASE, "6TL"),
    run_names=["out_TL_l1d_l2PML_CO_BA", "out_TL_l1d_l2PML_CO_F"],
    labels=[r"Uniform Stretch – Air Like PML 1.5 m", r"Varying Stretch – Earth-Like PML 0.75 m"],
    num_endfire=48,
)
runs2 = [l.endfire() for l in loaders2]

ReceiverLinePlot([evert_tl] + runs2).plot(
    suptitle="Endfire – Two-Layer CO",
    output_dir=POST2, fname="tl_co_comparison.png",
)
ReceiverLineErrorPlot(runs2, reference=evert_tl).plot(
    suptitle="Errors vs Evert – Two-Layer CO",
    output_dir=POST2, fname="tl_co_error.png",
)
ErrorHistogramPlot(runs2, reference=evert_tl).plot(
    suptitle="Error Histogram – Two-Layer CO",
    output_dir=POST2, fname="tl_co_hist.png",
)


# ===========================================================================
# Section 3 – PML thickness convergence  (plotting_thesis_final)
# ===========================================================================

POST3 = os.path.join(BASE, "Fin", "postprocess")

evert_fin = AnalyticalLoader(
    os.path.join(ANALYTICAL, "Exx_single_freq_4_9_100MHz_NR.csv"), label="Evert"
).endfire()

wave    = 3.0   # wavelength in second layer (m)
denoms  = [10, 12.5, 15, 17.5, 20, 22.5, 25]
loaders3 = load_elfe_batch(
    base_folder=os.path.join(BASE, "Fin"),
    run_names=[f"run_pml_{str(d).replace('.', 'p')}" for d in denoms],
    labels=[rf"$\lambda/{d}$" for d in denoms],
    num_endfire=48,
)
runs3 = [l.endfire() for l in loaders3]

xtick_labels = [rf"$\dfrac{{\lambda}}{{{d}}}$" for d in denoms]

ErrorStatPlot(
    runs3, reference=evert_fin,
    param_values=[wave / d for d in denoms],
    xtick_labels=xtick_labels,
).plot(
    suptitle="Error Statistics vs PML Thickness",
    output_dir=POST3, fname="pml_convergence.png",
)
