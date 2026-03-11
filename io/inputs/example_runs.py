"""
example_runs.py
---------------
Usage examples for GPRSurvey showing all combinations:

  Case 1 — Air only, no anomaly       (num_layers=0, anomaly=None)
  Case 2 — Layered earth, no anomaly  (num_layers=2, anomaly=None)
  Case 3 — Layered earth + anomaly    (num_layers=2, anomaly=BoxAnomaly)

"""

import numpy as np
from survey import GPRSurvey

# Shared frequency / wavelength setup

f = 100e6                   # 100 MHz central frequency
wave = 3e8 / f              # free-space wavelength in air = 3.0 m


# =============================================================================
# Case 1: Air only, no anomaly
# =============================================================================

survey_air = GPRSurvey.build(
    experiment_name="air_only",

    # Domain
    x_e=[-wave/10, 1 + wave/10],
    y_e=[-wave/10, wave/10],
    z_e=[-wave/10, wave/10],

    # Materials — air only, no earth layers
    air_eps_r=1.0,
    air_sigma=1e-16,
    layer_thicknesses=[], # no earth layers
    layer_eps_r=[],
    layer_sigma=[],

    # No anomaly — omit anomaly_x/y/z entirely

    # Source
    ricker_central_f=f,
    num_points_per_range=1,
    antenna_position=[0.0, 0.0, 0.025],
    source_type=6,
    current_direction=1,
    num_segments=1,
    s_f=250,
    bh_f=1.0,
    box_present=False,
    box_x=[-1.0, 1.0],

    # Receivers
    num_receivers_inline=48,
    num_receivers_endfire=0,
    num_receivers_oblique=0,

    # Solver
    solver_type=2,
    max_ref_steps=0,
    max_unknowns=5_000_000,
    accuracy_tol=3e-5,
    output_fields_vtk=1,

    # PML
    num_pml_layers=1,
    pml_layer_thickness=wave/10,
    pml_type="lin",
    pml_decay_type=1,

    # Mesh sizing
    least_samples_per_wavelength=20,
)

survey_air.generate()


# =============================================================================
# Case 2: Two earth layers, no anomaly
# =============================================================================

survey_layered = GPRSurvey.build(
    experiment_name="layered_no_anomaly",

    # Domain
    x_e=[-wave/10, 1 + wave/10],
    y_e=[-wave/10, wave/10],
    z_e=[-1.0 - wave/10/3, wave/10],

    # Materials — air + 2 earth layers
    air_eps_r=1.0,
    air_sigma=1e-16,
    layer_thicknesses=[1.0, wave/10/3],
    layer_eps_r=[4.0, 9.0],
    layer_sigma=[1e-4, 1e-3],
    layer_mu_r=[1.0, 1.0],
    layer_sigma_m=[0.0, 0.0],

    # No anomaly — omit anomaly_x/y/z entirely

    # Source
    ricker_central_f=f,
    num_points_per_range=1,
    antenna_position=[0.0, 0.0, 0.025],
    source_type=6,
    current_direction=1,
    num_segments=1,
    s_f=250,
    bh_f=1.0,
    box_present=False,
    box_x=[-1 + 0.75, 1 + 0.375],

    # Receivers
    num_receivers_inline=48,
    num_receivers_endfire=0,
    num_receivers_oblique=0,

    # Solver
    solver_type=2,
    max_ref_steps=0,
    max_unknowns=5_000_000,
    accuracy_tol=3e-5,
    output_fields_vtk=1,

    # PML
    num_pml_layers=1,
    pml_layer_thickness=wave/10,
    pml_type="lin",
    pml_decay_type=1,

    least_samples_per_wavelength=20,
)

survey_layered.generate()


# ==================================
# Case 3: Two earth layers + anomaly
# ==================================

survey_anomaly = GPRSurvey.build(
    experiment_name="AnAir",

    # Domain
    x_e=[-wave/10, 1 + wave/10],
    y_e=[-wave/10, wave/10],
    z_e=[-1.0 - wave/10/3, wave/10],

    # Materials — air + 2 earth layers
    air_eps_r=1.0,
    air_sigma=1e-16,
    layer_thicknesses=[1.0, wave/10/3],
    layer_eps_r=[4.0, 9.0],
    layer_sigma=[1e-4, 1e-3],
    layer_mu_r=[1.0, 1.0],
    layer_sigma_m=[0.0, 0.0],

    # Anomaly — all three coordinate tuples + properties required together
    anomaly_x=(0, wave/8),
    anomaly_y=(-wave/20, wave/20),
    anomaly_z=(-0.9, -0.5),
    anomaly_properties=(20, 1e-4, 1.0, 0.0),   # (eps_r, sigma, mu_r, sigma_m)

    # Source
    ricker_central_f=f,
    num_points_per_range=1,
    antenna_position=[0.0, 0.0, 0.025],
    source_type=6,
    current_direction=1,
    num_segments=1,
    s_f=250,
    bh_f=1.0,
    box_present=False,
    box_x=[-1 + 0.75, 1 + 0.375],
    m=5,

    # Receivers
    num_receivers_inline=48,
    num_receivers_endfire=0,
    num_receivers_oblique=0,

    # Solver
    solver_type=2,
    max_ref_steps=0,
    max_unknowns=5_000_000,
    accuracy_tol=3e-5,
    output_fields_vtk=1,

    # PML
    num_pml_layers=1,
    pml_layer_thickness=wave/10,
    pml_type="lin",
    pml_decay_type=1,

    least_samples_per_wavelength=20,
)

survey_anomaly.generate()
