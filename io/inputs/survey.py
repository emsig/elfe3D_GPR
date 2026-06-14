"""
survey.py
---------
Top-level entry point for elfe3D GPR input generation.

GPRSurvey composes all domain objects and runs the full input generation
pipeline: mesh sizing → assembly → .poly file → FEM input files.

It handles all combinations of:
  - num_layers = 0           : air-only domain (no earth)
  - num_layers >= 1          : layered earth
  - anomalies = []           : no anomaly
  - anomalies = [BoxAnomaly] : one or more rectangular bodies
  - anomalies = [SphereAnomaly] : one or more spherical bodies
  - mixed list               : any combination of the above
  - box_present = True/False : source refinement box

Two ways to supply anomalies via GPRSurvey.build():

  Option A — pass a ready-built list (preferred for multiple anomalies):
      anomalies=[
          BoxAnomaly(x=(...), y=(...), z=(...), properties=(...)),
          SphereAnomaly(center=(...), radius=..., properties=(...)),
      ]

  Option B — legacy flat parameters (single BoxAnomaly, fully backward-compat):
      anomaly_x=(-0.3, 0.3), anomaly_y=(...), anomaly_z=(...),
      anomaly_properties=(eps_r, sigma, mu_r, sigma_m)

Both options can be combined; Option B appends to whatever is in Option A.

Usage
-----
survey = GPRSurvey.build(...)
survey.generate()

Or step by step:
    survey = GPRSurvey.build(...)
    survey.validate()
    survey.generate()
"""

import numpy as np
from pathlib import Path
from dataclasses import dataclass, field

from .materials import Material
from .geolayers import GeoLayer, LayerStack
from .domain import ModelDomain
from .anomalies import BoxAnomaly, SphereAnomaly
from .sources import SourceAntenna
from .receivers import ReceiverArray
from .solver import SolverConfig
from .pml import PMLConfig
from .writetetgenpoly import PolyAssembler
from .writeinputfiles import FEMInputWriter


# =============================================================================
# IO path configuration
# =============================================================================

@dataclass
class IOConfig:
    """
    All file paths for one simulation run.

    Parameters
    ----------
    experiment_name : short label used to name all files and folders
    base_dir        : root directory where input/output folders are created
    """

    experiment_name: str
    base_dir: Path | str = "."

    def __post_init__(self) -> None:
        self.base_dir = Path(self.base_dir)

    @property
    def input_dir(self) -> Path:
        return self.base_dir / f"in_{self.experiment_name}"

    @property
    def output_dir(self) -> Path:
        return self.base_dir / f"out_{self.experiment_name}"

    @property
    def poly_file(self) -> Path:
        return self.input_dir / f"GPR_model_{self.experiment_name}.poly"

    @property
    def model_file_stem(self) -> str:
        """Path written into elfe3D_input.txt for the mesh file (no extension)."""
        return str(f"GPR_model_{self.experiment_name}.")

    @property
    def output_E_file(self) -> str:
        return str(f"../out_{self.experiment_name}/electric_fields")

    @property
    def output_H_file(self) -> str:
        return str(f"../out_{self.experiment_name}/magnetic_fields")

    def ensure_dirs(self) -> None:
        self.input_dir.mkdir(parents=True, exist_ok=True)
        self.output_dir.mkdir(parents=True, exist_ok=True)


# =============================================================================
# GPRSurvey — the single object you build and call .generate() on
# =============================================================================

@dataclass
class GPRSurvey:
    """
    Complete specification of one elfe3D GPR forward simulation.

    Do not instantiate directly — use GPRSurvey.build() which handles
    all conditional setup depending on num_layers and anomaly presence.

    Fields
    ------
    domain    : ModelDomain
    layers    : LayerStack  (air + any earth layers)
    source    : SourceAntenna
    receivers : ReceiverArray
    solver    : SolverConfig
    pml       : PMLConfig
    io        : IOConfig
    anomalies : list of BoxAnomaly / SphereAnomaly objects (may be empty)
    """

    domain:    ModelDomain
    layers:    LayerStack
    source:    SourceAntenna
    receivers: ReceiverArray
    solver:    SolverConfig
    pml:       PMLConfig
    io:        IOConfig
    anomalies: list = field(default_factory=list)

    # ------------------------------------------------------------------
    # Named constructor — the main entry point
    # ------------------------------------------------------------------

    @classmethod
    def build(
        cls,
        # ── Experiment identity ──────────────────────────────────────
        experiment_name: str,
        base_dir: str | Path = ".",

        # ── Domain extents ───────────────────────────────────────────
        x_e: list[float] = None,
        y_e: list[float] = None,
        z_e: list[float] = None,
        num_air_wavelengths: float = 0.0,

        # ── Materials ────────────────────────────────────────────────
        # air is always required
        air_eps_r: float = 1.0,
        air_sigma: float = 1e-16,
        # earth layers: parallel lists, one entry per layer
        layer_thicknesses: list[float] = None,    # [] means no earth layers
        layer_eps_r:       list[float] = None,    # one per earth layer
        layer_sigma:       list[float] = None,
        layer_mu_r:        list[float] = None,
        layer_sigma_m:     list[float] = None,

        # ── Anomalies ────────────────────────────────────────────────
        # Option A: pass a ready-built list of BoxAnomaly / SphereAnomaly objects
        anomalies: list | None = None,
        # Option B: legacy single-box flat parameters (still fully supported)
        anomaly_x: tuple | None = None,           # (x_min, x_max) or None
        anomaly_y: tuple | None = None,
        anomaly_z: tuple | None = None,
        anomaly_properties: tuple | None = None,  # (eps_r, sigma, mu_r, sigma_m)

        # ── Source ───────────────────────────────────────────────────
        source_type: int = 6,
        antenna_position: list[float] = None,
        current_direction: int = 1,
        current: float = 1.0,
        source_moment: int = 1,
        num_segments: int = 1,
        f_list: list[float] | None = None,
        ricker_central_f: float = 100e6,
        num_points_per_range: int = 1,
        m: int = 5,
        box_present: bool = False,
        s_f: float = 500.0,
        bh_f: float = 1.0,
        box_x: list[float] = None,

        # ── Receivers ────────────────────────────────────────────────
        num_receivers_inline:  int = 0,
        num_receivers_endfire: int = 0,
        num_receivers_oblique: int = 0,

        # ── Solver ───────────────────────────────────────────────────
        solver_type: int = 2,
        max_ref_steps: int = 0,
        max_unknowns: int = 5_000_000,
        beta_ref: float = 0.85,
        accuracy_tol: float = 3e-5,
        vtk: int = 0,
        error_est_method: int = 3,
        ref_strategy: int = 1,
        output_fields_vtk: int = 1,

        # ── PML ──────────────────────────────────────────────────────
        num_pml_layers: int = 1,
        pml_layer_thickness: float = 0.3,
        pml_type: str = "lin",
        pml_theory: int = 0,
        pml_decay_type: int = 1,

        # ── Mesh sizing ──────────────────────────────────────────────
        least_samples_per_wavelength: int = 20,

    ) -> "GPRSurvey":
        """
        Build a GPRSurvey from flat parameter lists.

        Handles:
          - num_layers = 0  (air only, no earth layers)
          - num_layers >= 1 (layered earth)
          - anomalies = []  (no anomaly)
          - anomalies present: any mix of BoxAnomaly / SphereAnomaly
        """

        # ── Defaults for mutable arguments ──────────────────────────
        layer_thicknesses = layer_thicknesses or []
        layer_eps_r       = layer_eps_r       or []
        layer_sigma       = layer_sigma       or []
        layer_mu_r        = layer_mu_r        or [1.0] * len(layer_thicknesses)
        layer_sigma_m     = layer_sigma_m     or [0.0] * len(layer_thicknesses)
        antenna_position  = antenna_position  or [0.0, 0.0, 0.025]
        box_x             = box_x             or [-1.0, 1.0]

        num_layers = len(layer_thicknesses)

        # ── Validate layer parameter consistency ─────────────────────
        for name, lst in [
            ("layer_eps_r",   layer_eps_r),
            ("layer_sigma",   layer_sigma),
            ("layer_mu_r",    layer_mu_r),
            ("layer_sigma_m", layer_sigma_m),
        ]:
            if len(lst) != num_layers:
                raise ValueError(
                    f"{name} has {len(lst)} entries but layer_thicknesses "
                    f"has {num_layers}. They must match."
                )

        # ── Build anomalies list ──────────────────────────────────────
        # Start from Option A list (or empty), then append any Option B box.
        anomalies_list: list = list(anomalies) if anomalies else []

        anomaly_coords = [anomaly_x, anomaly_y, anomaly_z]
        if any(c is not None for c in anomaly_coords):
            if not all(c is not None for c in anomaly_coords):
                raise ValueError(
                    "anomaly_x, anomaly_y, anomaly_z must all be provided "
                    "together, or all left as None."
                )
            if anomaly_properties is None:
                raise ValueError(
                    "anomaly_properties (eps_r, sigma, mu_r, sigma_m) "
                    "must be provided when anomaly_x/y/z are set."
                )
            anomalies_list.append(
                BoxAnomaly(
                    x=anomaly_x,
                    y=anomaly_y,
                    z=anomaly_z,
                    properties=anomaly_properties,
                )
            )

        # ── Build materials ──────────────────────────────────────────
        air = Material(
            name="air", eps_r=air_eps_r, sigma=air_sigma
        )

        earth_layers = [
            GeoLayer(
                thickness=layer_thicknesses[i],
                material=Material(
                    name=f"layer_{i+1}",
                    eps_r=layer_eps_r[i],
                    sigma=layer_sigma[i],
                    mu_r=layer_mu_r[i],
                    sigma_m=layer_sigma_m[i],
                ),
            )
            for i in range(num_layers)
        ]
        layers = LayerStack(air=air, layers=earth_layers)

        # ── Build domain ─────────────────────────────────────────────
        domain = ModelDomain(
            x_e=x_e or [-1.0, 1.0],
            y_e=y_e or [-0.1, 0.1],
            z_e=z_e or [-1.0, 0.1],
            num_air_wavelengths=num_air_wavelengths,
        )

        # ── PML ──────────────────────────────────────────────────────
        pml = PMLConfig(
            num_layers=num_pml_layers,
            layer_thickness=pml_layer_thickness,
            pml_type=pml_type,
            pml_theory=pml_theory,
            pml_decay_type=pml_decay_type,
        )

        # ── Solver ───────────────────────────────────────────────────
        solver = SolverConfig(
            solver=solver_type,
            max_ref_steps=max_ref_steps,
            max_unknowns=max_unknowns,
            beta_ref=beta_ref,
            accuracy_tol=accuracy_tol,
            vtk=vtk,
            error_est_method=error_est_method,
            ref_strategy=ref_strategy,
            output_fields_vtk=output_fields_vtk,
            pml_decay_type=pml_decay_type,
        )

        # ── Frequency list and mesh sizing ───────────────────────────
        if f_list is None:
            if ricker_central_f != 100e6 or num_points_per_range != 1:
                raise ValueError(
                    "Legacy ricker_central_f/num_points_per_range behavior is deprecated. "
                    "Use f_list=[...frequency values in Hz...] instead."
                )
            f_list = [100e6]
        if len(f_list) == 0:
            raise ValueError("f_list must contain at least one frequency value.")
        f_low = f_list[0]

        layers.compute_all_mesh_sizing(f_low, least_samples_per_wavelength)

        # Print mesh sizing info (matches original print statements)
        print(f"odepths: {[e/4.0 for e in layers.max_element_edges]}")

        # ── Source ───────────────────────────────────────────────────
        source = SourceAntenna(
            source_type=source_type,
            antenna_position=antenna_position,
            current_direction=current_direction,
            current=current,
            source_moment=source_moment,
            num_segments=num_segments,
            f_list=f_list,
            m=m,
            box_present=box_present,
            s_f=s_f,
            bh_f=bh_f,
            box_x_extents=box_x,
            num_layers=num_layers,
            largest_wavelengths=layers.max_element_edges,
            # NOTE: largest_wavelengths here uses max_element_edges as a proxy
            # for wavelength-based sizing, matching original logic.
            # The original used regions['largest_wavelengths'] which is
            # wavelength_from_f(f_low, eps_r) — replace with below:
            # largest_wavelengths=[m.wavelength_at(f_low) for m in layers.all_materials]
        )

        print(f"Source antenna length: {source.length} m")
        if box_present:
            print(f"Source antenna box extents: {source.box_x} {source.box_y} {source.box_z}")

        # ── Receivers ────────────────────────────────────────────────
        receivers = ReceiverArray(
            num_inline=num_receivers_inline,
            num_endfire=num_receivers_endfire,
            num_oblique=num_receivers_oblique,
            source_length=source.length,
            pml_type=pml_type,
            pml_thickness=pml_layer_thickness,
            num_pml_layers=num_pml_layers,
        )
        print(f"Receiver antenna depth: {receivers.depth} m")

        # ── IO paths ─────────────────────────────────────────────────
        io = IOConfig(
            experiment_name=experiment_name,
            base_dir=base_dir,
        )

        return cls(
            domain=domain,
            layers=layers,
            source=source,
            receivers=receivers,
            solver=solver,
            pml=pml,
            io=io,
            anomalies=anomalies_list,
        )

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def validate(self) -> None:
        """
        Cross-component validation. Raises ValueError for any
        configuration that would produce a bad TetGen input.
        """
        # Source must be inside domain
        sx, sy, sz = self.source.antenna_position
        if not self.domain.contains(sx, sy, sz):
            raise ValueError(
                f"Source position {self.source.antenna_position} is outside "
                f"the domain {self.domain}."
            )

        # All receivers must be inside domain
        for i, (rx, ry, rz) in enumerate(
            zip(self.receivers.x, self.receivers.y, self.receivers.z)
        ):
            if not self.domain.contains(rx, ry, rz):
                raise ValueError(
                    f"Receiver {i} at ({rx}, {ry}, {rz}) is outside the domain."
                )

        # Every anomaly centroid must be inside domain
        for i, anomaly in enumerate(self.anomalies):
            cx, cy, cz = anomaly.centroid
            if not self.domain.contains(cx, cy, cz):
                raise ValueError(
                    f"Anomaly {i} ({type(anomaly).__name__}) centroid "
                    f"{anomaly.centroid} is outside the domain."
                )

        # Mesh sizing must have been computed
        for mat in self.layers.all_materials:
            if mat.max_element_volume is None:
                raise ValueError(
                    f"Mesh sizing not computed for material '{mat.name}'. "
                    f"Call layers.compute_all_mesh_sizing() first."
                )

    # ------------------------------------------------------------------
    # Main generation method
    # ------------------------------------------------------------------

    def generate(self) -> None:
        """
        Run the full input generation pipeline:
          1. Validate
          2. Ensure output directories exist
          3. Assemble .poly geometry (nodes → regions → facets)
          4. Write .poly file
          5. Write FEM input files
        """
        self.validate()
        self.io.ensure_dirs()

        # ── Assemble .poly file ──────────────────────────────────────
        assembler = PolyAssembler(
            domain=self.domain,
            layers=self.layers,
            source=self.source,
            receivers=self.receivers,
            anomalies=self.anomalies,
            pml=self.pml,
        )
        assembler.evaluate_all_input_data()
        assembler.write(self.io.poly_file)
        # print(f"Written: {self.io.poly_file}")

        # ── Write FEM input files ────────────────────────────────────
        writer = FEMInputWriter(
            domain=self.domain,
            layers=self.layers,
            source=self.source,
            receivers=self.receivers,
            solver_cfg=self.solver,
            pml=self.pml,
            output_dir=self.io.input_dir,
            model_file_path=self.io.model_file_stem,
            output_E_file=self.io.output_E_file,
            output_H_file=self.io.output_H_file,
        )
        writer.write_all(regions=assembler.regions)
        # print(f"Written: {self.io.input_dir / 'elfe3D_input.txt'}")
        # print(f"Written: {self.io.input_dir / 'source.txt'}")
        # print(f"Written: {self.io.input_dir / 'regionparameters.txt'}")
        print("Input generation complete.")
