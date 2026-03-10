"""
sources.py
----------
Source antenna definition and derived geometry.

SourceAntenna holds the user-specified source parameters and computes
all derived quantities used by the assembler:
  - frequency list
  - physical antenna extents (x/y/z)
  - discretised segment positions
  - source refinement box extents

All computations preserve the original logic exactly; they are just
moved from the monolithic __init__ of elfe3DGPRTestDesign into
__post_init__ here.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import ClassVar


@dataclass
class SourceAntenna:
    """
    Dipole source antenna for elfe3D GPR.

    Parameters (user-specified)
    ---------------------------
    source_type        : source type integer (6 = Hertzian dipole)
    antenna_position   : [x, y, z] position of the antenna centre [m]
    current_direction  : 1=x, 2=y, 3=z
    current            : source current [A]
    source_moment      : source moment flag
    num_segments       : number of discretisation segments along antenna
    ricker_central_f   : central Ricker wavelet frequency [Hz]
    num_points_per_range : number of frequency points in the sweep
    m                  : mesh refinement factor around source box
    box_present        : whether a source refinement box is included
    s_f                : antenna length divisor (length = wavelength / s_f)
    bh_f               : box half-size factor (box half-extent = wavelength / (4*bh_f))
    box_x_extents      : user-specified x-extents for source box [m]
    num_layers         : number of earth layers (affects which wavelength
                         is used for antenna sizing)
    largest_wavelengths: largest wavelength per material region at f_low
                         (set from LayerStack after mesh sizing is computed)
    """

    # Required — source physics
    source_type: int
    antenna_position: list[float]
    current_direction: int
    current: float
    source_moment: int
    num_segments: int
    ricker_central_f: float

    # Required — frequency sweep
    num_points_per_range: int

    # Required — sizing
    m: int
    box_present: bool
    s_f: float
    bh_f: float
    box_x_extents: list[float]  # user-specified x range for refinement box

    # Required — context from LayerStack
    num_layers: int
    largest_wavelengths: list[float]  # one per region, computed externally

    # ------------------------------------------------------------------
    # Derived fields (computed in __post_init__, not set by user)
    # ------------------------------------------------------------------
    f_list: list[float] = field(default_factory=list, init=False)
    length: float = field(default=0.0, init=False)
    x_extents: list[float] = field(default_factory=list, init=False)
    y_extents: list[float] = field(default_factory=list, init=False)
    z_extents: list[float] = field(default_factory=list, init=False)
    x_disc: list[float] = field(default_factory=list, init=False)
    y_disc: list[float] = field(default_factory=list, init=False)
    z_disc: list[float] = field(default_factory=list, init=False)
    box_x: list[list[float]] = field(default_factory=list, init=False)
    box_y: list[list[float]] = field(default_factory=list, init=False)
    box_z: list[list[float]] = field(default_factory=list, init=False)

    def __post_init__(self) -> None:
        self._compute_frequency_list()
        self._compute_antenna_geometry()
        if self.box_present:
            self._compute_source_box()

    # ------------------------------------------------------------------
    # Frequency list — matches original logic exactly
    # ------------------------------------------------------------------

    def _compute_frequency_list(self) -> None:
        if self.num_points_per_range > 1:
            self.f_list = [
                i * self.ricker_central_f * (6.0 / self.num_points_per_range)
                for i in range(1, self.num_points_per_range + 1)
            ]
        else:
            self.f_list = [self.ricker_central_f]

    @property
    def f_low(self) -> float:
        return self.f_list[0]

    @property
    def f_high(self) -> float:
        return self.f_list[-1]

    # ------------------------------------------------------------------
    # Antenna physical geometry — matches original logic exactly
    # ------------------------------------------------------------------

    def _compute_antenna_geometry(self) -> None:
        # Antenna length: sized by wavelength of air (num_layers==0) or
        # first earth layer (num_layers>0), divided by s_f
        ref_wavelength = (
            self.largest_wavelengths[1]
            if self.num_layers > 0
            else self.largest_wavelengths[0]
        )
        self.length = ref_wavelength / self.s_f

        px, py, pz = self.antenna_position
        half = self.length / 2.0

        self.x_extents = [px - half, px + half]
        self.y_extents = [py, py]
        self.z_extents = [pz, pz]

        # Discretised segment positions along antenna
        self.x_disc = list(np.linspace(self.x_extents[0], self.x_extents[1], self.num_segments + 1))
        self.y_disc = list(np.linspace(self.y_extents[0], self.y_extents[1], self.num_segments + 1))
        self.z_disc = list(np.linspace(self.z_extents[0], self.z_extents[1], self.num_segments + 1))

    # ------------------------------------------------------------------
    # Source refinement box — matches original logic exactly
    # ------------------------------------------------------------------

    def _compute_source_box(self) -> None:
        # TODO: It might be better to remove it (as based on previous tests 
        # it didn't work well)
        px, py, pz = self.antenna_position
        wl_air = self.largest_wavelengths[0]

        # Air-side box (above z=0)
        self.box_x = [[px - wl_air / (4.0 * self.bh_f),
                        self.box_x_extents[0] if self.num_layers == 0 else self.box_x_extents[0]]]
        self.box_y = [[py - wl_air / (4.0 * self.bh_f),
                        py + wl_air / (4.0 * self.bh_f)]]
        self.box_z = [[pz + wl_air / (4.0 * self.bh_f), 0.0]]

        # Earth-side box (below z=0), only if there are layers
        if self.num_layers > 0:
            wl_earth = self.largest_wavelengths[1]
            self.box_x.append([px - wl_earth / (4.0 * self.bh_f),
                                self.box_x_extents[1]])
            self.box_y.append([py - wl_earth / (4.0 * self.bh_f),
                                py + wl_earth / (4.0 * self.bh_f)])
            self.box_z.append([0.0, pz - wl_earth / (4.0 * self.bh_f)])
        else:
            wl_air = self.largest_wavelengths[0]
            self.box_x.append([px - wl_air / (4.0 * self.bh_f),
                                self.box_x_extents[1] if isinstance(self.box_x_extents, list) else self.box_x_extents])
            self.box_y.append([py - wl_air / (4.0 * self.bh_f),
                                py + wl_air / (4.0 * self.bh_f)])
            self.box_z.append([0.0, pz - wl_air / (4.0 * self.bh_f)])

    # ------------------------------------------------------------------
    # Display
    # ------------------------------------------------------------------

    def __str__(self) -> str:
        return (
            f"SourceAntenna("
            f"type={self.source_type}, "
            f"pos={self.antenna_position}, "
            f"f_central={self.ricker_central_f/1e6:.1f} MHz, "
            f"length={self.length*1e3:.2f} mm)"
        )
