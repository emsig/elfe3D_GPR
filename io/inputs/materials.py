"""
materials.py
------------
Electromagnetic material properties for a single simulation region.

Each material holds its constitutive parameters. Resistivity (rho) is a
computed property. Mesh sizing (max_element_edge, max_element_volume) is
derived from frequency and must be set by calling compute_mesh_sizing()
once the frequency range is known.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import ClassVar


@dataclass
class Material:
    """
    Electromagnetic material for one simulation region.

    Parameters
    ----------
    name        : human-readable label (e.g. "air", "clay", "anomaly")
    eps_r       : relative permittivity [-]
    sigma       : electrical conductivity [S/m]
    mu_r        : relative permeability [-], default 1.0
    sigma_m     : magnetic conductivity [Ohm/m], default 0.0

    Mesh sizing fields (set by compute_mesh_sizing):
    max_element_edge   : maximum element edge length [m]
    max_element_volume : maximum element volume [m^3], stored as formatted string
                         to match TetGen .poly region field format
    """

    name: str
    eps_r: float
    sigma: float
    mu_r: float = 1.0
    sigma_m: float = 0.0

    max_element_edge: float | None = field(default=None, repr=False)
    max_element_volume: str | None = field(default=None, repr=False)

    _C: ClassVar[float] = 3e8  # speed of light [m/s]

    def __post_init__(self) -> None:
        if self.eps_r <= 0:
            raise ValueError(f"eps_r must be positive, got {self.eps_r} for '{self.name}'")
        if self.sigma <= 0:
            raise ValueError(f"sigma must be positive, got {self.sigma} for '{self.name}'")
        if self.mu_r <= 0:
            raise ValueError(f"mu_r must be positive, got {self.mu_r} for '{self.name}'")

    # ------------------------------------------------------------------
    # Derived electromagnetic properties
    # ------------------------------------------------------------------

    @property
    def rho(self) -> float:
        """Electrical resistivity [Ohm.m] = 1 / sigma."""
        return 1.0 / self.sigma

    def wavelength_at(self, frequency: float) -> float:
        """Free-space wavelength in this medium at a given frequency [m]."""
        return self._C / (frequency * np.sqrt(self.eps_r))

    # ------------------------------------------------------------------
    # Mesh sizing
    # ------------------------------------------------------------------

    def compute_mesh_sizing(self, f_low: float, samples_per_wavelength: int) -> None:
        """
        Compute max element edge and volume from the lowest simulation frequency.

        Uses the largest wavelength (lowest frequency) in this medium to set
        a mesh density of `samples_per_wavelength` elements per wavelength.
        Volume formula: edge^3 / (6 * sqrt(2)) : regular tetrahedron inscribed
        in a cube of that edge length.

        Must be called before the material is used in assembly.
        """
        largest_wavelength = self.wavelength_at(f_low)
        self.max_element_edge = largest_wavelength / samples_per_wavelength
        vol = self.max_element_edge ** 3 / (6.0 * np.sqrt(2.0))
        self.max_element_volume = f"{vol:.4e}"

    # ------------------------------------------------------------------
    # Display
    # ------------------------------------------------------------------

    def __str__(self) -> str:
        return (
            f"Material('{self.name}', eps_r={self.eps_r}, "
            f"sigma={self.sigma:.2e} S/m, mu_r={self.mu_r})"
        )
