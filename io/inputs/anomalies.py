"""
anomalies.py
------------
Anomalous body definitions embedded in the simulation domain.

Currently supports BoxAnomaly — a rectangular prismatic body defined by
its x, y, z coordinate ranges and a 4-tuple of material properties.

The anomaly geometry feeds directly into the TetGen .poly file as an
additional set of nodes + facets + a region marker.
"""

import numpy as np
from dataclasses import dataclass


@dataclass
class BoxAnomaly:
    """
    A rectangular box-shaped anomaly embedded in the model.

    Parameters
    ----------
    x : (x_min, x_max) — x-extent of the anomaly [m]
    y : (y_min, y_max) — y-extent of the anomaly [m]
    z : (z_min, z_max) — z-extent of the anomaly [m]
    properties : (eps_r, sigma, mu_r, sigma_m) — material properties tuple,
                 matching the original anomaly_properties convention
    """

    x: tuple[float, float]
    y: tuple[float, float]
    z: tuple[float, float]
    properties: tuple[float, float, float, float]

    def __post_init__(self) -> None:
        for label, extent in [("x", self.x), ("y", self.y), ("z", self.z)]:
            if len(extent) != 2:
                raise ValueError(f"Anomaly {label} must be (min, max), got {extent}")
            if extent[0] >= extent[1]:
                raise ValueError(
                    f"Anomaly {label}[0] must be < {label}[1], got {extent}"
                )

    # ------------------------------------------------------------------
    # Named extent accessors (match original anomaly['x'][0] style)
    # ------------------------------------------------------------------

    @property
    def x_min(self) -> float:
        return self.x[0]

    @property
    def x_max(self) -> float:
        return self.x[1]

    @property
    def y_min(self) -> float:
        return self.y[0]

    @property
    def y_max(self) -> float:
        return self.y[1]

    @property
    def z_min(self) -> float:
        return self.z[0]

    @property
    def z_max(self) -> float:
        return self.z[1]

    # ------------------------------------------------------------------
    # Material property accessors (match original a_mat indexing)
    # ------------------------------------------------------------------

    @property
    def eps_r(self) -> float:
        return self.properties[0]

    @property
    def sigma(self) -> float:
        return self.properties[1]

    @property
    def mu_r(self) -> float:
        return self.properties[2]

    @property
    def sigma_m(self) -> float:
        return self.properties[3]

    @property
    def rho(self) -> float:
        return 1.0 / self.sigma

    # ------------------------------------------------------------------
    # Geometry helpers
    # ------------------------------------------------------------------

    @property
    def centroid(self) -> tuple[float, float, float]:
        """Centre point of the anomaly box — used as TetGen region marker."""
        return (
            (self.x_min + self.x_max) / 2.0,
            (self.y_min + self.y_max) / 2.0,
            (self.z_min + self.z_max) / 2.0,
        )

    @property
    def z_midpoint(self) -> float:
        """Z-midpoint, used as the region seed z-coordinate in original code."""
        return (self.z_min + self.z_max) / 2.0

    def compute_max_element_volume(self, f_low: float) -> str:
        """
        Compute the max element volume for the anomaly region.
        Matches original: wavelength in anomaly medium at f_low, divided by 20.
        Returns formatted string for TetGen region field.
        """
        C = 3e8
        wavelength = C / (f_low * np.sqrt(self.eps_r))
        edge = wavelength / 20.0
        vol = edge ** 3 / (6.0 * np.sqrt(2.0))
        return f"{vol:.4e}"

    def __str__(self) -> str:
        return (
            f"BoxAnomaly("
            f"x={self.x}, y={self.y}, z={self.z}, "
            f"eps_r={self.eps_r}, sigma={self.sigma:.2e})"
        )
