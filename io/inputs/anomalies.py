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


@dataclass
class SphereAnomaly:
    """A spherical anomaly defined by center, radius, and material properties."""

    center: tuple[float, float, float]
    radius: float
    properties: tuple[float, float, float, float]
    subdivision_levels: int = 4

    def __post_init__(self) -> None:
        if len(self.center) != 3:
            raise ValueError(f"SphereAnomaly center must be 3D, got {self.center}")
        if self.radius <= 0:
            raise ValueError(f"SphereAnomaly radius must be > 0, got {self.radius}")
        if len(self.properties) != 4:
            raise ValueError(
                f"SphereAnomaly properties must be (eps_r, sigma, mu_r, sigma_m), got {self.properties}"
            )

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

    @property
    def centroid(self) -> tuple[float, float, float]:
        return self.center

    @property
    def z_midpoint(self) -> float:
        return self.center[2]

    def compute_max_element_volume(self, f_low: float) -> str:
        C = 3e8
        wavelength = C / (f_low * np.sqrt(self.eps_r))
        edge = wavelength / 20.0
        vol = edge ** 3 / (6.0 * np.sqrt(2.0))
        return f"{vol:.4e}"

    def surface_mesh(self) -> tuple[np.ndarray, np.ndarray]:
        """Return sphere surface vertices and triangular faces (0-based indices)."""
        # Icosahedron base
        t = (1.0 + np.sqrt(5.0)) / 2.0
        raw_vertices = np.array([
            (-1,  t,  0), ( 1,  t,  0), (-1, -t,  0), ( 1, -t,  0),
            ( 0, -1,  t), ( 0,  1,  t), ( 0, -1, -t), ( 0,  1, -t),
            ( t,  0, -1), ( t,  0,  1), (-t,  0, -1), (-t,  0,  1),
        ], dtype=float)

        faces = np.array([
            (0,11,5),(0,5,1),(0,1,7),(0,7,10),(0,10,11),
            (1,5,9),(5,11,4),(11,10,2),(10,7,6),(7,1,8),
            (3,9,4),(3,4,2),(3,2,6),(3,6,8),(3,8,9),
            (4,9,5),(2,4,11),(6,2,10),(8,6,7),(9,8,1),
        ], dtype=int)

        vertices = []
        for v in raw_vertices:
            norm = v / np.linalg.norm(v)
            vertices.append(norm)
        vertices = np.array(vertices, dtype=float)

        midpoint_cache = {}

        def get_midpoint(i1: int, i2: int) -> int:
            nonlocal vertices
            key = tuple(sorted((i1, i2)))
            if key in midpoint_cache:
                return midpoint_cache[key]
            v1 = vertices[i1]
            v2 = vertices[i2]
            midpoint = (v1 + v2) / 2.0
            midpoint = midpoint / np.linalg.norm(midpoint)
            vertices = np.vstack([vertices, midpoint])
            idx = len(vertices) - 1
            midpoint_cache[key] = idx
            return idx

        for _ in range(self.subdivision_levels):
            new_faces = []
            midpoint_cache = {}
            for tri in faces:
                v1, v2, v3 = tri
                a = get_midpoint(v1, v2)
                b = get_midpoint(v2, v3)
                c = get_midpoint(v3, v1)
                new_faces.extend([(v1, a, c), (v2, b, a), (v3, c, b), (a, b, c)])
            faces = np.array(new_faces, dtype=int)

        vertices = vertices * self.radius + np.array(self.center, dtype=float)
        return vertices, faces

    def __str__(self) -> str:
        return (
            f"SphereAnomaly(center={self.center}, radius={self.radius}, "
            f"eps_r={self.eps_r}, sigma={self.sigma:.2e})"
        )
