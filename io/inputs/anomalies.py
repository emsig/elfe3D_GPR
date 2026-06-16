"""
anomalies.py
------------
Anomalous body definitions embedded in the simulation domain.

Supports:
  BoxAnomaly    : rectangular prismatic body defined by x/y/z extents
  SphereAnomaly : spherical body defined by centre + radius, tessellated
                  as a subdivided icosphere

Both classes expose the same interface:
  .eps_r, .sigma, .mu_r, .sigma_m, .rho
  .centroid                     → (x, y, z) region seed point
  .compute_max_element_volume() → TetGen volume constraint string

PolyAssembler therefore handles a heterogeneous list of anomalies with
isinstance() dispatch only where the geometry differs (nodes / facets).
"""

import numpy as np
from dataclasses import dataclass


# =============================================================================
# Icosphere helper
# =============================================================================

def _build_icosphere(
    radius: float,
    center: np.ndarray,
    levels: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Build a subdivided icosphere and return world-space geometry.

    Parameters
    ----------
    radius  : sphere radius [m]
    center  : (3,) array — sphere centre in world coordinates
    levels  : subdivision iterations  (4 → ~2 500 triangles / ~1 300 verts)

    Returns
    -------
    vertices : (N, 3) float64 — world-space vertex positions
    faces    : (M, 3) int64  — 0-indexed vertex index triples
    """
    t = (1.0 + np.sqrt(5.0)) / 2.0
    raw = np.array([
        (-1,  t,  0), ( 1,  t,  0), (-1, -t,  0), ( 1, -t,  0),
        ( 0, -1,  t), ( 0,  1,  t), ( 0, -1, -t), ( 0,  1, -t),
        ( t,  0, -1), ( t,  0,  1), (-t,  0, -1), (-t,  0,  1),
    ], dtype=np.float64)

    # Project onto unit sphere, then scale and translate
    norms = np.linalg.norm(raw, axis=1, keepdims=True)
    verts = list((raw / norms) * radius + center)

    faces = [
        (0,11,5),(0,5,1),(0,1,7),(0,7,10),(0,10,11),
        (1,5,9),(5,11,4),(11,10,2),(10,7,6),(7,1,8),
        (3,9,4),(3,4,2),(3,2,6),(3,6,8),(3,8,9),
        (4,9,5),(2,4,11),(6,2,10),(8,6,7),(9,8,1),
    ]

    midpoint_cache: dict[tuple[int, int], int] = {}

    def _midpoint(i1: int, i2: int) -> int:
        key = (min(i1, i2), max(i1, i2))
        if key in midpoint_cache:
            return midpoint_cache[key]
        mid = (np.array(verts[i1]) + np.array(verts[i2])) / 2.0
        # Re-project back onto sphere surface
        local = mid - center
        local = local / np.linalg.norm(local) * radius
        new_v = local + center
        idx = len(verts)
        verts.append(new_v)
        midpoint_cache[key] = idx
        return idx

    for _ in range(levels):
        new_faces = []
        midpoint_cache.clear()
        for v1, v2, v3 in faces:
            a = _midpoint(v1, v2)
            b = _midpoint(v2, v3)
            c = _midpoint(v3, v1)
            new_faces += [(v1, a, c), (v2, b, a), (v3, c, b), (a, b, c)]
        faces = new_faces

    return np.array(verts, dtype=np.float64), np.array(faces, dtype=np.int64)


# =============================================================================
# BoxAnomaly
# =============================================================================

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


# =============================================================================
# SphereAnomaly
# =============================================================================

@dataclass
class SphereAnomaly:
    """
    A spherical anomaly embedded in the model, tessellated as an icosphere.

    Parameters
    ----------
    center            : (cx, cy, cz) — sphere centre [m]
    radius            : sphere radius [m]
    properties        : (eps_r, sigma, mu_r, sigma_m) — material properties
    subdivision_levels: icosphere subdivision depth
                        (default 4 → ~2 500 triangles, ~1 300 vertices)

    The icosphere mesh is generated once at construction time.
    PolyAssembler reads .vertices (N×3) and .faces (M×3) to emit TetGen
    nodes and triangular facets directly.
    """

    center: tuple[float, float, float]
    radius: float
    properties: tuple[float, float, float, float]
    subdivision_levels: int = 4

    def __post_init__(self) -> None:
        if self.radius <= 0:
            raise ValueError(f"SphereAnomaly radius must be > 0, got {self.radius}")
        if len(self.properties) != 4:
            raise ValueError(
                "SphereAnomaly properties must be (eps_r, sigma, mu_r, sigma_m), "
                f"got {self.properties}"
            )
        self._vertices, self._faces = _build_icosphere(
            self.radius,
            np.array(self.center, dtype=np.float64),
            self.subdivision_levels,
        )

    # ------------------------------------------------------------------
    # Mesh accessors
    # ------------------------------------------------------------------

    @property
    def vertices(self) -> np.ndarray:
        """(N, 3) float64 — icosphere vertex positions in world space."""
        return self._vertices

    @property
    def faces(self) -> np.ndarray:
        """(M, 3) int64 — 0-indexed vertex index triples (triangles)."""
        return self._faces

    # ------------------------------------------------------------------
    # Material property accessors  (same interface as BoxAnomaly)
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
    # Geometry helpers  (same interface as BoxAnomaly)
    # ------------------------------------------------------------------

    @property
    def centroid(self) -> tuple[float, float, float]:
        """Sphere centre — used as TetGen region seed."""
        return (float(self.center[0]), float(self.center[1]), float(self.center[2]))

    def compute_max_element_volume(self, f_low: float) -> str:
        """
        Max element volume for TetGen: same formula as BoxAnomaly —
        wavelength in medium at f_low divided by 20, cubed, / 6√2.
        """
        C = 3e8
        wavelength = C / (f_low * np.sqrt(self.eps_r))
        edge = wavelength / 20.0
        vol = edge ** 3 / (6.0 * np.sqrt(2.0))
        return f"{vol:.4e}"

    def __str__(self) -> str:
        return (
            f"SphereAnomaly("
            f"center={self.center}, radius={self.radius}, "
            f"eps_r={self.eps_r}, sigma={self.sigma:.2e}, "
            f"verts={len(self._vertices)}, faces={len(self._faces)})"
        )
