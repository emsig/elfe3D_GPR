"""
domain.py
---------
Model domain (the physical simulation box, excluding PML).

ModelDomain holds the x/y/z extents of the FEM domain that are specified
by the user (x_e, y_e, z_e). All other geometry objects (layers, PML)
reference these bounds.
"""

from dataclasses import dataclass
from pathlib import Path


@dataclass
class ModelDomain:
    """
    Axis-aligned bounding box of the FEM simulation domain, excluding PML.

    Parameters
    ----------
    x_e : [x_min, x_max]  — domain x-extent [m]
    y_e : [y_min, y_max]  — domain y-extent [m]
    z_e : [z_min, z_max]  — domain z-extent [m], z_min < 0 (into earth),
                            z_max > 0 (into air)
    num_air_wavelengths : number of air wavelengths above surface (informational)
    """

    x_e: list[float]
    y_e: list[float]
    z_e: list[float]
    num_air_wavelengths: float = 0.0

    def __post_init__(self) -> None:
        for label, ext in [("x_e", self.x_e), ("y_e", self.y_e), ("z_e", self.z_e)]:
            if len(ext) != 2:
                raise ValueError(f"{label} must be [min, max], got {ext}")
            if ext[0] >= ext[1]:
                raise ValueError(f"{label}[0] must be < {label}[1], got {ext}")

    # ------------------------------------------------------------------
    # Named extent accessors
    # ------------------------------------------------------------------

    @property
    def x_min(self) -> float:
        return self.x_e[0]

    @property
    def x_max(self) -> float:
        return self.x_e[1]

    @property
    def y_min(self) -> float:
        return self.y_e[0]

    @property
    def y_max(self) -> float:
        return self.y_e[1]

    @property
    def z_min(self) -> float:
        return self.z_e[0]

    @property
    def z_max(self) -> float:
        return self.z_e[1]

    @property
    def x_size(self) -> float:
        return self.x_max - self.x_min

    @property
    def y_size(self) -> float:
        return self.y_max - self.y_min

    @property
    def z_size(self) -> float:
        return self.z_max - self.z_min

    def contains(self, x: float, y: float, z: float) -> bool:
        """Return True if the point (x, y, z) is inside the domain."""
        return (
            self.x_min <= x <= self.x_max
            and self.y_min <= y <= self.y_max
            and self.z_min <= z <= self.z_max
        )

    def __str__(self) -> str:
        return (
            f"ModelDomain("
            f"x=[{self.x_min}, {self.x_max}], "
            f"y=[{self.y_min}, {self.y_max}], "
            f"z=[{self.z_min}, {self.z_max}])"
        )
