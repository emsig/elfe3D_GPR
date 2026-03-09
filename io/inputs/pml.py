"""
pml.py
------
Perfectly Matched Layer (PML) configuration.

PMLConfig holds the numerical PML parameters. The geometry of the PML
shells (nodes, facets, regions) is computed by the assembler in
writetetgenpoly.py, which uses these parameters directly.

The PML material properties are inherited from the LayerStack — each PML
region uses the same electromagnetic properties as the domain region at
the same z-depth. This is handled during region assembly.
"""

from dataclasses import dataclass
from enum import IntEnum


class PMLType(IntEnum):
    """Spatial distribution of PML layers."""
    LINEAR      = 0   # "lin" — uniform thickness per layer
    LOGARITHMIC = 1   # "log" — single thick slab


class PMLTheory(IntEnum):
    """PML formulation used in the FEM."""
    ANISOTROPIC          = 0
    DIFFERENTIAL_OPERATOR = 1


@dataclass
class PMLConfig:
    """
    Perfectly Matched Layer parameters.

    Parameters
    ----------
    num_layers      : number of PML shells (each adds one prism-width outward)
    layer_thickness : thickness of one PML shell [m]
    pml_type        : "lin" or "log" — spatial layout of shells
    pml_theory      : 0=anisotropic, 1=differential operator formulation
    pml_decay_type  : conductivity profile (passed to SolverConfig too)
                      1=simple reciprocal, 2=log, 3=exp, 4=polynomial
    """

    num_layers: int = 1
    layer_thickness: float = 0.3
    pml_type: str = "lin"
    pml_theory: int = PMLTheory.ANISOTROPIC
    pml_decay_type: int = 1

    def __post_init__(self) -> None:
        if self.num_layers < 0:
            raise ValueError(f"num_layers must be >= 0, got {self.num_layers}")
        if self.layer_thickness <= 0:
            raise ValueError(
                f"layer_thickness must be positive, got {self.layer_thickness}"
            )
        if self.pml_type not in ("lin", "log"):
            raise ValueError(f"pml_type must be 'lin' or 'log', got '{self.pml_type}'")

    @property
    def present(self) -> bool:
        """True if PML is active (num_layers > 0)."""
        return self.num_layers > 0

    @property
    def total_thickness(self) -> float:
        """
        Total PML thickness outward from domain boundary [m].

        For 'lin': num_layers * layer_thickness (each layer adds one shell)
        For 'log': layer_thickness (single slab regardless of num_layers)
        """
        if self.pml_type == "log":
            return self.layer_thickness
        return self.num_layers * self.layer_thickness

    def __str__(self) -> str:
        return (
            f"PMLConfig("
            f"n={self.num_layers}, "
            f"thickness={self.layer_thickness} m, "
            f"type='{self.pml_type}')"
        )
