"""
solver.py
---------
FEM solver configuration for elfe3D.

SolverConfig holds all parameters that control the adaptive FEM solver:
refinement strategy, convergence tolerances, output flags, and PML settings.
These map directly to the elfe3D_input.txt file fields.
"""

from dataclasses import dataclass
from enum import IntEnum


class SolverType(IntEnum):
    ITERATIVE = 1
    DIRECT    = 2


class ErrorEstMethod(IntEnum):
    ZIENKIEWICZ_ZHU = 1
    RESIDUAL        = 2
    RECOVERY        = 3
    IMPLICIT        = 4


class RefStrategy(IntEnum):
    UNIFORM    = 0
    ADAPTIVE   = 1


class PMLDecayType(IntEnum):
    SIMPLE_RECIPROCAL = 1
    LOGARITHMIC       = 2
    EXPONENTIAL       = 3
    POLYNOMIAL        = 4


@dataclass
class SolverConfig:
    """
    elfe3D FEM solver parameters.

    Parameters
    ----------
    solver           : solver type (1=iterative, 2=direct)
    max_ref_steps    : maximum adaptive refinement steps (0 = no refinement)
    max_unknowns     : maximum number of DOFs before stopping refinement
    beta_ref         : refinement fraction (fraction of elements to refine)
    accuracy_tol     : convergence tolerance for iterative solver
    vtk              : write VTK output during refinement (0/1)
    error_est_method : error estimator selection
    ref_strategy     : refinement strategy (0=uniform, 1=adaptive)
    output_fields_vtk: write final field output as VTK (0/1)
    pml_decay_type   : PML conductivity profile type
    """

    solver: int = SolverType.DIRECT
    max_ref_steps: int = 0
    max_unknowns: int = 5_000_000
    beta_ref: float = 0.85
    accuracy_tol: float = 3e-5
    vtk: int = 0
    error_est_method: int = ErrorEstMethod.RECOVERY
    ref_strategy: int = RefStrategy.ADAPTIVE
    output_fields_vtk: int = 1
    pml_decay_type: int = PMLDecayType.SIMPLE_RECIPROCAL

    def __post_init__(self) -> None:
        if not (0.0 < self.beta_ref <= 1.0):
            raise ValueError(f"beta_ref must be in (0, 1], got {self.beta_ref}")
        if self.accuracy_tol <= 0:
            raise ValueError(f"accuracy_tol must be positive, got {self.accuracy_tol}")
        if self.max_unknowns <= 0:
            raise ValueError(f"max_unknowns must be positive, got {self.max_unknowns}")

    def __str__(self) -> str:
        return (
            f"SolverConfig("
            f"type={self.solver}, "
            f"max_ref={self.max_ref_steps}, "
            f"tol={self.accuracy_tol:.1e})"
        )
