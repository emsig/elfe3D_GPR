"""
receivers.py
------------
Receiver antenna array definition and tetrahedra node geometry.

ReceiverArray computes inline, endfire, and oblique receiver positions
along with the three-node equilateral triangle tetrahedra vertices used
to embed receivers in the TetGen mesh.
"""

import numpy as np
from dataclasses import dataclass, field


@dataclass
class ReceiverArray:
    """
    Array of GPR receiver antennas in three orientations.

    Parameters (user-specified)
    ---------------------------
    num_inline   : number of inline receivers (along x-axis)
    num_endfire  : number of endfire receivers (along y-axis)
    num_oblique  : number of oblique receivers (45-degree, combined x+y)
    source_length : antenna length [m] — receivers match source length
    pml_type      : "log" or "lin" — affects PML buffer for placement
    pml_thickness : single-layer PML thickness [m]
    num_pml_layers : number of PML layers

    Derived fields (computed in __post_init__)
    -----------------------------------------
    depth   : receiver depth below surface [m] (= -length / 4)
    length  : same as source_length (receivers match)
    height  : equilateral triangle height for tetrahedra = length * sqrt(3) / 2

    Receiver positions (x, y, z) and tetrahedra node coordinates
    (x_tet, x_0_tet, x_1_tet, etc.) are computed and stored as flat
    combined lists across all three orientations, matching the original.
    """

    num_inline: int
    num_endfire: int
    num_oblique: int
    source_length: float
    pml_type: str
    pml_thickness: float
    num_pml_layers: int

    # ------------------------------------------------------------------
    # Derived geometry fields
    # ------------------------------------------------------------------
    depth: float = field(default=0.0, init=False)
    length: float = field(default=0.0, init=False)
    height: float = field(default=0.0, init=False)

    # Combined receiver positions (inline + endfire + oblique)
    x: list[float] = field(default_factory=list, init=False)
    y: list[float] = field(default_factory=list, init=False)
    z: list[float] = field(default_factory=list, init=False)

    # Combined tetrahedra node positions
    x_tet:   list[float] = field(default_factory=list, init=False)
    x_0_tet: list[float] = field(default_factory=list, init=False)
    x_1_tet: list[float] = field(default_factory=list, init=False)
    y_tet:   list[float] = field(default_factory=list, init=False)
    y_0_tet: list[float] = field(default_factory=list, init=False)
    y_1_tet: list[float] = field(default_factory=list, init=False)
    z_tet:   list[float] = field(default_factory=list, init=False)

    def __post_init__(self) -> None:
        self.length = self.source_length
        self.depth = -1.0 * self.length / 4.0
        self.height = (self.length * np.sqrt(3.0)) / 2.0

        inline   = self._compute_inline()
        endfire  = self._compute_endfire()
        oblique  = self._compute_oblique()

        # Combine all orientations into flat lists (matches original)
        self.x = inline["x"] + endfire["x"] + oblique["x"]
        self.y = inline["y"] + endfire["y"] + oblique["y"]
        self.z = inline["z"] + endfire["z"] + oblique["z"]

        self.x_tet   = inline["x_tet"]   + endfire["x_tet"]   + oblique["x_tet"]
        self.x_0_tet = inline["x_0_tet"] + endfire["x_0_tet"] + oblique["x_0_tet"]
        self.x_1_tet = inline["x_1_tet"] + endfire["x_1_tet"] + oblique["x_1_tet"]
        self.y_tet   = inline["y_tet"]   + endfire["y_tet"]   + oblique["y_tet"]
        self.y_0_tet = inline["y_0_tet"] + endfire["y_0_tet"] + oblique["y_0_tet"]
        self.y_1_tet = inline["y_1_tet"] + endfire["y_1_tet"] + oblique["y_1_tet"]
        self.z_tet   = inline["z_tet"]   + endfire["z_tet"]   + oblique["z_tet"]

    @property
    def num_receivers(self) -> int:
        return len(self.x)

    # ------------------------------------------------------------------
    # PML buffer (used for placement bounds)
    # ------------------------------------------------------------------

    @property
    def _pml_buffer(self) -> float:
        if self.pml_type == "log":
            return self.pml_thickness
        return self.pml_thickness * self.num_pml_layers

    # ------------------------------------------------------------------
    # Per-orientation position computation
    # ------------------------------------------------------------------

    def _compute_inline(self) -> dict:
        """Inline receivers along the x-axis, evenly spaced from 0.1 to 1.0."""
        n = self.num_inline
        if n == 0:
            return {k: [] for k in ["x","y","z","x_tet","x_0_tet","x_1_tet",
                                     "y_tet","y_0_tet","y_1_tet","z_tet"]}
        x_pos = [round(0.1 + i * (1.0 - 0.1) / (n - 1), 5) for i in range(n)]
        y_pos = [0.0] * n
        z_pos = [self.depth] * n

        return {
            "x": x_pos, "y": y_pos, "z": z_pos,
            "x_tet":   [x - self.length / 2.0 for x in x_pos],
            "x_0_tet": x_pos,
            "x_1_tet": [x + self.length / 2.0 for x in x_pos],
            "y_tet":   [-self.height / 2.0] * n,
            "y_0_tet": [ self.height / 2.0] * n,
            "y_1_tet": [-self.height / 2.0] * n,
            "z_tet":   [0.0] * n,
        }

    def _compute_endfire(self) -> dict:
        """Endfire receivers along the y-axis, evenly spaced from 0.1 to 1.0."""
        n = self.num_endfire
        if n == 0:
            return {k: [] for k in ["x","y","z","x_tet","x_0_tet","x_1_tet",
                                     "y_tet","y_0_tet","y_1_tet","z_tet"]}
        y_pos = [round(0.1 + i * (1.0 - 0.1) / (n - 1), 5) for i in range(n)]
        x_pos = [0.0] * n
        z_pos = [self.depth] * n

        return {
            "x": x_pos, "y": y_pos, "z": z_pos,
            "x_tet":   [ self.length / 2.0] * n,
            "x_0_tet": [-self.length / 2.0] * n,
            "x_1_tet": [0.0] * n,
            "y_tet":   [y - self.height / 2.0 for y in y_pos],
            "y_0_tet": [y - self.height / 2.0 for y in y_pos],
            "y_1_tet": [y + self.height / 2.0 for y in y_pos],
            "z_tet":   [0.0] * n,
        }

    def _compute_oblique(self) -> dict:
        """Oblique receivers at 45 degrees from origin, matching inline x-spacing."""
        n = self.num_oblique
        if n == 0:
            return {k: [] for k in ["x","y","z","x_tet","x_0_tet","x_1_tet",
                                     "y_tet","y_0_tet","y_1_tet","z_tet"]}

        # Oblique positions derived from inline x spacing, rotated 45 degrees
        inline_x = [round(0.1 + i * (1.0 - 0.1) / (n - 1), 5) for i in range(n)]
        x_pos = [round(x * np.cos(np.pi / 4.0), 5) for x in inline_x]
        y_pos = [round(x * np.sin(np.pi / 4.0), 5) for x in inline_x]
        z_pos = [self.depth] * n

        return {
            "x": x_pos, "y": y_pos, "z": z_pos,
            "x_tet":   [x - self.height / 2.0 for x in x_pos],
            "x_0_tet": [x - self.height / 2.0 for x in x_pos],
            "x_1_tet": [x + self.height / 2.0 for x in x_pos],
            "y_tet":   [y + self.length / 2.0 for y in y_pos],
            "y_0_tet": [y - self.length / 2.0 for y in y_pos],
            "y_1_tet": y_pos,
            "z_tet":   [0.0] * n,
        }

    # ------------------------------------------------------------------
    # Display
    # ------------------------------------------------------------------

    def __str__(self) -> str:
        return (
            f"ReceiverArray("
            f"inline={self.num_inline}, "
            f"endfire={self.num_endfire}, "
            f"oblique={self.num_oblique}, "
            f"depth={self.depth*100:.2f} cm)"
        )
