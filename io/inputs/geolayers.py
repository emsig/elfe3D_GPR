"""
geolayers.py
------------
Geological layer definitions and the layer stack that manages them.

GeoLayer    : a single horizontal layer with thickness and material.
LayerStack  : ordered collection (air + earth layers), computes z-interfaces.

The z-interface computation is based on the following logic:
    z[0] = 0.0  (air-earth surface)
    z[i] = -cumulative sum of thicknesses up to layer i
"""

import numpy as np
from dataclasses import dataclass, field
from .materials import Material


@dataclass
class GeoLayer:
    """
    A single horizontal geological layer.

    Parameters
    ----------
    thickness : layer thickness [m], positive downward
    material  : electromagnetic properties of this layer
    """

    thickness: float
    material: Material

    def __post_init__(self) -> None:
        if self.thickness <= 0:
            raise ValueError(
                f"Layer thickness must be positive, got {self.thickness} m "
                f"for material '{self.material.name}'"
            )

    def __str__(self) -> str:
        return f"GeoLayer({self.thickness} m | {self.material.name})"


@dataclass
class LayerStack:
    """
    Ordered stack of geological layers, surface downward.

    The stack always includes an air half-space above the first earth layer.
    Index 0 in `all_materials` is air; indices 1..n are earth layers.

    Parameters
    ----------
    air    : material for the air half-space above z = 0
    layers : list of GeoLayer objects in order from surface downward

    The x/y extents of layer interfaces are inherited from the ModelDomain
    at assembly time — they are not stored here.
    """

    air: Material
    layers: list[GeoLayer] = field(default_factory=list)

    # ------------------------------------------------------------------
    # Basic counts
    # ------------------------------------------------------------------

    @property
    def num_layers(self) -> int:
        """Number of earth layers (not counting air)."""
        return len(self.layers)

    @property
    def num_regions(self) -> int:
        """Total number of material regions: air + all earth layers."""
        return self.num_layers + 1

    # ------------------------------------------------------------------
    # Material access
    # ------------------------------------------------------------------

    @property
    def all_materials(self) -> list[Material]:
        """
        All materials in region order: air first, then earth layers.
        Index matches the region index used throughout the assembly.
        """
        return [self.air] + [layer.material for layer in self.layers]

    # ------------------------------------------------------------------
    # Z-interface positions
    # ------------------------------------------------------------------

    @property
    def z_interfaces(self) -> list[float]:
        """
        Z-coordinates of all layer interfaces, including the air-earth
        boundary at z = 0.

        Logic:
            z[0] = 0.0
            z[i] = -sum(thicknesses[0:i])   for i = 1 .. num_layers-1

        The last layer is assumed to extend to domain z_min (no bottom
        interface stored — the domain boundary closes it).

        Returns
        -------
        list of float, length = num_layers
            (One interface per layer boundary, starting at the surface.)
        """
        if self.num_layers == 0:
            return []

        thicknesses = [layer.thickness for layer in self.layers]
        interfaces = [0.0]
        for i in range(1, self.num_layers):
            depth = float(np.sum(thicknesses[:i], dtype=np.float64)) * -1.0
            interfaces.append(depth)
        return interfaces

    # ------------------------------------------------------------------
    # Mesh sizing delegation
    # ------------------------------------------------------------------

    def compute_all_mesh_sizing(
        self, f_low: float, samples_per_wavelength: int
    ) -> None:
        """
        Compute mesh sizing for all materials in this stack.
        Call once after the frequency range is finalised.
        """
        for material in self.all_materials:
            material.compute_mesh_sizing(f_low, samples_per_wavelength)

    # ------------------------------------------------------------------
    # Convenience accessors (match original regions dict structure)
    # ------------------------------------------------------------------

    @property
    def eps_r_list(self) -> list[float]:
        return [m.eps_r for m in self.all_materials]

    @property
    def sigma_list(self) -> list[float]:
        return [m.sigma for m in self.all_materials]

    @property
    def mu_r_list(self) -> list[float]:
        return [m.mu_r for m in self.all_materials]

    @property
    def sigma_m_list(self) -> list[float]:
        return [m.sigma_m for m in self.all_materials]

    @property
    def rho_list(self) -> list[float]:
        return [m.rho for m in self.all_materials]

    @property
    def max_element_edges(self) -> list[float]:
        return [m.max_element_edge for m in self.all_materials]

    @property
    def max_element_volumes(self) -> list[str]:
        return [m.max_element_volume for m in self.all_materials]

    # ------------------------------------------------------------------
    # Display
    # ------------------------------------------------------------------

    def __str__(self) -> str:
        lines = [f"LayerStack ({self.num_layers} earth layers):"]
        lines.append(f"  Air: {self.air}")
        for i, layer in enumerate(self.layers):
            lines.append(f"  Layer {i+1}: {layer}")
        return "\n".join(lines)
