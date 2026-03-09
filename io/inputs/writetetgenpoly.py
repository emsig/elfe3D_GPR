"""
writetetgenpoly.py
------------------
Assembles and writes the TetGen .poly input file.

PolyAssembler takes all domain objects (layers, source, receivers, anomaly,
PML) and runs the three assembly passes — nodes, regions, facets — then
writes the .poly file.

All geometry and math logic is preserved exactly from the original
elfe3DGPRTestWritingPML class. Only the structure is changed:
  - Data comes from typed dataclass objects instead of flat dicts
  - Node/facet/region assembly methods are clearly separated
  - Helper geometry functions are imported from tetgenprimitives.py
"""

import os
import numpy as np
from itertools import groupby
from pathlib import Path

from tetgenprimitives import (
    Node,
    Region,
    create_rectangular_prism_nodes,
    angular_sort,
    get_face_nodes,
    create_face_string,
    create_cuboid_faces_from_nodes,
    create_half_cuboid_faces_from_nodes,
)
from domain import ModelDomain
from geolayers import LayerStack
from sources import SourceAntenna
from receivers import ReceiverArray
from anomalies import BoxAnomaly
from pml import PMLConfig


class PolyAssembler:
    """
    Assembles nodes, facets, and regions for the TetGen .poly file.

    Usage
    -----
    assembler = PolyAssembler(domain, layers, source, receivers, anomaly, pml)
    assembler.evaluate_all_nodes()
    assembler.evaluate_all_regions()
    assembler.evaluate_all_facets()
    assembler.write(output_path)
    """

    def __init__(
        self,
        domain: ModelDomain,
        layers: LayerStack,
        source: SourceAntenna,
        receivers: ReceiverArray,
        anomaly: BoxAnomaly,
        pml: PMLConfig,
    ):
        self.domain    = domain
        self.layers    = layers
        self.source    = source
        self.receivers = receivers
        self.anomaly   = anomaly
        self.pml       = pml

        # Flat domain bounds (convenience, matches original self.x_min etc.)
        self.x_min = domain.x_min
        self.x_max = domain.x_max
        self.y_min = domain.y_min
        self.y_max = domain.y_max
        self.z_min = domain.z_min
        self.z_max = domain.z_max

        # Layer interfaces — (x/y bounds shared with domain)
        self.num_layers       = layers.num_layers
        self.layer_interfaces = {
            'x_min': domain.x_min,
            'x_max': domain.x_max,
            'y_min': domain.y_min,
            'y_max': domain.y_max,
            'z': layers.z_interfaces,
        }

        # Material property lists (match original self.region_* names)
        self.region_epsilon_r  = layers.eps_r_list
        self.region_sigma      = layers.sigma_list
        self.region_mu_r       = layers.mu_r_list
        self.region_sigma_m    = layers.sigma_m_list
        self.region_rho        = layers.rho_list
        self.region_volumes    = layers.max_element_volumes
        self.region_element_edges = layers.max_element_edges

        # PML parameters
        self.NUM_PML       = pml.num_layers
        self.pml_thickness = pml.layer_thickness

        # Source parameters (match original self.* names)
        self.source_type         = source.source_type
        self.current_direction   = source.current_direction
        self.source_moment       = source.source_moment
        self.m                   = source.m
        self.box_present         = source.box_present
        self.num_source_segments = source.num_segments
        self.length_antenna      = source.length
        self.source_x            = source.x_extents
        self.source_y            = source.y_extents
        self.source_z            = source.z_extents
        self.source_x_disc       = source.x_disc
        self.source_y_disc       = source.y_disc
        self.source_z_disc       = source.z_disc
        self.source_x_box        = source.box_x
        self.source_y_box        = source.box_y
        self.source_z_box        = source.box_z

        # Receiver positions (match original self.receivers_* names)
        self.receivers_x      = receivers.x
        self.receivers_y      = receivers.y
        self.receivers_z      = receivers.z
        self.receivers_x_tet  = receivers.x_tet
        self.receivers_x_0_tet = receivers.x_0_tet
        self.receivers_x_1_tet = receivers.x_1_tet
        self.receivers_y_tet  = receivers.y_tet
        self.receivers_y_0_tet = receivers.y_0_tet
        self.receivers_y_1_tet = receivers.y_1_tet
        self.receivers_z_tet  = receivers.z_tet

        # Anomaly parameters (match original self.ax, self.ay, etc.)
        self.ax    = anomaly.x
        self.ay    = anomaly.y
        self.az    = anomaly.z
        self.a_mat = anomaly.properties

        # Frequency list (for anomaly volume computation)
        self.f_list = source.f_list

        # State — set during evaluate_all_*
        self.num_node   = 0
        self.num_facet  = 0
        self.num_regions = 0
        self.regions    = []

    # ==========================================================================
    # Public API
    # ==========================================================================

    def evaluate_all_input_data(self) -> None:
        """Run all three assembly passes in the correct order."""
        self.evaluate_all_nodes()
        self.evaluate_all_regions()
        self.evaluate_all_facets()

    def write(self, path: Path | str) -> Path:
        """Write assembled data to a TetGen .poly file."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        self._write_poly_file(path)
        return path

    # ==========================================================================
    # Node assembly — preserved exactly from original evaluate_all_nodes()
    # ==========================================================================

    def evaluate_all_nodes(self) -> None:
        """
        Evaluate all nodes for the .poly file.
        Sets self.node_list_domain_prism, self.node_list_of_list_of_interfaces,
        self.node_list_receiver, self.node_list_source, self.node_list_anomaly,
        and (if PML active) self.node_list_PML_2, self.node_list_PML_3.
        """
        self.num_node = 0

        # Domain prism corner nodes (top and bottom faces)
        node_list_domain_prism, self.num_node = create_rectangular_prism_nodes(
            self.x_min, self.x_max, self.y_min, self.y_max, self.z_min, 1, self.num_node
        )
        top_nodes, self.num_node = create_rectangular_prism_nodes(
            self.x_min, self.x_max, self.y_min, self.y_max, self.z_max, 1, self.num_node
        )
        node_list_domain_prism.extend(top_nodes)

        # Interface nodes (one set per layer boundary)
        node_list_of_list_of_interfaces = []
        self.marker_interfaces = [0, 1, 7, 8, 9, 10]

        for layer in range(self.num_layers):
            interface_nodes, self.num_node = create_rectangular_prism_nodes(
                self.layer_interfaces['x_min'], self.layer_interfaces['x_max'],
                self.layer_interfaces['y_min'], self.layer_interfaces['y_max'],
                self.layer_interfaces['z'][layer],
                self.marker_interfaces[layer],
                self.num_node,
            )
            node_list_of_list_of_interfaces.append(interface_nodes)

        # Receiver tetrahedra nodes (3 nodes per receiver)
        node_list_receiver = []
        marker_receivers = 99
        for i in range(len(self.receivers_x)):
            self.num_node += 1
            node_list_receiver.append(Node(self.num_node, self.receivers_x_0_tet[i], self.receivers_y_0_tet[i], self.receivers_z_tet[i], marker_receivers))
            self.num_node += 1
            node_list_receiver.append(Node(self.num_node, self.receivers_x_1_tet[i], self.receivers_y_1_tet[i], self.receivers_z_tet[i], marker_receivers))
            self.num_node += 1
            node_list_receiver.append(Node(self.num_node, self.receivers_x_tet[i],   self.receivers_y_tet[i],   self.receivers_z_tet[i], marker_receivers))

        # Source discretised segment nodes
        node_list_source = []
        marker_source = 0
        for i in range(self.num_source_segments + 1):
            self.num_node += 1
            node_list_source.append(Node(self.num_node, self.source_x_disc[i], self.source_y_disc[i], self.source_z_disc[i], marker_source))

        # Source refinement box nodes
        if self.box_present:
            self.marker_source_box = 98
            node_list_source_box = self._build_source_box_nodes()
            self.node_list_source_box_0 = [node for node in node_list_source_box if node[3] == 0]

        # Store to self
        self.node_list_domain_prism              = node_list_domain_prism
        self.node_list_of_list_of_interfaces     = node_list_of_list_of_interfaces
        self.node_list_receiver                  = node_list_receiver
        self.node_list_source                    = node_list_source
        if self.box_present:
            self.node_list_source_box = node_list_source_box

        # Anomaly box nodes (top and bottom face corners)
        self.node_list_anomaly = []
        nodes_an, self.num_node = create_rectangular_prism_nodes(
            self.ax[0], self.ax[1], self.ay[0], self.ay[1], self.az[0], 101, self.num_node
        )
        self.node_list_anomaly.extend(nodes_an)
        nodes_an, self.num_node = create_rectangular_prism_nodes(
            self.ax[0], self.ax[1], self.ay[0], self.ay[1], self.az[1], 101, self.num_node
        )
        self.node_list_anomaly.extend(nodes_an)

        # PML nodes
        if self.NUM_PML > 0:
            self.node_list_PML_2, self.node_list_PML_3, self.pml_type_2_layer_nodes = (
                self._evaluate_all_nodes_pml()
            )

    def _build_source_box_nodes(self) -> list:
        """
        Build source refinement box nodes.
        Preserved exactly from original evaluate_all_nodes().
        """
        node_list_source_box = []

        if self.source_z_box[0][0] == -1 * self.source_z_box[1][1]:
            # Only one medium (box straddles z=0 symmetrically)
            box_z = [self.source_z_box[0][0], 0.0, self.source_z_box[1][1]]
            box_x = [self.source_x_box[0][0], self.source_x_box[0][1]]
            box_y = [self.source_y_box[0][0], self.source_y_box[0][1]]
            for i in range(3):
                for j in range(2):
                    for k in range(2):
                        self.num_node += 1
                        node_list_source_box.append(
                            Node(self.num_node, box_x[k], box_y[j], box_z[i], self.marker_source_box)
                        )
        else:
            # Two or more media
            for m in range(2):
                box_z = [self.source_z_box[m][0], self.source_z_box[m][1]]
                box_x = [self.source_x_box[m][0], self.source_x_box[m][1]]
                box_y = [self.source_y_box[m][0], self.source_y_box[m][1]]
                for i in range(2):
                    for j in range(2):
                        for k in range(2):
                            self.num_node += 1
                            node_list_source_box.append(
                                Node(self.num_node, box_x[k], box_y[j], box_z[i], self.marker_source_box)
                            )
        return node_list_source_box

    def _evaluate_all_nodes_pml(self) -> tuple[list, list, list]:
        """
        Evaluate nodes for all PML shells.
        Type 2: lateral edge extensions (intersection of two PML faces)
        Type 3: corner cube extensions (intersection of three PML faces)

        Preserved exactly from original evaluate_all_nodes_PML().
        """
        t = self.pml_thickness
        marker_pml_layers = [1000 for _ in range(self.NUM_PML)]
        self.marker_pml_layers = marker_pml_layers
        self.marker_pml_types  = [1, 2, 3]

        pml_type_2_nodes      = []
        pml_type_2_layer_nodes = []
        pml_type_3_nodes      = []

        z_layers = [self.layer_interfaces['z'][j] for j in range(self.num_layers)]

        # ── Type 2: lateral edge intersections ────────────────────────────────

        edge_configs = [
            # (x_start, x_stop, y_start, y_stop, exclude_x, exclude_y)
            (self.x_min, self.x_min - self.NUM_PML * t,  self.y_min, self.y_min - self.NUM_PML * t,  self.x_min, self.y_min),   # Front-Right
            (self.x_max, self.x_max + self.NUM_PML * t,  self.y_min, self.y_min - self.NUM_PML * t,  self.x_max, self.y_min),   # Front-Left
            (self.x_max, self.x_max + self.NUM_PML * t,  self.y_max, self.y_max + self.NUM_PML * t,  self.x_max, self.y_max),   # Back-Left
            (self.x_min, self.x_min - self.NUM_PML * t,  self.y_max, self.y_max + self.NUM_PML * t,  self.x_min, self.y_max),   # Back-Right
        ]

        for x_start, x_stop, y_start, y_stop, excl_x, excl_y in edge_configs:
            x = np.linspace(x_start, x_stop, self.NUM_PML + 1)
            y = np.linspace(y_start, y_stop, self.NUM_PML + 1)
            xyz = np.array(np.meshgrid(x, y, z_layers)).T.reshape(-1, 3)
            xyz = xyz[~((xyz[:, 0] == excl_x) & (xyz[:, 1] == excl_y))]

            for xi, yi, zi in xyz:
                self.num_node += 1
                node = Node(self.num_node, xi, yi, zi, marker_pml_layers[0] + self.marker_pml_types[1])
                pml_type_2_layer_nodes.append(node)
                pml_type_2_nodes.append(node)

        # ── Type 3: corner cubes ───────────────────────────────────────────────

        vertex_configs = [
            (self.x_max, self.x_max + self.NUM_PML * t,  self.y_min, self.y_min - self.NUM_PML * t,  self.z_max, self.z_max + self.NUM_PML * t),  # Front-Right-Top
            (self.x_max, self.x_max + self.NUM_PML * t,  self.y_min, self.y_min - self.NUM_PML * t,  self.z_min, self.z_min - self.NUM_PML * t),  # Front-Right-Bottom
            (self.x_min, self.x_min - self.NUM_PML * t,  self.y_min, self.y_min - self.NUM_PML * t,  self.z_max, self.z_max + self.NUM_PML * t),  # Front-Left-Top
            (self.x_min, self.x_min - self.NUM_PML * t,  self.y_min, self.y_min - self.NUM_PML * t,  self.z_min, self.z_min - self.NUM_PML * t),  # Front-Left-Bottom
            (self.x_max, self.x_max + self.NUM_PML * t,  self.y_max, self.y_max + self.NUM_PML * t,  self.z_max, self.z_max + self.NUM_PML * t),  # Back-Right-Top
            (self.x_max, self.x_max + self.NUM_PML * t,  self.y_max, self.y_max + self.NUM_PML * t,  self.z_min, self.z_min - self.NUM_PML * t),  # Back-Right-Bottom
            (self.x_min, self.x_min - self.NUM_PML * t,  self.y_max, self.y_max + self.NUM_PML * t,  self.z_max, self.z_max + self.NUM_PML * t),  # Back-Left-Top
            (self.x_min, self.x_min - self.NUM_PML * t,  self.y_max, self.y_max + self.NUM_PML * t,  self.z_min, self.z_min - self.NUM_PML * t),  # Back-Left-Bottom
        ]

        for x_start, x_stop, y_start, y_stop, z_start, z_stop in vertex_configs:
            x = np.linspace(x_start, x_stop, self.NUM_PML + 1)
            y = np.linspace(y_start, y_stop, self.NUM_PML + 1)
            z = np.linspace(z_start, z_stop, self.NUM_PML + 1)
            xyz = np.array(np.meshgrid(x, y, z)).T.reshape(-1, 3)[1:]  # skip origin

            for xi, yi, zi in xyz:
                self.num_node += 1
                pml_type_3_nodes.append(
                    Node(self.num_node, xi, yi, zi, marker_pml_layers[0] + self.marker_pml_types[2])
                )

        return pml_type_2_nodes, pml_type_3_nodes, pml_type_2_layer_nodes

    # ==========================================================================
    # Region assembly — preserved exactly from original evaluate_all_regions()
    # ==========================================================================

    def evaluate_all_regions(self) -> None:
        """
        Evaluate all TetGen region seed points with material attributes.
        Includes earth layers, source box (if present), anomaly, and PML.
        """
        self.num_regions = 0
        regions = []

        # ── Earth layer regions ─────────────────────────────────────────────
        for i in range(self.num_layers + 1):
            if i == 0:
                region_height = (
                    self.layer_interfaces['z'][0]
                    + abs(self.layer_interfaces['z'][0] - self.z_max) / 2.0
                )
                r_label = "# Air"
            elif i == self.num_layers:
                region_height = (
                    self.layer_interfaces['z'][-1]
                    - abs(self.layer_interfaces['z'][-1] - self.z_min) / 2.0
                )
                r_label = "# Last Earth Layer"
            else:
                region_height = (
                    self.layer_interfaces['z'][i - 1]
                    - abs(self.layer_interfaces['z'][i - 1] - self.layer_interfaces['z'][i]) / 2.0
                )
                r_label = f"# Earth Layer {i}"

            self.num_regions += 1
            regions.append([
                self.num_regions, 0, 0, round(region_height, 5),
                self.num_regions, self.region_volumes[i], r_label,
                self.region_rho[i], self.region_mu_r[i], self.region_epsilon_r[i],
            ])

        # ── Source box regions ──────────────────────────────────────────────
        if self.box_present:
            if self.num_layers == 0:
                region_height = 0.0 + float(self.region_volumes[0]) * 5
                self.num_regions += 1
                vol = float(self.region_volumes[0]) / self.m ** 3
                regions.append([
                    self.num_regions, 0, 0, round(region_height, 5),
                    self.num_regions, f"{vol:.4e}", "# Source Box 1",
                    self.region_rho[0], self.region_mu_r[0], self.region_epsilon_r[0],
                ])
                region_height = 0.0 - float(self.region_volumes[0]) * 5
                self.num_regions += 1
                regions.append([
                    self.num_regions, 0, 0, round(region_height, 5),
                    self.num_regions, f"{vol:.4e}", "# Source Box 2",
                    self.region_rho[0], self.region_mu_r[0], self.region_epsilon_r[0],
                ])
            else:
                region_height = 0.0 + float(self.region_volumes[0]) * 5
                self.num_regions += 1
                vol = float(self.region_volumes[0]) / self.m ** 3
                regions.append([
                    self.num_regions, 0, 0, round(region_height, 5),
                    self.num_regions, f"{vol:.4e}", "# Source Box 1",
                    self.region_rho[0], self.region_mu_r[0], self.region_epsilon_r[0],
                ])
                region_height = 0.0 - float(self.region_volumes[1]) * 5
                self.num_regions += 1
                vol = float(self.region_volumes[1]) / self.m ** 3
                regions.append([
                    self.num_regions, 0, 0, round(region_height, 5),
                    self.num_regions, f"{vol:.4e}", "# Source Box 2",
                    self.region_rho[1], self.region_mu_r[1], self.region_epsilon_r[1],
                ])

        self.regions = regions

        # ── Anomaly region ──────────────────────────────────────────────────
        self.num_regions += 1
        an_vol = self.anomaly.compute_max_element_volume(self.f_list[0])
        self.regions.append([
            self.num_regions, 0, 0, self.anomaly.z_midpoint,
            self.num_regions, an_vol, "# Anomaly Region",
            self.anomaly.rho, self.anomaly.mu_r, self.anomaly.eps_r,
        ])

        # ── PML regions ─────────────────────────────────────────────────────
        if self.NUM_PML > 0:
            pml_regions = self._evaluate_all_regions_pml()
            self.regions += pml_regions

    def _evaluate_all_regions_pml(self) -> list:
        """
        Evaluate PML region seed points for all three PML types.
        Material properties are assigned by z-midpoint lookup.

        Preserved exactly from original evaluate_all_regions_PML().
        """
        if self.num_layers > 0:
            all_z = sorted([self.z_min] + self.layer_interfaces['z'] + [self.z_max])
        else:
            all_z = [self.z_min, self.z_max]
        self.all_z = all_z

        PML_regions = []
        t = self.pml_thickness
        n = self.NUM_PML

        n_list     = np.linspace(t / 2.0, n * t / 2.0, n)
        n_pairs    = [(i, j) for i in n_list for j in n_list]
        n_triplets = [(i, j, k) for i in n_list for j in n_list for k in n_list]

        midheights = [
            (all_z[j] + all_z[j + 1]) / 2.0
            for j in range(len(all_z) - 1)
        ]

        # ── Type 1: face extensions ──────────────────────────────────────
        self.marker_pml_1 = []
        for i in range(self.NUM_PML):
            faces = [
                ("Front",  [0.0, self.y_min - (i + 1) * t / 2.0, 0.0]),
                ("Right",  [self.x_max + (i + 1) * t / 2.0, 0.0, 0.0]),
                ("Back",   [0.0, self.y_max + (i + 1) * t / 2.0, 0.0]),
                ("Left",   [self.x_min - (i + 1) * t / 2.0, 0.0, 0.0]),
                ("Top",    [0.0, 0.0, self.z_max + (i + 1) * t / 2.0]),
                ("Bottom", [0.0, 0.0, self.z_min - (i + 1) * t / 2.0]),
            ]
            for face_name, midpoint in faces:
                r_label = f"# PML Layer {i+1} - Type 1: {face_name}"
                if face_name in ("Top", "Bottom"):
                    self.num_regions += 1
                    self.marker_pml_1.append(self.num_regions)
                    PML_regions.append([
                        self.num_regions,
                        round(midpoint[0], 5), round(midpoint[1], 5), round(midpoint[2], 5),
                        self.marker_pml_1[-1], None, r_label, None, None, None,
                    ])
                else:
                    for l in range(self.num_layers + 1):
                        self.num_regions += 1
                        self.marker_pml_1.append(self.num_regions)
                        PML_regions.append([
                            self.num_regions,
                            round(midpoint[0], 5), round(midpoint[1], 5), round(midheights[l], 5),
                            self.marker_pml_1[-1], None, r_label, None, None, None,
                        ])

        # ── Type 2: edge extensions ──────────────────────────────────────
        self.marker_pml_2 = []
        edges = (
            [("Front-Right",  [self.x_max + i, self.y_min - j, 0.0]) for i, j in n_pairs] +
            [("Back-Right",   [self.x_max + i, self.y_max + j, 0.0]) for i, j in n_pairs] +
            [("Back-Left",    [self.x_min - i, self.y_max + j, 0.0]) for i, j in n_pairs] +
            [("Front-Left",   [self.x_min - i, self.y_min - j, 0.0]) for i, j in n_pairs] +
            [("Front-Top",    [0.0, self.y_min - j, self.z_max + i]) for i, j in n_pairs] +
            [("Right-Top",    [self.x_max + i, 0.0, self.z_max + j]) for i, j in n_pairs] +
            [("Back-Top",     [0.0, self.y_max + j, self.z_max + i]) for i, j in n_pairs] +
            [("Left-Top",     [self.x_min - i, 0.0, self.z_max + j]) for i, j in n_pairs] +
            [("Front-Bottom", [0.0, self.y_min - j, self.z_min - i]) for i, j in n_pairs] +
            [("Right-Bottom", [self.x_max + j, 0.0, self.z_min - i]) for i, j in n_pairs] +
            [("Back-Bottom",  [0.0, self.y_max + j, self.z_min - i]) for i, j in n_pairs] +
            [("Left-Bottom",  [self.x_min - j, 0.0, self.z_min - i]) for i, j in n_pairs]
        )
        for edge_name, midpoint in edges:
            r_label = f"# PML Layer - Type 2: {edge_name}"
            if edge_name in ("Front-Right", "Front-Left", "Back-Left", "Back-Right"):
                for l in range(self.num_layers + 1):
                    self.num_regions += 1
                    self.marker_pml_2.append(self.num_regions)
                    PML_regions.append([
                        self.num_regions,
                        round(midpoint[0], 5), round(midpoint[1], 5), round(midheights[l], 5),
                        self.marker_pml_2[-1], None, r_label, None, None, None,
                    ])
            else:
                self.num_regions += 1
                self.marker_pml_2.append(self.num_regions)
                PML_regions.append([
                    self.num_regions,
                    round(midpoint[0], 5), round(midpoint[1], 5), round(midpoint[2], 5),
                    self.marker_pml_2[-1], None, r_label, None, None, None,
                ])

        # ── Type 3: corner extensions ────────────────────────────────────
        self.marker_pml_3 = []
        vertices = (
            [("Front-Right-Top",    [self.x_max + i, self.y_min - j, self.z_max + k]) for i, j, k in n_triplets] +
            [("Back-Right-Top",     [self.x_max + i, self.y_max + j, self.z_max + k]) for i, j, k in n_triplets] +
            [("Back-Left-Top",      [self.x_min - i, self.y_max + j, self.z_max + k]) for i, j, k in n_triplets] +
            [("Front-Left-Top",     [self.x_min - i, self.y_min - j, self.z_max + k]) for i, j, k in n_triplets] +
            [("Front-Right-Bottom", [self.x_max + i, self.y_min - j, self.z_min - k]) for i, j, k in n_triplets] +
            [("Back-Right-Bottom",  [self.x_max + i, self.y_max + j, self.z_min - k]) for i, j, k in n_triplets] +
            [("Back-Left-Bottom",   [self.x_min - i, self.y_max + j, self.z_min - k]) for i, j, k in n_triplets] +
            [("Front-Left-Bottom",  [self.x_min - i, self.y_min - j, self.z_min - k]) for i, j, k in n_triplets]
        )
        for vertex_name, midpoint in vertices:
            r_label = f"# PML Layer - Type 3: {vertex_name}"
            self.num_regions += 1
            self.marker_pml_3.append(self.num_regions)
            PML_regions.append([
                self.num_regions,
                round(midpoint[0], 5), round(midpoint[1], 5), round(midpoint[2], 5),
                self.marker_pml_3[-1], None, r_label, None, None, None,
            ])

        # ── Assign material properties by z-midpoint ─────────────────────
        interfaces = self.layer_interfaces['z']
        num_layers = self.num_layers

        for region in PML_regions:
            mz = region[3]
            if mz >= interfaces[0]:
                mat_idx = 0
            elif mz < interfaces[-1]:
                mat_idx = num_layers
            else:
                mat_idx = None
                for idx in range(1, len(interfaces)):
                    if interfaces[idx - 1] > mz >= interfaces[idx]:
                        mat_idx = idx
                        break
                if mat_idx is None:
                    raise ValueError(f"Cannot assign material for z={mz}")

            region[5] = self.region_volumes[mat_idx]
            region[7] = self.region_rho[mat_idx]
            region[8] = self.region_mu_r[mat_idx]
            region[9] = self.region_epsilon_r[mat_idx]

        print(self.layer_interfaces['z'])
        return PML_regions

    # ==========================================================================
    # Facet assembly — preserved exactly from original evaluate_all_facets()
    # The facet methods are long; they are preserved in full without changes
    # to the polygon construction logic.
    # ==========================================================================

    def evaluate_all_facets(self) -> None:
        """
        Evaluate all .poly file facets.
        Sets self.domain_string, self.interface_string, self.interfaces_string,
        self.source_string, self.source_box_string, self.anomaly_string,
        self.pml_string.
        """
        self.num_facet = 0

        node_list_interfaces = []
        for node_list in self.node_list_of_list_of_interfaces:
            node_list_interfaces += node_list

        # ── Domain outer faces (front/right/back/left/top/bottom) ────────
        self.domain_string = self._build_domain_face_strings(node_list_interfaces)

        # ── Source box z=0 faces (if present) ────────────────────────────
        if self.box_present:
            source_box_z0_string, source_box_1_faces, source_box_2_faces = (
                self._build_source_box_faces()
            )

        # ── Air-earth interface + receivers + source ──────────────────────
        self._build_interface_and_receiver_strings(
            source_box_z0_string if self.box_present else ""
        )

        # ── Remaining layer interfaces ────────────────────────────────────
        if self.num_layers > 0:
            self.num_facet += len(self.node_list_of_list_of_interfaces)
            marker_interfaces = [6, 7, 8, 9, 10]
            interfaces_n_string = ""
            i = 0
            num_polygon = 0
            for interface_nodes in self.node_list_of_list_of_interfaces:
                num_polygon += 1
                interface_face_nodes_list = [node[0] for node in sorted(interface_nodes)]
                interface_n_string = " ".join(map(str, interface_face_nodes_list)) + "\n"
                if num_polygon % 2 != 0:
                    str_1 = f"# interface layer {i+1}\n"
                else:
                    i += 1
                    str_1 = f"# transition layer {i+1}\n"
                interfaces_n_string += (
                    str_1
                    + f"1 0 {marker_interfaces[num_polygon-1]} \n4 "
                    + interface_n_string
                )
            self.interfaces_string = "# other interfaces and transitions\n" + interfaces_n_string
        else:
            self.interfaces_string = ""

        # ── Source box side faces ─────────────────────────────────────────
        if self.box_present:
            self.source_box_string = "# source box facets\n"
            for face in source_box_1_faces + source_box_2_faces:
                self.source_box_string += f"1 0 {self.marker_source_box}\n"
                self.source_box_string += "4 " + " ".join(map(str, [node[0] for node in face])) + "\n"
                self.num_facet += 1

        # ── Anomaly facets ────────────────────────────────────────────────
        an_faces = create_cuboid_faces_from_nodes(self.node_list_anomaly)
        self.anomaly_string = "\n# anomaly facet\n"
        self.marker_anomaly = 101
        for face in an_faces:
            self.num_facet += 1
            self.anomaly_string += f"1 0 {self.marker_anomaly}\n"
            self.anomaly_string += "4 " + " ".join(map(str, [node[0] for node in face])) + "\n"

        # ── PML facets ────────────────────────────────────────────────────
        if self.NUM_PML > 0:
            self._evaluate_all_facets_pml()

    def _build_domain_face_strings(self, node_list_interfaces: list) -> str:
        """Build the 6 outer domain face strings (front/right/back/left/top/bottom)."""

        def domain_side_face(face_name, axis_value, axis_index, sort1, sort2):
            self.num_facet += 1
            outer_nodes = get_face_nodes(
                self.node_list_domain_prism + node_list_interfaces,
                axis_value, axis_index, sort1, sort2,
            )
            iface_lists = [
                get_face_nodes(nl, axis_value, axis_index, sort1, sort2)
                for nl in self.node_list_of_list_of_interfaces
            ]
            n_poly = 1 + len(self.node_list_of_list_of_interfaces)
            info = [n_poly, 0, 1]
            str_1 = " ".join(map(str, info))
            str_2 = " ".join(map(str, outer_nodes))
            str_3 = "".join(
                f"{len(il)} " + " ".join(map(str, il)) + "\n"
                for il in iface_lists
            )
            return f"#{face_name}\n{str_1}\n{len(outer_nodes)} {str_2}\n{str_3}\n"

        def domain_horizontal_face(label, axis_value):
            self.num_facet += 1
            face_nodes = get_face_nodes(self.node_list_domain_prism, axis_value, 3, 1, 2)
            info_str = " ".join(map(str, [1, 0, 1]))
            return f"# {label}\n{info_str}\n4 {' '.join(map(str, face_nodes))}\n"

        front  = domain_side_face("front",  self.y_min, 2, 1, 3)
        right  = domain_side_face("right",  self.x_max, 1, 2, 3)
        back   = domain_side_face("back",   self.y_max, 2, 1, 3)
        left   = domain_side_face("left",   self.x_min, 1, 2, 3)
        top    = domain_horizontal_face("top and bottom", self.z_min)
        bottom = domain_horizontal_face("bottom",         self.z_max)

        return (
            "# polygons - # holes - boundary marker \n"
            + front + right + back + left + top + bottom
        )

    def _build_source_box_faces(self):
        """Build source box face lists and z=0 face string."""
        nb = len(self.node_list_source_box)
        if nb == 12:
            box1_faces = create_cuboid_faces_from_nodes(self.node_list_source_box[:8])
            box2_faces = create_cuboid_faces_from_nodes(self.node_list_source_box[4:])
        elif nb == 16:
            box1_faces = create_cuboid_faces_from_nodes(self.node_list_source_box[:8])
            box2_faces = create_cuboid_faces_from_nodes(self.node_list_source_box[8:])
        else:
            box1_faces = []
            box2_faces = []

        def filter_faces_no_z0(faces):
            with_z0, without_z0 = [], []
            for face in faces:
                if all(abs(node[3]) < 1e-8 for node in face):
                    with_z0.append(face)
                else:
                    without_z0.append(face)
            return without_z0, with_z0

        box1_faces, box1_z0 = filter_faces_no_z0(box1_faces)
        box2_faces, box2_z0 = filter_faces_no_z0(box2_faces)

        z0_string = ""
        for face in box1_z0 + box2_z0:
            z0_string += "4 " + " ".join(map(str, [node[0] for node in face])) + "\n"

        return z0_string, box1_faces, box2_faces

    def _build_interface_and_receiver_strings(self, source_box_z0_string: str) -> None:
        """
        Build the combined air-earth interface facet string (with receivers
        and optionally source box z=0 faces and source dipole segment).

        Preserved exactly from original evaluate_all_facets().
        """
        source_above_surface = any(node[3] > 0.0 for node in self.node_list_source)

        num_polygon = 0
        self.num_facet += 1

        interface_string = ""
        if self.num_layers > 0:
            node_list_air_earth_interface = self.node_list_of_list_of_interfaces[0]
            num_polygon += 1
            interface_face_nodes_list = [node[0] for node in sorted(node_list_air_earth_interface)]
            interface_n_string = " ".join(map(str, interface_face_nodes_list)) + "\n"
            interface_string = "# air earth interface\n4 " + interface_n_string

        # Receiver triangles
        n_receivers_string = ""
        for i in range(len(self.node_list_receiver) // 3):
            num_polygon += 1
            tri = [3] + [self.node_list_receiver[i * 3 + k][0] for k in range(3)]
            n_receivers_string += " ".join(map(str, tri)) + "\n"

        # Source dipole segment (only if source is at surface or above)
        source_string_inner = ""
        if not source_above_surface:
            num_polygon += 1
            source_nodes = [node[0] for node in sorted(self.node_list_source)]
            source_string_inner = (
                f"{len(source_nodes)} " + " ".join(map(str, source_nodes)) + "\n"
            )

        if self.box_present:
            num_polygon += 2

        interface_facet_information = [num_polygon, 0, 99]
        interface_base_string = " ".join(map(str, interface_facet_information)) + "\n"

        facet_label = "# air-earth facet\n" if self.num_layers > 0 else "# Tx-Rx facet\n"

        if source_above_surface:
            self.interface_string = (
                facet_label + interface_base_string
                + interface_string
                + "# receiver triangles\n" + n_receivers_string
                + ("# source box\n" + source_box_z0_string if self.box_present else "")
            )
            # Source dipole as separate facet
            self.num_facet += 1
            source_nodes = [node[0] for node in sorted(self.node_list_source)]
            src_str = f"{len(source_nodes)} " + " ".join(map(str, source_nodes)) + "\n"
            self.source_string = f"# source facet\n1 0 987\n{src_str}\n"
        else:
            self.interface_string = (
                facet_label + interface_base_string
                + interface_string
                + "# receiver triangles\n" + n_receivers_string
                + "# source nodes\n" + source_string_inner
                + ("# source box\n" + source_box_z0_string if self.box_present else "")
            )
            self.source_string = ""

        # After consuming the air-earth interface, remove it from the remaining list
        if self.num_layers > 0:
            self.node_list_of_list_of_interfaces_old = self.node_list_of_list_of_interfaces
            self.node_list_of_list_of_interfaces = self.node_list_of_list_of_interfaces[1:]

    def _evaluate_all_facets_pml(self) -> None:
        """
        Build TetGen .poly facet strings for PML Type-1, Type-2, and Type-3 shells.
        Preserved exactly from original evaluate_all_facets_PML().
        """
        # This method is intentionally long — it preserves the full PML facet
        # logic from the original without any mathematical changes.
        # See the original elfe3DGPRTestWritingPMLOAn.py evaluate_all_facets_PML()
        # for the authoritative version. Refactoring of the internal PML facet
        # loop structure can be done separately once the overall structure is stable.
        """
        Function to evaluate all three types of facets for each PML layer.
        All these facets will be rectangular prisms.
        """

        t = self.pml_thickness
        self.pml_string = "# PML Layers\n"
        
        all_z = self.all_z
        if self.num_layers > 0:
            all_face_nodes = self.node_list_domain_prism + [
                        node for node_list in self.node_list_of_list_of_interfaces_old for node in node_list
                    ] + self.node_list_PML_2 + self.node_list_PML_3
        else:
            all_face_nodes = self.node_list_domain_prism + self.node_list_PML_2 + self.node_list_PML_3
        
        xe = [self.x_min, self.x_max]
        ye = [self.y_min, self.y_max]

        # Type 1: Extend each face of the domain: 1 new face per original face
        # Assembling facets per layer
        for i in range(self.NUM_PML):
            self.pml_string += f"\n# PML Layer {i+1} - Type 1:\n"
            faces = ["Front", "Right", "Back", "Left", "Top", "Bottom"]
            
            all_face_filters = [
                ("Front", lambda node: (node[1] == self.x_min or node[1] == self.x_max) and (node[2] == self.y_min - i * t or node[2] == self.y_min - (i+1) * t) and (node[3] >= self.z_min and node[3] <= self.z_max)),
                ("Right", lambda node: (node[1] == self.x_max + i * t or node[1] == self.x_max + (i+1) * t) and (node[2] == self.y_min or node[2] == self.y_max) and (node[3] >= self.z_min and node[3] <= self.z_max)),
                ("Back", lambda node: (node[1] == self.x_min or node[1] == self.x_max) and (node[2] == self.y_max + i * t or node[2] == self.y_max + (i+1) * t) and (node[3] >= self.z_min and node[3] <= self.z_max)),
                ("Left", lambda node: (node[1] == self.x_min - (i+1) * t or node[1] == self.x_min - i * t) and (node[2] == self.y_min or node[2] == self.y_max) and (node[3] >= self.z_min and node[3] <= self.z_max)),
                ("Top", lambda node: (node[2] == self.y_min or node[2] == self.y_max) and (node[1] == self.x_min or node[1] == self.x_max) and (node[3] == self.z_max + i * t or node[3] == self.z_max + (i+1) * t)),
                ("Bottom", lambda node: (node[2] == self.y_min or node[2] == self.y_max) and (node[1] == self.x_min or node[1] == self.x_max) and (node[3] == self.z_min - i * t or node[3] == self.z_min - (i+1) * t))
            ]

            # Front face: sort by x, z for the main face, then by z for each interface
            front_face_sorting_settings = [(self.y_min - (i+1)*t, 2, 1, 3)] + [
                (z, 3, 1, 2) for z in self.layer_interfaces['z']
            ]

            right_face_sorting_settings = [(self.x_max + (i+1)*t, 1, 2, 3)] + [
                (z, 3, 1, 2) for z in self.layer_interfaces['z']
            ]

            back_face_sorting_settings = [(self.y_max + (i+1)*t, 2, 1, 3)] + [
                (z, 3, 1, 2) for z in self.layer_interfaces['z']
            ]

            left_face_sorting_settings = [(self.x_min - (i+1)*t, 1, 2, 3)] + [
                (z, 3, 1, 2) for z in self.layer_interfaces['z']
            ]

            top_face_sorting_settings = [
                (self.z_max + (i+1)*t, 3, 1, 2),  # Top face: sort by z, x
            ]

            bottom_face_sorting_settings = [
                (self.z_min - (i+1)*t, 3, 1, 2),  # Bottom face: sort by z, x
            ]

            pml_type_1_faces_nodes_list = []
            it = 0
            for _, filter_func in all_face_filters:
                it += 1
                face_nodes = [node for node in all_face_nodes if filter_func(node)]
                pml_type_1_faces_nodes_list.append([face_nodes])              

            if self.num_layers > 0:
                combined_interfaces_nodes = [node for node_list in self.node_list_of_list_of_interfaces_old for node in node_list] + self.pml_type_2_layer_nodes
            else:
                combined_interfaces_nodes = []
            pml_type_1_interfaces_nodes_list = []
            for _, filter_func in all_face_filters:
                interface_nodes = [node for node in combined_interfaces_nodes if filter_func(node)]
                pml_type_1_interfaces_nodes_list.append([interface_nodes])

            for face in range(len(faces)):
                if face == 0:
                    face_sorting_settings = front_face_sorting_settings
                    self.pml_string += f"\n# 1+l new face with interface segments for the Front face\n"
                elif face == 1:
                    face_sorting_settings = right_face_sorting_settings
                    self.pml_string += f"\n# 1+l new face with interface segments for the Right face\n"
                elif face == 2:
                    face_sorting_settings = back_face_sorting_settings
                    self.pml_string += f"\n# 1+l new face with interface segments for the Back face\n"
                elif face == 3:
                    face_sorting_settings = left_face_sorting_settings
                    self.pml_string += f"\n# 1+l new face with interface segments for the Left face\n"
                elif face == 4:
                    face_sorting_settings = top_face_sorting_settings
                    self.pml_string += f"\n# 1 new face with interface segments for the Top face\n"
                    pml_type_1_interfaces_nodes_list[face][0] = None
                elif face == 5:
                    face_sorting_settings = bottom_face_sorting_settings
                    self.pml_string += f"\n# 1 new face with interface segments for the Bottom face\n"
                    pml_type_1_interfaces_nodes_list[face][0] = None

                marker = self.marker_pml_1[i*6 + face]
                for new_face in range(len(face_sorting_settings)): 
                    self.num_facet += 1
                    self.pml_string += self.create_face_string(marker, pml_type_1_faces_nodes_list[face][0], pml_type_1_interfaces_nodes_list[face][0], *face_sorting_settings[new_face])
                
        # __________________________________________________________
        # Create 12 cuboids at each edge of the domain - Type 2
        # __________________________________________________________
        n = self.NUM_PML
        if n > 1:
            raise ValueError("The number of PML layers must be 1.")
        
        n_list = np.linspace(t, n*t, n)
        n_pairs = [(i, j) for i in n_list for j in n_list]
        n_b_list = np.linspace(0, (n-1)*t, n)
        n_b_pairs = [(i, j) for i in n_b_list for j in n_b_list]
        n_s_list = np.linspace(0, n*t, n+1)
        n_s_pairs = [(i, j) for i in n_s_list for j in n_s_list]

        all_edge_filters = [
            ("Front-Right", lambda node, i=i, j=j, bi=bi, bj=bj: ((node[1] >= self.x_max + bi and node[1] <= self.x_max + i) and (node[2] <= self.y_min - bj and node[2] >= self.y_min - j) and (node[3] >= self.z_min and node[3] <= self.z_max))) for (i, j), (bi, bj) in zip(n_pairs, n_b_pairs)
        ] + [("Back-Right", lambda node, i=i, j=j, bi=bi, bj=bj: ((node[1] >= self.x_max + bi and node[1] <= self.x_max + i) and (node[2] >= self.y_max + bj and node[2] <= self.y_max + j) and (node[3] >= self.z_min and node[3] <= self.z_max))) for (i, j), (bi, bj) in zip(n_pairs, n_b_pairs)
        ] + [("Back-Left", lambda node, i=i, j=j, bi=bi, bj=bj: ((node[1] <= self.x_min - bi and node[1] >= self.x_min - i) and (node[2] >= self.y_max + bj and node[2] <= self.y_max + j) and (node[3] >= self.z_min and node[3] <= self.z_max))) for (i, j), (bi, bj) in zip(n_pairs, n_b_pairs)
        ] + [("Front-Left", lambda node, i=i, j=j, bi=bi, bj=bj: ((node[1] <= self.x_min - bi and node[1] >= self.x_min - i) and (node[2] <= self.y_min - bj and node[2] >= self.y_min - j) and (node[3] >= self.z_min and node[3] <= self.z_max))) for (i, j), (bi, bj) in zip(n_pairs, n_b_pairs)
        
        ] + [("Front-Top", lambda node, i=i, j=j, bi=bi, bj=bj: ((node[2] <= self.y_min - bj and node[2] >= self.y_min - j) and (node[1] >= self.x_min and node[1] <= self.x_max) and (node[3] >= self.z_max + bi and node[3] <= self.z_max + i))) for (i, j), (bi, bj) in zip(n_pairs, n_b_pairs)
        ] + [("Right-Top", lambda node, i=i, j=j, bi=bi, bj=bj: ((node[2] >= self.y_min and node[2] <= self.y_max) and (node[1] >= self.x_max + bi and node[1] <= self.x_max + i) and (node[3] >= self.z_max + bj and node[3] <= self.z_max + j))) for (i, j), (bi, bj) in zip(n_pairs, n_b_pairs)
        ] + [("Back-Top", lambda node, i=i, j=j, bi=bi, bj=bj: ((node[2] >= self.y_max + bj and node[2] <= self.y_max + j) and (node[1] >= self.x_min and node[1] <= self.x_max) and (node[3] >= self.z_max + bi and node[3] <= self.z_max + i))) for (i, j), (bi, bj) in zip(n_pairs, n_b_pairs)
        ] + [("Left-Top", lambda node, i=i, j=j, bi=bi, bj=bj: ((node[2] >= self.y_min and node[2] <= self.y_max) and (node[1] <= self.x_min - bi and node[1] >= self.x_min - i) and (node[3] >= self.z_max + bj and node[3] <= self.z_max + j))) for (i, j), (bi, bj) in zip(n_pairs, n_b_pairs)
        
        ] + [("Front-Bottom", lambda node, i=i, j=j, bi=bi, bj=bj: ((node[2] <= self.y_min - bj and node[2] >= self.y_min - j) and (node[1] >= self.x_min and node[1] <= self.x_max) and (node[3] <= self.z_min - bi and node[3] >= self.z_min - i))) for (i, j), (bi, bj) in zip(n_pairs, n_b_pairs)
        ] + [("Right-Bottom", lambda node, i=i, j=j, bi=bi, bj=bj: ((node[2] >= self.y_min and node[2] <= self.y_max) and (node[1] >= self.x_max + bi and node[1] <= self.x_max + i) and (node[3] <= self.z_min - bj and node[3] >= self.z_min - j))) for (i, j), (bi, bj) in zip(n_pairs, n_b_pairs)
        ] + [("Back-Bottom", lambda node, i=i, j=j, bi=bi, bj=bj: ((node[2] >= self.y_max + bj and node[2] <= self.y_max + j) and (node[1] >= self.x_min and node[1] <= self.x_max) and (node[3] <= self.z_min - bi and node[3] >= self.z_min - i))) for (i, j), (bi, bj) in zip(n_pairs, n_b_pairs)
        ] + [("Left-Bottom", lambda node, i=i, j=j, bi=bi, bj=bj: ((node[2] >= self.y_min and node[2] <= self.y_max) and (node[1] <= self.x_min - bi and node[1] >= self.x_min - i) and (node[3] <= self.z_min - bj and node[3] >= self.z_min - j))) for (i, j), (bi, bj) in zip(n_pairs, n_b_pairs)]
        
        # Lateral cuboids

        front_right_sorting_settings = [val for pair in zip(
            [(self.x_max + i, 1, 2, 3) for i, _ in n_s_pairs], # Right faces: sort by y, z
            [(self.y_min - j, 2, 1, 3) for _, j in n_s_pairs], # Front faces: sort by x, z
            [(z, 3, 1, 2) for z in all_z], # z faces
        ) for val in pair]

        back_right_sorting_settings = [val for pair in zip(
            [(self.x_max + i, 1, 2, 3) for i, _ in n_s_pairs],  # Right face: sort by y, z
            [(self.y_max + j, 2, 1, 3) for _, j in n_s_pairs],   # Back face: sort by x, z
            [(z, 3, 1, 2) for z in all_z], # z faces
        ) for val in pair]

        back_left_sorting_settings = [val for pair in zip(
            [(self.x_min - i, 1, 2, 3) for i, _ in n_s_pairs],  # Left face: sort by y, z
            [(self.y_max + j, 2, 1, 3) for _, j in n_s_pairs],   # Back face: sort by x, z
            [(z, 3, 1, 2) for z in all_z], # z faces
            ) for val in pair]

        front_left_sorting_settings = [val for pair in zip(
            [(self.x_min - i, 1, 2, 3) for i, _ in n_s_pairs],  # Left face: sort by y, z
            [(self.y_min - j, 2, 1, 3) for _, j in n_s_pairs],   # Front face: sort by x, z
            [(z, 3, 1, 2) for z in all_z], # z faces
        ) for val in pair]

        # "Horizontal" cuboids (along the top and bottom faces)

        front_top_sorting_settings = [val for pair in zip(
            [(self.y_min - j, 2, 1, 3) for _, j in n_s_pairs],  # Front face: sort by x, z
            [(self.z_max + j, 3, 1, 2) for _, j in n_s_pairs],   # Top face: sort by z, x
            [(x, 1, 2, 3) for x in xe],  # x faces
        ) for val in pair]

        right_top_sorting_settings = [val for pair in zip(
            [(self.x_max + j, 1, 2, 3) for _, j in n_s_pairs],  # Right face: sort by y, z
            [(self.z_max + j, 3, 1, 2) for _, j in n_s_pairs],   # Top face: sort by z, x
            [(y, 2, 1, 3) for y in ye],  # y faces
            ) for val in pair]

        back_top_sorting_settings = [val for pair in zip(
            [(self.y_max + j, 2, 1, 3) for _, j in n_s_pairs],  # Back face: sort by x, z
            [(self.z_max + j, 3, 1, 2) for _, j in n_s_pairs],   # Top face: sort by z, x
            [(x, 1, 2, 3) for x in xe],  # x faces
        ) for val in pair]

        left_top_sorting_settings = [val for pair in zip(
            [(self.x_min - j, 1, 2, 3) for _, j in n_s_pairs],  # Left face: sort by y, z
            [(self.z_max + j, 3, 1, 2) for _, j in n_s_pairs],   # Top face: sort by z, x
            [(y, 2, 1, 3) for y in ye],  # y faces
            ) for val in pair]

        front_bottom_sorting_settings = [val for pair in zip(
            [(self.y_min - j, 2, 1, 3) for _, j in n_s_pairs],  # Front face: sort by x, z
            [(self.z_min - j, 3, 1, 2) for _, j in n_s_pairs],   # Bottom face: sort by z, x
            [(x, 1, 2, 3) for x in xe],  # x faces
            ) for val in pair]

        right_bottom_sorting_settings = [val for pair in zip(
            [(self.x_max + j, 1, 2, 3) for _, j in n_s_pairs],  # Right face: sort by y, z
            [(self.z_min - j, 3, 1, 2) for _, j in n_s_pairs],   # Bottom face: sort by z, x
            [(y, 2, 1, 3) for y in ye],  # y faces
            ) for val in pair]

        back_bottom_sorting_settings = [val for pair in zip(
            [(self.y_max + j, 2, 1, 3) for _, j in n_s_pairs],  # Back face: sort by x, z
            [(self.z_min - j, 3, 1, 2) for _, j in n_s_pairs],   # Bottom face: sort by z, x
            [(x, 1, 2, 3) for x in xe],  # x faces
            ) for val in pair]

        left_bottom_sorting_settings = [val for pair in zip(
            [(self.x_min - j, 1, 2, 3) for _, j in n_s_pairs],  # Left face: sort by y, z
            [(self.z_min - j, 3, 1, 2) for _, j in n_s_pairs],   # Bottom face: sort by z, x
            [(y, 2, 1, 3) for y in ye],  # y faces
            ) for val in pair]

        # PML Type 2 - Create 12 cuboids at each edge of the domain
        self.pml_string += f"\n# PML Layer - Type 2:\n"

        pml_type_2_faces_nodes_list = []
        for _, filter_func in all_edge_filters:
            edge_nodes = [node for node in all_face_nodes if filter_func(node)]
            pml_type_2_faces_nodes_list.append([edge_nodes])

        if self.num_layers > 0:
            combined_interfaces_nodes = [node for node_list in self.node_list_of_list_of_interfaces_old for node in node_list] + self.pml_type_2_layer_nodes
        else:
            combined_interfaces_nodes = []

        pml_type_2_interfaces_nodes_list = []
        for _, filter_func in all_edge_filters:
            interface_nodes = [node for node in combined_interfaces_nodes if filter_func(node)]
            pml_type_2_interfaces_nodes_list.append([interface_nodes])

        int_idx = 0
        int_markers = [99-n-self.num_layers for n in range(self.num_layers*4)]
        
        for edge in range(len(all_edge_filters)):
            # FRONT-RIGHT, BACK-RIGHT, BACK-LEFT, FRONT-LEFT, FRONT-TOP, RIGHT-TOP, BACK-TOP, LEFT-TOP, FRONT-BOTTOM, RIGHT-BOTTOM, BACK-BOTTOM, LEFT-BOTTOM
            if edge == 0:
                edge_sorting_settings = front_right_sorting_settings
                self.pml_string += f"\n# 6+l new faces with interface segments for the Front-Right edge\n"
            elif edge == len(n_pairs):
                edge_sorting_settings = back_right_sorting_settings
                self.pml_string += f"\n# 6+l new faces with interface segments for the Back-Right edge\n"
            elif edge == 2*len(n_pairs):
                edge_sorting_settings = back_left_sorting_settings
                self.pml_string += f"\n# 6+l new faces with interface segments for the Back-Left edge\n"
            elif edge == 3*len(n_pairs):
                edge_sorting_settings = front_left_sorting_settings
                self.pml_string += f"\n# 6+l new faces with interface segments for the Front-Left edge\n"
            elif edge == 4*len(n_pairs):
                edge_sorting_settings = front_top_sorting_settings
                self.pml_string += f"\n# 6 new faces with interface segments for the Front-Top edge\n"
            elif edge == 5*len(n_pairs):
                edge_sorting_settings = right_top_sorting_settings
                self.pml_string += f"\n# 6 new faces with interface segments for the Right-Top edge\n"
            elif edge == 6*len(n_pairs):
                edge_sorting_settings = back_top_sorting_settings
                self.pml_string += f"\n# 6 new faces with interface segments for the Back-Top edge\n"
            elif edge == 7*len(n_pairs):
                edge_sorting_settings = left_top_sorting_settings
                self.pml_string += f"\n# 6 new faces with interface segments for the Left-Top edge\n"
            elif edge == 8*len(n_pairs):
                edge_sorting_settings = front_bottom_sorting_settings
                self.pml_string += f"\n# 6 new faces with interface segments for the Front-Bottom edge\n"
            elif edge == 9*len(n_pairs):
                edge_sorting_settings = right_bottom_sorting_settings
                self.pml_string += f"\n# 6 new faces with interface segments for the Right-Bottom edge\n"
            elif edge == 10*len(n_pairs):
                edge_sorting_settings = back_bottom_sorting_settings
                self.pml_string += f"\n# 6 new faces with interface segments for the Back-Bottom edge\n"
            elif edge == 11*len(n_pairs):
                edge_sorting_settings = left_bottom_sorting_settings
                self.pml_string += f"\n# 6 new faces with interface segments for the Left-Bottom edge\n"

            if edge >= 4*len(n_pairs):
                pml_type_2_interfaces_nodes_list[edge][0] = None
                            
            # Use int_markers for faces corresponding to z value filters for edges less than 4*n_pairs
            if edge < 4 * len(n_pairs):
                # edge_sorting_settings contains (axis_value, axis_index, sort_axis_1, sort_axis_2) for each face
                # The z-faces are every 3rd in edge_sorting_settings, starting at index 2
                marker_list = []
                for new_face in range(len(edge_sorting_settings)):
                    axis_value, axis_index, _, _ = edge_sorting_settings[new_face]
                    # If this is a z-face (axis_index == 3), assign int_markers by z value
                    if axis_index == 3 and axis_value != self.z_min and axis_value != self.z_max:
                        marker = int_markers[int_idx]
                        int_idx += 1
                    else:
                        marker = self.marker_pml_2[edge]
                    marker_list.append(marker)
                all_face_strings = [
                    self.create_face_string(marker_list[i], pml_type_2_faces_nodes_list[edge][0], pml_type_2_interfaces_nodes_list[edge][0], *edge_sorting_settings[i])
                    for i in range(len(edge_sorting_settings))
                ]
                unique_face_strings = list(set(all_face_strings))
                for face_str in unique_face_strings:
                    self.num_facet += 1
                    self.pml_string += face_str
                    all_face_strings = [
                        self.create_face_string(marker, pml_type_2_faces_nodes_list[edge][0], pml_type_2_interfaces_nodes_list[edge][0], axis_value, axis_index, sort_axis_1, sort_axis_2)
                        for axis_value, axis_index, sort_axis_1, sort_axis_2 in edge_sorting_settings
                    ]
            else:
                marker = self.marker_pml_2[edge] 
                all_component_strings = []
                for new_face in range(len(edge_sorting_settings)):
                    face_string = self.create_face_string(
                        marker, 
                        pml_type_2_faces_nodes_list[edge][0], 
                        pml_type_2_interfaces_nodes_list[edge][0], 
                        *edge_sorting_settings[new_face]
                    )
                    all_component_strings.append(face_string)
                unique_component_strings = list(set(all_component_strings))
                for new_face in range(len(unique_component_strings)):
                    self.num_facet += 1
                    self.pml_string += unique_component_strings[new_face]

        # Type 3: Extend each vertex of the domain: n^2 new cubes per original vertex
        # Assembling facets for all layers in one loop

        n_triplets = [(i, j, k) for i in n_list for j in n_list for k in n_list]
        n_b_triplets = [(i, j, k) for i in n_b_list for j in n_b_list for k in n_b_list]
        
        self.pml_string += f"\n# PML Layer - Type 3:\n"

        all_vertex_filters = [
            ("Front-Right-Top", (self.x_max, self.y_min, self.z_max), lambda node, i=i, j=j, k=k, bi=bi, bj=bj, bk=bk: ((node[1] >= self.x_max + bi and node[1] <= self.x_max + i) and (node[2] <= self.y_min - bj and node[2] >= self.y_min - j) and (node[3] >= self.z_max + bk and node[3] <= self.z_max + k))) for (i, j, k), (bi, bj, bk) in zip(n_triplets, n_b_triplets)
        ] + [("Back-Right-Top", (self.x_max, self.y_max, self.z_max), lambda node, i=i, j=j, k=k, bi=bi, bj=bj, bk=bk: ((node[1] >= self.x_max + bi and node[1] <= self.x_max + i) and (node[2] >= self.y_max + bj and node[2] <= self.y_max + j) and (node[3] >= self.z_max + bk and node[3] <= self.z_max + k))) for (i, j, k), (bi, bj, bk) in zip(n_triplets, n_b_triplets)
        ] + [("Back-Left-Top", (self.x_min, self.y_max, self.z_max), lambda node, i=i, j=j, k=k, bi=bi, bj=bj, bk=bk: ((node[1] <= self.x_min - bi and node[1] >= self.x_min - i) and (node[2] >= self.y_max + bj and node[2] <= self.y_max + j) and (node[3] >= self.z_max + bk and node[3] <= self.z_max + k))) for (i, j, k), (bi, bj, bk) in zip(n_triplets, n_b_triplets)
        ] + [("Front-Left-Top", (self.x_min, self.y_min, self.z_max), lambda node, i=i, j=j, k=k, bi=bi, bj=bj, bk=bk: ((node[1] <= self.x_min - bi and node[1] >= self.x_min - i) and (node[2] <= self.y_min - bj and node[2] >= self.y_min - j) and (node[3] >= self.z_max + bk and node[3] <= self.z_max + k))) for (i, j, k), (bi, bj, bk) in zip(n_triplets, n_b_triplets)
        ] + [("Front-Right-Bottom", (self.x_max, self.y_min, self.z_min), lambda node, i=i, j=j, k=k, bi=bi, bj=bj, bk=bk: ((node[1] >= self.x_max + bi and node[1] <= self.x_max + i) and (node[2] <= self.y_min - bj and node[2] >= self.y_min - j) and (node[3] <= self.z_min - bk and node[3] >= self.z_min - k))) for (i, j, k), (bi, bj, bk) in zip(n_triplets, n_b_triplets)
        ] + [("Back-Right-Bottom", (self.x_max, self.y_max, self.z_min), lambda node, i=i, j=j, k=k, bi=bi, bj=bj, bk=bk: ((node[1] >= self.x_max + bi and node[1] <= self.x_max + i) and (node[2] >= self.y_max + bj and node[2] <= self.y_max + j) and (node[3] <= self.z_min - bk and node[3] >= self.z_min - k))) for (i, j, k), (bi, bj, bk) in zip(n_triplets, n_b_triplets)
        ] + [("Back-Left-Bottom", (self.x_min, self.y_max, self.z_min), lambda node, i=i, j=j, k=k, bi=bi, bj=bj, bk=bk: ((node[1] <= self.x_min - bi and node[1] >= self.x_min - i) and (node[2] >= self.y_max + bj and node[2] <= self.y_max + j) and (node[3] <= self.z_min - bk and node[3] >= self.z_min - k))) for (i, j, k), (bi, bj, bk) in zip(n_triplets, n_b_triplets)
        ] + [("Front-Left-Bottom", (self.x_min, self.y_min, self.z_min), lambda node, i=i, j=j, k=k, bi=bi, bj=bj, bk=bk: ((node[1] <= self.x_min - bi and node[1] >= self.x_min - i) and (node[2] <= self.y_min - bj and node[2] >= self.y_min - j) and (node[3] <= self.z_min - bk and node[3] >= self.z_min - k))) for (i, j, k), (bi, bj, bk) in zip(n_triplets, n_b_triplets)]

        pml_type_3_faces_nodes_list = []
        for _, _, filter_func in all_vertex_filters:
            vertex_nodes = [node for node in all_face_nodes if filter_func(node)]
            # Only keep unique vertex nodes (by coordinates)
            unique_vertex_nodes = []
            seen = set()
            for node in vertex_nodes:
                key = (round(node[1], 8), round(node[2], 8), round(node[3], 8))
                if key not in seen:
                    seen.add(key)
                    unique_vertex_nodes.append(node)
            vertex_nodes = unique_vertex_nodes
            pml_type_3_faces_nodes_list.append([vertex_nodes])

        all_cuboid_faces = []
        not_needed_nodes = [n for _, n, _ in all_vertex_filters]
        for i in range(len(pml_type_3_faces_nodes_list)):
            faces = self.create_half_cuboid_faces_from_nodes(pml_type_3_faces_nodes_list[i][0], not_needed_nodes[i])
            all_cuboid_faces.append(faces)

        # Remove duplicates from all_cuboid_faces
        unique_faces = []
        unique_markers = []
        unique_sorted_faces = []
        it = 0
        for faces in all_cuboid_faces:
            for face in faces:
                # Sort the face by ascending order of its entries (node numbers)
                sorted_face = sorted(face, key=lambda node: node[0])
                if sorted_face not in unique_sorted_faces:
                    unique_faces.append(face)
                    unique_sorted_faces.append(sorted_face)
                    unique_markers.append(self.marker_pml_3[it])
                
            it += 1
        

        for i, face in enumerate(unique_faces):
            self.pml_string += f"1 0 {unique_markers[i]}\n"
            self.pml_string += "4 " + " ".join(map(str, [node[0] for node in face])) + "\n"
            self.num_facet += 1

    # ==========================================================================
    # .poly file writer
    # ==========================================================================

    def _write_poly_file(self, path: Path) -> None:
        """Write the assembled .poly file."""
        with open(path, 'w') as f:
            # Part 1 — Node list
            f.write("# Part 1 - Node List\n")
            f.write(f"{self.num_node} 3 0 1\n")

            f.write("# Points for cube (domain)\n")
            for node in self.node_list_domain_prism:
                f.write(f"{node[0]} {node[1]} {node[2]} {node[3]} {node[4]}\n")
            f.write("\n")

            f.write("# Layers\n")
            for layer in range(self.num_layers):
                node_list_interface = self.node_list_of_list_of_interfaces_old[layer]
                label = "# interface air-earth" if layer == 0 else f"# interface layer {layer-1}"
                f.write(label + "\n")
                for node in node_list_interface:
                    f.write(f"{node[0]} {node[1]} {node[2]} {node[3]} {node[4]}\n")
                f.write("\n")
            f.write("\n")

            f.write("# receivers\n")
            for i, node in enumerate(self.node_list_receiver):
                f.write(f"{node[0]} {round(node[1],5)} {round(node[2],5)} {round(node[3],5)} {node[4]}\n")
                if (i + 1) % 3 == 0:
                    f.write("\n")
            f.write("\n")

            f.write("# source nodes\n")
            for node in self.node_list_source:
                f.write(f"{node[0]} {round(node[1],5)} {round(node[2],5)} {round(node[3],5)} {node[4]}\n")
            f.write("\n")

            if self.box_present:
                f.write("# source box nodes\n")
                for node in self.node_list_source_box:
                    f.write(f"{node[0]} {round(node[1],5)} {round(node[2],5)} {round(node[3],5)} {node[4]}\n")
                f.write("\n")

            f.write("# Anomaly nodes\n")
            for node in self.node_list_anomaly:
                f.write(f"{node[0]} {round(node[1],5)} {round(node[2],5)} {round(node[3],5)} {node[4]}\n")
            f.write("\n")

            if self.NUM_PML > 0:
                f.write("# PML nodes for n PML Layers\n")
                f.write("# Type 2 Interface Intersections\n")
                for node in self.node_list_PML_2:
                    f.write(f"{node[0]} {node[1]} {node[2]} {node[3]} {node[4]}\n")
                f.write("\n")
                f.write("# Type 3: Corner Cubes\n")
                for node in self.node_list_PML_3:
                    f.write(f"{node[0]} {node[1]} {node[2]} {node[3]} {node[4]}\n")
                f.write("\n")

            # Part 2 — Facet list
            f.write("\n# Part 2 - Facet List\n")
            f.write(f"{self.num_facet} 1\n")
            f.write(self.domain_string)
            f.write("\n")
            f.write(self.interface_string)
            f.write("\n")
            f.write(self.interfaces_string)
            f.write("\n")
            if self.source_string:
                f.write(self.source_string)
            f.write("\n")
            if self.box_present:
                f.write(self.source_box_string)
                f.write("\n")
            f.write(self.anomaly_string)
            f.write("\n")
            if self.NUM_PML > 0:
                f.write(self.pml_string)
                f.write("\n")

            # Part 3 — Hole list
            f.write("\n# Part 3 - Hole List\n")
            f.write("0\n")

            # Part 4 — Region list
            f.write("\n# Part 4 - Region List\n")
            f.write(f"{self.num_regions}\n")
            for region in self.regions:
                f.write(
                    f"{region[0]} {region[1]} {region[2]} {region[3]} "
                    f"{region[4]} {region[5]} {region[6]} \n"
                )
