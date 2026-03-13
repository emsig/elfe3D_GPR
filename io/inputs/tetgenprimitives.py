"""
tetgenprimitives.py
-------------------
Low-level TetGen data structures and geometry helper functions.

Contains:
  - Node    : a 3D mesh node with index and boundary marker
  - Facet   : a planar polygon face (list of node indices + marker)
  - Region  : a TetGen region seed point (attribute, volume, material params)
  - Geometry helpers: all helper functions used during assembly

Node implements __getitem__ and __iter__ so existing code that uses
list-style node[0], node[1] etc. continues to work unchanged, while
new code can use node.x, node.y, node.z for clarity.

Node format convention (matches original throughout):
    node[0] = index
    node[1] = x
    node[2] = y
    node[3] = z
    node[4] = marker
"""

import math
import numpy as np
from dataclasses import dataclass, field
from itertools import groupby


# =============================================================================
# Core data structures
# =============================================================================

@dataclass
class Node:
    """
    A mesh node in 3D space with a TetGen boundary marker.

    Supports both attribute access (node.x) and legacy list-style
    index access (node[1]) for compatibility with existing geometry code.
    """

    index: int
    x: float
    y: float
    z: float
    marker: int = 0

    # ------------------------------------------------------------------
    # Backward-compatible list-style access
    # node[0] = index, node[1] = x, node[2] = y, node[3] = z, node[4] = marker
    # ------------------------------------------------------------------

    def __getitem__(self, i: int):
        return (self.index, self.x, self.y, self.z, self.marker)[i]

    def __iter__(self):
        return iter((self.index, self.x, self.y, self.z, self.marker))

    def __len__(self) -> int:
        return 5

    def __lt__(self, other: "Node") -> bool:
        """Needed for sorted() calls on node lists (sorts by index)."""
        return self.index < other.index

    def __str__(self) -> str:
        return f"Node({self.index}: {self.x:.4f}, {self.y:.4f}, {self.z:.4f} m={self.marker})"


@dataclass
class Region:
    """
    A TetGen region seed: a point inside a volumetric region that
    carries a region attribute (material ID), max volume constraint,
    and material electromagnetic parameters for the regionparameters file.

    The region list format used throughout the assembler is:
        [label, x, y, z, attribute, max_volume, comment, rho, mu_r, eps_r]
    Indices 7, 8, 9 are written to regionparameters.txt only.
    """
    label: int              # region attribute / material ID
    x: float
    y: float
    z: float
    attribute: int          # TetGen region attribute (same as label here)
    max_volume: str         # formatted string, e.g. "1.2345e-06"
    comment: str = ""       # human-readable label written as # comment
    rho: float = 0.0
    mu_r: float = 1.0
    eps_r: float = 1.0

    def to_list(self) -> list:
        """Return the 10-element list format used in the original assembler."""
        return [self.label, self.x, self.y, self.z, self.attribute,
                self.max_volume, self.comment, self.rho, self.mu_r, self.eps_r]


# =============================================================================
# Node creation helpers
# =============================================================================

def create_rectangular_prism_nodes(
    x_min: float, x_max: float,
    y_min: float, y_max: float,
    z: float,
    marker: int,
    num_node: int
) -> tuple[list[Node], int]:
    """
    Create 4 corner nodes of a rectangular face at a given z-level.

    Nodes are created in counterclockwise order (indices 0,1,3,2 of the
    corner grid) to match the original convention.

    Returns
    -------
    nodes    : list of 4 Node objects
    num_node : updated global node counter
    """
    nodes = []
    for i in [0, 1, 3, 2]:  # counterclockwise order — DO NOT CHANGE
        x = x_max if i & 1 else x_min
        y = y_max if i & 2 else y_min
        num_node += 1
        nodes.append(Node(index=num_node, x=x, y=y, z=z, marker=marker))
    return nodes, num_node


# =============================================================================
# Angular sorting helpers (used for face polygon generation)
# =============================================================================

def angular_sort(nodes: list, center_axis_1: int, center_axis_2: int) -> list[int]:
    """
    Sort nodes angularly based on their relative positions in a 2D plane
    defined by center_axis_1 and center_axis_2 (1=x, 2=y, 3=z in node indexing).

    Returns a list of node indices (node[0]) in angular order, starting
    from the node with the smallest index.
    """
    center_x = sum(node[center_axis_1] for node in nodes) / len(nodes)
    center_z = sum(node[center_axis_2] for node in nodes) / len(nodes)

    def angle_from_center(node):
        x = node[center_axis_1]
        z = node[center_axis_2]
        return -math.atan2(z - center_z, x - center_x)  # negated to reverse order

    sorted_nodes = sorted(nodes, key=angle_from_center)

    start_node = min(sorted_nodes, key=lambda node: node[0])
    start_index = sorted_nodes.index(start_node)
    sorted_nodes = sorted_nodes[start_index:] + sorted_nodes[:start_index]

    return [node[0] for node in sorted_nodes]


def get_face_nodes(
    node_list: list,
    axis_value: float,
    axis_index: int,
    sort_axis_1: int,
    sort_axis_2: int
) -> list[int]:
    """
    Filter nodes to those on a given axis plane, then sort angularly.

    Parameters
    ----------
    node_list   : list of nodes (Node objects or 5-element lists)
    axis_value  : the coordinate value to filter on
    axis_index  : which axis to filter (1=x, 2=y, 3=z in node indexing)
    sort_axis_1 : first sort axis
    sort_axis_2 : second sort axis

    Returns list of node indices (ints).
    """
    filtered = [node for node in node_list if node[axis_index] == axis_value]
    return angular_sort(filtered, center_axis_1=sort_axis_1, center_axis_2=sort_axis_2)


# =============================================================================
# Face string builders (used during poly file assembly)
# =============================================================================

def create_face_string(
    marker: int,
    nodes: list,
    interfaces_nodes,
    axis_value: float,
    axis_index: int,
    sort_axis_1: int,
    sort_axis_2: int,
) -> str:
    """
    Build the TetGen .poly facet string for one face.

    Used for PML Type-1 and Type-2 faces. The facet may include
    polygon holes (interface lines) if interfaces_nodes is provided.
    """
    face_nodes_list = get_face_nodes(
        node_list=nodes,
        axis_value=axis_value,
        axis_index=axis_index,
        sort_axis_1=sort_axis_1,
        sort_axis_2=sort_axis_2,
    )
    num_face_nodes = len(face_nodes_list)

    if axis_index == 3:
        interfaces_nodes_list = []
    elif interfaces_nodes is None:
        interfaces_nodes_list = []
    else:
        grouped = [
            list(g)
            for _, g in groupby(
                sorted(interfaces_nodes, key=lambda n: n[3]),
                key=lambda n: n[3],
            )
        ]
        interfaces_nodes_list = [
            get_face_nodes(
                node_list=iface_nodes,
                axis_value=axis_value,
                axis_index=axis_index,
                sort_axis_1=sort_axis_1,
                sort_axis_2=sort_axis_2,
            )
            for iface_nodes in grouped
        ]

    if axis_index == 3 or interfaces_nodes is None:
        num_polygon = 1
    else:
        num_polygon = 1 + len(interfaces_nodes_list)

    face_information = [num_polygon, 0, marker]

    str_1 = " ".join(map(str, face_information))
    str_2 = " ".join(map(str, face_nodes_list))
    str_3 = "\n".join(
        f"{len(iface_nodes)} " + " ".join(map(str, iface_nodes))
        for iface_nodes in interfaces_nodes_list
    )

    if interfaces_nodes is None:
        return f"{str_1}\n{num_face_nodes} {str_2}\n"
    else:
        return f"\n{str_1}\n{num_face_nodes} {str_2}\n{str_3}\n"


# =============================================================================
# Cuboid face construction (used for anomaly, source box, PML type-3)
# =============================================================================

def create_cuboid_faces_from_nodes(nodes: list, tol: float = 1e-6) -> list[list]:
    """
    Given 8 nodes at the vertices of a cuboid, return the 6 faces,
    each as a list of 4 nodes in proper order.

    Nodes are grouped by y-coordinate into top and bottom faces, then
    lateral faces are built from corresponding edges.
    """
    groups = {}
    for node in nodes:
        y = node[2]
        found_key = None
        for key in groups:
            if abs(key - y) < tol:
                found_key = key
                break
        if found_key is None:
            groups[y] = []
            found_key = y
        groups[found_key].append(node)

    if len(groups) != 2:
        raise ValueError("Expected exactly 2 distinct y values for a cuboid.")

    ys = sorted(groups.keys())
    bottom_nodes = groups[ys[0]]
    top_nodes    = groups[ys[1]]

    if len(bottom_nodes) != 4 or len(top_nodes) != 4:
        raise ValueError("Expected exactly 4 nodes per face (bottom and top).")

    def order_face(face_nodes):
        cx = sum(n[1] for n in face_nodes) / len(face_nodes)
        cz = sum(n[3] for n in face_nodes) / len(face_nodes)
        return sorted(face_nodes, key=lambda n: math.atan2(n[3] - cz, n[1] - cx))

    bottom_ordered = order_face(bottom_nodes)
    top_ordered    = order_face(top_nodes)

    lateral_faces = []
    for i in range(4):
        next_i = (i + 1) % 4
        face = [
            bottom_ordered[i], bottom_ordered[next_i],
            top_ordered[next_i], top_ordered[i],
        ]
        lateral_faces.append(face)

    return lateral_faces + [top_ordered[:], bottom_ordered[:]]


def create_half_cuboid_faces_from_nodes(
    nodes: list,
    not_needed_node: tuple,
    tol: float = 1e-6
) -> list[list]:
    """
    Like create_cuboid_faces_from_nodes but removes the face that
    contains the given node coordinate (used for PML corner cubes that
    share a face with the main domain).

    not_needed_node : (x, y, z) tuple of the corner to exclude.
    """
    faces = create_cuboid_faces_from_nodes(nodes, tol=tol)

    needed_faces = []
    for face in faces:
        if not any(
            float(node[1]) == float(not_needed_node[0])
            and float(node[2]) == float(not_needed_node[1])
            and float(node[3]) == float(not_needed_node[2])
            for node in face
        ):
            needed_faces.append(face)
    return needed_faces
