"""Create the frozen 2 x 2 x 4 mm Akantu extrinsic cohesive contract mesh.

The interface is a declared physical surface used only to test the insertion
API. It is not a candidate natural crack path.
"""

from __future__ import annotations

from collections import Counter
from itertools import permutations
from math import sqrt
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "input"
OUT.mkdir(parents=True, exist_ok=True)

LENGTHS_M = (2.0e-3, 2.0e-3, 4.0e-3)
CELLS = (2, 2, 4)
INTERFACE_Z_M = 2.0e-3


def node_index(i: int, j: int, k: int) -> int:
    nx, ny, _ = CELLS
    return (k * (ny + 1) + j) * (nx + 1) + i


nodes: list[tuple[float, float, float]] = []
for k in range(CELLS[2] + 1):
    for j in range(CELLS[1] + 1):
        for i in range(CELLS[0] + 1):
            nodes.append(
                (
                    -LENGTHS_M[0] / 2 + LENGTHS_M[0] * i / CELLS[0],
                    -LENGTHS_M[1] / 2 + LENGTHS_M[1] * j / CELLS[1],
                    LENGTHS_M[2] * k / CELLS[2],
                )
            )


def determinant(a, b, c, d) -> float:
    u = [b[q] - a[q] for q in range(3)]
    v = [c[q] - a[q] for q in range(3)]
    w = [d[q] - a[q] for q in range(3)]
    return (
        u[0] * (v[1] * w[2] - v[2] * w[1])
        - u[1] * (v[0] * w[2] - v[2] * w[0])
        + u[2] * (v[0] * w[1] - v[1] * w[0])
    )


tetrahedra: list[tuple[int, int, int, int]] = []
for k in range(CELLS[2]):
    for j in range(CELLS[1]):
        for i in range(CELLS[0]):
            origin = (i, j, k)
            for order in permutations(range(3)):
                p1 = list(origin)
                p1[order[0]] += 1
                p2 = p1.copy()
                p2[order[1]] += 1
                p3 = [origin[q] + 1 for q in range(3)]
                tet = [
                    node_index(*origin),
                    node_index(*p1),
                    node_index(*p2),
                    node_index(*p3),
                ]
                # Akantu expects positive tetrahedral orientation.
                if determinant(*(nodes[n] for n in tet)) < 0:
                    tet[1], tet[2] = tet[2], tet[1]
                tetrahedra.append(tuple(tet))


faces: Counter[tuple[int, int, int]] = Counter()
for tet in tetrahedra:
    a, b, c, d = tet
    for face in ((a, b, c), (a, b, d), (a, c, d), (b, c, d)):
        faces[tuple(sorted(face))] += 1


def on_z(face: tuple[int, int, int], z: float) -> bool:
    return all(abs(nodes[n][2] - z) < 1.0e-15 for n in face)


interface = sorted(face for face, n in faces.items() if n == 2 and on_z(face, INTERFACE_Z_M))
top = sorted(face for face, n in faces.items() if n == 1 and on_z(face, LENGTHS_M[2]))
bottom = sorted(face for face, n in faces.items() if n == 1 and on_z(face, 0.0))


def triangle_area(face: tuple[int, int, int]) -> float:
    a, b, c = (nodes[n] for n in face)
    u = [b[q] - a[q] for q in range(3)]
    v = [c[q] - a[q] for q in range(3)]
    cross = (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )
    return 0.5 * sqrt(sum(q * q for q in cross))


volume = sum(abs(determinant(*(nodes[n] for n in tet))) / 6 for tet in tetrahedra)
interface_area = sum(triangle_area(face) for face in interface)
assert len(tetrahedra) == 96
assert len(interface) == 8
assert abs(volume - 16.0e-9) < 1.0e-20
assert abs(interface_area - 4.0e-6) < 1.0e-18
assert len(top) == len(bottom) == 8


elements: list[tuple[int, int, tuple[int, ...]]] = []
for face in interface:
    elements.append((2, 2, tuple(n + 1 for n in face)))
for face in top:
    elements.append((3, 3, tuple(n + 1 for n in face)))
for face in bottom:
    elements.append((4, 4, tuple(n + 1 for n in face)))
for tet in tetrahedra:
    elements.append((1, 1, tuple(n + 1 for n in tet)))

mesh_path = OUT / "block_2x2x4mm.msh"
with mesh_path.open("w", encoding="ascii", newline="\n") as stream:
    stream.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
    stream.write(
        '$PhysicalNames\n4\n2 2 "interface"\n2 3 "top"\n'
        '2 4 "bottom"\n3 1 "body"\n$EndPhysicalNames\n'
    )
    stream.write(f"$Nodes\n{len(nodes)}\n")
    for idx, xyz in enumerate(nodes, start=1):
        stream.write(f"{idx} {xyz[0]:.16g} {xyz[1]:.16g} {xyz[2]:.16g}\n")
    stream.write("$EndNodes\n")
    stream.write(f"$Elements\n{len(elements)}\n")
    for idx, (physical, geometric, connectivity) in enumerate(elements, start=1):
        element_type = 2 if len(connectivity) == 3 else 4
        stream.write(
            f"{idx} {element_type} 2 {physical} {geometric} "
            + " ".join(str(n) for n in connectivity)
            + "\n"
        )
    stream.write("$EndElements\n")


gc = 7.6  # J/m^2; illustrative official example value, not ceramic calibration.
sigma_c = 70.0e6  # Pa; illustrative official example value, not ceramic calibration.
delta_c = 2 * gc / sigma_c
material_path = OUT / "material.dat"
material_path.write_text(
    """model solid_mechanics_model_cohesive [
  material elastic [
    name = body
    rho = 2500
    E = 70e9
    nu = 0.22
    finite_deformation = false
  ]
  material cohesive_linear [
    name = insertion
    beta = 3
    G_c = 7.6
    sigma_c = 70e6
    delta_c = 2.1714285714285714e-7
  ]
]
""",
    encoding="ascii",
)

metadata = {
    "units": "SI",
    "dimensions_m": list(LENGTHS_M),
    "cells": list(CELLS),
    "bulk_tetrahedra": len(tetrahedra),
    "nodes_before_insertion": len(nodes),
    "interface_facets": len(interface),
    "interface_area_m2": interface_area,
    "bulk_volume_m3": volume,
    "toy_material_unvalidated": True,
    "material": {
        "rho_kg_m3": 2500.0,
        "E_Pa": 70.0e9,
        "nu": 0.22,
        "Gc_J_m2": gc,
        "sigma_c_Pa": sigma_c,
        "beta": 3.0,
        "delta_c_m": delta_c,
    },
}
(OUT / "generated_geometry.json").write_text(
    json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps(metadata, indent=2))
