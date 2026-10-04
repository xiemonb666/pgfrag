from __future__ import annotations

import base64
import re
import struct
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

FIELD_COMPONENTS = {"opening": 3, "tractions": 3, "damage": 1, "delta max": 1}


def _array(node: ET.Element, name: str, components: int, tuples: int) -> np.ndarray:
    if node.get("Name") != name or node.get("type") != "Float64":
        raise ValueError(f"unexpected VTK field metadata for {name}: {node.attrib}")
    if node.get("format") != "binary":
        raise ValueError(f"{name}: expected inline binary VTK data")
    declared_components = int(node.get("NumberOfComponents", "1"))
    if declared_components != components:
        raise ValueError(f"{name}: components={declared_components}, expected {components}")
    raw = base64.b64decode("".join((node.text or "").split()), validate=True)
    if len(raw) < 4:
        raise ValueError(f"{name}: truncated UInt32 length header")
    declared_bytes = struct.unpack("<I", raw[:4])[0]
    if declared_bytes != len(raw):
        raise ValueError(f"{name}: declared block bytes {declared_bytes} != decoded bytes {len(raw)}")
    payload = raw[4:]
    if len(payload) % 8 or len(payload) != tuples * components * 8:
        raise ValueError(f"{name}: payload length {len(payload)} is inconsistent with tuples/components")
    return np.frombuffer(payload, dtype="<f8").copy().reshape(tuples, components)


def parse_pvtu(path: Path, output_root: Path, expected_cells: int = 8,
               expected_points: int = 54) -> dict:
    path = Path(path).resolve()
    output_root = Path(output_root).resolve()
    if not path.is_relative_to(output_root):
        raise ValueError("PVTU path escaped run output")
    pv_root = ET.parse(path).getroot()
    if pv_root.get("type") != "PUnstructuredGrid" or pv_root.get("compressor"):
        raise ValueError("unsupported PVTU type/compression")
    declarations = {}
    for node in pv_root.findall(".//PCellData/PDataArray"):
        name = node.get("Name")
        if name not in FIELD_COMPONENTS or name in declarations:
            raise ValueError(f"unexpected or duplicate PVTU field {name!r}")
        declarations[name] = node
    if set(declarations) != set(FIELD_COMPONENTS):
        raise ValueError(f"native cohesive fields differ: {sorted(declarations)}")
    pieces = pv_root.findall(".//Piece")
    if len(pieces) != 1 or not pieces[0].get("Source"):
        raise ValueError("expected one serial PVTU Piece")
    piece_path = (path.parent / pieces[0].get("Source")).resolve()
    if not piece_path.is_relative_to(output_root) or not piece_path.is_file():
        raise ValueError("PVTU Piece missing or escaped run output")
    root = ET.parse(piece_path).getroot()
    if root.get("type") != "UnstructuredGrid" or root.get("compressor"):
        raise ValueError("unsupported VTU type/compression")
    piece = root.find(".//Piece")
    if piece is None:
        raise ValueError("VTU Piece missing")
    cells, points = int(piece.get("NumberOfCells", "-1")), int(piece.get("NumberOfPoints", "-1"))
    if (cells, points) != (expected_cells, expected_points):
        raise ValueError(f"piece counts {(cells, points)} != {(expected_cells, expected_points)}")
    data = {}
    nodes = {node.get("Name"): node for node in piece.findall("./CellData/DataArray")}
    if set(nodes) != set(FIELD_COMPONENTS):
        raise ValueError(f"native CellData fields differ: {sorted(nodes)}")
    for name, components in FIELD_COMPONENTS.items():
        schema = declarations[name]
        if schema.get("type") != "Float64" or int(schema.get("NumberOfComponents", "1")) != components:
            raise ValueError(f"invalid PVTU schema for {name}")
        data[name] = _array(nodes[name], name, components, cells)
    match = re.search(r"_(\d+)\.pvtu$", path.name)
    return {"pvtu": str(path), "vtu": str(piece_path),
            "step_index": int(match.group(1)) if match else None,
            "cells": cells, "points": points, "fields": data}


def new_pvtu_files(output_root: Path, previous: set[Path]) -> set[Path]:
    current = {p.resolve() for p in Path(output_root).rglob("*.pvtu")}
    return current - {p.resolve() for p in previous}