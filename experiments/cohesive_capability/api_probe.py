"""Record relevant public Python symbols for ancestry/event API capability."""

from __future__ import annotations

import importlib.metadata as metadata
import json
from pathlib import Path

import akantu as aka


symbols = sorted(dir(aka))
relevant = [
    name
    for name in symbols
    if any(key in name.lower() for key in ("event", "listener", "parent", "origin", "lineage"))
]
result = {
    "akantu_version": metadata.version("akantu"),
    "module_path": aka.__file__,
    "relevant_public_symbols": relevant,
    "new_cohesive_nodes_event_exported": hasattr(aka, "NewCohesiveNodesEvent"),
    "element_inserter_api": dir(aka.CohesiveElementInserter),
    "mesh_connectivity_api": [name for name in dir(aka.Mesh) if "connect" in name.lower()],
    "finding": "No dedicated Python event/parent-lineage symbol found; connectivity-based mapping must pass its runtime consistency checks or remain UNKNOWN.",
}
print(json.dumps(result, indent=2))
