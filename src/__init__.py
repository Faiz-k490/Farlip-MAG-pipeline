"""FaRLiP Ingestion Pipeline.

A config-driven Biopython + ENA ingestion layer that locks in the canonical
21-gene FaRLiP cluster organism (Leptolyngbya sp. JSC-1) and a polar
cyanobacterial MAG corpus (ENA PRJEB59431) for downstream comparative analysis.

Iteration 1 scope: data acquisition only. See ``DATA_SOURCES.md`` for the lock-in
spec and ``DEVLOG.md`` for change history.
"""

from __future__ import annotations

__version__ = "0.1.0"
