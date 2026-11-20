"""Shared, session-scoped registry for sequence records.

This module provides a small compatibility layer between the application pages.
Biopython SeqRecord objects are kept in Streamlit session state so that records
retrieved or selected in one module can be reused by another module without
re-downloading them.
"""
from __future__ import annotations

from datetime import datetime, timezone
from typing import Any, Iterable

import streamlit as st


REGISTRY_KEY = "shared_transcript_registry"
ACTIVE_ACCESSION_KEY = "shared_active_accession"


def _session_state():
    """Return Streamlit's session state (isolated here for easier maintenance)."""
    return st.session_state


def _registry() -> dict[str, dict[str, Any]]:
    """Get or initialize the registry stored in the current user session."""
    state = _session_state()
    registry = state.get(REGISTRY_KEY)
    if not isinstance(registry, dict):
        registry = {}
        state[REGISTRY_KEY] = registry
    return registry


def _record_accession(record: Any) -> str:
    """Return a stable accession-like identifier for a sequence record."""
    accession = getattr(record, "id", None)
    if accession is None:
        return ""
    accession = str(accession).strip()
    if not accession or accession == "<unknown id>":
        return ""
    return accession


def _safe_length(record: Any) -> int | None:
    try:
        return len(record.seq)
    except (AttributeError, TypeError):
        return None


def register_records(
    records: Iterable[Any] | Any,
    *,
    source: str = "Unknown",
    gene: str = "",
    organism: str = "",
) -> int:
    """Add or update records in the shared registry.

    Records are keyed by ``record.id``. Registering the same accession again
    updates its metadata and record object instead of creating a duplicate.
    Returns the number of valid records processed.
    """
    if records is None:
        return 0

    # Accept a single SeqRecord as well as an iterable of SeqRecords.
    if _record_accession(records):
        records = [records]
    else:
        try:
            records = iter(records)
        except TypeError:
            return 0

    registry = _registry()
    registered = 0
    timestamp = datetime.now(timezone.utc).isoformat(timespec="seconds")

    for record in records:
        accession = _record_accession(record)
        if not accession:
            continue

        previous = registry.get(accession, {})
        registry[accession] = {
            "record": record,
            "source": str(source or previous.get("source") or "Unknown"),
            "gene": str(gene or previous.get("gene") or ""),
            "organism": str(organism or previous.get("organism") or ""),
            "accession": accession,
            "length": _safe_length(record),
            "record_type": type(record).__name__,
            "registered_at": timestamp,
        }
        registered += 1

    return registered


def register_active_record(
    record: Any,
    *,
    gene: str = "",
    organism: str = "",
    source: str = "Gene & Transcript",
) -> bool:
    """Register a record and mark its accession as the active sequence."""
    accession = _record_accession(record) if record is not None else ""
    if not accession:
        return False

    register_records([record], source=source, gene=gene, organism=organism)
    _session_state()[ACTIVE_ACCESSION_KEY] = accession
    return True


def get_registered_records() -> list[Any]:
    """Return all registered SeqRecord objects in insertion order."""
    return [
        item["record"]
        for item in _registry().values()
        if isinstance(item, dict) and item.get("record") is not None
    ]


def get_registry_summary() -> list[dict[str, Any]]:
    """Return registry metadata without the non-displayable SeqRecord objects."""
    summary: list[dict[str, Any]] = []
    for accession, item in _registry().items():
        if not isinstance(item, dict):
            continue
        summary.append({
            key: value
            for key, value in item.items()
            if key != "record"
        })
    return summary


def get_active_accession() -> str | None:
    """Return the accession currently marked active, if it is registered."""
    accession = _session_state().get(ACTIVE_ACCESSION_KEY)
    if accession and accession in _registry():
        return str(accession)
    return None


def get_active_record() -> Any | None:
    """Return the active SeqRecord object, if available."""
    accession = get_active_accession()
    if accession is None:
        return None
    return _registry()[accession].get("record")


def get_record_by_accession(accession: str) -> Any | None:
    """Retrieve one registered record by its exact accession/record ID."""
    key = str(accession or "").strip()
    item = _registry().get(key)
    return item.get("record") if isinstance(item, dict) else None


def remove_record(accession: str) -> bool:
    """Remove a record from the current session registry."""
    key = str(accession or "").strip()
    registry = _registry()
    if key not in registry:
        return False
    del registry[key]
    if _session_state().get(ACTIVE_ACCESSION_KEY) == key:
        _session_state()[ACTIVE_ACCESSION_KEY] = None
    return True


def clear_registry() -> None:
    """Clear all shared sequence records for the current session."""
    _session_state()[REGISTRY_KEY] = {}
    _session_state()[ACTIVE_ACCESSION_KEY] = None


def set_active_accession(accession: str) -> bool:
    """Mark an already registered accession as active."""
    key = str(accession or "").strip()
    if key not in _registry():
        return False
    _session_state()[ACTIVE_ACCESSION_KEY] = key
    return True


# Backward-compatible aliases for pages that may use older naming.
get_shared_records = get_registered_records
get_shared_registry_summary = get_registry_summary
