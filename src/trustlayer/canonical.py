"""
Canonical JSON hashing utility (fixes the Day-0 Check A hash-canonicalization bug).

BUG BEING FIXED: earlier receipts hashed a `sort_keys=True` string but wrote the file
with `sort_keys=False` plus the hash field appended, so re-hashing the on-disk file did
NOT reproduce the stored hash -- breaking the tamper-evidence claim.

CONTRACT: write_hashed_json(path, payload, hash_key) writes the file such that:
  1. the bytes on disk are exactly the canonical (sorted, indent=2) serialization of
     `payload` WITH the hash field included, and
  2. verify_hashed_json(path, hash_key) recomputes the hash over the on-disk bytes minus
     the hash field and confirms it matches.

The hash is computed over the canonical serialization of the payload with the hash field
set to an empty-string PLACEHOLDER, then the real hash is inserted and the file is written
in that same canonical order. Verification reverses this exactly.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


def _canonical(payload: dict) -> str:
    return json.dumps(payload, indent=2, sort_keys=True)


def compute_hash(payload: dict, hash_key: str) -> str:
    """Hash the canonical form of `payload` with hash_key set to a fixed placeholder."""
    tmp = dict(payload)
    tmp[hash_key] = ""  # placeholder so the hash never depends on itself
    return hashlib.sha256(_canonical(tmp).encode()).hexdigest()


def write_hashed_json(path: str | Path, payload: dict, hash_key: str) -> str:
    """Write `payload` to `path` in canonical order with a self-consistent hash.

    Returns the hash. The bytes on disk ARE the canonical serialization including the
    hash field, so verify_hashed_json() will reproduce the hash exactly.
    """
    h = compute_hash(payload, hash_key)
    out = dict(payload)
    out[hash_key] = h
    Path(path).write_text(_canonical(out))
    return h


def verify_hashed_json(path: str | Path, hash_key: str) -> tuple[bool, str, str]:
    """Re-derive the hash from the on-disk file and compare to the stored value.

    Returns (ok, stored_hash, recomputed_hash).
    """
    data = json.loads(Path(path).read_text())
    stored = data.get(hash_key, "")
    recomputed = compute_hash(data, hash_key)
    return (stored == recomputed, stored, recomputed)
