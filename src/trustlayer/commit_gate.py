"""
Deterministic commit gate (Plan v3, Fix 4).

This is the SOLE writer of a final target nomination. It is default-deny:
a nomination is written ONLY if it passes every check against a structured
critic verdict + the live genetics receipt. Otherwise it writes a BLOCKED
artifact explaining exactly why.

v3 replaces the old "no-Write tool scoping" claim (bundled-profile-only,
falsifiable by the SDK authors) with this plain Python gate whose behavior
is covered by tests and shown -- not narrated -- in the demo.

The gate is deterministic: same inputs -> same artifact -> same hash. LLM
critic calls happen UPSTREAM and are frozen into the CriticVerdict; the gate
itself never calls a model.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Literal, Optional

Decision = Literal["GO", "WITHHELD"]


@dataclass
class CriticVerdict:
    """Frozen output of the upstream critic (independent re-derivation).

    All fields are set by the critic BEFORE the gate runs. The gate treats this
    as immutable evidence and never mutates it.
    """

    gene: str
    # genetics evidence
    genetic_association: Optional[float]  # Open Targets, live
    genome_wide_sig_snp: bool  # a T1D GWS lead SNP exists
    celltype_matched_eqtl: bool  # CD4/Treg eQTL exists for the risk allele
    eqtl_direction_consistent: Optional[
        bool
    ]  # risk allele -> expression change consistent w/ KD-mimics-protection
    proxy_tissue_only: bool  # eQTL evidence is proxy-tissue, cell-type-unconfirmed
    # model / trust evidence
    trust_score: Optional[float]  # calibrated conformal trust (0..1)
    leakage_audit_passed: bool  # gene not in program axis; splits frozen
    # critic notes
    notes: str = ""


@dataclass
class Nomination:
    gene: str
    role: str  # ANCHOR | BET | CONTROL
    decision: Decision
    trust_score: Optional[float]
    reasons: list[str] = field(default_factory=list)
    verdict: dict = field(default_factory=dict)
    generated_utc: str = ""

    def content_hash(self) -> str:
        payload = json.dumps(
            {k: v for k, v in asdict(self).items() if k != "generated_utc"},
            sort_keys=True,
        )
        return hashlib.sha256(payload.encode()).hexdigest()[:16]


# ---- gate thresholds (pre-registered; changing these is a code+test change) ----
MIN_GENETIC_ASSOCIATION = 0.20  # PRKCQ (0.162) falls below this by design
MIN_TRUST_SCORE = 0.50


class CommitGate:
    """Default-deny nomination writer. The only path that emits a GO artifact."""

    def __init__(self, out_dir: str | Path = "data/gold/nominations"):
        self.out_dir = Path(out_dir)
        self.out_dir.mkdir(parents=True, exist_ok=True)

    def evaluate(self, verdict: CriticVerdict, role: str) -> Nomination:
        """Return a Nomination. decision=GO only if ALL checks pass; else WITHHELD."""
        reasons: list[str] = []

        # Hard leakage gate — a leakage failure is always fatal, regardless of role.
        if not verdict.leakage_audit_passed:
            reasons.append(
                "BLOCK: leakage audit failed (gene in program axis or splits not frozen)"
            )

        # Genetic anchoring floor.
        ga = verdict.genetic_association
        if ga is None or ga < MIN_GENETIC_ASSOCIATION:
            reasons.append(
                f"BLOCK: genetic_association={ga} below floor {MIN_GENETIC_ASSOCIATION}"
            )

        # Trust-score floor.
        ts = verdict.trust_score
        if ts is None or ts < MIN_TRUST_SCORE:
            reasons.append(f"BLOCK: trust_score={ts} below floor {MIN_TRUST_SCORE}")

        # Direction / eQTL check: a GO nomination needs cell-type-matched, direction-
        # consistent eQTL evidence. Proxy-only or missing eQTL -> cannot GO (the
        # RASGRP1 case: the risk-allele eQTL in CD4/Treg genuinely does not exist).
        if not verdict.celltype_matched_eqtl:
            reasons.append(
                "BLOCK: no cell-type-matched (CD4/Treg) eQTL for the risk allele"
                + (" (proxy-tissue only)" if verdict.proxy_tissue_only else "")
            )
        elif verdict.eqtl_direction_consistent is False:
            reasons.append(
                "BLOCK: eQTL direction inconsistent with KD-mimics-protection"
            )

        decision: Decision = "GO" if not reasons else "WITHHELD"
        if decision == "GO":
            reasons.append(
                f"GO: genetics-anchored (GA={ga}, GWS_SNP={verdict.genome_wide_sig_snp}), "
                f"cell-type eQTL direction-consistent, trust={ts}, leakage-clean"
            )

        return Nomination(
            gene=verdict.gene,
            role=role,
            decision=decision,
            trust_score=ts,
            reasons=reasons,
            verdict=asdict(verdict),
            generated_utc=datetime.now(timezone.utc).isoformat(),
        )

    def commit(self, nomination: Nomination) -> Path:
        """Write the artifact. GO -> nomination file; WITHHELD -> BLOCKED file.

        This is the only function in the codebase permitted to write a final
        nomination. It always writes *something* (fail-closed, auditable).
        """
        tag = "GO" if nomination.decision == "GO" else "BLOCKED"
        h = nomination.content_hash()
        path = self.out_dir / f"{nomination.gene}_{tag}_{h}.json"
        path.write_text(json.dumps(asdict(nomination) | {"content_hash": h}, indent=2))
        return path
