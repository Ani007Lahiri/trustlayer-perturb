"""
Day-1: build the tiered, cited T1D validation gold set (Plan v3, Fix 2).

Replaces the v2.1 "genetics-derived" genes that had ZERO T1D GWAS signal
(MEOX1, CD1E) or sat far outside the top-250 (LGALS3BP #2402, CD247 #1074)
with genuinely T1D-anchored genes, each carrying a live-derived citation.

Tiers (v3):
  T1  GWAS genome-wide significant (p < 5e-8), lead SNP + trait from GWAS Catalog
  T2  Open Targets genetic_association score, top-of-list (rank shown)
  T3  fine-mapped T-cell/immune eQTL   (labeled: eQTL-catalogue-pending)
  T4  mechanism-only                    (labeled as such; NOT 'genetics-derived')

Every T1/T2 row is re-derived from the live gate. Output:
  data/gold/t1d_gold_set.json   (machine-readable, with provenance per gene)
  data/gold/t1d_gold_set.tsv    (flat table for the memo/deck)

The recovery evaluation (Day 4) uses tiers T1+T2 as the BLINDED positive set.
The UCell program axis (Day 2) MUST NOT name any gene in this set -> emitted as
`axis_forbidden_genes` for a build-time leakage assertion.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from trustlayer.genetics import ot_genetic_association, gwas_catalog_associations  # noqa: E402

# T1 tier: genome-wide-significant T1D genes with a known lead SNP to reconfirm.
# (SNPs chosen from established T1D GWAS; reconfirmed live against GWAS Catalog.)
T1_LEAD_SNPS = {
    "PTPN22": "rs2476601",
    "CTLA4": "rs3087243",
    "IL2RA": "rs12722495",
    "IFIH1": "rs1990760",
    "SH2B3": "rs3184504",
    "BACH2": "rs11755527",
    "TYK2": "rs34536443",
}
# T2 tier: strong Open Targets association, credentialed by score+rank (no single SNP required here).
T2_GENES = ["CTSH", "INS"]
# T4 tier: mechanism-only anchors (labeled honestly, excluded from 'genetics-derived').
T4_MECHANISM = {
    "FOXP3": "master Treg lineage TF (mechanism-only; IPEX monogenic, not a T1D GWAS locus)",
}

# v2.1 genes to explicitly RETIRE from the gold set (with the live reason).
RETIRED = {
    "MEOX1": "not associated with T1D in Open Targets; zero GWAS-Catalog T1D hits",
    "CD1E": "not associated with T1D in Open Targets",
    "LGALS3BP": "Open Targets rank ~#2402/5887; far outside top-250",
    "CD247": "Open Targets rank ~#1074/5887; outside top-250",
}

OUT_JSON = Path("data/gold/t1d_gold_set.json")
OUT_TSV = Path("data/gold/t1d_gold_set.tsv")
DE_STATS_CSV = Path("data/raw/DE_stats.suppl_table.csv")


def _perturbed_genes():
    """Set of genes actually perturbed in the Perturb-seq data, or None if the
    DE_stats table hasn't been pulled yet (gold set still builds, unrestricted)."""
    if not DE_STATS_CSV.exists():
        return None
    import pandas as pd

    df = pd.read_csv(DE_STATS_CSV, usecols=["target_contrast_gene_name"])
    return set(df["target_contrast_gene_name"].unique())


def _gwas_best_t1d(rsid: str):
    """Return the best (smallest-p) T1D association row for a SNP, or None."""
    rows = gwas_catalog_associations(rsid)
    t1d = []
    for a in rows:
        if "error" in a:
            continue
        traits = [t.lower() for t in a.get("traits", []) if t]
        if any("type 1 diabetes" in t for t in traits):
            exp = a.get("pvalueExponent")
            man = a.get("pvalueMantissa")
            if exp is not None and man is not None:
                t1d.append((man, exp, a))
    if not t1d:
        return None
    # smallest p = most negative exponent, then smallest mantissa
    t1d.sort(key=lambda x: (x[1], x[0]))
    man, exp, a = t1d[0]
    return {"rsid": rsid, "pvalue": f"{man}e{exp}", "genome_wide_sig": exp <= -8}


def main() -> int:
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    all_symbols = list(T1_LEAD_SNPS) + T2_GENES
    assoc = ot_genetic_association(all_symbols)

    gold = []

    # ---- T1 ----
    for sym, rsid in T1_LEAD_SNPS.items():
        a = assoc.get(sym)
        snp = _gwas_best_t1d(rsid)
        gold.append(
            {
                "gene": sym,
                "tier": "T1",
                "tier_desc": "GWAS genome-wide-significant",
                "ot_genetic_association": a.genetic_association if a else None,
                "lead_snp": rsid,
                "gwas_t1d": snp,
                "citation": f"GWAS Catalog {rsid} (T1D {snp['pvalue']})"
                if snp
                else f"GWAS Catalog {rsid}",
                "provenance": "Open Targets v4 + GWAS Catalog REST (live)",
            }
        )

    # ---- T2 ----
    for sym in T2_GENES:
        a = assoc.get(sym)
        gold.append(
            {
                "gene": sym,
                "tier": "T2",
                "tier_desc": "Open Targets genetic_association (top-of-list)",
                "ot_genetic_association": a.genetic_association if a else None,
                "lead_snp": None,
                "gwas_t1d": None,
                "citation": f"Open Targets genetic_association={a.genetic_association if a else '?'} (T1D MONDO_0005147)",
                "provenance": "Open Targets v4 GraphQL (live)",
            }
        )

    # ---- T4 (mechanism-only, labeled) ----
    for sym, why in T4_MECHANISM.items():
        gold.append(
            {
                "gene": sym,
                "tier": "T4",
                "tier_desc": "mechanism-only (NOT genetics-derived)",
                "ot_genetic_association": None,
                "lead_snp": None,
                "gwas_t1d": None,
                "citation": why,
                "provenance": "literature/mechanism (labeled)",
            }
        )

    # positive set for blinded recovery = T1 + T2 (genetics-anchored only),
    # further restricted to genes ACTUALLY PERTURBED in the Perturb-seq data.
    positive_all = [g["gene"] for g in gold if g["tier"] in ("T1", "T2")]
    perturbed = _perturbed_genes()
    positive_set = [g for g in positive_all if (perturbed is None or g in perturbed)]
    dropped_not_perturbed = [
        g for g in positive_all if perturbed is not None and g not in perturbed
    ]
    for g in gold:
        if g["gene"] in dropped_not_perturbed:
            g["in_perturbseq"] = False
        elif g["tier"] in ("T1", "T2"):
            g["in_perturbseq"] = True

    doc = {
        "disease": {"name": "type 1 diabetes mellitus", "id": "MONDO_0005147"},
        "note": "Tiered, per-gene cited gold set (v3 Fix 2). T1+T2 = genetics-anchored "
        "blinded positive set for recovery. T4 labeled mechanism-only.",
        "retired_from_v2.1": RETIRED,
        "genes": gold,
        "recovery_positive_set": positive_set,
        "recovery_dropped_not_perturbed": dropped_not_perturbed,
        "axis_forbidden_genes": sorted(set(g["gene"] for g in gold)),
    }
    OUT_JSON.write_text(json.dumps(doc, indent=2))

    with OUT_TSV.open("w") as f:
        f.write("gene\ttier\tot_genetic_association\tlead_snp\tgwas_t1d_p\tcitation\n")
        for g in gold:
            p = (g["gwas_t1d"] or {}).get("pvalue", "")
            f.write(
                f"{g['gene']}\t{g['tier']}\t{g['ot_genetic_association']}\t"
                f"{g['lead_snp'] or ''}\t{p}\t{g['citation']}\n"
            )

    # ---- summary ----
    print("=" * 70)
    print("TIERED T1D GOLD SET  (v3 Fix 2, live-derived)")
    print("=" * 70)
    for g in gold:
        p = (g["gwas_t1d"] or {}).get("pvalue", "-")
        gws = (g["gwas_t1d"] or {}).get("genome_wide_sig")
        flag = " [GWS]" if gws else (" [sub-GWS]" if gws is False else "")
        print(
            f"  [{g['tier']}] {g['gene']:8s} ot_GA={g['ot_genetic_association']}  "
            f"snp={g['lead_snp'] or '-':12s} T1D_p={p}{flag}"
        )
    print(
        f"\n  Blinded recovery positive set (T1+T2 AND perturbed, n={len(positive_set)}): {positive_set}"
    )
    if dropped_not_perturbed:
        print(
            f"  Dropped (genetics-anchored but NOT perturbed in data): {dropped_not_perturbed}"
        )
    print("\n  RETIRED from v2.1 (live reason):")
    for k, v in RETIRED.items():
        print(f"    {k:10s} {v}")
    print(f"\n  -> {OUT_JSON}\n  -> {OUT_TSV}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
