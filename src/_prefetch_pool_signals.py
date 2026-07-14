"""Prefetch all structured signals for the enriched pool, with a slow-API cache to /tmp.
Writes data/gold/_checkA_pool_signals_cache.json so the frozen-set builder is fast + offline-reproducible."""

from __future__ import annotations
import json, re, sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import requests

CAND = [
    "SIRPG",
    "CD28",
    "ICOS",
    "IL2RB",
    "CCR7",
    "ITK",
    "CD5",
    "LCK",
    "ZAP70",
    "STAT4",
    "IKZF3",
    "CD6",
    "IL10",
    "GATA3",
    "RUNX1",
    "PTPN2",
    "TNFAIP3",
    "SUOX",
    "TCF7L2",
    "APOBR",
    "HNF1A",
    "FTO",
    "CDKAL1",
    "KCNQ1",
    "PGM1",
]
OUT = Path("data/gold/_checkA_pool_signals_cache.json")
MYGENE = "https://mygene.info/v3/query"
OT = "https://api.platform.opentargets.org/api/v4/graphql"
GW = "https://www.ebi.ac.uk/gwas/rest/api"
IM = re.compile(
    "|".join(
        [
            r"\bt cell\b",
            r"\bt-cell\b",
            r"\blymphocyte\b",
            r"\bleukocyte\b",
            r"immune (response|system|process)",
            r"\binterleukin\b",
            r"antigen (processing|presentation|receptor)",
            r"\bcytokine\b",
            r"inflammat",
            r"adaptive immun",
            r"\bthymus\b",
            r"regulatory t",
        ]
    )
)


def mygene(g):
    r = requests.get(
        MYGENE,
        params={
            "q": f"symbol:{g}",
            "species": "human",
            "fields": "go.BP.term,go.MF.term,ensembl.gene",
        },
        timeout=45,
    ).json()
    h = r.get("hits", [])
    hit = h[0] if h else {}
    go = hit.get("go", {})
    terms = []
    for a in ("BP", "MF"):
        it = go.get(a, [])
        it = [it] if isinstance(it, dict) else it
        terms += [i.get("term", "") for i in it]
    imm = sorted({t for t in terms if IM.search(t.lower())})
    e = hit.get("ensembl", {})
    e = e[0] if isinstance(e, list) and e else e
    return {"immune": imm, "has_immune": len(imm) > 0, "ensembl": (e or {}).get("gene")}


def otdiff(ens):
    if not ens:
        return {"t1d": 0.0, "t2d": 0.0}
    q = "query T($id:String!){ target(ensemblId:$id){ associatedDiseases(page:{index:0,size:30}){ rows{ disease{ name } score } } } }"
    rows = []
    for _ in range(3):
        r = requests.post(
            OT, json={"query": q, "variables": {"id": ens}}, timeout=45
        ).json()
        d = (r or {}).get("data") or {}
        tgt = d.get("target") or {}
        ad = tgt.get("associatedDiseases") or {}
        rows = ad.get("rows")
        if rows is not None:
            break
        time.sleep(1.0)
    rows = rows or []
    t1 = t2 = 0.0
    for row in rows:
        n = row["disease"]["name"].lower()
        if "type 1 diabetes" in n:
            t1 = max(t1, round(row["score"], 3))
        elif "type 2 diabetes" in n:
            t2 = max(t2, round(row["score"], 3))
    return {"t1d": t1, "t2d": t2}


def gwas(g):
    tags = set()
    try:
        r = requests.get(
            f"{GW}/singleNucleotidePolymorphisms/search/findByGene",
            params={"geneName": g},
            timeout=45,
        )
        snps = r.json().get("_embedded", {}).get("singleNucleotidePolymorphisms", [])
        for snp in snps[:6]:
            try:
                a = requests.get(
                    f"{GW}/singleNucleotidePolymorphisms/{snp['rsId']}/associations",
                    params={"projection": "associationBySnp"},
                    timeout=15,
                ).json()
                for assoc in a.get("_embedded", {}).get("associations", []):
                    for t in assoc.get("efoTraits", []):
                        tr = (t.get("trait") or "").lower()
                        if "diabet" in tr:
                            tags.add(tr)
            except:
                pass
    except:
        pass
    return {"t2d_tag": any("type 2" in t for t in tags), "traits": sorted(tags)}


def one(g):
    mg = mygene(g)
    ot = otdiff(mg["ensembl"])
    gw = gwas(g)
    print(
        f"  {g:8s} immune={mg['has_immune']!s:5s} T1D={ot['t1d']:.3f} T2D={ot['t2d']:.3f} t2d_tag={gw['t2d_tag']}",
        flush=True,
    )
    return g, {
        "immune_go_terms": mg["immune"],
        "has_immune": mg["has_immune"],
        "ensembl": mg["ensembl"],
        "ot_t1d": ot["t1d"],
        "ot_t2d": ot["t2d"],
        "gwas_t2d_tag": gw["t2d_tag"],
        "gwas_traits": gw["traits"],
    }


if __name__ == "__main__":
    print("prefetching pool signals (parallel)...")
    res = {}
    with ThreadPoolExecutor(max_workers=6) as ex:
        for g, d in ex.map(one, CAND):
            res[g] = d
    OUT.write_text(json.dumps(res, indent=2))
    print(f"wrote {OUT}  ({len(res)} genes)")
