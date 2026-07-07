# Research State: Chaining beats leaping — layered genotype→phenotype prediction

## Current Stage
SCOPE

## Research Question
Does predicting a cellular phenotype by *chaining* intermediate molecular layers
(DNA sequence → regulatory activity → gene expression → protein-network context)
outperform *leaping* directly from genotype to phenotype — and does the omnigenic
core/peripheral distinction predict *which* genes chain cleanly?

## Context (Hackathon)
- Event: Built with Claude: Life Sciences (Cerebral Valley x Gladstone), Research track.
- Team size ≤2. Prize: $100k API/usage credits. 1 week. $200 API credits + Claude Max 20x.
- Three provided Gladstone datasets (each = one molecular layer):
  1. T-cell Perturb-seq (Marson lab + Pritchard, Stanford) — gene KO → transcriptome.
  2. DNA→regulatory-activity + MPRA single-base effects (Pollard lab) — sequence → regulatory.
  3. Protein interaction networks (Krogan lab) — protein-protein context.

## Core Thesis
The genotype→phenotype map resists single-leap prediction but becomes more tractable
when decomposed into chained sub-mappings (the "chain the layers, don't leap the gap"
idea). Pritchard's omnigenic model supplies the theoretical spine: a few "core" genes
act directly; many "peripheral" genes act indirectly through the network.

## Key Decisions
- SCOPE: Chose integrative design (Option A) over single-layer projects — highest novelty,
  uses ≥2 of the 3 datasets, directly tests a named theory (omnigenic).

## Experiment Log
| Attempt | Method | Result | Status |
|---------|--------|--------|--------|

## Critique History
- Pre-COMPUTE: pending
- Post-COMPUTE: pending

## What Worked
-

## What Didn't Work
-

## Open Questions
- Which specific public Perturb-seq/MPRA datasets are downloadable pre-hackathon (Replogle
  et al. genome-wide Perturb-seq; Marson T-cell CRISPR screens; Pollard MPRA)?
- Does a small fine-tuned sequence model (Enformer/DNABERT/HyenaDNA head) predict MPRA
  activity well enough to be a usable "layer"?
- How to define "chaining accuracy" rigorously vs. a direct-leap baseline.

## Artifacts
- literature-review.md: missing
- reasoning.md: missing
- methodology.md: missing
- figures/: empty
- experiments.tsv: missing
