# Zenodo deposition — instant citable DOI (you own the Zenodo account)

## Steps
1. Log in at https://zenodo.org (GitHub SSO works). New Upload.
2. Upload the release tarball (release/trustlayer_release.tar.gz) OR connect the GitHub
   repo (Zenodo -> GitHub -> flip the repo switch ON, then re-tag a release; Zenodo mints
   a DOI automatically on each GitHub release).
3. Metadata (paste):

   Title:        TrustLayer: a calibrated-trust layer for perturbation-biology prediction
   Authors:      Lahiri, Anirudh
   Upload type:  Software
   Description:  A calibrated-trust layer for perturbation-effect prediction: split-conformal
                 coverage, a default-deny commit gate, and pre-registered shift-robustness
                 diagnostics. Across 5 architectures x a 4-level shift ladder on public
                 perturb-seq data, conformal coverage@0.90 is model-independent up to
                 cross-donor shift and breaks under cross-dataset shift; a small labeled
                 anchor recalibrates coverage under shift. Includes hashed receipts and a
                 self-red-team audit trail.
   License:      MIT
   Keywords:     conformal prediction; perturbation biology; calibration; single-cell; uncertainty
   Related:      Ahlmann-Eltze, Huber & Anders, Nat Methods 2025, doi:10.1038/s41592-025-02772-6 (references)

4. Publish -> Zenodo returns a DOI like 10.5281/zenodo.XXXXXXX immediately.
5. Put the DOI badge in the README top-fold and cite it in the submission + CITATION.cff.

## bioRxiv (optional, indexing lags past 9PM but the DOI above suffices for the deadline)
Submit release/PREPRINT_ABSTRACT.md expanded to a short methods note at
https://www.biorxiv.org/submit — category: Bioinformatics. Not required for the deadline
if the Zenodo DOI is live.
