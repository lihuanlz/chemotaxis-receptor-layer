# Erratum — v2.0.0 (2026-11-14)

## What happened

Release v1.0.x of this repository (manuscript: *"Amplitude and midpoint statistics of the same
chemotaxis cells demand incompatible receptor models"*, Zenodo doi:10.5281/zenodo.23402857)
claimed that no parameter set of the standard MWC-plus-perfect-adaptation model class can jointly
reproduce the amplitude table and the midpoint distribution of the Moore et al. (2024) single-cell
FRET archive.

**That claim is retracted.** It was produced by two definition mismatches in our own analysis
code (`code54_K12_joint_fit_rupture_quantification.py`), found on 2026-10-19 by an internal
zero-memory recheck and re-verified from scratch on 2026-11-14:

| | Observation side | Model side (v1 `K12_model`) |
|---|---|---|
| **E1 concentration axis** | total ligand T = B + F | returned foreground increment F |
| **E2 response level** | absolute threshold da = 0.5 (censoring 45–68% of cells) | relative half-maximum (censoring 2/268) |

E1 alone fabricates up to 1.21 dex of apparent mismatch (at B = 100 µM: model 7.3 µM vs measured
119.6 µM). A third, independent defect: the v1 λ sweep jumped from 3 to 10 and skipped the
λ ≈ 3.5 point that enters its own acceptance region.

## The corrected result

With consistent definitions (relative half-maximum, total-ligand axis), a single parameter set of
the standard model reproduces **both** statistics:

- amplitude table R² = 0.980, mean midpoint error **0.023 dex** (v1 reported 0.622 dex);
- every weighting λ ∈ [0, 30] lands inside the v1 acceptance region;
- closed-form vs numerical model midpoint agree to 1.1×10⁻¹².

Corrected manuscript: `manuscript/main.tex` — *"A midpoint-definition artefact reconciles the
standard MWC model with single-cell chemotaxis dose responses."* The refuted v1 manuscript and SI
are preserved in `manuscript/_archive_refuted_v1/` for the record.

## Claim-by-claim disposition

| v1 claim | Status |
|---|---|
| No parameter set meets the acceptance region | **Retracted** (entire corrected front is inside) |
| Midpoint error 0.622 dex at amplitude-optimal fit | **Retracted** (0.023 dex with consistent definitions) |
| N≈2–4 (midpoints) vs N≈12–39 (amplitude) contradiction | **Retracted** (a single N fits both; the two statistics are not independent constraints on N) |
| Eight repair hypotheses all fail | **Void** (premise collapsed) |
| Weber line K₁/₂ = 1.17·(1.95 + B), max 0.10 dex | Stands (absolute-threshold convention); the relative convention gives 0.91 + 1.06·B, max 0.05 dex |
| Prospective mixed-ligand hits +0.014/+0.038/+0.031 dex | Stands as numbers, **downgraded**: frozen band belonged to the absolute-threshold statistic; not used for model-class discrimination; the third reading was a control arm marked "not scored" |
| Adaptation relay: baseline restored within 2.2% after 5 s washout (2,217 cells, 77,595 events) | **Stands** (independent of the rupture analysis) |
| N = 39 as a fitted physical parameter | **Retracted**: data give only N ≳ 8, with a strong N–K_A ridge |
| "Weak-responder subpopulation" | **Retracted**: threshold artefact (threshold sat at the response-distribution median) |

New in v2: N-unidentifiability profile (`results/corrected_v2/deep.json`) and the residual
single-cell heterogeneity analysis (dispersion 0.260 → 0.022 dex from low to high background;
sensitivity–responsiveness Spearman ρ = +0.205, p = 6.6×10⁻⁷, n = 579).

## Provenance

- Internal recheck reports (Chinese, working record):
  `05_主线纲领与设计/_chem_recheck/复核报告_趋化_2026-10-19.md`,
  `深挖C_趋化_2026-10-19.md` in the project workspace.
- Re-verification 2026-11-14: `code/recheck_chemotaxis.py` rerun inside this repository
  (repo-relative paths) reproduces protocol A (0.980/0.622) and protocol D (0.980/0.023) exactly;
  output in `results/corrected_v2/recheck_results.json`.
