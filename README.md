# A midpoint-definition artefact reconciles the standard MWC model with single-cell chemotaxis dose responses

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23402856.svg)](https://doi.org/10.5281/zenodo.23402856)

> ## ⚠️ ERRATUM (v2.0.0, 2026-11-14)
>
> **The v1.0.x claim of this repository — that no standard MWC parameter set jointly reproduces the amplitude table and the midpoint distribution — is retracted.** It resulted from two definition mismatches in our own code (concentration axis F vs B+F; absolute vs relative half-maximum threshold). With consistent definitions a single parameter set reproduces both statistics (R² = 0.980, midpoint error 0.023 dex, not 0.622 dex). Full claim-by-claim record: [`ERRATUM.md`](ERRATUM.md). The refuted manuscript is preserved in `manuscript/_archive_refuted_v1/`; Zenodo record 10.5281/zenodo.23402857 hosts the superseded version.

Companion code, results and manuscript source for the corrected paper:

> H. Li. *A midpoint-definition artefact reconciles the standard MWC model with single-cell chemotaxis dose responses.* Companion manuscript (2026), split from the umbrella project [cell-communication-audit](https://github.com/lihuanlz/cell-communication-audit) (archived at Zenodo, doi:10.5281/zenodo.23008824).

## Layout

- `manuscript/` — corrected manuscript (main.tex/pdf, si.tex), figures fig1_v2.png / fig2_v3.png; the refuted v1 sources are archived in `manuscript/_archive_refuted_v1/`.
- `code/` — **corrected pipeline (entry points)**: `recheck_chemotaxis.py` (2×2 protocol crossing and Pareto fronts), `sensitivity_chemotaxis.py`, `fine_sweep.py`, `deep_N_diversity.py`, `deep_confound.py`, `fig_v2.py`, `fig2_v3.py` — all repo-relative, runnable as-is. The legacy v1 scripts (code49–code70 series) are kept for the historical record; eight of them carry hard-coded absolute paths from the original run environment.
- `results/corrected_v2/` — outputs of the corrected pipeline (`recheck_results.json`, `sens.json`, `deep.json`, `fig1_v2_data.json`). Other `results/` content is v1 lineage, superseded where it concerns the rupture analysis; the P7 adjudication tables and the adaptation-relay audit are unaffected.
- `data_external/moore2024_fret/` — the single-cell FRET archive of Moore et al. (2024), redistributed under its CC0 dedication; see `data_external/PROVENANCE.md`.

## Running

Python 3 (numpy/scipy/matplotlib). Pipeline entry points read the `.mat` files in `data_external/moore2024_fret/` directly.

## Ownership

Code, derived results and manuscript: the author's work product (MIT for code, CC BY 4.0 for text). The FRET archive is third-party data under CC0, mirrored here for self-containment with full attribution.

## Companion papers

- *The p53 pulse encoder counts by topological necessity* → [p53-pulse-encoder](https://github.com/lihuanlz/p53-pulse-encoder)
- *The quantitative currency of GPCR bias signalling fails its own algebraic identities* → [gpcr-field-audit](https://github.com/lihuanlz/gpcr-field-audit)
- *The price of absolute scale* (framework) → [scale-free-framework](https://github.com/lihuanlz/scale-free-framework)
