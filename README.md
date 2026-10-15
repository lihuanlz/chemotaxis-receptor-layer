# Amplitude and midpoint statistics of the same chemotaxis cells demand incompatible receptor models

Companion code, results and manuscript source for the paper:

> H. Li. *Amplitude and midpoint statistics of the same chemotaxis cells demand incompatible receptor models.* Companion manuscript (2026), split from the umbrella project [cell-communication-audit](https://github.com/lihuanlz/cell-communication-audit) (archived at Zenodo, doi:10.5281/zenodo.23008824).

## Layout

- `manuscript/` — LaTeX source, compiled PDF, figure PNG, figure-generation script.
- `code/` — analysis scripts (ours): P7 blind adjudication and erratum (code49/49b), receptor-layer self-check (code51chemo), protocol-artifact test (code52), CheA node mechanism (code53), joint-fit rupture quantification (code54), multi-species/local-adaptation tests (code55), FRET readout nonlinearity (code56), TCS depletion robustness (code57), P8 independent arms (code58), residual waterfall (code59), censoring test (code60), row-wise Hill steepness (code61), same-source statistics verdict (code62), prospective mixed-ligand arms (code70), and the adaptation-relay baseline audit.
- `results/` — derived outputs (ours), including the P7 adjudication tables and the prospective-arm CSVs.
- `data_external/moore2024_fret/` — the single-cell FRET archive of Moore et al. (2024), redistributed under its CC0 dedication; see `data_external/PROVENANCE.md`.

## Running

Python 3 (numpy/scipy/matplotlib). Pipeline entry points read the `.mat` files in `data_external/moore2024_fret/` directly.

## Ownership

Code, derived results and manuscript: the author's work product (MIT for code, CC BY 4.0 for text). The FRET archive is third-party data under CC0, mirrored here for self-containment with full attribution.

## Companion papers

- *The p53 pulse encoder counts by topological necessity* → [p53-pulse-encoder](https://github.com/lihuanlz/p53-pulse-encoder)
- *The quantitative currency of GPCR bias signalling fails its own algebraic identities* → [gpcr-field-audit](https://github.com/lihuanlz/gpcr-field-audit)
- *The price of absolute scale* (framework) → [scale-free-framework](https://github.com/lihuanlz/scale-free-framework)
