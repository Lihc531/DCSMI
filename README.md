# DCSMI evidence-calibration repository v1.4.1

Code, frozen outputs, prespecification documents and reproducibility metadata for the manuscript:

**Structural calibration and evidence boundary analysis reveal convergent molecular space with residual syndrome differentiation in treatment knowledge derived inference**

## Repository status

This package is designed to be **accurate about what is and is not reproducible from the supplied files**. It does not relabel scaffolds as completed endpoint code and it does not silently replace frozen inputs.

| Manuscript layer | Supplied status |
|---|---|
| Treatment / target projection | Frozen code + inputs + outputs; clean-rerun key-table equivalence checked |
| Exact herb-target topology null | Frozen code + prespecification + outputs present |
| GO/KEGG functional calibration | Executed code + outputs present; third-party annotation snapshots not redistributed here |
| CAP disease correspondence | Frozen reconstruction/matched-null code + derived cohort inputs + outputs present |
| Psoriasis DCSMI validation | Code-frozen endpoint runner + outputs present; executed primary disease background = GSE147339 Tier-B |
| URTI DCSMI validation | Frozen outputs + exact replay/audit reconstruction + full functional reimplementation are provided. The lost historical generator and exact GSE63990/patient-level raw inputs remain unavailable, so historical null/bootstrap byte-for-byte regeneration is not claimed |
| Cross-layer contraction | Executed code/output present; portable adapter explicitly restores the functional-eligibility dependency |

## Important execution facts

- Master structural/null seed: `20260922` where specified by the frozen analyses.
- Phase 12A: 10,000 Curveball null graphs, 4 chains, exact row/column degree preservation.
- Functional calibration: 10,000 within-stratum structural permutations; 374-target universe; 19 matched strata.
- Psoriasis DCSMI: 10,000 structural permutations; 2,000 patient bootstrap resamples.
- STRING RWR secondary analysis: primary restart `r=0.5`; sensitivity `0.3/0.7`; `max_iter=10000`; L1 tolerance `1e-12`.
- URTI P3 is `NOT_IDENTIFIABLE / NOT_TESTED` under a frozen pre-result amendment.

## Original vs portable code

Exact frozen/executed source bytes are preserved under `provenance/original_frozen_scripts/`. `scripts/*PORTABLE.py` are convenience adaptations that change file locations (and, for cross-layer analysis, make an implicit functional-eligibility dependency explicit) without changing the statistical definitions. Every adaptation is listed with original/runnable SHA256 values in `provenance/PORTABILITY_ADAPTATIONS.csv`.

## Data redistribution

Author-derived frozen/derived inputs needed for audit are bundled under `data/frozen_inputs/`. Third-party database snapshots and licensed resources are not redistributed in this GitHub package; exact filenames and known hashes are recorded in `data/manifests/EXPECTED_EXTERNAL_INPUTS.csv`. The complete internal basic archive preserves the supplied snapshots separately.

## Reproducibility boundary

The repository supports audit/reproduction of the current frozen respiratory/CAP/psoriasis workflow to the extent permitted by the supplied source files and external-resource availability. The exact historical URTI endpoint-generating source remains unavailable. v1.4 adds clearly labeled reconstructed code: an exact replay/audit script that reproduces all independently identifiable endpoint statistics from frozen assets, a full functional reimplementation for use if the missing raw inputs are recovered, and a forensic partial-D reconstruction. None is presented as the lost historical source.

See `RUNNING_THE_REPOSITORY.md`, `KNOWN_LIMITATIONS.md`, and `MANUSCRIPT_ALIGNMENT_NOTES.md` before public release.
