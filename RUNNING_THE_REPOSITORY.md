# Running the repository

## 1. Create an environment

```bash
conda env create -f environment.yml
conda activate dcsmi-evidence-calibration
```

The environment file is a portable reconstruction; the exact historical package lockfile was not supplied.

## 2. Set the repository root

```bash
export DCSMI_PROJECT_ROOT=/absolute/path/to/DCSMI_GitHub_Repository_v1_4
```

On Windows PowerShell:

```powershell
$env:DCSMI_PROJECT_ROOT = "C:\path\to\DCSMI_GitHub_Repository_v1_4"
```

## 3. Verify packaged files

```bash
python scripts/utilities/check_repository_integrity.py
```

## 4. Target projection

```bash
python scripts/01_target_projection/weighted_herb_target_projection_v1_1_PORTABLE.py
```

During package construction, the bundled inputs reproduced the frozen key target-projection tables.

## 5. Phase 12A topology null

Point the topology runner to either the newly generated projection output or the bundled frozen output:

```bash
export PHASE12A_INPUT_DIR="$DCSMI_PROJECT_ROOT/results/full_frozen_outputs/01_target_projection"
export PHASE12A_OUTPUT_DIR="$DCSMI_PROJECT_ROOT/outputs/02_topology_null"
python scripts/02_topology_null/run_degree_preserving_herb_target_convergence_null_v1_0_PORTABLE.py
```

## 6. Functional and cross-layer analyses

Supply the external GO/KEGG snapshots listed in `data/manifests/EXPECTED_EXTERNAL_INPUTS.csv`, either under `external_data/` or by setting the corresponding environment variables. Then run the portable functional and cross-layer scripts.

## 7. CAP

The bundled GSE103119/GSE196399 derived gene-level tables support CAP meta-Z reconstruction and downstream matched-null analyses. The secondary STRING RWR additionally requires the external/bundled-in-internal-archive STRING edge file listed in the input manifest.

## 8. Psoriasis

The exact endpoint runner requires the original GSE192867 series matrix and frozen Hallmark universe in addition to the bundled prior, annotation, sample manifest, strata and GSE147339 disease-background statistics. Use the runner's command-line arguments; do not substitute the two-cohort GSE147339+GSE61281 stress-test background for the executed primary input.

## 9. URTI

The original supplied `DCSMI_URTI_endpoint_runner_EXECUTION_v1_0_SCAFFOLD.py` remains a scaffold and must not be relabeled as the historical endpoint generator. v1.4 adds reconstructed code under `scripts/06_URTI_validation/reconstructed/`.

For an exact audit/replay of the supplied endpoint statistics, run `reconstruct_DCSMI_URTI_endpoints_REPLAY_v1_0.py` with the frozen Mtilde/O/P1-null/P2-bootstrap assets and the frozen Hallmark resource. This reproduces every independently reconstructable endpoint field within 1e-10 of the frozen table.

`run_DCSMI_URTI_endpoints_FUNCTIONAL_RECONSTRUCTION_v1_0.py` is a full functional reimplementation of the prespecified P1/P2 workflow. Exact historical regeneration still requires the missing frozen GSE63990 disease input and patient-level expression matrix. Do not claim that this reconstructed runner is the lost historical script.

URTI P3 remains `NOT_IDENTIFIABLE / NOT_TESTED`.

Example replay command (after placing the exact Hallmark TSV listed in `data/manifests/EXPECTED_EXTERNAL_INPUTS.csv` under `external_data/`):

```bash
python scripts/06_URTI_validation/reconstructed/reconstruct_DCSMI_URTI_endpoints_REPLAY_v1_0.py \
  --prior data/frozen_inputs/URTI/WindHeat_WindCold_target_prior_374_FROZEN.csv \
  --strata data/frozen_inputs/treatment/matched_null_target_universe_and_strata_v1_0_1.csv \
  --hallmark external_data/Hallmark_2026_1_Hs_FROZEN.tsv \
  --mtilde results/full_frozen_outputs/07_URTI_external_validation/DCSMI_URTI_Mtilde.csv \
  --observed-o results/full_frozen_outputs/07_URTI_external_validation/DCSMI_URTI_O.csv \
  --p1-null results/full_frozen_outputs/07_URTI_external_validation/DCSMI_URTI_P1_null_distribution.csv \
  --p2-bootstrap results/full_frozen_outputs/07_URTI_external_validation/DCSMI_URTI_P2_bootstrap_distribution.csv \
  --reference-endpoints results/full_frozen_outputs/07_URTI_external_validation/DCSMI_URTI_primary_endpoints_rawP.csv \
  --outdir outputs/06_URTI_reconstructed_replay
```

## 10. Global multiplicity reconstruction

The transparent utility below reconstructs the frozen global P1 Holm family across psoriasis + URTI and the psoriasis-only P3 family after the URTI P3 amendment:

```bash
python scripts/utilities/reconstruct_global_holm_from_frozen_raw_p.py \
  --psoriasis results/full_frozen_outputs/06_psoriasis_external_validation/DCSMI_psoriasis_primary_endpoints.csv \
  --urti results/full_frozen_outputs/07_URTI_external_validation/DCSMI_URTI_primary_endpoints_rawP.csv \
  --out outputs/multiplicity
```

This utility is reconstructed from the prespecified multiplicity rule and is clearly separated from the historical frozen endpoint code.
