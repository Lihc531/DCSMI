# Reconstructed URTI endpoint code

These files are **functional reconstructions**, not recovered historical source code.

The supplied archive contains frozen URTI endpoint outputs but not the original script that generated the final GSE63990-conditioned P1 null distribution and patient-level P2 bootstrap distribution. The reconstruction is therefore separated from `provenance/original_frozen_scripts/` and is never represented as the historical execution source.

## 1. Exact replay/audit reconstruction

`reconstruct_DCSMI_URTI_endpoints_REPLAY_v1_0.py`

Uses the frozen WH/WC target prior, 374-target/19-stratum null, Hallmark definitions, frozen `Mtilde`, frozen observed `O`, frozen 10,000-row P1 null, and frozen 2,000-per-syndrome P2 bootstrap distribution. It independently reconstructs A+, the syndrome-only baselines, P1 observed statistics/empirical P values, P2 deltas/CIs, and applies the frozen P3 amendment.

Internal validation against the supplied frozen endpoint table reproduces every independently reconstructable numeric field at absolute error <= 1e-10. See `provenance/reconstructed_URTI_validation/RECONSTRUCTION_QA_vs_frozen_endpoints.csv`.

## 2. Full functional reimplementation

`run_DCSMI_URTI_endpoints_FUNCTIONAL_RECONSTRUCTION_v1_0.py`

Implements the DCSMI v1.0 equations/null/bootstrap logic using the frozen prespecification and the supplied psoriasis endpoint runner as the shared reference implementation. It can generate `Mtilde`, a 10,000-replicate structural P1 null, and the 2,000-replicate patient bootstrap when the exact frozen GSE63990 disease input and syndrome-labelled patient expression matrix are available.

The exact historical GSE63990 input and original patient-level source matrix are not present in the supplied archive, so this full mode has been syntax/smoke tested but cannot currently be claimed to reproduce the historical null/bootstrap bytes.

## 3. Forensic partial disease-vector reconstruction

`reconstruct_URTI_partial_disease_D_from_frozen_outputs_v1_0.py`

Uses frozen `Mtilde`, reconstructed A+, and frozen compatibility R to identify the disease Hallmark vector `D` wherever at least one observed syndrome has A+ > 0. It recovers 34/50 Hallmarks exactly across WH/WC consistency checks; 16/50 Hallmarks remain mathematically unidentifiable from the supplied outputs because both observed A+ profiles are zero there. This is why the historical P1-null distribution cannot be regenerated exactly from the endpoint outputs alone.

## P3

URTI P3 remains **NOT_IDENTIFIABLE / NOT_TESTED** under the frozen pre-result amendment. The reconstructed code intentionally does not replace it with a new confirmatory endpoint.
