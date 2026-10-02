# 06 URTI DCSMI external validation

Final frozen URTI P1/P2 outputs are included under `results/full_frozen_outputs/07_URTI_external_validation`. The exact historical source that generated those final endpoint distributions was not present in any supplied archive. `DCSMI_URTI_endpoint_runner_EXECUTION_v1_0_SCAFFOLD.py` remains only an input-binding/hash-check scaffold.

v1.4 adds transparent reconstructions under `reconstructed/`:

- `reconstruct_DCSMI_URTI_endpoints_REPLAY_v1_0.py` — exact replay/audit from frozen components; validated against the frozen endpoint table.
- `run_DCSMI_URTI_endpoints_FUNCTIONAL_RECONSTRUCTION_v1_0.py` — full functional reimplementation based on the frozen DCSMI prespecification and psoriasis reference runner; historical equivalence is not claimed until the missing raw inputs are recovered.
- `reconstruct_URTI_partial_disease_D_from_frozen_outputs_v1_0.py` — forensic recovery of the identifiable 34/50 Hallmark disease-vector components.

URTI P3 is **NOT_IDENTIFIABLE / NOT_TESTED** under the frozen pre-result amendment.
