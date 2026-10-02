# URTI reconstruction status

## What is now reconstructable exactly from supplied frozen assets

The replay reconstruction independently reproduces, to absolute error <= 1e-10:

- eligible Hallmark count (50);
- observed A+ nonzero counts (WH 32, WC 24);
- P1 observed Spearman correlation;
- P1 empirical upper-tail P value from the supplied 10,000-row frozen null;
- syndrome-only Spearman baseline;
- P2 disease-conditioning delta;
- P2 95% percentile bootstrap CI from the supplied 2,000-per-syndrome bootstrap distribution.

URTI P3 is not reconstructed as a numeric test because the frozen pre-result amendment designates it NOT_IDENTIFIABLE / NOT_TESTED.

## What cannot be independently regenerated from currently supplied raw inputs

The lost historical generator and two source inputs are absent:

1. the exact frozen GSE63990 gene-level disease-background input used to construct the full 50-Hallmark disease vector D;
2. the original syndrome-labelled patient-level expression/count matrix used to generate the 2,000 patient bootstrap replicates (the observed-O freeze manifest records the source `data.zip` SHA256, but the bytes were not supplied).

From frozen Mtilde, A+, and compatibility R, 34/50 components of D are identifiable exactly. Sixteen Hallmarks have A+=0 for both WH and WC and therefore cannot be recovered from these outputs. Because structural-null permutations can produce A+>0 in those Hallmarks, the historical P1 null distribution cannot be regenerated exactly without the missing full D vector.

## Code classification

- `reconstruct_DCSMI_URTI_endpoints_REPLAY_v1_0.py`: validated audit/replay reconstruction.
- `run_DCSMI_URTI_endpoints_FUNCTIONAL_RECONSTRUCTION_v1_0.py`: full algorithmic reimplementation for future use if missing inputs are recovered; not historical source.
- `reconstruct_URTI_partial_disease_D_from_frozen_outputs_v1_0.py`: forensic identifiability utility; not historical source.
