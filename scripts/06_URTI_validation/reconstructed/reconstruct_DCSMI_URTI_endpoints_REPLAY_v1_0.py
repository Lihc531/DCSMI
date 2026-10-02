#!/usr/bin/env python3
"""
Functional reconstruction / audit replay for the missing URTI endpoint generator.

STATUS
------
RECONSTRUCTED_FROM_FROZEN_INPUTS_AND_OUTPUTS; NOT THE HISTORICAL EXECUTION SCRIPT.

This script independently reconstructs all URTI endpoint quantities that are
identifiable from the supplied frozen assets:
  * structural-null A+ from the frozen WH/WC 374-target priors;
  * A+ nonzero counts;
  * syndrome-only Spearman baselines;
  * P1 observed Spearman correlations from frozen Mtilde and O;
  * P1 empirical upper-tail P values from the frozen 10,000-null distribution;
  * P2 delta and 95% percentile bootstrap CIs from the frozen 2,000-per-syndrome
    bootstrap distribution;
  * P3 status = NOT_IDENTIFIABLE / NOT_TESTED under the frozen pre-result
    two-syndrome amendment.

It does NOT regenerate the historical GSE63990-derived disease vector, the
historical P1 null distribution, or the historical patient bootstrap draws.
Those require source inputs/code that are not present in the supplied archive.

The structural A+ reconstruction follows DCSMI v1.0 exactly: the same donor
permutation is applied simultaneously to WH and WC in each replicate, within
the frozen 19 matched strata, N=10,000, seed=20260922.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

SEED = 20260922
N_NULL = 10000
SYNDROMES = ["WH", "WC"]


def sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_hallmarks(path: str | Path) -> dict[str, set[str]]:
    df = pd.read_csv(path, sep="\t")
    required = {"hallmark", "genes"}
    if not required.issubset(df.columns):
        raise ValueError(f"Hallmark file must contain {sorted(required)}")
    return {
        str(r.hallmark): set(str(r.genes).split(";"))
        for _, r in df.iterrows()
    }


def reconstruct_aplus(prior_path: str | Path, strata_path: str | Path,
                      hallmark_path: str | Path, hallmark_names: list[str]):
    prior = pd.read_csv(prior_path)
    strata = pd.read_csv(strata_path).reset_index(drop=True)
    if len(strata) != 374 or strata["gene"].nunique() != 374:
        raise ValueError("Frozen target universe must contain exactly 374 unique genes")
    if strata["matched_stratum"].nunique() != 19:
        raise ValueError("Frozen structural null must contain exactly 19 strata")

    genes = strata["gene"].astype(str).tolist()
    hallmarks = load_hallmarks(hallmark_path)
    missing_h = [h for h in hallmark_names if h not in hallmarks]
    if missing_h:
        raise ValueError(f"Hallmark definitions missing {missing_h[:5]}")

    inc = np.zeros((374, len(hallmark_names)), dtype=float)
    for j, h in enumerate(hallmark_names):
        hs = hallmarks[h]
        inc[:, j] = [g in hs for g in genes]

    T = np.zeros((2, 374), dtype=float)
    for i, s in enumerate(SYNDROMES):
        x = prior.loc[prior["syndrome_id"].astype(str) == s]
        if x.empty:
            raise ValueError(f"Missing syndrome {s} in prior")
        z = x.set_index("target_id")["target_weight"]
        T[i] = np.array([z.get(g, 0.0) for g in genes], dtype=float)
        if not np.isclose(T[i].sum(), 1.0, atol=1e-8):
            raise ValueError(f"Target-prior mass conservation failed for {s}")

    groups = [
        np.where(strata["matched_stratum"].values == s)[0]
        for s in sorted(strata["matched_stratum"].unique())
    ]
    qobs = T @ inc
    qnull = np.empty((N_NULL, 2, len(hallmark_names)), dtype=float)
    rng = np.random.default_rng(SEED)
    for b in range(N_NULL):
        donor = np.arange(374)
        for ix in groups:
            donor[ix] = rng.permutation(ix)
        qnull[b] = T[:, donor] @ inc

    mu = qnull.mean(axis=0)
    sd = qnull.std(axis=0, ddof=1)
    A = np.divide(qobs - mu, sd, out=np.zeros_like(qobs), where=sd > 0)
    Aplus = np.maximum(A, 0.0)
    Anorm = np.divide(
        Aplus,
        Aplus.sum(axis=1, keepdims=True),
        out=np.zeros_like(Aplus),
        where=Aplus.sum(axis=1, keepdims=True) > 0,
    )
    return Aplus, Anorm


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--prior", required=True)
    p.add_argument("--strata", required=True)
    p.add_argument("--hallmark", required=True)
    p.add_argument("--mtilde", required=True)
    p.add_argument("--observed-o", required=True)
    p.add_argument("--p1-null", required=True)
    p.add_argument("--p2-bootstrap", required=True)
    p.add_argument("--reference-endpoints", default=None,
                   help="Optional frozen raw endpoint table used only for QA comparison; values are never used to compute reconstructed statistics.")
    p.add_argument("--outdir", required=True)
    a = p.parse_args()

    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    M = pd.read_csv(a.mtilde, index_col=0)
    O = pd.read_csv(a.observed_o, index_col=0)
    if list(M.index.astype(str)) != SYNDROMES or list(O.index.astype(str)) != SYNDROMES:
        raise ValueError("Expected WH/WC rows in Mtilde and O")
    if list(M.columns) != list(O.columns):
        raise ValueError("Mtilde and O Hallmark columns do not match")
    hallmark_names = list(M.columns)

    p1null = pd.read_csv(a.p1_null)
    if len(p1null) != N_NULL or not set(SYNDROMES).issubset(p1null.columns):
        raise ValueError("Frozen P1 null must contain 10,000 rows and WH/WC columns")
    boot = pd.read_csv(a.p2_bootstrap)
    if not {"syndrome", "delta_boot"}.issubset(boot.columns):
        raise ValueError("P2 bootstrap file must contain syndrome, delta_boot")
    for s in SYNDROMES:
        if (boot["syndrome"].astype(str) == s).sum() != 2000:
            raise ValueError(f"Expected exactly 2,000 P2 bootstrap replicates for {s}")

    Aplus, Anorm = reconstruct_aplus(a.prior, a.strata, a.hallmark, hallmark_names)

    rows = []
    for i, s in enumerate(SYNDROMES):
        p1_obs = float(spearmanr(M.loc[s].values, O.loc[s].values).statistic)
        p1_emp = float((1 + (p1null[s].to_numpy(float) >= p1_obs).sum()) / (N_NULL + 1))
        syndrome_only = float(spearmanr(Anorm[i], O.loc[s].values).statistic)
        delta = p1_obs - syndrome_only
        b = boot.loc[boot["syndrome"].astype(str) == s, "delta_boot"].to_numpy(float)
        ci_lo, ci_hi = np.quantile(b, [0.025, 0.975])
        rows.append({
            "syndrome": s,
            "n_hallmarks": len(hallmark_names),
            "Aplus_nonzero": int((Aplus[i] > 0).sum()),
            "P1_spearman": p1_obs,
            "P1_empirical_p": p1_emp,
            "syndrome_only_spearman": syndrome_only,
            "P2_delta": float(delta),
            "P2_boot_CI_low": float(ci_lo),
            "P2_boot_CI_high": float(ci_hi),
            "P2_support": bool(ci_lo > 0),
            "P3_status": "NOT_IDENTIFIABLE_NOT_TESTED",
            "P3_selectivity": np.nan,
            "P3_empirical_p": np.nan,
        })

    rec = pd.DataFrame(rows)
    rec.to_csv(out / "DCSMI_URTI_primary_endpoints_RECONSTRUCTED_REPLAY.csv", index=False)
    pd.DataFrame(Aplus, index=SYNDROMES, columns=hallmark_names).to_csv(
        out / "DCSMI_URTI_Aplus_RECONSTRUCTED.csv"
    )

    qa = []
    if a.reference_endpoints:
        ref = pd.read_csv(a.reference_endpoints).set_index("syndrome")
        for _, r in rec.iterrows():
            s = r["syndrome"]
            for col in [
                "n_hallmarks", "Aplus_nonzero", "P1_spearman", "P1_empirical_p",
                "syndrome_only_spearman", "P2_delta", "P2_boot_CI_low",
                "P2_boot_CI_high"
            ]:
                rv = float(r[col])
                fv = float(ref.loc[s, col])
                qa.append({
                    "syndrome": s,
                    "field": col,
                    "reconstructed": rv,
                    "frozen_reference": fv,
                    "absolute_difference": abs(rv - fv),
                    "pass_1e-10": abs(rv - fv) <= 1e-10,
                })
        pd.DataFrame(qa).to_csv(out / "RECONSTRUCTION_QA_vs_frozen_endpoints.csv", index=False)

    inputs = {
        "prior": a.prior, "strata": a.strata, "hallmark": a.hallmark,
        "mtilde": a.mtilde, "observed_o": a.observed_o,
        "p1_null": a.p1_null, "p2_bootstrap": a.p2_bootstrap,
    }
    if a.reference_endpoints:
        inputs["reference_endpoints_QA_only"] = a.reference_endpoints
    manifest = {
        "status": "RECONSTRUCTED_FROM_FROZEN_INPUTS_AND_OUTPUTS_NOT_HISTORICAL_EXECUTION_SCRIPT",
        "seed": SEED,
        "N_structural_null_for_Aplus_reconstruction": N_NULL,
        "P1_null_source": "frozen distribution supplied as input; not regenerated",
        "P2_bootstrap_source": "frozen distribution supplied as input; not regenerated",
        "P3": "NOT_IDENTIFIABLE_NOT_TESTED per frozen pre-result amendment",
        "inputs_sha256": {k: sha256(v) for k, v in inputs.items()},
        "qa_all_pass_1e-10": bool(all(x["pass_1e-10"] for x in qa)) if qa else None,
    }
    (out / "RECONSTRUCTION_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(rec.to_string(index=False))
    if qa:
        print("QA all pass at 1e-10:", manifest["qa_all_pass_1e-10"])


if __name__ == "__main__":
    main()
