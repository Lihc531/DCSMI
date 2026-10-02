#!/usr/bin/env python3
"""
Functional reconstruction of the missing DCSMI URTI endpoint generator.

STATUS
------
RECONSTRUCTED_FROM_FROZEN_PRESPECIFICATION_AND_PSORIASIS_REFERENCE_RUNNER.
THIS IS NOT THE LOST HISTORICAL EXECUTION SCRIPT.

Purpose
-------
Implement the prespecified URTI P1/P2 workflow so it can be rerun if the exact
frozen GSE63990 disease input and the syndrome-labelled patient expression
matrix are recovered. The equations/nulls follow DCSMI v1.0 and mirror the
supplied psoriasis endpoint runner where the design is shared.

Outputs
-------
* DCSMI_URTI_Mtilde_RECONSTRUCTED.csv
* DCSMI_URTI_O_RECONSTRUCTED.csv
* DCSMI_URTI_P1_null_distribution_RECONSTRUCTED.csv
* DCSMI_URTI_P2_bootstrap_distribution_RECONSTRUCTED.csv
* DCSMI_URTI_primary_endpoints_RECONSTRUCTED.csv
* RECONSTRUCTED_ENDPOINT_RUN_MANIFEST.json

P3 is deliberately NOT calculated. Under the frozen pre-result amendment,
P3 is NOT_IDENTIFIABLE / NOT_TESTED for the two-syndrome direction-agnostic
|Z| observed representation.

Input assumptions
-----------------
1) prior: syndrome_id,target_id,target_weight for WH/WC, 374 targets each.
2) strata: frozen 374-target table with columns gene, matched_stratum.
3) hallmark: frozen Hallmark TSV with columns hallmark, genes.
4) disease: gene-level GSE63990 disease-vs-comparator statistics. Supported:
   gene + signed_Z, OR gene + logFC + pvalue. Gene identifiers must match
   Hallmark symbols unless --disease-gene-map is supplied.
5) For O/P2, either:
   a) --observed-o plus --p2-bootstrap-replay (audit/replay of frozen patient
      evaluation), OR
   b) --patient-expression plus --sample-manifest to recompute O and 2,000
      stratified patient bootstraps.

Because the exact historical GSE63990 input and raw patient expression file are
not present in the supplied archive, exact historical regeneration cannot be
validated until those inputs are recovered.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import stdtr
from scipy.stats import norm, rankdata, spearmanr

SEED = 20260922
N_NULL = 10000
N_BOOT = 2000
SYND = ["WH", "WC"]


def sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def signed_z(fc, p):
    p = np.clip(np.asarray(p, float), 1e-300, 1.0)
    return np.sign(np.asarray(fc, float)) * norm.ppf(1 - p / 2)


def load_hallmarks(path):
    df = pd.read_csv(path, sep="\t")
    if not {"hallmark", "genes"}.issubset(df.columns):
        raise ValueError("Hallmark TSV requires hallmark, genes columns")
    return {str(r.hallmark): set(str(r.genes).split(";")) for _, r in df.iterrows()}


def load_gene_map(path):
    if path is None:
        return None
    df = pd.read_csv(path)
    if df.shape[1] < 2:
        raise ValueError("Gene map needs at least two columns")
    return df.iloc[:, :2].rename(columns={df.columns[0]: "source_id", df.columns[1]: "gene"})


def map_gene_table(df, gene_col, gene_map=None):
    x = df.copy()
    x[gene_col] = x[gene_col].astype(str).str.replace(r"\.\d+$", "", regex=True)
    if gene_map is None:
        x["gene_symbol"] = x[gene_col]
    else:
        mp = gene_map.copy()
        mp["source_id"] = mp["source_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)
        x = x.merge(mp, left_on=gene_col, right_on="source_id", how="left")
        x["gene_symbol"] = x["gene"].astype("string")
        x = x[x["gene_symbol"].notna() & (x["gene_symbol"].str.len() > 0)].copy()
    return x


def percentile_module_from_z(z, genes, hall, eligible_names=None):
    ok = np.isfinite(z)
    zz = np.asarray(z, float)[ok]
    gg = np.asarray(genes, dtype=str)[ok]
    # Collapse duplicate symbols before ranking, using the largest |Z| as a
    # deterministic representative. For exact historical replication, the
    # recovered frozen disease input should ideally already contain unique
    # symbols, making this branch irrelevant.
    tmp = pd.DataFrame({"gene": gg, "z": zz})
    tmp["absz"] = tmp["z"].abs()
    tmp = tmp.sort_values(["gene", "absz"], ascending=[True, False]).drop_duplicates("gene")
    gg = tmp["gene"].to_numpy(str)
    zz = tmp["z"].to_numpy(float)
    rr = rankdata(np.abs(zz), method="average") / (len(zz) + 1)
    mp = dict(zip(gg, rr))
    names = list(hall) if eligible_names is None else list(eligible_names)
    raw, n = [], []
    for k in names:
        vals = [mp[g] for g in hall[k] if g in mp]
        raw.append(np.mean(vals) if vals else np.nan)
        n.append(len(vals))
    raw = np.asarray(raw, float)
    n = np.asarray(n, int)
    elig = n >= 15
    if eligible_names is None:
        if elig.sum() < 40:
            raise ValueError("fewer than 40 eligible Hallmarks")
        names = np.asarray(names)[elig]
        raw = raw[elig]
        n = n[elig]
    else:
        if not np.all(elig):
            bad = np.asarray(names)[~elig]
            raise ValueError(f"Frozen eligible Hallmark fell below 15 measured genes: {bad[:5]}")
        names = np.asarray(names)
    v = raw / np.sqrt(np.sum(raw ** 2))
    return names, v, n


def welch_fast(A, B):
    n1, n2 = A.shape[1], B.shape[1]
    m1, m2 = np.nanmean(A, axis=1), np.nanmean(B, axis=1)
    v1, v2 = np.nanvar(A, axis=1, ddof=1), np.nanvar(B, axis=1, ddof=1)
    a, b = v1 / n1, v2 / n2
    se = np.sqrt(a + b)
    t = np.divide(m1 - m2, se, out=np.zeros_like(m1), where=se > 0)
    den = a * a / (n1 - 1) + b * b / (n2 - 1)
    df = np.divide((a + b) ** 2, den, out=np.full_like(m1, np.inf), where=den > 0)
    p = 2 * stdtr(df, -np.abs(t))
    p = np.where(se == 0, 1.0, p)
    return m1 - m2, p


def read_expression(path, sep=None):
    if sep is None:
        sep = "\t" if str(path).lower().endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(path, sep=sep)
    if df.shape[1] < 3:
        raise ValueError("Patient expression matrix must contain gene column plus sample columns")
    return df


def normalize_counts_log2cpm(X):
    lib = np.nansum(X, axis=0)
    if np.any(lib <= 0):
        raise ValueError("Nonpositive library size in raw counts")
    cpm = X / lib[None, :] * 1e6
    return np.log2(cpm + 1.0)


def compute_observed_and_bootstrap(expr_df, sample_manifest, hall, hallmark_names,
                                   gene_map, expression_scale, Mtilde, Anorm, rng):
    gene_col = expr_df.columns[0]
    x = map_gene_table(expr_df, gene_col, gene_map)
    sample_cols = [c for c in expr_df.columns[1:] if c in set(sample_manifest["sample_id"].astype(str))]
    if len(sample_cols) != len(sample_manifest):
        missing = sorted(set(sample_manifest["sample_id"].astype(str)) - set(sample_cols))
        raise ValueError(f"Expression/sample-manifest mismatch; missing samples: {missing[:5]}")
    vals = x[sample_cols].apply(pd.to_numeric, errors="coerce")
    x2 = pd.concat([x[["gene_symbol"]].reset_index(drop=True), vals.reset_index(drop=True)], axis=1)
    # collapse duplicate mapped genes by median, as in the supplied psoriasis runner
    gene = x2.groupby("gene_symbol").median(numeric_only=True)
    gene = gene.loc[:, sample_manifest["sample_id"].astype(str).tolist()]
    G = gene.to_numpy(float)
    if expression_scale == "raw_counts":
        G = normalize_counts_log2cpm(G)
        gene = pd.DataFrame(G, index=gene.index, columns=gene.columns)
    elif expression_scale != "log2cpm":
        raise ValueError("--expression-scale must be raw_counts or log2cpm")

    groups = sample_manifest.set_index("sample_id").loc[gene.columns, "group"].astype(str)
    if set(groups.unique()) != {"WH", "WC"}:
        raise ValueError("URTI patient manifest must contain exactly WH and WC groups")
    pos = {s: np.where(groups.values == s)[0] for s in SYND}

    # Observed O: WH-vs-WC direction-agnostic profile. The inverse WC-vs-WH
    # yields the identical |Z| profile by construction.
    fc, pv = welch_fast(G[:, pos["WH"]], G[:, pos["WC"]])
    _, oval, _ = percentile_module_from_z(signed_z(fc, pv), gene.index.astype(str).values,
                                          hall, hallmark_names)
    O = np.vstack([oval, oval])

    # Freeze bootstrap draws for each group before endpoint calculations.
    draws = {}
    for s in SYND:
        n = len(pos[s])
        C = np.zeros((N_BOOT, n), dtype=np.int16)
        for b in range(N_BOOT):
            ix = rng.choice(n, size=n, replace=True)
            C[b] = np.bincount(ix, minlength=n)
        draws[s] = C

    ng = G.shape[0]
    Hmat = np.zeros((ng, len(hallmark_names)), float)
    for j, h in enumerate(hallmark_names):
        hs = hall[h]
        Hmat[:, j] = [g in hs for g in gene.index.astype(str)]
    Hcount = Hmat.sum(axis=0)
    if np.any(Hcount < 15):
        raise ValueError("bootstrap Hallmark eligibility invariant failed")

    rM = np.array([rankdata(Mtilde[i], method="average") for i in range(2)])
    rA = np.array([rankdata(Anorm[i], method="average") for i in range(2)])

    def rowcorr_fixed(rfixed, Y):
        ry = rankdata(Y, axis=1, method="average")
        xc = rfixed - rfixed.mean()
        yc = ry - ry.mean(axis=1, keepdims=True)
        return (yc @ xc) / (np.sqrt((yc * yc).sum(axis=1)) * np.sqrt((xc * xc).sum()))

    boot = np.empty((N_BOOT, 2), float)
    batch = 50
    # Same resampled WH/WC patient sets define O for both syndrome endpoints.
    for lo in range(0, N_BOOT, batch):
        hi = min(N_BOOT, lo + batch)
        Cw, Cc = draws["WH"][lo:hi], draws["WC"][lo:hi]
        n1, n2 = len(pos["WH"]), len(pos["WC"])
        S1 = Cw @ G[:, pos["WH"]].T
        SS1 = Cw @ (G[:, pos["WH"]] ** 2).T
        S2 = Cc @ G[:, pos["WC"]].T
        SS2 = Cc @ (G[:, pos["WC"]] ** 2).T
        m1, m2 = S1 / n1, S2 / n2
        v1 = np.maximum((SS1 - S1 * S1 / n1) / (n1 - 1), 0)
        v2 = np.maximum((SS2 - S2 * S2 / n2) / (n2 - 1), 0)
        aa, bb = v1 / n1, v2 / n2
        se = np.sqrt(aa + bb)
        tt = np.divide(m1 - m2, se, out=np.zeros_like(m1), where=se > 0)
        den = aa * aa / (n1 - 1) + bb * bb / (n2 - 1)
        df = np.divide((aa + bb) ** 2, den, out=np.full_like(m1, np.inf), where=den > 0)
        pv = 2 * stdtr(df, -np.abs(tt))
        pv = np.where(se == 0, 1.0, pv)
        zz = signed_z(m1 - m2, pv)
        rr = rankdata(np.abs(zz), axis=1, method="average") / (ng + 1)
        raw = (rr @ Hmat) / Hcount
        ob = raw / np.sqrt((raw * raw).sum(axis=1, keepdims=True))
        for i in range(2):
            boot[lo:hi, i] = rowcorr_fixed(rM[i], ob) - rowcorr_fixed(rA[i], ob)
    return O, boot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prior", required=True)
    ap.add_argument("--strata", required=True)
    ap.add_argument("--hallmark", required=True)
    ap.add_argument("--disease", required=True)
    ap.add_argument("--disease-gene-map", default=None)
    ap.add_argument("--observed-o", default=None)
    ap.add_argument("--patient-expression", default=None)
    ap.add_argument("--patient-gene-map", default=None)
    ap.add_argument("--sample-manifest", default=None)
    ap.add_argument("--expression-scale", choices=["raw_counts", "log2cpm"], default="log2cpm")
    ap.add_argument("--p2-bootstrap-replay", default=None,
                    help="Frozen P2 distribution can be supplied when patient-level expression is unavailable.")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    hall = load_hallmarks(a.hallmark)
    strata = pd.read_csv(a.strata).reset_index(drop=True)
    if len(strata) != 374 or strata["gene"].nunique() != 374 or strata["matched_stratum"].nunique() != 19:
        raise ValueError("Frozen 374-target / 19-stratum invariant failed")
    genes374 = strata["gene"].astype(str).tolist()
    stratum_groups = [np.where(strata["matched_stratum"].values == s)[0]
                      for s in sorted(strata["matched_stratum"].unique())]

    prior = pd.read_csv(a.prior)
    T = np.zeros((2, 374), float)
    for i, s in enumerate(SYND):
        z = prior.loc[prior["syndrome_id"].astype(str) == s].set_index("target_id")["target_weight"]
        T[i] = [z.get(g, 0.0) for g in genes374]
        if not np.isclose(T[i].sum(), 1.0, atol=1e-8):
            raise ValueError(f"prior weight conservation failed for {s}")

    # Disease vector D from the recovered/frozen GSE63990 gene-level input.
    dis = pd.read_csv(a.disease)
    gene_col = "gene" if "gene" in dis.columns else ("gene_id" if "gene_id" in dis.columns else dis.columns[0])
    dmap = load_gene_map(a.disease_gene_map)
    dis = map_gene_table(dis, gene_col, dmap)
    if "signed_Z" in dis.columns:
        zd = pd.to_numeric(dis["signed_Z"], errors="coerce").to_numpy(float)
    elif {"logFC", "pvalue"}.issubset(dis.columns):
        zd = signed_z(dis["logFC"], dis["pvalue"])
    else:
        raise ValueError("Disease input requires signed_Z or logFC+pvalue")
    dnames, D, _ = percentile_module_from_z(zd, dis["gene_symbol"].astype(str).values, hall)
    hallmark_names = list(dnames)
    K = len(hallmark_names)

    # 374-target x eligible-Hallmark incidence.
    Inc = np.zeros((374, K), float)
    for j, h in enumerate(hallmark_names):
        hs = hall[h]
        Inc[:, j] = [g in hs for g in genes374]

    qobs = T @ Inc
    qnull = np.empty((N_NULL, 2, K), float)
    for b in range(N_NULL):
        donor = np.arange(374)
        for ix in stratum_groups:
            donor[ix] = rng.permutation(ix)
        qnull[b] = T[:, donor] @ Inc
    mu, sd = qnull.mean(axis=0), qnull.std(axis=0, ddof=1)
    A = np.divide(qobs - mu, sd, out=np.zeros_like(qobs), where=sd > 0)
    Aplus = np.maximum(A, 0.0)
    nonzero = (Aplus > 0).sum(axis=1)
    if np.any(nonzero < 10):
        raise SystemExit("STOP: Aplus has fewer than 10 nonzero eligible Hallmarks")
    Anorm = Aplus / Aplus.sum(axis=1, keepdims=True)
    M = D[None, :] * Aplus
    Mtilde = M / M.sum(axis=1, keepdims=True)

    # Observed O and/or P2 bootstrap.
    boot = None
    if a.patient_expression:
        if not a.sample_manifest:
            raise ValueError("--sample-manifest is required with --patient-expression")
        expr = read_expression(a.patient_expression)
        sm = pd.read_csv(a.sample_manifest)
        if not {"sample_id", "group"}.issubset(sm.columns):
            raise ValueError("Sample manifest requires sample_id, group")
        pmap = load_gene_map(a.patient_gene_map)
        O, boot = compute_observed_and_bootstrap(expr, sm, hall, hallmark_names,
                                                 pmap, a.expression_scale, Mtilde, Anorm, rng)
    elif a.observed_o:
        o_df = pd.read_csv(a.observed_o, index_col=0)
        if list(o_df.columns) != hallmark_names:
            raise ValueError("Observed O Hallmarks do not match disease-eligible Hallmarks")
        O = o_df.loc[SYND].to_numpy(float)
        if a.p2_bootstrap_replay:
            br = pd.read_csv(a.p2_bootstrap_replay)
            if not {"syndrome", "delta_boot"}.issubset(br.columns):
                raise ValueError("P2 replay requires syndrome,delta_boot")
            boot = np.column_stack([
                br.loc[br["syndrome"].astype(str) == s, "delta_boot"].to_numpy(float)
                for s in SYND
            ])
            if boot.shape != (N_BOOT, 2):
                raise ValueError("P2 replay must contain exactly 2,000 values per syndrome")
    else:
        raise ValueError("Provide either patient expression or frozen observed O")

    p1obs = np.array([spearmanr(Mtilde[i], O[i]).statistic for i in range(2)])
    p1null = np.empty((N_NULL, 2), float)
    for b in range(N_NULL):
        Ab = np.divide(qnull[b] - mu, sd, out=np.zeros((2, K)), where=sd > 0)
        Ap = np.maximum(Ab, 0.0)
        Mb = D[None, :] * Ap
        for i in range(2):
            if Mb[i].sum() <= 0:
                p1null[b, i] = -1.0
            else:
                mt = Mb[i] / Mb[i].sum()
                p1null[b, i] = spearmanr(mt, O[i]).statistic
    p1p = (1 + (p1null >= p1obs).sum(axis=0)) / (N_NULL + 1)

    syndonly = np.array([spearmanr(Anorm[i], O[i]).statistic for i in range(2)])
    delta = p1obs - syndonly
    ci = np.full((2, 2), np.nan)
    if boot is not None:
        ci = np.quantile(boot, [0.025, 0.975], axis=0)

    rows = []
    for i, s in enumerate(SYND):
        rows.append({
            "syndrome": s,
            "n_hallmarks": K,
            "Aplus_nonzero": int(nonzero[i]),
            "compatibility_R": float(np.dot(D, Aplus[i]) / (np.linalg.norm(D) * np.linalg.norm(Aplus[i]))),
            "P1_spearman": float(p1obs[i]),
            "P1_empirical_p": float(p1p[i]),
            "syndrome_only_spearman": float(syndonly[i]),
            "P2_delta": float(delta[i]),
            "P2_boot_CI_low": float(ci[0, i]),
            "P2_boot_CI_high": float(ci[1, i]),
            "P2_support": bool(ci[0, i] > 0) if np.isfinite(ci[0, i]) else False,
            "P3_status": "NOT_IDENTIFIABLE_NOT_TESTED",
            "P3_selectivity": np.nan,
            "P3_empirical_p": np.nan,
        })

    pd.DataFrame(rows).to_csv(out / "DCSMI_URTI_primary_endpoints_RECONSTRUCTED.csv", index=False)
    pd.DataFrame(Mtilde, index=SYND, columns=hallmark_names).to_csv(out / "DCSMI_URTI_Mtilde_RECONSTRUCTED.csv")
    pd.DataFrame(O, index=SYND, columns=hallmark_names).to_csv(out / "DCSMI_URTI_O_RECONSTRUCTED.csv")
    pd.DataFrame(p1null, columns=SYND).to_csv(out / "DCSMI_URTI_P1_null_distribution_RECONSTRUCTED.csv", index=False)
    if boot is not None:
        pd.DataFrame({"syndrome": np.repeat(SYND, N_BOOT), "delta_boot": boot.T.reshape(-1)}).to_csv(
            out / "DCSMI_URTI_P2_bootstrap_distribution_RECONSTRUCTED.csv", index=False
        )

    inputs = {"prior": a.prior, "strata": a.strata, "hallmark": a.hallmark, "disease": a.disease}
    for key, val in {
        "disease_gene_map": a.disease_gene_map,
        "observed_o": a.observed_o,
        "patient_expression": a.patient_expression,
        "patient_gene_map": a.patient_gene_map,
        "sample_manifest": a.sample_manifest,
        "p2_bootstrap_replay": a.p2_bootstrap_replay,
    }.items():
        if val:
            inputs[key] = val
    manifest = {
        "status": "FUNCTIONAL_RECONSTRUCTION_NOT_HISTORICAL_EXECUTION_SCRIPT",
        "basis": "DCSMI_v1.0 prespecification + supplied psoriasis v1.0.1 endpoint runner + frozen URTI P3 amendment",
        "seed": SEED,
        "N_null": N_NULL,
        "N_bootstrap": N_BOOT,
        "n_targets": 374,
        "n_strata": 19,
        "n_hallmarks": K,
        "P3": "NOT_IDENTIFIABLE_NOT_TESTED",
        "P2_mode": "recomputed_patient_bootstrap" if a.patient_expression else ("frozen_bootstrap_replay" if a.p2_bootstrap_replay else "not_available"),
        "inputs_sha256": {k: sha256(v) for k, v in inputs.items()},
        "historical_equivalence_claim": False,
    }
    (out / "RECONSTRUCTED_ENDPOINT_RUN_MANIFEST.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
