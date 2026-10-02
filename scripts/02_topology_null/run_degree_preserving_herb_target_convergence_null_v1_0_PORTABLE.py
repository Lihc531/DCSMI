# -*- coding: utf-8 -*-
"""
Phase 12A — Degree-preserving herb–target convergence null v1.0

Purpose
-------
Test whether the three frozen LOPO syndrome target profiles are MORE CONVERGENT
than expected solely from the exact bipartite degree structure of the frozen
TCM-ID herb-target graph.

IMPORTANT
---------
This analysis is CAP-INDEPENDENT and was prespecified before the full null run.
Do not change the primary statistic, null model, number of nulls, burn-in,
thinning, seed, target mapping, or projection formula after inspecting results.
Any methodological change requires a new version.

Primary test
------------
Primary representation: frozen LOPO ubiquity-corrected target profile.
Primary statistic: generalized Jensen-Shannon divergence (GJSD, log base 2)
among PDOL, PHOL, and WHIL. Smaller GJSD = stronger convergence.
Primary empirical P: LOWER tail under an exact degree-preserving bipartite null.

Null randomization
------------------
Curveball trades on herb rows preserve exactly:
  * every herb degree,
  * every target degree,
  * total edge count,
  * simple bipartite graph status.

Four independent chains are used for mixing diagnostics.
"""

from __future__ import annotations

from pathlib import Path
import hashlib
import json
import math
import os
from datetime import datetime, timezone
from itertools import combinations

import numpy as np
import pandas as pd


# =============================================================================
# 1. USER SETTINGS — EDIT PATHS HERE
# =============================================================================

BASE_DIR = Path(os.environ.get("DCSMI_PROJECT_ROOT", Path(__file__).resolve().parents[2]))

# The formal frozen v1.1 target-projection directory. If your files are inside
# a subfolder named "frozen_final", point directly to that folder instead.
TARGET_PROJECTION_DIR = BASE_DIR / "05_weighted_herb_target_projection_v1_1_FROZEN"

OUTPUT_DIR = BASE_DIR / "12A_degree_preserving_herb_target_convergence_null_v1_0"

# Optional environment overrides used only for reproducible testing / another machine.
# Normal PyCharm use does not require them.
if os.environ.get("PHASE12A_INPUT_DIR"):
    TARGET_PROJECTION_DIR = Path(os.environ["PHASE12A_INPUT_DIR"])
if os.environ.get("PHASE12A_OUTPUT_DIR"):
    OUTPUT_DIR = Path(os.environ["PHASE12A_OUTPUT_DIR"])

CANONICAL_PAIRS_FILE = TARGET_PROJECTION_DIR / "canonical_herb_target_pairs_v1_1.csv"
MAPPING_FILE = TARGET_PROJECTION_DIR / "FROZEN_herb_to_target_mapping_dictionary_v1_1.csv"
FROZEN_PRIMARY_FILE = TARGET_PROJECTION_DIR / "LOPO_FROZEN_primary_target_weights_for_RWR_v1_1.csv"
TARGET_UBIQUITY_FILE = TARGET_PROJECTION_DIR / "target_database_target_ubiquity_v1_1.csv"

# Optional: place the prespecification file beside this script or set a full path.
PRESPECIFICATION_FILE = Path(__file__).with_name("PHASE12A_PRESPECIFICATION_FREEZE_v1_0.txt")

# Strictly verify the four frozen upstream files.
STRICT_HASH_CHECK = True
EXPECTED_SHA256 = {
    "canonical_pairs": "a6b283fe945885046b78141a94b6474d1acdef41a4ace8e941096fee95dff983",
    "mapping": "709fac26761667242c5754e4bc7c67be047f9398aea010542e25d195eea95e66",
    "frozen_primary": "9256cd3fe7ae142a91c7e081d97ba8a9225aa5aea3434c5f5673d53e5a10c9ef",
    "target_ubiquity": "7a714cbfadd44033aabcc7178d4ff045bbc826f6cc5606705acc76fcad002bba",
}

# =============================================================================
# 2. FROZEN ANALYSIS SETTINGS — DO NOT TUNE AFTER RESULTS
# =============================================================================

SYNDROME_CODES = ["PDOL", "PHOL", "WHIL"]
REPRESENTATION = "LOPO"

N_NULL = 10_000
N_CHAINS = 4
NULLS_PER_CHAIN = N_NULL // N_CHAINS
MASTER_SEED = 20260922

# Defined as multiples of the number of target-database herbs.
BURN_IN_TRADES_PER_HERB = 100
THIN_TRADES_PER_HERB = 5

# QA thresholds from the frozen prespecification.
PROFILE_RECONSTRUCTION_TOL = 1e-12
RHAT_PASS_MAX = 1.01
ABS_LAG1_ACF_PASS_MAX = 0.10

# Save degree QA every this many retained null graphs per chain.
DEGREE_QA_EVERY = 250

# If True, also save the randomized graph edge list for the FIRST retained
# null graph of each chain. This is audit-only and does not affect inference.
SAVE_FIRST_NULL_GRAPH_PER_CHAIN = True

# =============================================================================
# 3. HELPERS
# =============================================================================


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk_size)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def require_file(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Required input not found: {path}")


def require_columns(df: pd.DataFrame, cols: list[str], label: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"{label} missing columns {missing}; available={list(df.columns)}")


def verify_hash(path: Path, expected: str, label: str) -> str:
    actual = sha256_file(path)
    if STRICT_HASH_CHECK and actual != expected:
        raise RuntimeError(
            f"Frozen input hash mismatch for {label}:\n"
            f"  expected={expected}\n  actual  ={actual}\n  file={path}"
        )
    return actual


def safe_sum1(x: np.ndarray, axis: int = 1) -> np.ndarray:
    denom = x.sum(axis=axis, keepdims=True)
    if np.any(denom <= 0):
        raise RuntimeError("Cannot normalize a target profile with non-positive total mass.")
    return x / denom


def generalized_jsd_base2(P: np.ndarray) -> float:
    """Generalized JSD = mean_s KL_2(P_s || mean_s P_s)."""
    M = P.mean(axis=0)
    vals = []
    for p in P:
        mask = p > 0
        vals.append(float(np.sum(p[mask] * np.log2(p[mask] / M[mask]))))
    return float(np.mean(vals))


def pairwise_jsd_base2(p: np.ndarray, q: np.ndarray) -> float:
    m = 0.5 * (p + q)
    out = 0.0
    mp = p > 0
    mq = q > 0
    out += 0.5 * float(np.sum(p[mp] * np.log2(p[mp] / m[mp])))
    out += 0.5 * float(np.sum(q[mq] * np.log2(q[mq] / m[mq])))
    return out


def cosine_similarity(p: np.ndarray, q: np.ndarray) -> float:
    denom = float(np.linalg.norm(p) * np.linalg.norm(q))
    return float(np.dot(p, q) / denom) if denom > 0 else np.nan


def convergence_metrics(P: np.ndarray) -> dict[str, float]:
    """P shape = (3 syndromes, n_targets), each row sum=1."""
    out: dict[str, float] = {}
    out["gjsd"] = generalized_jsd_base2(P)
    out["common_mass"] = float(np.min(P, axis=0).sum())

    pair_jsd = []
    pair_cos = []
    for i, j in combinations(range(len(SYNDROME_CODES)), 2):
        a, b = SYNDROME_CODES[i], SYNDROME_CODES[j]
        js = pairwise_jsd_base2(P[i], P[j])
        cs = cosine_similarity(P[i], P[j])
        out[f"jsd_{a}_{b}"] = js
        out[f"cosine_{a}_{b}"] = cs
        pair_jsd.append(js)
        pair_cos.append(cs)

    out["mean_pairwise_jsd"] = float(np.mean(pair_jsd))
    out["mean_pairwise_cosine"] = float(np.mean(pair_cos))
    return out


def curveball_trade(neighbors: list[set[int]], rng: np.random.Generator) -> bool:
    """
    One Curveball trade between two herb rows.
    Exact herb degrees and target degrees are preserved.
    Returns True if a non-trivial trade occurred.
    """
    n = len(neighbors)
    a = int(rng.integers(n))
    b = int(rng.integers(n - 1))
    if b >= a:
        b += 1

    A = neighbors[a]
    B = neighbors[b]
    only_a = A - B
    only_b = B - A
    if not only_a or not only_b:
        return False

    pool = np.fromiter(only_a | only_b, dtype=np.int32)
    rng.shuffle(pool)
    n_a = len(only_a)
    common = A & B

    neighbors[a] = common | set(pool[:n_a].tolist())
    neighbors[b] = common | set(pool[n_a:].tolist())
    return True


def run_trades(neighbors: list[set[int]], rng: np.random.Generator, n_attempts: int) -> int:
    effective = 0
    for _ in range(n_attempts):
        effective += int(curveball_trade(neighbors, rng))
    return effective


def degree_vectors(neighbors: list[set[int]], n_targets: int) -> tuple[np.ndarray, np.ndarray]:
    row = np.fromiter((len(s) for s in neighbors), dtype=np.int32, count=len(neighbors))
    col = np.zeros(n_targets, dtype=np.int32)
    for s in neighbors:
        if s:
            idx = np.fromiter(s, dtype=np.int32)
            np.add.at(col, idx, 1)
    return row, col


def profiles_from_graph(
    neighbors: list[set[int]],
    syndrome_herb_coeff: np.ndarray,
    target_idf: np.ndarray,
    n_targets: int,
) -> np.ndarray:
    """
    Frozen ubiquity-corrected projection.
    syndrome_herb_coeff[s,h] already equals aggregated W_hs / degree_h.
    """
    scores = np.zeros((len(SYNDROME_CODES), n_targets), dtype=np.float64)
    active = np.any(syndrome_herb_coeff != 0, axis=0)

    for h_idx, target_set in enumerate(neighbors):
        if not active[h_idx] or not target_set:
            continue
        tidx = np.fromiter(target_set, dtype=np.int32)
        scores[:, tidx] += syndrome_herb_coeff[:, h_idx, None]

    scores *= target_idf[None, :]
    return safe_sum1(scores, axis=1)


def empirical_p_lower(obs: float, null: np.ndarray) -> float:
    return float((1 + np.sum(null <= obs)) / (len(null) + 1))


def empirical_p_upper(obs: float, null: np.ndarray) -> float:
    return float((1 + np.sum(null >= obs)) / (len(null) + 1))


def empirical_p_two_sided(obs: float, null: np.ndarray) -> float:
    return float(min(1.0, 2.0 * min(empirical_p_lower(obs, null), empirical_p_upper(obs, null))))


def autocorr(x: np.ndarray, lag: int) -> float:
    x = np.asarray(x, dtype=float)
    if len(x) <= lag or np.std(x[:-lag]) == 0 or np.std(x[lag:]) == 0:
        return np.nan
    return float(np.corrcoef(x[:-lag], x[lag:])[0, 1])


def gelman_rhat(chains: list[np.ndarray]) -> float:
    """Classic scalar Gelman-Rubin R-hat for equal-length chains."""
    if len(chains) < 2:
        return np.nan
    n = min(len(c) for c in chains)
    if n < 2:
        return np.nan
    arr = np.vstack([np.asarray(c[:n], dtype=float) for c in chains])
    m = arr.shape[0]
    means = arr.mean(axis=1)
    vars_ = arr.var(axis=1, ddof=1)
    W = float(vars_.mean())
    B = float(n * means.var(ddof=1))
    if W <= 0:
        return np.nan
    var_hat = ((n - 1) / n) * W + B / n
    return float(math.sqrt(var_hat / W))


def save_graph_edges(
    neighbors: list[set[int]], herbs: list[str], targets: list[str], path: Path
) -> None:
    rows = []
    for h_idx, s in enumerate(neighbors):
        for t_idx in sorted(s):
            rows.append((herbs[h_idx], targets[t_idx]))
    pd.DataFrame(rows, columns=["canonical_name", "target_id"]).to_csv(
        path, index=False, encoding="utf-8-sig"
    )


def make_manifest(output_dir: Path) -> pd.DataFrame:
    rows = []
    for p in sorted(output_dir.iterdir()):
        if p.is_file() and p.name != "SHA256_MANIFEST_v1_0.csv":
            rows.append({"file": p.name, "sha256": sha256_file(p), "bytes": p.stat().st_size})
    return pd.DataFrame(rows)


# =============================================================================
# 4. MAIN
# =============================================================================


def main() -> None:
    if N_NULL % N_CHAINS != 0:
        raise ValueError("N_NULL must be divisible by N_CHAINS.")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    for p in [CANONICAL_PAIRS_FILE, MAPPING_FILE, FROZEN_PRIMARY_FILE, TARGET_UBIQUITY_FILE]:
        require_file(p)

    print("[1/9] Verifying frozen upstream hashes...")
    input_hashes = {
        "canonical_pairs": verify_hash(CANONICAL_PAIRS_FILE, EXPECTED_SHA256["canonical_pairs"], "canonical_pairs"),
        "mapping": verify_hash(MAPPING_FILE, EXPECTED_SHA256["mapping"], "mapping"),
        "frozen_primary": verify_hash(FROZEN_PRIMARY_FILE, EXPECTED_SHA256["frozen_primary"], "frozen_primary"),
        "target_ubiquity": verify_hash(TARGET_UBIQUITY_FILE, EXPECTED_SHA256["target_ubiquity"], "target_ubiquity"),
    }

    pairs = pd.read_csv(CANONICAL_PAIRS_FILE)
    mapping = pd.read_csv(MAPPING_FILE)
    frozen_primary = pd.read_csv(FROZEN_PRIMARY_FILE)
    ubiq = pd.read_csv(TARGET_UBIQUITY_FILE)

    require_columns(pairs, ["canonical_name", "target_id"], "canonical pairs")
    require_columns(
        mapping,
        ["representation", "syndrome_code", "canonical_name", "target_mapping_name",
         "primary_weight_W", "mapped_to_target_db"],
        "mapping dictionary",
    )
    require_columns(
        frozen_primary,
        ["syndrome_code", "target_id", "ubiquity_corrected_sum1"],
        "frozen primary",
    )
    require_columns(ubiq, ["target_id", "target_herb_degree", "target_idf"], "target ubiquity")

    # Enforce simple unique graph.
    pairs = pairs[["canonical_name", "target_id"]].dropna().drop_duplicates().copy()
    herbs = sorted(pairs["canonical_name"].astype(str).unique().tolist())
    targets = sorted(pairs["target_id"].astype(str).unique().tolist())
    herb_to_i = {h: i for i, h in enumerate(herbs)}
    target_to_i = {t: i for i, t in enumerate(targets)}

    if len(herbs) != 172 or len(targets) != 374 or len(pairs) != 13285:
        raise RuntimeError(
            f"Unexpected frozen graph size: herbs={len(herbs)}, targets={len(targets)}, edges={len(pairs)}"
        )

    print("[2/9] Building observed frozen bipartite graph...")
    observed_neighbors: list[set[int]] = [set() for _ in herbs]
    for h, t in pairs.itertuples(index=False):
        observed_neighbors[herb_to_i[str(h)]].add(target_to_i[str(t)])

    observed_row_degree, observed_col_degree = degree_vectors(observed_neighbors, len(targets))
    if int(observed_row_degree.sum()) != len(pairs) or int(observed_col_degree.sum()) != len(pairs):
        raise RuntimeError("Observed degree totals do not equal frozen edge count.")

    # Cross-check target degrees and IDF against frozen target-ubiquity table.
    ubiq_idx = ubiq.set_index("target_id").reindex(targets)
    if ubiq_idx["target_herb_degree"].isna().any() or ubiq_idx["target_idf"].isna().any():
        raise RuntimeError("Frozen target ubiquity table does not cover the 374-target universe.")
    max_target_degree_error = int(np.max(np.abs(
        observed_col_degree - ubiq_idx["target_herb_degree"].to_numpy(dtype=int)
    )))
    if max_target_degree_error != 0:
        raise RuntimeError(f"Target-degree crosscheck failed: max error={max_target_degree_error}")
    target_idf = ubiq_idx["target_idf"].to_numpy(dtype=float)

    print("[3/9] Building frozen LOPO syndrome-herb coefficients...")
    lopo = mapping[
        mapping["representation"].astype(str).eq(REPRESENTATION)
        & mapping["syndrome_code"].astype(str).isin(SYNDROME_CODES)
        & pd.to_numeric(mapping["mapped_to_target_db"], errors="coerce").fillna(0).astype(int).eq(1)
    ].copy()
    lopo["primary_weight_W"] = pd.to_numeric(lopo["primary_weight_W"], errors="raise")
    lopo["target_mapping_name"] = lopo["target_mapping_name"].astype(str)

    syndrome_herb_coeff = np.zeros((len(SYNDROME_CODES), len(herbs)), dtype=np.float64)
    mapping_rows = []
    for row in lopo.itertuples(index=False):
        s = str(row.syndrome_code)
        h = str(row.target_mapping_name)
        if h not in herb_to_i:
            raise RuntimeError(f"Mapped herb '{h}' is absent from frozen target graph.")
        s_idx = SYNDROME_CODES.index(s)
        h_idx = herb_to_i[h]
        coeff = float(row.primary_weight_W) / float(observed_row_degree[h_idx])
        syndrome_herb_coeff[s_idx, h_idx] += coeff
        mapping_rows.append({
            "syndrome_code": s,
            "canonical_name": str(row.canonical_name),
            "target_mapping_name": h,
            "primary_weight_W": float(row.primary_weight_W),
            "frozen_herb_degree": int(observed_row_degree[h_idx]),
            "degree_normalized_coefficient": coeff,
        })
    pd.DataFrame(mapping_rows).to_csv(
        OUTPUT_DIR / "phase12A_frozen_LOPO_herb_coefficients_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )

    print("[4/9] Reconstructing observed frozen target profiles...")
    observed_profiles = profiles_from_graph(
        observed_neighbors, syndrome_herb_coeff, target_idf, len(targets)
    )

    frozen_pivot = (
        frozen_primary.pivot(index="target_id", columns="syndrome_code", values="ubiquity_corrected_sum1")
        .reindex(targets)
        .fillna(0.0)
    )
    frozen_matrix = np.vstack([frozen_pivot[s].to_numpy(dtype=float) for s in SYNDROME_CODES])
    max_profile_error = float(np.max(np.abs(observed_profiles - frozen_matrix)))
    reconstruction_pass = bool(max_profile_error <= PROFILE_RECONSTRUCTION_TOL)
    if not reconstruction_pass:
        raise RuntimeError(
            f"Frozen primary reconstruction failed: max_abs_error={max_profile_error:.3e}"
        )

    recon_rows = []
    for s_idx, s in enumerate(SYNDROME_CODES):
        recon_rows.append({
            "syndrome_code": s,
            "reconstructed_sum": float(observed_profiles[s_idx].sum()),
            "frozen_sum": float(frozen_matrix[s_idx].sum()),
            "max_abs_target_error": float(np.max(np.abs(observed_profiles[s_idx] - frozen_matrix[s_idx]))),
            "pass": bool(np.max(np.abs(observed_profiles[s_idx] - frozen_matrix[s_idx])) <= PROFILE_RECONSTRUCTION_TOL),
        })
    pd.DataFrame(recon_rows).to_csv(
        OUTPUT_DIR / "phase12A_frozen_primary_reconstruction_QA_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )

    observed_metrics = convergence_metrics(observed_profiles)
    observed_df = pd.DataFrame([{"representation": "LOPO_ubiquity_corrected_sum1", **observed_metrics}])
    observed_df.to_csv(
        OUTPUT_DIR / "phase12A_observed_convergence_metrics_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )

    # Save observed profiles for audit.
    profile_rows = []
    for s_idx, s in enumerate(SYNDROME_CODES):
        for t_idx, t in enumerate(targets):
            profile_rows.append({
                "syndrome_code": s,
                "target_id": t,
                "reconstructed_ubiquity_corrected_sum1": observed_profiles[s_idx, t_idx],
                "frozen_ubiquity_corrected_sum1": frozen_matrix[s_idx, t_idx],
                "absolute_error": abs(observed_profiles[s_idx, t_idx] - frozen_matrix[s_idx, t_idx]),
            })
    pd.DataFrame(profile_rows).to_csv(
        OUTPUT_DIR / "phase12A_observed_target_profiles_crosscheck_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )

    print("[5/9] Running four-chain degree-preserving Curveball null...")
    burn_in = BURN_IN_TRADES_PER_HERB * len(herbs)
    thinning = THIN_TRADES_PER_HERB * len(herbs)
    seed_sequence = np.random.SeedSequence(MASTER_SEED)
    child_seeds = seed_sequence.spawn(N_CHAINS)

    null_rows: list[dict] = []
    degree_qa_rows: list[dict] = []
    chain_gjsd: list[np.ndarray] = []
    chain_trade_rows: list[dict] = []

    for chain_id, child_seed in enumerate(child_seeds, start=1):
        rng = np.random.default_rng(child_seed)
        neighbors = [set(s) for s in observed_neighbors]

        effective_burn = run_trades(neighbors, rng, burn_in)
        row_deg, col_deg = degree_vectors(neighbors, len(targets))
        burn_row_err = int(np.max(np.abs(row_deg - observed_row_degree)))
        burn_col_err = int(np.max(np.abs(col_deg - observed_col_degree)))
        if burn_row_err != 0 or burn_col_err != 0:
            raise RuntimeError(f"Degree preservation failed after burn-in in chain {chain_id}.")

        chain_values = []
        effective_between_total = 0
        print(f"    Chain {chain_id}/{N_CHAINS}: burn-in effective trades={effective_burn}/{burn_in}")

        for iteration in range(1, NULLS_PER_CHAIN + 1):
            effective_between = run_trades(neighbors, rng, thinning)
            effective_between_total += effective_between

            P = profiles_from_graph(neighbors, syndrome_herb_coeff, target_idf, len(targets))
            metrics = convergence_metrics(P)
            chain_values.append(metrics["gjsd"])
            null_rows.append({
                "chain_id": chain_id,
                "iteration_within_chain": iteration,
                "global_null_index": (chain_id - 1) * NULLS_PER_CHAIN + iteration,
                **metrics,
            })

            if iteration == 1 and SAVE_FIRST_NULL_GRAPH_PER_CHAIN:
                save_graph_edges(
                    neighbors, herbs, targets,
                    OUTPUT_DIR / f"phase12A_chain{chain_id}_first_retained_null_graph_v1_0.csv"
                )

            if iteration == 1 or iteration % DEGREE_QA_EVERY == 0 or iteration == NULLS_PER_CHAIN:
                row_deg, col_deg = degree_vectors(neighbors, len(targets))
                row_err = int(np.max(np.abs(row_deg - observed_row_degree)))
                col_err = int(np.max(np.abs(col_deg - observed_col_degree)))
                degree_qa_rows.append({
                    "chain_id": chain_id,
                    "iteration_within_chain": iteration,
                    "max_abs_herb_degree_error": row_err,
                    "max_abs_target_degree_error": col_err,
                    "edge_count": int(row_deg.sum()),
                    "edge_count_error": int(row_deg.sum() - len(pairs)),
                    "pass": bool(row_err == 0 and col_err == 0 and int(row_deg.sum()) == len(pairs)),
                })
                if row_err != 0 or col_err != 0 or int(row_deg.sum()) != len(pairs):
                    raise RuntimeError(
                        f"Degree preservation failed in chain {chain_id}, iteration {iteration}."
                    )

            if iteration % 500 == 0:
                print(f"      retained {iteration}/{NULLS_PER_CHAIN}")

        chain_arr = np.asarray(chain_values, dtype=float)
        chain_gjsd.append(chain_arr)
        chain_trade_rows.append({
            "chain_id": chain_id,
            "burn_in_trade_attempts": burn_in,
            "burn_in_effective_trades": effective_burn,
            "burn_in_effective_fraction": effective_burn / burn_in,
            "between_sample_trade_attempts_total": thinning * NULLS_PER_CHAIN,
            "between_sample_effective_trades_total": effective_between_total,
            "between_sample_effective_fraction": effective_between_total / (thinning * NULLS_PER_CHAIN),
        })

    null_df = pd.DataFrame(null_rows)
    if len(null_df) != N_NULL:
        raise RuntimeError(f"Expected {N_NULL} null rows; got {len(null_df)}")
    null_df.to_csv(
        OUTPUT_DIR / "phase12A_degree_preserving_null_distribution_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )
    pd.DataFrame(degree_qa_rows).to_csv(
        OUTPUT_DIR / "phase12A_degree_preservation_QA_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )
    pd.DataFrame(chain_trade_rows).to_csv(
        OUTPUT_DIR / "phase12A_curveball_trade_QA_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )

    print("[6/9] Computing frozen primary inference and secondary summaries...")
    summary_rows = []
    metric_directions = {
        "gjsd": "lower_is_more_convergent_PRIMARY",
        "common_mass": "higher_is_more_convergent",
        "mean_pairwise_jsd": "lower_is_more_convergent",
        "mean_pairwise_cosine": "higher_is_more_convergent",
        "jsd_PDOL_PHOL": "lower_is_more_convergent",
        "jsd_PDOL_WHIL": "lower_is_more_convergent",
        "jsd_PHOL_WHIL": "lower_is_more_convergent",
        "cosine_PDOL_PHOL": "higher_is_more_convergent",
        "cosine_PDOL_WHIL": "higher_is_more_convergent",
        "cosine_PHOL_WHIL": "higher_is_more_convergent",
    }

    for metric, direction in metric_directions.items():
        obs = float(observed_metrics[metric])
        null = null_df[metric].to_numpy(dtype=float)
        null_mean = float(np.mean(null))
        null_sd = float(np.std(null, ddof=1))
        z = float((obs - null_mean) / null_sd) if null_sd > 0 else np.nan
        summary_rows.append({
            "metric": metric,
            "role": "PRIMARY" if metric == "gjsd" else "SECONDARY_DIAGNOSTIC",
            "convergence_direction": direction,
            "observed": obs,
            "null_mean": null_mean,
            "null_sd": null_sd,
            "z_vs_null": z,
            "null_q025": float(np.quantile(null, 0.025)),
            "null_q050": float(np.quantile(null, 0.50)),
            "null_q975": float(np.quantile(null, 0.975)),
            "empirical_p_lower": empirical_p_lower(obs, null),
            "empirical_p_upper": empirical_p_upper(obs, null),
            "empirical_p_two_sided_2min": empirical_p_two_sided(obs, null),
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(
        OUTPUT_DIR / "phase12A_convergence_null_summary_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )

    print("[7/9] Computing chain mixing QA...")
    rhat = gelman_rhat(chain_gjsd)
    mixing_rows = []
    for chain_id, arr in enumerate(chain_gjsd, start=1):
        row = {
            "chain_id": chain_id,
            "n_retained": len(arr),
            "gjsd_mean": float(np.mean(arr)),
            "gjsd_sd": float(np.std(arr, ddof=1)),
            "gjsd_acf_lag1": autocorr(arr, 1),
            "gjsd_acf_lag5": autocorr(arr, 5),
            "gjsd_acf_lag10": autocorr(arr, 10),
            "pooled_gjsd_rhat": rhat,
        }
        row["lag1_pass"] = bool(abs(row["gjsd_acf_lag1"]) <= ABS_LAG1_ACF_PASS_MAX)
        row["rhat_pass"] = bool(np.isfinite(rhat) and rhat <= RHAT_PASS_MAX)
        row["mixing_pass"] = bool(row["lag1_pass"] and row["rhat_pass"])
        mixing_rows.append(row)

    mixing_df = pd.DataFrame(mixing_rows)
    mixing_df.to_csv(
        OUTPUT_DIR / "phase12A_chain_mixing_QA_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )

    all_degree_pass = bool(pd.DataFrame(degree_qa_rows)["pass"].all())
    all_mixing_pass = bool(mixing_df["mixing_pass"].all())
    technical_pass = bool(reconstruction_pass and all_degree_pass and all_mixing_pass)

    primary_row = summary_df.loc[summary_df["metric"].eq("gjsd")].iloc[0]
    primary_p = float(primary_row["empirical_p_lower"])
    opposite_p = float(primary_row["empirical_p_upper"])

    if not technical_pass:
        interpretation = (
            "TECHNICAL_QA_FAILED: Phase 12A v1.0 must not be biologically interpreted. "
            "A new preregistered version is required after diagnosing mixing/QA."
        )
    elif primary_p < 0.05:
        interpretation = (
            "PRIMARY_SUPPORT: observed LOPO target profiles are more convergent than expected "
            "under the exact herb- and target-degree-preserving null. This supports structured "
            "shared target architecture beyond degree structure alone, but does not establish CAP relevance."
        )
    elif opposite_p < 0.05:
        interpretation = (
            "PRIMARY_NOT_SUPPORTED_OPPOSITE_DEPARTURE: the prespecified convergence hypothesis is not supported; "
            "instead, observed profiles are more divergent than expected under the exact degree-preserving null. "
            "This must not be relabeled as support for convergence."
        )
    else:
        interpretation = (
            "PRIMARY_NOT_SUPPORTED: observed convergence does not differ in the prespecified lower-tail direction "
            "from the exact degree-preserving herb-target null."
        )

    print("[8/9] Writing metadata and freeze-candidate statement...")
    metadata = {
        "analysis_name": "Degree-preserving herb-target convergence null v1.0",
        "phase": "12A",
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "status": "FULL_NULL_COMPLETED_AWAITING_FINAL_AUDIT",
        "CAP_transcriptomics_used": False,
        "RWR_used": False,
        "pathway_enrichment_used": False,
        "frozen_representation": "LOPO ubiquity_corrected_sum1",
        "input_files": {
            "canonical_pairs": str(CANONICAL_PAIRS_FILE),
            "mapping": str(MAPPING_FILE),
            "frozen_primary": str(FROZEN_PRIMARY_FILE),
            "target_ubiquity": str(TARGET_UBIQUITY_FILE),
        },
        "input_sha256": input_hashes,
        "graph": {
            "n_herbs": len(herbs),
            "n_targets": len(targets),
            "n_edges": len(pairs),
            "max_target_degree_crosscheck_error": max_target_degree_error,
        },
        "null": {
            "algorithm": "Curveball row trades on simple bipartite herb-target graph",
            "n_null": N_NULL,
            "n_chains": N_CHAINS,
            "nulls_per_chain": NULLS_PER_CHAIN,
            "master_seed": MASTER_SEED,
            "burn_in_trade_attempts_per_chain": burn_in,
            "thinning_trade_attempts_between_retained_graphs": thinning,
            "preserves_exact_herb_degrees": True,
            "preserves_exact_target_degrees": True,
        },
        "primary": {
            "statistic": "generalized Jensen-Shannon divergence, log base 2",
            "direction": "lower tail = more convergence",
            "empirical_p_definition": "(1 + count(null <= observed)) / (N + 1)",
            "alpha": 0.05,
        },
        "QA": {
            "profile_reconstruction_max_abs_error": max_profile_error,
            "profile_reconstruction_tolerance": PROFILE_RECONSTRUCTION_TOL,
            "all_degree_checks_pass": all_degree_pass,
            "gjsd_rhat": rhat,
            "rhat_pass_threshold": RHAT_PASS_MAX,
            "abs_lag1_acf_pass_threshold": ABS_LAG1_ACF_PASS_MAX,
            "all_mixing_checks_pass": all_mixing_pass,
            "technical_pass": technical_pass,
        },
        "primary_result": {
            "observed_gjsd": float(primary_row["observed"]),
            "null_mean_gjsd": float(primary_row["null_mean"]),
            "z_vs_null": float(primary_row["z_vs_null"]),
            "empirical_p_lower_PRIMARY": primary_p,
            "empirical_p_upper_diagnostic": opposite_p,
            "empirical_p_two_sided_diagnostic": float(primary_row["empirical_p_two_sided_2min"]),
        },
        "frozen_interpretation_rule_applied": interpretation,
    }

    with open(OUTPUT_DIR / "phase12A_run_metadata_v1_0.json", "w", encoding="utf-8") as f:
        json.dump(metadata, f, ensure_ascii=False, indent=2)

    statement = f"""PHASE 12A — DEGREE-PRESERVING HERB–TARGET CONVERGENCE NULL v1.0\nRESULT FREEZE CANDIDATE — AWAITING INDEPENDENT AUDIT\n\nTechnical QA pass: {technical_pass}\nObserved GJSD: {float(primary_row['observed']):.12g}\nNull mean GJSD: {float(primary_row['null_mean']):.12g}\nNull SD: {float(primary_row['null_sd']):.12g}\nZ vs null: {float(primary_row['z_vs_null']):.6f}\nPrimary lower-tail empirical P: {primary_p:.12g}\nDiagnostic upper-tail empirical P: {opposite_p:.12g}\nDiagnostic two-sided empirical P: {float(primary_row['empirical_p_two_sided_2min']):.12g}\nPooled four-chain R-hat: {rhat:.8f}\nMax frozen-profile reconstruction error: {max_profile_error:.3e}\n\nFrozen interpretation:\n{interpretation}\n\nNo CAP transcriptomic data were used in Phase 12A.\nDo not alter the null model/statistic after this result. Any change requires a new version.\n"""
    with open(OUTPUT_DIR / "PHASE12A_RESULT_FREEZE_CANDIDATE_v1_0.txt", "w", encoding="utf-8") as f:
        f.write(statement)

    # Copy prespecification into output if available, preserving the pre-run document.
    if PRESPECIFICATION_FILE.exists():
        prespec_text = PRESPECIFICATION_FILE.read_text(encoding="utf-8")
        (OUTPUT_DIR / "PHASE12A_PRESPECIFICATION_FREEZE_v1_0.txt").write_text(
            prespec_text, encoding="utf-8"
        )

    print("[9/9] Writing SHA256 manifest...")
    manifest = make_manifest(OUTPUT_DIR)
    manifest.to_csv(
        OUTPUT_DIR / "SHA256_MANIFEST_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )

    print("\n=== PHASE 12A COMPLETE ===")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Technical QA pass: {technical_pass}")
    print(f"Observed GJSD: {float(primary_row['observed']):.8f}")
    print(f"Null mean GJSD: {float(primary_row['null_mean']):.8f}")
    print(f"Primary lower-tail P: {primary_p:.8g}")
    print(f"Diagnostic upper-tail P: {opposite_p:.8g}")
    print(f"R-hat: {rhat:.6f}")
    print(interpretation)


if __name__ == "__main__":
    main()
