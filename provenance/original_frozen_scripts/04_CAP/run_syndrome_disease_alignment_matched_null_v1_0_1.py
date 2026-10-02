# -*- coding: utf-8 -*-
"""
run_syndrome_disease_alignment_matched_null_v1_0_1.py

Primary purpose
---------------
Quantify disease-context alignment between the FROZEN LOPO syndrome target
representations and the FROZEN CAP continuous meta-Z program, with a matched
null that preserves the target-database/network structure.

This script is downstream-only. It does NOT alter:
- syndrome-herb weights;
- herb-target mapping;
- target ubiquity correction;
- STRING threshold/edge weights;
- RWR restart probability;
- CAP signed-Z/meta-Z construction.

Primary syndrome representation
-------------------------------
LOPO ubiquity-corrected target weights:
    ubiquity_corrected_sum1

Primary disease representation
------------------------------
All 12,083 shared genes:
    Z_meta_equal

Primary raw-target alignment
----------------------------
For syndrome S, with frozen target weights w_g and CAP meta-Z z_g:

    A_raw(S) =
        sum_g w_g * I_CAP(g) * z_g
        --------------------------------
        sum_g w_g * I_CAP(g)

Thus the statistic is the weighted mean CAP meta-Z among measured frozen
targets. The denominator is reported explicitly because only a subset of the
target-database universe is measured in both CAP cohorts.

Secondary RWR alignment
-----------------------
STRING RWR is already frozen as a SECONDARY network-localization layer.

For RWR stationary probabilities p_g:

    A_RWR(S) =
        sum_g p_g * I_CAP(g) * z_g
        --------------------------------
        sum_g p_g * I_CAP(g)

The script computes the observed value from the frozen RWR output and
independently reproduces it with an adjoint RWR calculation.

Matched null
------------
The null universe is the frozen TCM-ID target-database universe.

For each syndrome, the complete target-universe vector is formed by assigning:
- the frozen positive target weight to observed syndrome targets;
- zero to target-database genes not present in that syndrome profile.

Genes are stratified by:
1. STRING unweighted-degree quartile; and
2. TCM-ID target-herb-ubiquity quartile.

Genes absent from the frozen STRING graph use NO_STRING as their degree category,
but remain stratified by the same target-herb-ubiquity bins as the full target universe.

Within EACH stratum, the entire syndrome weight vector (including zeros) is
permuted. Therefore every null replicate preserves exactly:
- the global multiset of weights;
- the number of positive/zero weights;
- the positive/zero count within every matched stratum;
- the total weight within every matched stratum;
- the target-database universe.

The SAME random stratum permutation is applied to all three syndromes in each
replicate. This enables paired null comparisons of syndrome differences.

For RWR nulls, the permuted seed vector is propagated analytically via the
adjoint stationary RWR equation. This is mathematically equivalent to rerunning
RWR for every null replicate but is much faster.

No CAP result is used to choose the matching bins, STRING parameters, restart
probability, or syndrome representation.

Primary inferential quantity
----------------------------
For each method and syndrome:

    Z_align = (A_obs - mean(A_null)) / sd(A_null)

The script reports:
- empirical upper-tail P: positive alignment beyond matched background;
- empirical lower-tail P: negative departure;
- empirical two-sided P around the null mean.

All empirical P values use the +1 correction:
    (count + 1) / (N_null + 1)

Important interpretation
------------------------
Positive CAP meta-Z means higher expression in CAP than healthy control.
Negative CAP meta-Z means lower expression in CAP than healthy control.
This is disease-context localization, NOT therapeutic reversal.

Shared/specific syndrome decomposition is intentionally NOT reconstructed here.
It must be generated and frozen upstream from observed-only syndrome-herb
weights before it is added to a later alignment version.
"""

from __future__ import annotations

from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import math
import re
import warnings

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix, csr_matrix, diags


# =============================================================================
# 1. USER SETTINGS — EDIT PATHS HERE IF NEEDED
# =============================================================================

PROJECT_ROOT = Path(r"E:\00project\Tanre Yongfei\重构")

# Preferred locations. If a preferred path is missing, the script searches
# recursively under PROJECT_ROOT for the exact filename.
FROZEN_TARGET_SEED_FILE = (
    PROJECT_ROOT
    / "05_weighted_herb_target_projection_v1_1_FROZEN"
    / "LOPO_FROZEN_primary_target_weights_for_RWR_v1_1.csv"
)

TARGET_UBIQUITY_FILE = (
    PROJECT_ROOT
    / "05_weighted_herb_target_projection_v1_1_FROZEN"
    / "target_database_target_ubiquity_v1_1.csv"
)

CAP_META_Z_FILE = (
    PROJECT_ROOT
    / "07_CAP_transcriptomic_metaZ_v1_0_1"
    / "CAP_CONTINUOUS_PRIMARY_equal_weight_metaZ_v1_0_1.csv"
)

STRING_CLEAN_EDGE_FILE = (
    PROJECT_ROOT
    / "06_STRING_RWR_frozen_v1_0"
    / "STRING_clean_high_confidence_edges_v1_0.csv"
)

RWR_PRIMARY_SCORE_FILE = (
    PROJECT_ROOT
    / "06_STRING_RWR_frozen_v1_0"
    / "RWR_PRIMARY_r0p50_scores_v1_0.csv"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "08_syndrome_disease_alignment_matched_null_v1_0_1"
)

# Frozen input SHA256 values. These prevent accidental mixing of analysis
# branches. Set ENFORCE_EXPECTED_SHA256=False only for deliberate debugging,
# never to silently accept a different frozen input.
ENFORCE_EXPECTED_SHA256 = True

EXPECTED_SHA256 = {
    "frozen_target_seed": (
        "9256cd3fe7ae142a91c7e081d97ba8a9225aa5aea3434c5f5673d53e5a10c9ef"
    ),
    "target_ubiquity": (
        "7a714cbfadd44033aabcc7178d4ff045bbc826f6cc5606705acc76fcad002bba"
    ),
    "cap_meta_z": (
        "168d9a23c06df388b69f73d92b2e7811b1f5332dfb507893511367fd752be58e"
    ),
    "string_clean_edges": (
        "deb2d74bdf0bcb7f249883c39d18bfecab9c9e9f5567612e0fc1fdc6db5ac6ba"
    ),
    "rwr_primary_scores": (
        "b97e578a8dfcf95014e4d61c899434a02ab867edd7e08330f945da035d776978"
    ),
}

# -------------------------------------------------------------------------
# Frozen analysis settings
# -------------------------------------------------------------------------

SYNDROME_CODES = ["PDOL", "PHOL", "WHIL"]

TARGET_WEIGHT_COLUMN = "ubiquity_corrected_sum1"
CAP_PRIMARY_Z_COLUMN = "Z_meta_equal"
CAP_SENSITIVITY_Z_COLUMN = "Z_meta_neff_weighted"

PRIMARY_RWR_RESTART = 0.5

# Final null count. 10,000 gives a minimum empirical P of 1/10,001.
N_NULL = 10_000
RANDOM_SEED = 20260921

# Matched-null strata.
N_DEGREE_BINS = 4
N_UBIQUITY_BINS = 4

# Adjoint RWR convergence.
ADJOINT_TOL_MAXABS = 1e-12
ADJOINT_MAX_ITER = 10_000

# Numerical QA.
WEIGHT_SUM_TOL = 1e-8
RWR_CROSSCHECK_TOL = 1e-10
NULL_INVARIANT_TOL = 1e-12

# Minimum disease-measured weight/mass before warning.
MIN_RECOMMENDED_MEASURED_MASS = 0.50

# Save all 10,000 null values for reproducibility.
SAVE_FULL_NULL_DISTRIBUTION = True


# =============================================================================
# 2. BASIC HELPERS
# =============================================================================

def norm_gene(x) -> str:
    if pd.isna(x):
        return ""
    return re.sub(r"\s+", "", str(x).strip()).upper()


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk_size)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def require_columns(df: pd.DataFrame, cols: list[str], name: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"{name} missing required columns: {missing}\n"
            f"Available columns: {list(df.columns)}"
        )


def resolve_input(preferred: Path, label: str) -> Path:
    """
    Resolve a frozen input by preferred path or exact-filename recursive search.

    If multiple files with the same name exist:
    - if their SHA256 values are identical, use the first sorted path;
    - if hashes differ, stop rather than silently choosing an analysis branch.
    """
    if preferred.exists():
        return preferred

    if not PROJECT_ROOT.exists():
        raise FileNotFoundError(
            f"{label}: preferred file not found:\n{preferred}\n"
            f"PROJECT_ROOT also does not exist:\n{PROJECT_ROOT}"
        )

    matches = sorted(PROJECT_ROOT.rglob(preferred.name))
    if not matches:
        raise FileNotFoundError(
            f"{label}: could not find '{preferred.name}' under:\n"
            f"{PROJECT_ROOT}"
        )

    if len(matches) == 1:
        print(f"  Auto-resolved {label}:\n    {matches[0]}")
        return matches[0]

    hashes = {str(p): sha256_file(p) for p in matches}
    unique_hashes = set(hashes.values())
    if len(unique_hashes) == 1:
        chosen = matches[0]
        print(
            f"  Multiple identical copies found for {label}; "
            f"using:\n    {chosen}"
        )
        return chosen

    details = "\n".join(
        f"  {p}\n    SHA256={h}" for p, h in hashes.items()
    )
    raise RuntimeError(
        f"{label}: multiple non-identical files named '{preferred.name}' "
        f"were found. Resolve manually.\n{details}"
    )


def verify_hash(path: Path, key: str) -> str:
    observed = sha256_file(path)
    expected = EXPECTED_SHA256.get(key)

    if ENFORCE_EXPECTED_SHA256 and expected and observed != expected:
        raise RuntimeError(
            f"Frozen input SHA256 mismatch for {key}.\n"
            f"File: {path}\n"
            f"Expected: {expected}\n"
            f"Observed: {observed}\n"
            "Do not mix analysis branches. If this is an intentional new "
            "version, create a new downstream version rather than bypassing "
            "the frozen v1.0 inputs."
        )
    return observed


def empirical_p_values(
    null_values: np.ndarray,
    observed: float,
) -> dict:
    null_values = np.asarray(null_values, dtype=float)
    n = len(null_values)

    mu = float(np.mean(null_values))
    sd = float(np.std(null_values, ddof=1))

    z_align = (
        (float(observed) - mu) / sd
        if sd > 0 else np.nan
    )

    p_upper = (
        1 + int(np.sum(null_values >= observed))
    ) / (n + 1)

    p_lower = (
        1 + int(np.sum(null_values <= observed))
    ) / (n + 1)

    obs_dev = abs(float(observed) - mu)
    p_two = (
        1 + int(np.sum(np.abs(null_values - mu) >= obs_dev))
    ) / (n + 1)

    return {
        "null_mean": mu,
        "null_sd": sd,
        "Z_align": z_align,
        "empirical_p_upper": p_upper,
        "empirical_p_lower": p_lower,
        "empirical_p_two_sided": p_two,
        "null_q025": float(np.quantile(null_values, 0.025)),
        "null_q25": float(np.quantile(null_values, 0.25)),
        "null_median": float(np.quantile(null_values, 0.50)),
        "null_q75": float(np.quantile(null_values, 0.75)),
        "null_q975": float(np.quantile(null_values, 0.975)),
    }


# =============================================================================
# 3. LOAD FROZEN INPUTS
# =============================================================================

def load_frozen_inputs():
    seed_path = resolve_input(
        FROZEN_TARGET_SEED_FILE, "frozen target seed"
    )
    ubiq_path = resolve_input(
        TARGET_UBIQUITY_FILE, "target ubiquity"
    )
    cap_path = resolve_input(
        CAP_META_Z_FILE, "CAP meta-Z"
    )
    edge_path = resolve_input(
        STRING_CLEAN_EDGE_FILE, "STRING clean edges"
    )
    rwr_path = resolve_input(
        RWR_PRIMARY_SCORE_FILE, "RWR primary scores"
    )

    hashes = {
        "frozen_target_seed": verify_hash(
            seed_path, "frozen_target_seed"
        ),
        "target_ubiquity": verify_hash(
            ubiq_path, "target_ubiquity"
        ),
        "cap_meta_z": verify_hash(
            cap_path, "cap_meta_z"
        ),
        "string_clean_edges": verify_hash(
            edge_path, "string_clean_edges"
        ),
        "rwr_primary_scores": verify_hash(
            rwr_path, "rwr_primary_scores"
        ),
    }

    seeds = pd.read_csv(seed_path)
    require_columns(
        seeds,
        ["syndrome_code", "target_id", TARGET_WEIGHT_COLUMN],
        seed_path.name,
    )
    seeds = seeds.copy()
    seeds["syndrome_code"] = (
        seeds["syndrome_code"].astype(str).str.strip()
    )
    seeds["target_id"] = seeds["target_id"].map(norm_gene)
    seeds[TARGET_WEIGHT_COLUMN] = pd.to_numeric(
        seeds[TARGET_WEIGHT_COLUMN], errors="coerce"
    )
    seeds = seeds[
        seeds["syndrome_code"].isin(SYNDROME_CODES)
        & seeds["target_id"].ne("")
        & seeds[TARGET_WEIGHT_COLUMN].notna()
        & seeds[TARGET_WEIGHT_COLUMN].gt(0)
    ].copy()

    if seeds.duplicated(
        ["syndrome_code", "target_id"]
    ).any():
        raise RuntimeError(
            "Duplicate syndrome-target rows in frozen seed file."
        )

    sums = seeds.groupby("syndrome_code")[
        TARGET_WEIGHT_COLUMN
    ].sum()
    bad = sums[(sums - 1.0).abs() > WEIGHT_SUM_TOL]
    if len(bad):
        raise RuntimeError(
            "Frozen syndrome target weights do not sum to 1:\n"
            f"{bad.to_string()}"
        )

    ubiq = pd.read_csv(ubiq_path)
    require_columns(
        ubiq,
        [
            "target_id",
            "target_herb_degree",
            "N_database_herbs",
            "target_idf",
        ],
        ubiq_path.name,
    )
    ubiq = ubiq.copy()
    ubiq["target_id"] = ubiq["target_id"].map(norm_gene)
    ubiq = ubiq[ubiq["target_id"].ne("")].copy()

    if ubiq["target_id"].duplicated().any():
        raise RuntimeError(
            "Duplicate target IDs in target ubiquity file."
        )

    cap = pd.read_csv(cap_path)
    require_columns(
        cap,
        ["gene", CAP_PRIMARY_Z_COLUMN],
        cap_path.name,
    )
    cap = cap.copy()
    cap["gene"] = cap["gene"].map(norm_gene)
    cap[CAP_PRIMARY_Z_COLUMN] = pd.to_numeric(
        cap[CAP_PRIMARY_Z_COLUMN], errors="coerce"
    )

    if CAP_SENSITIVITY_Z_COLUMN in cap.columns:
        cap[CAP_SENSITIVITY_Z_COLUMN] = pd.to_numeric(
            cap[CAP_SENSITIVITY_Z_COLUMN], errors="coerce"
        )

    cap = cap[
        cap["gene"].ne("")
        & cap[CAP_PRIMARY_Z_COLUMN].notna()
    ].copy()

    if cap["gene"].duplicated().any():
        raise RuntimeError(
            "Duplicate genes in frozen CAP meta-Z file."
        )

    edges = pd.read_csv(edge_path)
    require_columns(
        edges, ["u", "v", "edge_weight"], edge_path.name
    )
    edges = edges.copy()
    edges["u"] = edges["u"].map(norm_gene)
    edges["v"] = edges["v"].map(norm_gene)
    edges["edge_weight"] = pd.to_numeric(
        edges["edge_weight"], errors="coerce"
    )
    edges = edges[
        edges["u"].ne("")
        & edges["v"].ne("")
        & edges["edge_weight"].notna()
        & edges["edge_weight"].gt(0)
        & edges["u"].ne(edges["v"])
    ].copy()

    rwr = pd.read_csv(rwr_path)
    require_columns(
        rwr,
        [
            "gene",
            "syndrome_code",
            "restart_probability",
            "rwr_score",
        ],
        rwr_path.name,
    )
    rwr = rwr.copy()
    rwr["gene"] = rwr["gene"].map(norm_gene)
    rwr["syndrome_code"] = (
        rwr["syndrome_code"].astype(str).str.strip()
    )
    rwr["restart_probability"] = pd.to_numeric(
        rwr["restart_probability"], errors="coerce"
    )
    rwr["rwr_score"] = pd.to_numeric(
        rwr["rwr_score"], errors="coerce"
    )
    rwr = rwr[
        rwr["syndrome_code"].isin(SYNDROME_CODES)
        & np.isclose(
            rwr["restart_probability"],
            PRIMARY_RWR_RESTART,
        )
        & rwr["gene"].ne("")
        & rwr["rwr_score"].notna()
    ].copy()

    paths = {
        "frozen_target_seed": seed_path,
        "target_ubiquity": ubiq_path,
        "cap_meta_z": cap_path,
        "string_clean_edges": edge_path,
        "rwr_primary_scores": rwr_path,
    }

    return seeds, ubiq, cap, edges, rwr, paths, hashes


# =============================================================================
# 4. BUILD FROZEN TARGET UNIVERSE + MATCHED STRATA
# =============================================================================

def build_string_network(edges: pd.DataFrame) -> dict:
    nodes = sorted(set(edges["u"]) | set(edges["v"]))
    index = {g: i for i, g in enumerate(nodes)}

    i = edges["u"].map(index).to_numpy(dtype=int)
    j = edges["v"].map(index).to_numpy(dtype=int)
    w = edges["edge_weight"].to_numpy(dtype=float)

    A = coo_matrix(
        (
            np.concatenate([w, w]),
            (
                np.concatenate([i, j]),
                np.concatenate([j, i]),
            ),
        ),
        shape=(len(nodes), len(nodes)),
        dtype=float,
    ).tocsr()

    degree_weighted = np.asarray(A.sum(axis=1)).ravel()
    degree_unweighted = np.diff(A.indptr).astype(int)

    if np.any(degree_weighted <= 0):
        raise RuntimeError(
            "Zero weighted-degree node in STRING network."
        )

    P = (
        diags(
            1.0 / degree_weighted,
            offsets=0,
            format="csr",
        )
        @ A
    ).tocsr()

    row_sums = np.asarray(P.sum(axis=1)).ravel()
    max_row_error = float(
        np.max(np.abs(row_sums - 1.0))
    )
    if max_row_error > 1e-12:
        raise RuntimeError(
            "STRING transition matrix is not row-stochastic."
        )

    return {
        "nodes": nodes,
        "index": index,
        "P": P,
        "degree_unweighted": degree_unweighted,
        "degree_weighted": degree_weighted,
        "max_row_sum_error": max_row_error,
    }


def quantile_bin_preserve_ties(
    values: pd.Series,
    n_bins: int,
    prefix: str,
) -> pd.Series:
    """
    Value-based quantile bins that keep identical numeric values together.

    v1.0 used rank(method="first"), which could split tied STRING degrees or
    tied target-ubiquity values across adjacent bins according to row order.
    v1.0.1 removes that arbitrary tie-breaking. Bin sizes may therefore be
    unequal, which is preferable to separating identical matching values.
    """
    x = pd.to_numeric(values, errors="coerce")
    if x.isna().any():
        raise ValueError("Quantile matching values contain missing data.")
    if len(x) < n_bins:
        raise ValueError(
            f"Not enough values ({len(x)}) for {n_bins} bins."
        )

    codes = pd.qcut(
        x,
        q=n_bins,
        labels=False,
        duplicates="drop",
    )
    if codes.isna().any():
        raise ValueError("Could not construct quantile matching bins.")

    return codes.astype(int).map(
        lambda k: f"{prefix}{k + 1}"
    )


def build_target_universe(
    ubiq: pd.DataFrame,
    network: dict,
    cap: pd.DataFrame,
) -> pd.DataFrame:
    u = ubiq.rename(columns={"target_id": "gene"}).copy()

    deg = pd.DataFrame({
        "gene": network["nodes"],
        "network_degree_unweighted":
            network["degree_unweighted"],
        "network_degree_weighted":
            network["degree_weighted"],
    })

    u = u.merge(deg, on="gene", how="left")
    u["in_STRING"] = (
        u["network_degree_unweighted"].notna()
    )

    cap_primary = cap[
        ["gene", CAP_PRIMARY_Z_COLUMN]
    ].copy()

    if CAP_SENSITIVITY_Z_COLUMN in cap.columns:
        cap_primary[CAP_SENSITIVITY_Z_COLUMN] = (
            cap[CAP_SENSITIVITY_Z_COLUMN]
        )

    u = u.merge(
        cap_primary,
        on="gene",
        how="left",
    )
    u["measured_in_CAP_primary"] = (
        u[CAP_PRIMARY_Z_COLUMN].notna()
    )

    # Degree bins are defined among target-database genes represented in
    # STRING. Identical degree values are never split across bins.
    u["degree_bin"] = "NO_STRING"
    present = u["in_STRING"]
    u.loc[present, "degree_bin"] = (
        quantile_bin_preserve_ties(
            u.loc[present, "network_degree_unweighted"],
            N_DEGREE_BINS,
            "D",
        )
        .to_numpy()
    )

    # Target-herb ubiquity is a target-database property, so its bins are
    # defined over the ENTIRE frozen target universe, including genes absent
    # from STRING. This also prevents the seven NO_STRING genes from being
    # pooled across very different ubiquity values.
    u["ubiquity_bin"] = quantile_bin_preserve_ties(
        u["target_herb_degree"],
        N_UBIQUITY_BINS,
        "U",
    )

    u["matched_stratum"] = (
        u["degree_bin"].astype(str)
        + "_"
        + u["ubiquity_bin"].astype(str)
    )

    u = u.sort_values("gene").reset_index(drop=True)

    # Every frozen syndrome target must belong to this universe.
    return u


def build_weight_matrix(
    seeds: pd.DataFrame,
    universe: pd.DataFrame,
) -> np.ndarray:
    genes = universe["gene"].tolist()
    W = np.zeros(
        (len(SYNDROME_CODES), len(genes)),
        dtype=float,
    )

    for s_idx, syndrome in enumerate(SYNDROME_CODES):
        d = dict(
            zip(
                seeds.loc[
                    seeds["syndrome_code"].eq(syndrome),
                    "target_id",
                ],
                seeds.loc[
                    seeds["syndrome_code"].eq(syndrome),
                    TARGET_WEIGHT_COLUMN,
                ],
            )
        )

        outside = set(d) - set(genes)
        if outside:
            raise RuntimeError(
                f"{syndrome}: frozen targets absent from the frozen "
                f"target-database universe: {sorted(outside)[:20]}"
            )

        W[s_idx, :] = np.array(
            [float(d.get(g, 0.0)) for g in genes],
            dtype=float,
        )

        if abs(float(W[s_idx].sum()) - 1.0) > WEIGHT_SUM_TOL:
            raise RuntimeError(
                f"{syndrome}: complete target-universe weight vector "
                f"does not sum to 1."
            )

    return W


# =============================================================================
# 5. OBSERVED RAW-TARGET ALIGNMENT
# =============================================================================

def weighted_measured_alignment(
    weights: np.ndarray,
    z: np.ndarray,
    measured: np.ndarray,
) -> tuple[float, float]:
    measured_weight = float(weights @ measured)

    if measured_weight <= 0:
        return np.nan, measured_weight

    numerator = float(weights @ z)
    return numerator / measured_weight, measured_weight


def observed_raw_alignment(
    universe: pd.DataFrame,
    W: np.ndarray,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    measured = (
        universe["measured_in_CAP_primary"]
        .astype(float)
        .to_numpy()
    )
    z = (
        universe[CAP_PRIMARY_Z_COLUMN]
        .fillna(0.0)
        .to_numpy(dtype=float)
    )

    z_sens = None
    if CAP_SENSITIVITY_Z_COLUMN in universe.columns:
        z_sens = (
            universe[CAP_SENSITIVITY_Z_COLUMN]
            .fillna(0.0)
            .to_numpy(dtype=float)
        )

    rows = []
    contribution_parts = []

    for s_idx, syndrome in enumerate(SYNDROME_CODES):
        w = W[s_idx]
        observed, measured_mass = weighted_measured_alignment(
            w, z, measured
        )

        sens = np.nan
        if z_sens is not None:
            sens, _ = weighted_measured_alignment(
                w, z_sens, measured
            )

        positive_mask = (z > 0) & (measured > 0)
        negative_mask = (z < 0) & (measured > 0)

        up_component = (
            float(w @ np.where(positive_mask, z, 0.0))
            / measured_mass
        )
        down_component_magnitude = (
            float(
                w
                @ np.where(
                    negative_mask,
                    -z,
                    0.0,
                )
            )
            / measured_mass
        )

        rows.append({
            "method": "RAW_FROZEN_TARGET_PRIMARY",
            "syndrome_code": syndrome,
            "observed_alignment": observed,
            "CAP_measured_weight_mass": measured_mass,
            "CAP_unmeasured_weight_mass": 1.0 - measured_mass,
            "positive_Z_component": up_component,
            "negative_Z_component_magnitude":
                down_component_magnitude,
            "observed_alignment_neff_weighted_sensitivity":
                sens,
            "target_universe_n": int(len(universe)),
            "positive_weight_targets_total":
                int(np.sum(w > 0)),
            "positive_weight_targets_measured_in_CAP":
                int(np.sum((w > 0) & (measured > 0))),
        })

        c = universe[
            [
                "gene",
                "target_herb_degree",
                "target_idf",
                "network_degree_unweighted",
                "network_degree_weighted",
                "in_STRING",
                "matched_stratum",
                CAP_PRIMARY_Z_COLUMN,
                "measured_in_CAP_primary",
            ]
        ].copy()
        c["syndrome_code"] = syndrome
        c["frozen_target_weight"] = w
        c["measured_weight_normalized"] = np.where(
            measured > 0,
            w / measured_mass,
            0.0,
        )
        c["alignment_contribution"] = (
            c["measured_weight_normalized"]
            * c[CAP_PRIMARY_Z_COLUMN].fillna(0.0)
        )
        contribution_parts.append(c)

    return (
        pd.DataFrame(rows),
        pd.concat(contribution_parts, ignore_index=True),
    )


# =============================================================================
# 6. ADJOINT RWR FOR FAST, EXACT NULL PROPAGATION
# =============================================================================

def solve_adjoint_rwr(
    P: csr_matrix,
    vector: np.ndarray,
    restart: float,
) -> tuple[np.ndarray, dict]:
    """
    If stationary RWR is:
        p = r*p0 + (1-r)*P.T@p

    then for any network vector x:
        x.T p = p0.T q

    where:
        q = r*x + (1-r)*P@q

    Thus q can be solved once and used for thousands of null seed vectors.
    """
    x = np.asarray(vector, dtype=float).reshape(-1)
    if len(x) != P.shape[0]:
        raise ValueError(
            "Adjoint vector length does not match STRING network."
        )

    q = np.zeros_like(x)
    converged = False
    last_delta = np.nan

    for iteration in range(1, ADJOINT_MAX_ITER + 1):
        q_new = (
            restart * x
            + (1.0 - restart) * (P @ q)
        )

        last_delta = float(
            np.max(np.abs(q_new - q))
        )
        q = q_new

        if last_delta < ADJOINT_TOL_MAXABS:
            converged = True
            break

    if not converged:
        raise RuntimeError(
            "Adjoint RWR failed to converge. "
            f"Last max-abs delta={last_delta:.3e}"
        )

    return q, {
        "iterations": int(iteration),
        "final_max_abs_delta": last_delta,
        "converged": True,
    }


def build_adjoint_disease_vectors(
    network: dict,
    cap: pd.DataFrame,
) -> tuple[np.ndarray, np.ndarray, dict]:
    cap_z = dict(
        zip(
            cap["gene"],
            cap[CAP_PRIMARY_Z_COLUMN],
        )
    )

    z = np.zeros(len(network["nodes"]), dtype=float)
    measured = np.zeros(
        len(network["nodes"]), dtype=float
    )

    for i, gene in enumerate(network["nodes"]):
        if gene in cap_z:
            z[i] = float(cap_z[gene])
            measured[i] = 1.0

    q_z, qa_z = solve_adjoint_rwr(
        network["P"],
        z,
        PRIMARY_RWR_RESTART,
    )
    q_m, qa_m = solve_adjoint_rwr(
        network["P"],
        measured,
        PRIMARY_RWR_RESTART,
    )

    qa = {
        "CAP_genes_in_STRING": int(measured.sum()),
        "adjoint_Z_iterations": qa_z["iterations"],
        "adjoint_Z_final_max_abs_delta":
            qa_z["final_max_abs_delta"],
        "adjoint_measurement_iterations":
            qa_m["iterations"],
        "adjoint_measurement_final_max_abs_delta":
            qa_m["final_max_abs_delta"],
    }

    return q_z, q_m, qa


def map_adjoint_to_target_universe(
    universe: pd.DataFrame,
    network: dict,
    q_z: np.ndarray,
    q_m: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    qzu = np.zeros(len(universe), dtype=float)
    qmu = np.zeros(len(universe), dtype=float)

    for i, gene in enumerate(universe["gene"]):
        j = network["index"].get(gene)
        if j is not None:
            qzu[i] = q_z[j]
            qmu[i] = q_m[j]

    return qzu, qmu


def observed_rwr_alignment_from_frozen_scores(
    rwr: pd.DataFrame,
    cap: pd.DataFrame,
) -> pd.DataFrame:
    cap_small = cap[
        ["gene", CAP_PRIMARY_Z_COLUMN]
    ].copy()

    x = rwr.merge(
        cap_small,
        on="gene",
        how="left",
    )
    x["measured"] = (
        x[CAP_PRIMARY_Z_COLUMN].notna()
    ).astype(float)
    x["z0"] = x[CAP_PRIMARY_Z_COLUMN].fillna(0.0)

    rows = []
    for syndrome in SYNDROME_CODES:
        s = x[x["syndrome_code"].eq(syndrome)]
        total_mass = float(s["rwr_score"].sum())
        measured_mass = float(
            np.sum(
                s["rwr_score"]
                * s["measured"]
            )
        )
        numerator = float(
            np.sum(
                s["rwr_score"]
                * s["z0"]
            )
        )
        alignment = numerator / measured_mass

        rows.append({
            "method": "STRING_RWR_SECONDARY",
            "syndrome_code": syndrome,
            "observed_alignment": alignment,
            "RWR_total_probability_mass": total_mass,
            "CAP_measured_RWR_mass": measured_mass,
            "CAP_unmeasured_RWR_mass":
                total_mass - measured_mass,
            "RWR_network_nodes": int(len(s)),
            "CAP_measured_network_nodes":
                int(s["measured"].sum()),
        })

    return pd.DataFrame(rows)


def observed_rwr_alignment_adjoint(
    W: np.ndarray,
    qzu: np.ndarray,
    qmu: np.ndarray,
) -> pd.DataFrame:
    rows = []

    for s_idx, syndrome in enumerate(SYNDROME_CODES):
        w = W[s_idx]
        numerator = float(w @ qzu)
        denominator = float(w @ qmu)
        alignment = numerator / denominator

        rows.append({
            "syndrome_code": syndrome,
            "adjoint_alignment": alignment,
            "adjoint_CAP_measured_RWR_mass_scaled":
                denominator,
        })

    return pd.DataFrame(rows)


# =============================================================================
# 7. MATCHED-NULL PERMUTATION
# =============================================================================

def stratum_indices(
    universe: pd.DataFrame,
) -> dict[str, np.ndarray]:
    result = {}
    strata = universe["matched_stratum"].astype(str).to_numpy()
    for stratum in sorted(set(strata)):
        result[stratum] = np.flatnonzero(
            strata == stratum
        )
    return result


def null_invariant_qa(
    universe: pd.DataFrame,
    W: np.ndarray,
    strata: dict[str, np.ndarray],
) -> pd.DataFrame:
    """
    Deterministic QA of what the within-stratum permutation preserves.
    """
    rows = []
    for s_idx, syndrome in enumerate(SYNDROME_CODES):
        w = W[s_idx]

        for name, idx in strata.items():
            rows.append({
                "syndrome_code": syndrome,
                "matched_stratum": name,
                "stratum_n": int(len(idx)),
                "positive_weight_n": int(
                    np.sum(w[idx] > 0)
                ),
                "zero_weight_n": int(
                    np.sum(w[idx] == 0)
                ),
                "stratum_weight_sum": float(
                    w[idx].sum()
                ),
                "mean_network_degree_unweighted": float(
                    universe.loc[
                        idx,
                        "network_degree_unweighted",
                    ].fillna(0).mean()
                ),
                "mean_target_herb_degree": float(
                    universe.loc[
                        idx,
                        "target_herb_degree",
                    ].mean()
                ),
            })

    return pd.DataFrame(rows)


def generate_matched_null(
    universe: pd.DataFrame,
    W: np.ndarray,
    qzu: np.ndarray,
    qmu: np.ndarray,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(RANDOM_SEED)
    strata = stratum_indices(universe)

    measured = (
        universe["measured_in_CAP_primary"]
        .astype(float)
        .to_numpy()
    )
    z = (
        universe[CAP_PRIMARY_Z_COLUMN]
        .fillna(0.0)
        .to_numpy(dtype=float)
    )

    z_sens = None
    if CAP_SENSITIVITY_Z_COLUMN in universe.columns:
        z_sens = (
            universe[CAP_SENSITIVITY_Z_COLUMN]
            .fillna(0.0)
            .to_numpy(dtype=float)
        )

    n_s = len(SYNDROME_CODES)
    raw_null = np.empty((N_NULL, n_s), dtype=float)
    rwr_null = np.empty((N_NULL, n_s), dtype=float)

    raw_sens_null = (
        np.empty((N_NULL, n_s), dtype=float)
        if z_sens is not None else None
    )

    raw_measured_mass_null = np.empty(
        (N_NULL, n_s), dtype=float
    )
    rwr_measured_mass_scaled_null = np.empty(
        (N_NULL, n_s), dtype=float
    )

    base_index = np.arange(len(universe), dtype=int)

    for k in range(N_NULL):
        perm = base_index.copy()

        # One paired permutation shared by all syndromes.
        for idx in strata.values():
            perm[idx] = rng.permutation(idx)

        # w_perm[gene i] = original weight from permuted donor gene.
        WP = W[:, perm]

        raw_den = WP @ measured
        if np.any(raw_den <= 0):
            raise RuntimeError(
                "A null replicate has zero CAP-measured raw-target mass."
            )
        raw_num = WP @ z
        raw_null[k, :] = raw_num / raw_den
        raw_measured_mass_null[k, :] = raw_den

        if z_sens is not None:
            raw_sens_null[k, :] = (
                (WP @ z_sens) / raw_den
            )

        rwr_den = WP @ qmu
        if np.any(rwr_den <= 0):
            raise RuntimeError(
                "A null replicate has zero CAP-measured RWR mass."
            )
        rwr_num = WP @ qzu
        rwr_null[k, :] = rwr_num / rwr_den
        rwr_measured_mass_scaled_null[k, :] = rwr_den

    parts = []
    for s_idx, syndrome in enumerate(SYNDROME_CODES):
        d = pd.DataFrame({
            "null_iteration": np.arange(1, N_NULL + 1),
            "syndrome_code": syndrome,
            "raw_alignment_null":
                raw_null[:, s_idx],
            "rwr_alignment_null":
                rwr_null[:, s_idx],
            "raw_CAP_measured_weight_mass_null":
                raw_measured_mass_null[:, s_idx],
            "rwr_CAP_measured_mass_scaled_null":
                rwr_measured_mass_scaled_null[:, s_idx],
        })
        if raw_sens_null is not None:
            d["raw_alignment_neff_sensitivity_null"] = (
                raw_sens_null[:, s_idx]
            )
        parts.append(d)

    null_long = pd.concat(parts, ignore_index=True)

    return null_long, null_invariant_qa(
        universe, W, strata
    )


# =============================================================================
# 8. SUMMARIZE OBSERVED VS NULL
# =============================================================================

def summarize_null_results(
    raw_obs: pd.DataFrame,
    rwr_obs: pd.DataFrame,
    null_long: pd.DataFrame,
) -> pd.DataFrame:
    rows = []

    for syndrome in SYNDROME_CODES:
        ro = raw_obs[
            raw_obs["syndrome_code"].eq(syndrome)
        ].iloc[0]

        rn = null_long[
            null_long["syndrome_code"].eq(syndrome)
        ]["raw_alignment_null"].to_numpy()

        stats = empirical_p_values(
            rn,
            float(ro["observed_alignment"]),
        )
        rows.append({
            "method": "RAW_FROZEN_TARGET_PRIMARY",
            "syndrome_code": syndrome,
            "observed_alignment":
                float(ro["observed_alignment"]),
            "observed_CAP_measured_mass":
                float(ro["CAP_measured_weight_mass"]),
            "n_null": N_NULL,
            "null_model":
                "within_degree_x_target_ubiquity_stratum_weight_permutation",
            **stats,
        })

        rw = rwr_obs[
            rwr_obs["syndrome_code"].eq(syndrome)
        ].iloc[0]
        wn = null_long[
            null_long["syndrome_code"].eq(syndrome)
        ]["rwr_alignment_null"].to_numpy()

        stats = empirical_p_values(
            wn,
            float(rw["observed_alignment"]),
        )
        rows.append({
            "method": "STRING_RWR_SECONDARY",
            "syndrome_code": syndrome,
            "observed_alignment":
                float(rw["observed_alignment"]),
            "observed_CAP_measured_mass":
                float(rw["CAP_measured_RWR_mass"]),
            "n_null": N_NULL,
            "null_model":
                "same_matched_seed_permutation_then_exact_adjoint_RWR",
            **stats,
        })

    return pd.DataFrame(rows)


def summarize_neff_sensitivity(
    raw_obs: pd.DataFrame,
    null_long: pd.DataFrame,
) -> pd.DataFrame:
    if (
        "raw_alignment_neff_sensitivity_null"
        not in null_long.columns
    ):
        return pd.DataFrame()

    rows = []
    for syndrome in SYNDROME_CODES:
        obs = raw_obs[
            raw_obs["syndrome_code"].eq(syndrome)
        ][
            "observed_alignment_neff_weighted_sensitivity"
        ].iloc[0]

        null = null_long[
            null_long["syndrome_code"].eq(syndrome)
        ][
            "raw_alignment_neff_sensitivity_null"
        ].to_numpy()

        rows.append({
            "method":
                "RAW_TARGET_CAP_NEFF_META_Z_SENSITIVITY",
            "syndrome_code": syndrome,
            "observed_alignment": float(obs),
            "n_null": N_NULL,
            **empirical_p_values(null, float(obs)),
        })

    return pd.DataFrame(rows)


def paired_syndrome_difference_null(
    raw_obs: pd.DataFrame,
    rwr_obs: pd.DataFrame,
    null_long: pd.DataFrame,
) -> pd.DataFrame:
    """
    Descriptive paired-null comparison. It is NOT used to rank syndromes.
    """
    pairs = [
        ("PDOL", "PHOL"),
        ("PDOL", "WHIL"),
        ("PHOL", "WHIL"),
    ]
    rows = []

    raw_obs_map = dict(
        zip(
            raw_obs["syndrome_code"],
            raw_obs["observed_alignment"],
        )
    )
    rwr_obs_map = dict(
        zip(
            rwr_obs["syndrome_code"],
            rwr_obs["observed_alignment"],
        )
    )

    null_wide_raw = null_long.pivot(
        index="null_iteration",
        columns="syndrome_code",
        values="raw_alignment_null",
    )
    null_wide_rwr = null_long.pivot(
        index="null_iteration",
        columns="syndrome_code",
        values="rwr_alignment_null",
    )

    for method, obs_map, nw in [
        (
            "RAW_FROZEN_TARGET_PRIMARY",
            raw_obs_map,
            null_wide_raw,
        ),
        (
            "STRING_RWR_SECONDARY",
            rwr_obs_map,
            null_wide_rwr,
        ),
    ]:
        for a, b in pairs:
            obs_diff = float(
                obs_map[a] - obs_map[b]
            )
            null_diff = (
                nw[a].to_numpy()
                - nw[b].to_numpy()
            )

            rows.append({
                "method": method,
                "syndrome_A": a,
                "syndrome_B": b,
                "observed_A_minus_B": obs_diff,
                "n_null": N_NULL,
                **empirical_p_values(
                    null_diff,
                    obs_diff,
                ),
            })

    return pd.DataFrame(rows)


# =============================================================================
# 9. QA TABLES
# =============================================================================

def target_universe_qa(
    universe: pd.DataFrame,
    W: np.ndarray,
) -> pd.DataFrame:
    rows = []

    for s_idx, syndrome in enumerate(SYNDROME_CODES):
        w = W[s_idx]
        measured = (
            universe["measured_in_CAP_primary"]
            .to_numpy(dtype=bool)
        )
        in_string = (
            universe["in_STRING"]
            .to_numpy(dtype=bool)
        )

        rows.append({
            "syndrome_code": syndrome,
            "target_database_universe_n":
                int(len(universe)),
            "positive_weight_targets":
                int(np.sum(w > 0)),
            "positive_weight_targets_in_STRING":
                int(np.sum((w > 0) & in_string)),
            "positive_weight_targets_measured_CAP":
                int(np.sum((w > 0) & measured)),
            "positive_weight_targets_in_STRING_and_CAP":
                int(
                    np.sum(
                        (w > 0)
                        & in_string
                        & measured
                    )
                ),
            "weight_sum": float(w.sum()),
            "weight_mass_in_STRING":
                float(w @ in_string.astype(float)),
            "weight_mass_measured_CAP":
                float(w @ measured.astype(float)),
            "weight_mass_in_STRING_and_CAP":
                float(
                    w
                    @ (
                        in_string
                        & measured
                    ).astype(float)
                ),
        })

    return pd.DataFrame(rows)


def matched_strata_qa(
    universe: pd.DataFrame,
) -> pd.DataFrame:
    rows = []

    for name, x in universe.groupby(
        "matched_stratum",
        sort=True,
    ):
        rows.append({
            "matched_stratum": name,
            "n_genes": int(len(x)),
            "n_CAP_measured": int(
                x["measured_in_CAP_primary"].sum()
            ),
            "n_in_STRING": int(
                x["in_STRING"].sum()
            ),
            "degree_min": (
                float(
                    x["network_degree_unweighted"].min()
                )
                if x["in_STRING"].any()
                else np.nan
            ),
            "degree_median": (
                float(
                    x["network_degree_unweighted"].median()
                )
                if x["in_STRING"].any()
                else np.nan
            ),
            "degree_max": (
                float(
                    x["network_degree_unweighted"].max()
                )
                if x["in_STRING"].any()
                else np.nan
            ),
            "target_herb_degree_min":
                float(x["target_herb_degree"].min()),
            "target_herb_degree_median":
                float(x["target_herb_degree"].median()),
            "target_herb_degree_max":
                float(x["target_herb_degree"].max()),
        })

    return pd.DataFrame(rows)


# =============================================================================
# 10. MAIN
# =============================================================================

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("SYNDROME–DISEASE ALIGNMENT + MATCHED NULL v1.0.1")
    print("Frozen raw target = PRIMARY; STRING RWR = SECONDARY")
    print("=" * 80)

    print("\n[1/9] Resolving and verifying frozen inputs...")
    (
        seeds,
        ubiq,
        cap,
        edges,
        rwr,
        input_paths,
        input_hashes,
    ) = load_frozen_inputs()

    print(
        f"  Frozen target database universe: "
        f"{len(ubiq):,} genes"
    )
    print(
        f"  Frozen CAP continuous universe: "
        f"{len(cap):,} genes"
    )
    print(
        f"  Frozen STRING edges: "
        f"{len(edges):,}"
    )

    print("\n[2/9] Building frozen STRING transition matrix...")
    network = build_string_network(edges)
    print(f"  STRING nodes: {len(network['nodes']):,}")
    print(
        f"  max transition row-sum error: "
        f"{network['max_row_sum_error']:.3e}"
    )

    print("\n[3/9] Building target universe and matched strata...")
    universe = build_target_universe(
        ubiq, network, cap
    )
    W = build_weight_matrix(seeds, universe)

    universe_qa = target_universe_qa(
        universe, W
    )
    strata_qa = matched_strata_qa(universe)

    print(universe_qa.to_string(index=False))
    print(
        f"\n  matched strata: "
        f"{universe['matched_stratum'].nunique()}"
    )
    print(
        f"  smallest stratum: "
        f"{strata_qa['n_genes'].min()} genes"
    )

    if (
        universe_qa["weight_mass_measured_CAP"].min()
        < MIN_RECOMMENDED_MEASURED_MASS
    ):
        warnings.warn(
            "At least one syndrome has <50% frozen target weight "
            "measured in the CAP shared-gene universe."
        )

    print("\n[4/9] Computing observed primary raw-target alignment...")
    raw_obs, raw_contrib = observed_raw_alignment(
        universe, W
    )
    print(
        raw_obs[
            [
                "syndrome_code",
                "observed_alignment",
                "CAP_measured_weight_mass",
                "positive_weight_targets_measured_in_CAP",
            ]
        ].to_string(index=False)
    )

    print("\n[5/9] Computing secondary RWR alignment + exact cross-check...")
    q_z, q_m, adjoint_qa = build_adjoint_disease_vectors(
        network, cap
    )
    qzu, qmu = map_adjoint_to_target_universe(
        universe, network, q_z, q_m
    )

    rwr_obs = observed_rwr_alignment_from_frozen_scores(
        rwr, cap
    )
    rwr_adjoint = observed_rwr_alignment_adjoint(
        W, qzu, qmu
    )

    rwr_crosscheck = rwr_obs.merge(
        rwr_adjoint,
        on="syndrome_code",
        how="inner",
    )
    rwr_crosscheck["absolute_alignment_difference"] = (
        rwr_crosscheck["observed_alignment"]
        - rwr_crosscheck["adjoint_alignment"]
    ).abs()

    max_crosscheck_error = float(
        rwr_crosscheck[
            "absolute_alignment_difference"
        ].max()
    )
    if max_crosscheck_error > RWR_CROSSCHECK_TOL:
        raise RuntimeError(
            "Frozen RWR score alignment and adjoint RWR alignment "
            f"do not agree. Max difference={max_crosscheck_error:.3e}"
        )

    print(
        rwr_crosscheck[
            [
                "syndrome_code",
                "observed_alignment",
                "adjoint_alignment",
                "CAP_measured_RWR_mass",
                "absolute_alignment_difference",
            ]
        ].to_string(index=False)
    )

    print(
        "\n[6/9] Generating matched null "
        f"(N={N_NULL:,}, seed={RANDOM_SEED})..."
    )
    null_long, invariant_qa = generate_matched_null(
        universe, W, qzu, qmu
    )

    print("\n[7/9] Summarizing observed vs matched null...")
    summary = summarize_null_results(
        raw_obs, rwr_obs, null_long
    )
    neff_sensitivity = summarize_neff_sensitivity(
        raw_obs, null_long
    )
    pairwise = paired_syndrome_difference_null(
        raw_obs, rwr_obs, null_long
    )

    print(
        summary[
            [
                "method",
                "syndrome_code",
                "observed_alignment",
                "null_mean",
                "null_sd",
                "Z_align",
                "empirical_p_upper",
                "empirical_p_two_sided",
            ]
        ].to_string(index=False)
    )

    print("\n[8/9] Writing QA and reproducibility outputs...")

    # Frozen observed results.
    raw_obs.to_csv(
        OUTPUT_DIR
        / "CAP_alignment_observed_RAW_primary_v1_0_1.csv",
        index=False,
        encoding="utf-8-sig",
    )
    rwr_obs.to_csv(
        OUTPUT_DIR
        / "CAP_alignment_observed_RWR_secondary_v1_0_1.csv",
        index=False,
        encoding="utf-8-sig",
    )
    summary.to_csv(
        OUTPUT_DIR
        / "CAP_alignment_matched_null_summary_v1_0_1.csv",
        index=False,
        encoding="utf-8-sig",
    )
    pairwise.to_csv(
        OUTPUT_DIR
        / "CAP_alignment_paired_syndrome_difference_null_v1_0_1.csv",
        index=False,
        encoding="utf-8-sig",
    )
    neff_sensitivity.to_csv(
        OUTPUT_DIR
        / "CAP_alignment_neff_metaZ_sensitivity_v1_0_1.csv",
        index=False,
        encoding="utf-8-sig",
    )

    # Contributions / QA.
    raw_contrib.to_csv(
        OUTPUT_DIR
        / "CAP_alignment_RAW_gene_contributions_v1_0_1.csv",
        index=False,
        encoding="utf-8-sig",
    )
    universe.to_csv(
        OUTPUT_DIR
        / "matched_null_target_universe_and_strata_v1_0_1.csv",
        index=False,
        encoding="utf-8-sig",
    )
    universe_qa.to_csv(
        OUTPUT_DIR
        / "alignment_target_coverage_QA_v1_0_1.csv",
        index=False,
        encoding="utf-8-sig",
    )
    strata_qa.to_csv(
        OUTPUT_DIR
        / "matched_null_strata_QA_v1_0_1.csv",
        index=False,
        encoding="utf-8-sig",
    )
    invariant_qa.to_csv(
        OUTPUT_DIR
        / "matched_null_weight_invariant_QA_v1_0_1.csv",
        index=False,
        encoding="utf-8-sig",
    )
    rwr_crosscheck.to_csv(
        OUTPUT_DIR
        / "RWR_alignment_adjoint_crosscheck_QA_v1_0_1.csv",
        index=False,
        encoding="utf-8-sig",
    )

    if SAVE_FULL_NULL_DISTRIBUTION:
        null_long.to_csv(
            OUTPUT_DIR
            / "CAP_alignment_full_matched_null_distribution_v1_0_1.csv.gz",
            index=False,
            compression="gzip",
        )

    print("\n[9/9] Writing metadata and freeze-candidate statement...")

    metadata = {
        "analysis_name":
            "syndrome-disease alignment + matched null v1.0.1",
        "created_utc":
            datetime.now(timezone.utc).isoformat(
                timespec="seconds"
            ),
        "status":
            "FREEZE_CANDIDATE_REQUIRES_OUTPUT_AUDIT",
        "CAP_transcriptomics_used_to_tune_upstream": False,
        "primary_syndrome_representation":
            "LOPO ubiquity_corrected_sum1 frozen target weights",
        "primary_disease_representation":
            "CAP equal-weight signed-Z Stouffer meta-Z; all 12083 shared genes",
        "primary_alignment_formula":
            "sum(w_g * I_CAP(g) * Z_CAP,g) / sum(w_g * I_CAP(g))",
        "secondary_RWR_alignment_formula":
            "sum(p_g * I_CAP(g) * Z_CAP,g) / sum(p_g * I_CAP(g))",
        "RWR_restart_probability": PRIMARY_RWR_RESTART,
        "null": {
            "N": N_NULL,
            "random_seed": RANDOM_SEED,
            "candidate_universe":
                "frozen TCM-ID target-database universe",
            "degree_matching":
                f"{N_DEGREE_BINS} value-based quantile bins of STRING unweighted degree; tied values kept together",
            "target_ubiquity_matching":
                f"{N_UBIQUITY_BINS} value-based quantile bins of target_herb_degree over the full frozen target universe; tied values kept together",
            "no_STRING_policy":
                "NO_STRING is a separate degree category but remains stratified by target-ubiquity bin",
            "permutation":
                "entire syndrome weight vector including zeros permuted within each matched stratum",
            "paired_across_syndromes": True,
            "preserves": [
                "global weight multiset",
                "global positive/zero target count",
                "positive/zero target count within every stratum",
                "total target weight within every stratum",
                "target-database universe",
            ],
            "empirical_p_correction":
                "(count + 1)/(N + 1)",
            "primary_directional_test":
                "upper tail = positive CAP alignment beyond matched target background",
            "also_reported":
                [
                    "lower-tail empirical P",
                    "two-sided empirical P around null mean",
                    "Z_align standardized against null",
                ],
        },
        "RWR_null_acceleration": {
            "method":
                "exact adjoint stationary RWR identity",
            "stationary_RWR":
                "p = r*p0 + (1-r)*P.T@p",
            "adjoint":
                "q = r*x + (1-r)*P@q",
            "identity":
                "x.T@p = p0.T@q",
            "observed_crosscheck_max_abs_error":
                max_crosscheck_error,
            **adjoint_qa,
        },
        "target_database_universe_n":
            int(len(universe)),
        "CAP_genes_n":
            int(len(cap)),
        "STRING_nodes_n":
            int(len(network["nodes"])),
        "STRING_edges_n":
            int(len(edges)),
        "matched_strata_n":
            int(
                universe["matched_stratum"].nunique()
            ),
        "input_paths": {
            k: str(v)
            for k, v in input_paths.items()
        },
        "input_sha256": input_hashes,
        "null_stratification_revision":
            "v1.0.1 keeps tied matching values together and retains target-ubiquity matching for NO_STRING genes; this revision was triggered by QA of the null construction, not by significance optimization",
        "shared_specific_decomposition":
            "NOT INCLUDED; must be frozen upstream from observed-only syndrome-herb weights before later alignment version",
        "interpretation_guardrail":
            "Positive CAP meta-Z is disease-up relative to healthy control; alignment is disease-context localization, not therapeutic reversal.",
    }

    with open(
        OUTPUT_DIR
        / "syndrome_disease_alignment_matched_null_run_metadata_v1_0_1.json",
        "w",
        encoding="utf-8",
    ) as f:
        json.dump(
            metadata,
            f,
            ensure_ascii=False,
            indent=2,
        )

    statement = f"""SYNDROME–DISEASE ALIGNMENT + MATCHED NULL v1.0.1
FREEZE CANDIDATE — REQUIRES OUTPUT AUDIT

Primary syndrome representation:
LOPO frozen ubiquity-corrected target weights.

Primary CAP representation:
all {len(cap):,} frozen shared genes, equal-weight signed-Z Stouffer meta-Z.

Primary alignment:
A_raw(S) =
sum_g w_g * I_CAP(g) * Z_CAP,g
--------------------------------
sum_g w_g * I_CAP(g)

Secondary network-localization alignment:
STRING RWR r={PRIMARY_RWR_RESTART}

Matched null:
N = {N_NULL:,}
random seed = {RANDOM_SEED}
universe = frozen TCM-ID target-database universe ({len(universe)} genes)
matching = value-based STRING unweighted-degree quartile x target-herb-ubiquity quartile; tied values remain together
permutation = entire syndrome weight vector including zeros within each stratum
same stratum permutation shared across syndromes in each null replicate

Primary positive-alignment empirical P:
P_upper = (1 + # null >= observed)/(N + 1)

Also report:
P_lower, two-sided P around null mean, and Z_align.

RWR null:
computed with the exact adjoint stationary-RWR identity and cross-checked
against the frozen observed RWR score vectors.
Maximum observed RWR alignment cross-check error:
{max_crosscheck_error:.3e}

Important:
1. Do not tune target weights, STRING, r, CAP meta-Z, or matching bins after
   inspecting these results.
2. Positive CAP meta-Z means higher expression in CAP than healthy control.
3. Alignment is disease-context localization, not therapeutic reversal.
4. Raw frozen target alignment is primary.
5. RWR is secondary/sensitivity only.
6. Shared/specific decomposition is not reconstructed here and must not be
   invented from downstream CAP results.
7. Final freeze requires audit of coverage, null invariants, null variance,
   and observed-vs-null results.
"""
    (
        OUTPUT_DIR
        / "SYNDROME_DISEASE_ALIGNMENT_FREEZE_CANDIDATE_v1_0_1.txt"
    ).write_text(
        statement,
        encoding="utf-8",
    )

    print("\n" + "=" * 80)
    print("Alignment + matched-null computation completed.")
    print("STATUS: FREEZE CANDIDATE — audit before final freeze.")
    print("=" * 80)
    print(f"\nOutputs:\n{OUTPUT_DIR}")


if __name__ == "__main__":
    main()
