# -*- coding: utf-8 -*-
"""
run_projection_method_matched_null_ablation_v1_0.py

Projection-method × matched-null ablation for the frozen syndrome–CAP analysis.

PURPOSE
-------
This is a downstream, prespecified ablation. It asks whether the apparent CAP
alignment changes as the SAME frozen LOPO syndrome-herb evidence is projected
through increasingly bias-aware target representations:

    1) BINARY_UNION
       Every target linked to any mapped syndrome herb receives equal weight.

    2) WEIGHTED_RAW
       Target score = sum_h W_hs * A_ht.

    3) HERB_DEGREE_NORMALIZED
       Target score = sum_h W_hs * A_ht / |T_h|.

    4) UBIQUITY_CORRECTED
       Target score = sum_h W_hs * A_ht / |T_h| * IDF_t.
       This is the already frozen primary target projection.

The analysis does NOT select a projection method using CAP results. The primary
projection (LOPO UBIQUITY_CORRECTED) was frozen upstream before this ablation.
The other methods are diagnostic ablations only.

PRIMARY VS SENSITIVITY
----------------------
PRIMARY representation: LOPO (leave-pneumonia-out).
SENSITIVITY representation: FULL.

For comparability, every method-specific target vector is normalized to sum 1
before alignment. This does not change the weighted-mean alignment statistic;
it only makes CAP-measured mass directly comparable across methods.

CAP ALIGNMENT
-------------
For representation R, projection method M, syndrome S:

    A(R,M,S) = sum_g w_g * I_CAP(g) * Z_CAP,g
               --------------------------------
               sum_g w_g * I_CAP(g)

where Z_CAP is the frozen equal-weight CAP meta-Z.

MATCHED NULL
------------
The null reuses the EXACT 374-gene / 19-stratum target universe already frozen
in syndrome-disease alignment v1.0.1. Matching strata therefore are NOT rebuilt
or tuned here.

Within each null replicate, the SAME within-stratum target-label permutation is
applied simultaneously to:
    - all 4 projection methods;
    - all 3 syndromes;
    - both LOPO and FULL representations.

This preserves, for every representation × method × syndrome vector:
    - the complete 374-target universe;
    - the global weight multiset;
    - positive/zero target counts;
    - positive/zero counts within each frozen matched stratum;
    - total target weight within each frozen matched stratum.

Using one shared permutation also makes method-to-method contrasts paired.

INFERENCE
---------
For each representation × method × syndrome, report:
    Z_align = (A_obs - mean(A_null)) / sd(A_null)
    empirical upper/lower/two-sided P with +1 correction.

Prespecified paired method contrasts (AFTER - BEFORE):
    WEIGHTED_RAW - BINARY_UNION
    HERB_DEGREE_NORMALIZED - WEIGHTED_RAW
    UBIQUITY_CORRECTED - HERB_DEGREE_NORMALIZED
    UBIQUITY_CORRECTED - BINARY_UNION

These contrasts quantify how the alignment statistic changes across projection
choices relative to the SAME matched-background permutation. They are method
ablation diagnostics, not a post-hoc method-selection test.

MULTIPLICITY / INTERPRETATION
-----------------------------
This module is primarily an ablation/robustness analysis. Raw empirical P
values are reported for transparency. Holm-adjusted P values are additionally
reported across the four projection methods within each representation ×
syndrome family for the upper-tail tests, and across the four prespecified
method contrasts within each representation × syndrome family for two-sided
contrast tests.

Positive CAP meta-Z means higher expression in CAP than healthy control.
This analysis is disease-context localization, NOT therapeutic reversal and
NOT treatment efficacy.

HARD CROSS-CHECK
----------------
The LOPO UBIQUITY_CORRECTED null must reproduce the already frozen v1.0.1 RAW
complete-profile null exactly (within floating-point tolerance), because it
uses the same target vector, frozen strata, RNG seed, and permutation order.
If the frozen reference null is available, the script stops on mismatch.

No RWR is run here. RWR has already been frozen as a secondary network-
localization analysis and is not a projection-method ablation.
"""

from __future__ import annotations

from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import re
import warnings

import numpy as np
import pandas as pd


# =============================================================================
# 1. USER SETTINGS — EDIT PATHS HERE
# =============================================================================

PROJECT_ROOT = Path(r"E:\00project\Tanre Yongfei\重构")

# Frozen v1.1 target projection containing FULL and LOPO, all four methods.
PROJECTION_ALL_METHODS_FILE = (
    PROJECT_ROOT
    / "05_weighted_herb_target_projection_v1_1_FROZEN"
    / "syndrome_target_projection_all_methods_v1_1.csv"
)

# IMPORTANT: reuse the already frozen 374-target / 19-stratum universe.
# Do NOT rebuild matching bins in this script.
FROZEN_MATCHED_STRATA_FILE = (
    PROJECT_ROOT
    / "08_syndrome_disease_alignment_matched_null_v1_0_1"
    / "matched_null_target_universe_and_strata_v1_0_1.csv"
)

CAP_META_Z_FILE = (
    PROJECT_ROOT
    / "07_CAP_transcriptomic_metaZ_v1_0_1"
    / "CAP_CONTINUOUS_PRIMARY_equal_weight_metaZ_v1_0_1.csv"
)

# Reference null from the already frozen complete-profile analysis.
# Used ONLY as a hard reproducibility cross-check for LOPO UBIQUITY_CORRECTED.
FROZEN_COMPLETE_NULL_FILE = (
    PROJECT_ROOT
    / "08_syndrome_disease_alignment_matched_null_v1_0_1"
    / "CAP_alignment_full_matched_null_distribution_v1_0_1.csv.gz"
)

FROZEN_COMPLETE_SUMMARY_FILE = (
    PROJECT_ROOT
    / "08_syndrome_disease_alignment_matched_null_v1_0_1"
    / "CAP_alignment_matched_null_summary_v1_0_1.csv"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "11_projection_method_matched_null_ablation_v1_0"
)

# Frozen input hashes from the audited analysis branch.
ENFORCE_EXPECTED_SHA256 = True
EXPECTED_SHA256 = {
    "projection_all_methods":
        "ab5a99b6c51ea0518eb210843f40e048bf2055a4a7a76864a16f4fa744c2aeaa",
    "frozen_matched_strata":
        "a3871d4da3b3f4005833875f4e8d7c11e529214f6f93a95255f5acc5fffef297",
    "cap_meta_z":
        "168d9a23c06df388b69f73d92b2e7811b1f5332dfb507893511367fd752be58e",
    "frozen_complete_null":
        "7062555f8443581a46c3d328d77afa8d1a08095db411ea60e4e68f6cdf744f09",
    "frozen_complete_summary":
        "f2e67bc857e7d7d5055ac3b1a9062571b15e5cc07e028a6d0bd43644578a040a",
}

# -----------------------------------------------------------------------------
# Frozen / prespecified analysis settings
# -----------------------------------------------------------------------------
SYNDROME_CODES = ["PDOL", "PHOL", "WHIL"]
REPRESENTATIONS = ["LOPO", "FULL"]  # LOPO primary; FULL sensitivity.
PRIMARY_REPRESENTATION = "LOPO"
PRIMARY_METHOD = "UBIQUITY_CORRECTED"

METHOD_ORDER = [
    "BINARY_UNION",
    "WEIGHTED_RAW",
    "HERB_DEGREE_NORMALIZED",
    "UBIQUITY_CORRECTED",
]

# Columns already present in the frozen v1.1 projection file.
# Binary union is normalized to sum 1 inside this script.
METHOD_SOURCE_COLUMNS = {
    "BINARY_UNION": "binary_union",
    "WEIGHTED_RAW": "weighted_raw_sum1",
    "HERB_DEGREE_NORMALIZED": "herb_degree_normalized_sum1",
    "UBIQUITY_CORRECTED": "ubiquity_corrected_sum1",
}

# Prespecified staged ablation contrasts: AFTER - BEFORE.
METHOD_CONTRASTS = [
    ("WEIGHTED_RAW", "BINARY_UNION"),
    ("HERB_DEGREE_NORMALIZED", "WEIGHTED_RAW"),
    ("UBIQUITY_CORRECTED", "HERB_DEGREE_NORMALIZED"),
    ("UBIQUITY_CORRECTED", "BINARY_UNION"),
]

CAP_PRIMARY_Z_COLUMN = "Z_meta_equal"
CAP_SENSITIVITY_Z_COLUMN = "Z_meta_neff_weighted"

N_NULL = 10_000
RANDOM_SEED = 20260921

WEIGHT_SUM_TOL = 1e-10
NULL_CROSSCHECK_TOL = 1e-12
OBSERVED_CROSSCHECK_TOL = 1e-12

SAVE_FULL_NULL_DISTRIBUTION = True
RUN_EFFECTIVE_N_SENSITIVITY = True


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
    """Resolve preferred path or exact filename recursively under PROJECT_ROOT."""
    if preferred.exists():
        return preferred

    if not PROJECT_ROOT.exists():
        raise FileNotFoundError(
            f"{label}: preferred file not found:\n{preferred}\n"
            f"PROJECT_ROOT does not exist:\n{PROJECT_ROOT}"
        )

    matches = sorted(PROJECT_ROOT.rglob(preferred.name))
    if not matches:
        raise FileNotFoundError(
            f"{label}: could not find '{preferred.name}' under:\n{PROJECT_ROOT}"
        )

    if len(matches) == 1:
        print(f"  Auto-resolved {label}:\n    {matches[0]}")
        return matches[0]

    hashes = {str(p): sha256_file(p) for p in matches}
    unique_hashes = set(hashes.values())
    if len(unique_hashes) == 1:
        chosen = matches[0]
        print(f"  Multiple identical copies found for {label}; using:\n    {chosen}")
        return chosen

    details = "\n".join(f"  {p}\n    SHA256={h}" for p, h in hashes.items())
    raise RuntimeError(
        f"{label}: multiple non-identical files named '{preferred.name}' found.\n"
        f"Resolve manually; do not mix branches.\n{details}"
    )


def verify_hash(path: Path, key: str) -> str:
    observed = sha256_file(path)
    expected = EXPECTED_SHA256.get(key)
    if ENFORCE_EXPECTED_SHA256 and expected and observed != expected:
        raise RuntimeError(
            f"Frozen input SHA256 mismatch for {key}.\n"
            f"File: {path}\nExpected: {expected}\nObserved: {observed}\n"
            "Do not silently mix analysis branches. Create a new version if the "
            "input was intentionally changed."
        )
    return observed


def empirical_p_values(null_values: np.ndarray, observed: float) -> dict:
    null_values = np.asarray(null_values, dtype=float)
    if not np.isfinite(null_values).all():
        raise RuntimeError("Null distribution contains non-finite values.")

    n = len(null_values)
    mu = float(np.mean(null_values))
    sd = float(np.std(null_values, ddof=1))
    z_align = (float(observed) - mu) / sd if sd > 0 else np.nan

    p_upper = (1 + int(np.sum(null_values >= observed))) / (n + 1)
    p_lower = (1 + int(np.sum(null_values <= observed))) / (n + 1)
    obs_dev = abs(float(observed) - mu)
    p_two = (
        1 + int(np.sum(np.abs(null_values - mu) >= obs_dev))
    ) / (n + 1)

    return {
        "null_mean": mu,
        "null_sd": sd,
        "observed_minus_null_mean": float(observed) - mu,
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


def holm_adjust(p_values: pd.Series) -> pd.Series:
    """Holm step-down family-wise adjustment; preserves original index."""
    p = pd.to_numeric(p_values, errors="coerce")
    out = pd.Series(np.nan, index=p.index, dtype=float)
    valid = p.dropna()
    if valid.empty:
        return out

    order = valid.sort_values().index.tolist()
    m = len(order)
    running = 0.0
    for rank, idx in enumerate(order, start=1):
        adjusted = (m - rank + 1) * float(valid.loc[idx])
        running = max(running, adjusted)
        out.loc[idx] = min(running, 1.0)
    return out


# =============================================================================
# 3. LOAD FROZEN INPUTS
# =============================================================================

def load_inputs():
    projection_path = resolve_input(PROJECTION_ALL_METHODS_FILE, "projection all methods")
    strata_path = resolve_input(FROZEN_MATCHED_STRATA_FILE, "frozen matched strata")
    cap_path = resolve_input(CAP_META_Z_FILE, "CAP meta-Z")
    ref_null_path = resolve_input(FROZEN_COMPLETE_NULL_FILE, "frozen complete null")
    ref_summary_path = resolve_input(
        FROZEN_COMPLETE_SUMMARY_FILE, "frozen complete summary"
    )

    paths = {
        "projection_all_methods": projection_path,
        "frozen_matched_strata": strata_path,
        "cap_meta_z": cap_path,
        "frozen_complete_null": ref_null_path,
        "frozen_complete_summary": ref_summary_path,
    }
    hashes = {key: verify_hash(path, key) for key, path in paths.items()}

    projection = pd.read_csv(projection_path)
    require_columns(
        projection,
        [
            "representation", "syndrome_code", "target_id",
            "binary_union", "weighted_raw_sum1",
            "herb_degree_normalized_sum1", "ubiquity_corrected_sum1",
        ],
        projection_path.name,
    )
    projection = projection.copy()
    projection["representation"] = projection["representation"].astype(str).str.strip()
    projection["syndrome_code"] = projection["syndrome_code"].astype(str).str.strip()
    projection["target_id"] = projection["target_id"].map(norm_gene)

    for c in METHOD_SOURCE_COLUMNS.values():
        projection[c] = pd.to_numeric(projection[c], errors="coerce")

    projection = projection[
        projection["representation"].isin(REPRESENTATIONS)
        & projection["syndrome_code"].isin(SYNDROME_CODES)
        & projection["target_id"].ne("")
    ].copy()

    if projection.duplicated(["representation", "syndrome_code", "target_id"]).any():
        raise RuntimeError("Duplicate representation × syndrome × target rows in projection file.")

    strata = pd.read_csv(strata_path)
    require_columns(
        strata,
        [
            "gene", "target_herb_degree", "target_idf",
            "network_degree_unweighted", "in_STRING",
            "Z_meta_equal", "measured_in_CAP_primary",
            "degree_bin", "ubiquity_bin", "matched_stratum",
        ],
        strata_path.name,
    )
    strata = strata.copy()
    strata["gene"] = strata["gene"].map(norm_gene)
    strata = strata[strata["gene"].ne("")].copy()
    if strata["gene"].duplicated().any():
        raise RuntimeError("Duplicate genes in frozen matched-strata file.")
    if len(strata) != 374:
        raise RuntimeError(f"Frozen target universe must contain 374 genes; observed {len(strata)}.")
    if strata["matched_stratum"].nunique() != 19:
        raise RuntimeError(
            f"Frozen matched universe must contain 19 strata; observed "
            f"{strata['matched_stratum'].nunique()}."
        )

    # Sort by gene to reproduce the frozen v1.0.1 permutation order exactly.
    strata = strata.sort_values("gene").reset_index(drop=True)

    cap = pd.read_csv(cap_path)
    require_columns(cap, ["gene", CAP_PRIMARY_Z_COLUMN], cap_path.name)
    cap = cap.copy()
    cap["gene"] = cap["gene"].map(norm_gene)
    cap[CAP_PRIMARY_Z_COLUMN] = pd.to_numeric(cap[CAP_PRIMARY_Z_COLUMN], errors="coerce")
    if CAP_SENSITIVITY_Z_COLUMN in cap.columns:
        cap[CAP_SENSITIVITY_Z_COLUMN] = pd.to_numeric(
            cap[CAP_SENSITIVITY_Z_COLUMN], errors="coerce"
        )
    cap = cap[cap["gene"].ne("") & cap[CAP_PRIMARY_Z_COLUMN].notna()].copy()
    if cap["gene"].duplicated().any():
        raise RuntimeError("Duplicate genes in CAP meta-Z file.")

    # Hard cross-check: CAP values in frozen strata must equal frozen CAP file.
    cap_cols = ["gene", CAP_PRIMARY_Z_COLUMN]
    if CAP_SENSITIVITY_Z_COLUMN in cap.columns:
        cap_cols.append(CAP_SENSITIVITY_Z_COLUMN)
    check = strata[["gene", CAP_PRIMARY_Z_COLUMN]].merge(
        cap[cap_cols], on="gene", how="left", suffixes=("_strata", "_cap")
    )
    a = check[f"{CAP_PRIMARY_Z_COLUMN}_strata"].to_numpy(dtype=float)
    b = check[f"{CAP_PRIMARY_Z_COLUMN}_cap"].to_numpy(dtype=float)
    mask = np.isfinite(a) | np.isfinite(b)
    if not np.array_equal(np.isfinite(a), np.isfinite(b)):
        raise RuntimeError("CAP measurement mask differs between frozen strata and CAP file.")
    if mask.any() and np.nanmax(np.abs(a[mask] - b[mask])) > 1e-12:
        raise RuntimeError("CAP primary meta-Z values differ from frozen matched-strata branch.")

    ref_null = pd.read_csv(ref_null_path)
    require_columns(
        ref_null,
        ["null_iteration", "syndrome_code", "raw_alignment_null"],
        ref_null_path.name,
    )
    ref_null["syndrome_code"] = ref_null["syndrome_code"].astype(str).str.strip()

    ref_summary = pd.read_csv(ref_summary_path)
    require_columns(
        ref_summary,
        ["method", "syndrome_code", "observed_alignment", "null_mean", "Z_align"],
        ref_summary_path.name,
    )

    return projection, strata, cap, ref_null, ref_summary, paths, hashes


# =============================================================================
# 4. BUILD COMPLETE METHOD-SPECIFIC WEIGHT TENSORS
# =============================================================================

def build_weight_tensor(
    projection: pd.DataFrame,
    universe: pd.DataFrame,
) -> tuple[np.ndarray, pd.DataFrame]:
    """
    Return W with shape:
        [representation, method, syndrome, target]

    Every vector is normalized to sum 1.
    """
    genes = universe["gene"].tolist()
    gene_set = set(genes)
    W = np.zeros(
        (len(REPRESENTATIONS), len(METHOD_ORDER), len(SYNDROME_CODES), len(genes)),
        dtype=float,
    )

    qa_rows = []

    for r_idx, rep in enumerate(REPRESENTATIONS):
        for s_idx, syndrome in enumerate(SYNDROME_CODES):
            x = projection[
                projection["representation"].eq(rep)
                & projection["syndrome_code"].eq(syndrome)
            ].copy()

            outside = set(x["target_id"]) - gene_set
            if outside:
                raise RuntimeError(
                    f"{rep}/{syndrome}: projection contains targets outside frozen universe: "
                    f"{sorted(outside)[:20]}"
                )

            if set(x["target_id"]) != gene_set:
                missing = sorted(gene_set - set(x["target_id"]))
                raise RuntimeError(
                    f"{rep}/{syndrome}: projection does not contain complete 374-target universe. "
                    f"Missing examples: {missing[:20]}"
                )

            x = x.set_index("target_id").loc[genes]

            for m_idx, method in enumerate(METHOD_ORDER):
                source_col = METHOD_SOURCE_COLUMNS[method]
                values = x[source_col].to_numpy(dtype=float)
                if not np.isfinite(values).all():
                    raise RuntimeError(f"Non-finite weights in {rep}/{method}/{syndrome}.")
                if np.any(values < -1e-15):
                    raise RuntimeError(f"Negative target weights in {rep}/{method}/{syndrome}.")
                values = np.clip(values, 0.0, None)

                total_before = float(values.sum())
                if total_before <= 0:
                    raise RuntimeError(f"Zero total target weight in {rep}/{method}/{syndrome}.")
                values = values / total_before

                if abs(float(values.sum()) - 1.0) > WEIGHT_SUM_TOL:
                    raise RuntimeError(f"Weight normalization failed for {rep}/{method}/{syndrome}.")

                W[r_idx, m_idx, s_idx, :] = values
                qa_rows.append({
                    "representation": rep,
                    "method": method,
                    "syndrome_code": syndrome,
                    "source_column": source_col,
                    "source_weight_sum_before_normalization": total_before,
                    "normalized_weight_sum": float(values.sum()),
                    "positive_target_n": int(np.sum(values > 0)),
                    "zero_target_n": int(np.sum(values == 0)),
                })

    return W, pd.DataFrame(qa_rows)


def build_universe_weight_table(universe: pd.DataFrame, W: np.ndarray) -> pd.DataFrame:
    out = universe.copy()
    for r_idx, rep in enumerate(REPRESENTATIONS):
        for m_idx, method in enumerate(METHOD_ORDER):
            for s_idx, syndrome in enumerate(SYNDROME_CODES):
                out[f"W__{rep}__{method}__{syndrome}"] = W[r_idx, m_idx, s_idx, :]
    return out


# =============================================================================
# 5. OBSERVED ALIGNMENT
# =============================================================================

def disease_vectors(universe: pd.DataFrame, cap: pd.DataFrame):
    measured = universe["measured_in_CAP_primary"].astype(float).to_numpy()
    z_primary = universe[CAP_PRIMARY_Z_COLUMN].fillna(0.0).to_numpy(dtype=float)

    z_sens = None
    if RUN_EFFECTIVE_N_SENSITIVITY and CAP_SENSITIVITY_Z_COLUMN in cap.columns:
        sens_map = dict(zip(cap["gene"], cap[CAP_SENSITIVITY_Z_COLUMN]))
        z_sens = np.array(
            [float(sens_map.get(g, np.nan)) for g in universe["gene"]], dtype=float
        )
        sens_mask = np.isfinite(z_sens)
        if not np.array_equal(sens_mask, measured.astype(bool)):
            raise RuntimeError(
                "Effective-N CAP sensitivity measurement mask differs from primary CAP mask."
            )
        z_sens = np.nan_to_num(z_sens, nan=0.0)

    return measured, z_primary, z_sens


def compute_observed(W: np.ndarray, measured: np.ndarray, z: np.ndarray) -> pd.DataFrame:
    rows = []
    for r_idx, rep in enumerate(REPRESENTATIONS):
        for m_idx, method in enumerate(METHOD_ORDER):
            for s_idx, syndrome in enumerate(SYNDROME_CODES):
                w = W[r_idx, m_idx, s_idx, :]
                den = float(w @ measured)
                if den <= 0:
                    raise RuntimeError(f"Zero CAP-measured mass for {rep}/{method}/{syndrome}.")
                obs = float(w @ z) / den
                rows.append({
                    "representation": rep,
                    "analysis_role": "PRIMARY" if rep == PRIMARY_REPRESENTATION else "SENSITIVITY",
                    "method": method,
                    "method_order": METHOD_ORDER.index(method) + 1,
                    "syndrome_code": syndrome,
                    "observed_alignment": obs,
                    "CAP_measured_weight_mass": den,
                    "CAP_unmeasured_weight_mass": 1.0 - den,
                    "positive_target_n": int(np.sum(w > 0)),
                    "positive_target_measured_CAP_n": int(np.sum((w > 0) & (measured > 0))),
                })
    return pd.DataFrame(rows)


# =============================================================================
# 6. FROZEN MATCHED-NULL PERMUTATION
# =============================================================================

def stratum_indices(universe: pd.DataFrame) -> dict[str, np.ndarray]:
    strata = universe["matched_stratum"].astype(str).to_numpy()
    return {
        name: np.flatnonzero(strata == name)
        for name in sorted(set(strata))
    }


def invariant_qa(universe: pd.DataFrame, W: np.ndarray) -> pd.DataFrame:
    strata = stratum_indices(universe)
    rows = []
    for r_idx, rep in enumerate(REPRESENTATIONS):
        for m_idx, method in enumerate(METHOD_ORDER):
            for s_idx, syndrome in enumerate(SYNDROME_CODES):
                w = W[r_idx, m_idx, s_idx, :]
                for name, idx in strata.items():
                    rows.append({
                        "representation": rep,
                        "method": method,
                        "syndrome_code": syndrome,
                        "matched_stratum": name,
                        "stratum_n": int(len(idx)),
                        "positive_weight_n": int(np.sum(w[idx] > 0)),
                        "zero_weight_n": int(np.sum(w[idx] == 0)),
                        "stratum_weight_sum": float(w[idx].sum()),
                    })
    return pd.DataFrame(rows)


def generate_null(
    universe: pd.DataFrame,
    W: np.ndarray,
    measured: np.ndarray,
    z_primary: np.ndarray,
    z_sens: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray | None]:
    """
    Returns:
      null_primary shape [N_NULL, representation, method, syndrome]
      null_sens    same shape or None
    """
    rng = np.random.default_rng(RANDOM_SEED)
    strata = stratum_indices(universe)
    base_index = np.arange(len(universe), dtype=int)

    shape = (N_NULL, len(REPRESENTATIONS), len(METHOD_ORDER), len(SYNDROME_CODES))
    null_primary = np.empty(shape, dtype=float)
    null_sens = np.empty(shape, dtype=float) if z_sens is not None else None

    for k in range(N_NULL):
        perm = base_index.copy()
        # EXACT paired permutation shared by every method, syndrome, representation.
        for idx in strata.values():
            perm[idx] = rng.permutation(idx)

        WP = W[..., perm]  # [R, M, S, target]
        den = np.einsum("rmst,t->rms", WP, measured)
        if np.any(den <= 0):
            raise RuntimeError(f"Null replicate {k+1} has zero CAP-measured target mass.")

        num = np.einsum("rmst,t->rms", WP, z_primary)
        null_primary[k, ...] = num / den

        if z_sens is not None:
            num_sens = np.einsum("rmst,t->rms", WP, z_sens)
            null_sens[k, ...] = num_sens / den

        if (k + 1) % 1000 == 0:
            print(f"    completed null replicate {k+1:,}/{N_NULL:,}")

    return null_primary, null_sens


def null_to_long(null_array: np.ndarray, value_name: str) -> pd.DataFrame:
    parts = []
    iterations = np.arange(1, N_NULL + 1)
    for r_idx, rep in enumerate(REPRESENTATIONS):
        for m_idx, method in enumerate(METHOD_ORDER):
            for s_idx, syndrome in enumerate(SYNDROME_CODES):
                parts.append(pd.DataFrame({
                    "null_iteration": iterations,
                    "representation": rep,
                    "method": method,
                    "syndrome_code": syndrome,
                    value_name: null_array[:, r_idx, m_idx, s_idx],
                }))
    return pd.concat(parts, ignore_index=True)


# =============================================================================
# 7. SUMMARIES + PAIRED METHOD CONTRASTS
# =============================================================================

def summarize_observed_vs_null(observed: pd.DataFrame, null_array: np.ndarray) -> pd.DataFrame:
    rows = []
    for r_idx, rep in enumerate(REPRESENTATIONS):
        for m_idx, method in enumerate(METHOD_ORDER):
            for s_idx, syndrome in enumerate(SYNDROME_CODES):
                obs_row = observed[
                    observed["representation"].eq(rep)
                    & observed["method"].eq(method)
                    & observed["syndrome_code"].eq(syndrome)
                ].iloc[0]
                obs = float(obs_row["observed_alignment"])
                null = null_array[:, r_idx, m_idx, s_idx]
                rows.append({
                    "representation": rep,
                    "analysis_role": "PRIMARY" if rep == PRIMARY_REPRESENTATION else "SENSITIVITY",
                    "method": method,
                    "method_order": m_idx + 1,
                    "syndrome_code": syndrome,
                    "observed_alignment": obs,
                    "observed_CAP_measured_mass": float(obs_row["CAP_measured_weight_mass"]),
                    "n_null": N_NULL,
                    "null_model": "frozen_19_strata_same_target_label_permutation_all_methods_syndromes_representations",
                    **empirical_p_values(null, obs),
                })

    out = pd.DataFrame(rows)
    out["holm_p_upper_within_representation_syndrome"] = np.nan
    for (_, _), idx in out.groupby(["representation", "syndrome_code"]).groups.items():
        out.loc[idx, "holm_p_upper_within_representation_syndrome"] = holm_adjust(
            out.loc[idx, "empirical_p_upper"]
        )
    return out


def paired_method_contrasts(observed: pd.DataFrame, null_array: np.ndarray) -> pd.DataFrame:
    rows = []
    method_index = {m: i for i, m in enumerate(METHOD_ORDER)}

    for r_idx, rep in enumerate(REPRESENTATIONS):
        for s_idx, syndrome in enumerate(SYNDROME_CODES):
            obs_map = observed[
                observed["representation"].eq(rep)
                & observed["syndrome_code"].eq(syndrome)
            ].set_index("method")["observed_alignment"].to_dict()

            for after, before in METHOD_CONTRASTS:
                a_idx = method_index[after]
                b_idx = method_index[before]
                obs_diff = float(obs_map[after] - obs_map[before])
                null_diff = (
                    null_array[:, r_idx, a_idx, s_idx]
                    - null_array[:, r_idx, b_idx, s_idx]
                )
                rows.append({
                    "representation": rep,
                    "analysis_role": "PRIMARY" if rep == PRIMARY_REPRESENTATION else "SENSITIVITY",
                    "syndrome_code": syndrome,
                    "method_after": after,
                    "method_before": before,
                    "contrast": f"{after}_MINUS_{before}",
                    "observed_after_minus_before": obs_diff,
                    "n_null": N_NULL,
                    **empirical_p_values(null_diff, obs_diff),
                })

    out = pd.DataFrame(rows)
    out["holm_p_two_sided_within_representation_syndrome"] = np.nan
    for (_, _), idx in out.groupby(["representation", "syndrome_code"]).groups.items():
        out.loc[idx, "holm_p_two_sided_within_representation_syndrome"] = holm_adjust(
            out.loc[idx, "empirical_p_two_sided"]
        )
    return out


def make_primary_trajectory(summary: pd.DataFrame) -> pd.DataFrame:
    x = summary[summary["representation"].eq(PRIMARY_REPRESENTATION)].copy()
    x["method_order"] = pd.Categorical(
        x["method"], categories=METHOD_ORDER, ordered=True
    )
    x = x.sort_values(["syndrome_code", "method_order"]).reset_index(drop=True)
    x["method_order"] = x["method"].map({m: i + 1 for i, m in enumerate(METHOD_ORDER)})
    return x[
        [
            "representation", "syndrome_code", "method_order", "method",
            "observed_alignment", "null_mean", "observed_minus_null_mean",
            "Z_align", "empirical_p_upper", "empirical_p_two_sided",
            "holm_p_upper_within_representation_syndrome",
            "observed_CAP_measured_mass", "null_q025", "null_q975",
        ]
    ]


# =============================================================================
# 8. HARD CROSS-CHECK AGAINST FROZEN COMPLETE-PROFILE v1.0.1
# =============================================================================

def frozen_primary_crosscheck(
    observed: pd.DataFrame,
    null_primary: np.ndarray,
    ref_null: pd.DataFrame,
    ref_summary: pd.DataFrame,
) -> pd.DataFrame:
    r_idx = REPRESENTATIONS.index("LOPO")
    m_idx = METHOD_ORDER.index("UBIQUITY_CORRECTED")
    rows = []

    ref_summary_raw = ref_summary[
        ref_summary["method"].eq("RAW_FROZEN_TARGET_PRIMARY")
    ].copy()

    for s_idx, syndrome in enumerate(SYNDROME_CODES):
        new_obs = float(observed[
            observed["representation"].eq("LOPO")
            & observed["method"].eq("UBIQUITY_CORRECTED")
            & observed["syndrome_code"].eq(syndrome)
        ]["observed_alignment"].iloc[0])

        old_obs = float(ref_summary_raw[
            ref_summary_raw["syndrome_code"].eq(syndrome)
        ]["observed_alignment"].iloc[0])

        old_null = ref_null[
            ref_null["syndrome_code"].eq(syndrome)
        ].sort_values("null_iteration")["raw_alignment_null"].to_numpy(dtype=float)
        new_null = null_primary[:, r_idx, m_idx, s_idx]

        if len(old_null) != N_NULL:
            raise RuntimeError(
                f"Frozen reference null for {syndrome} has {len(old_null)} rows, expected {N_NULL}."
            )

        max_null_error = float(np.max(np.abs(new_null - old_null)))
        obs_error = abs(new_obs - old_obs)

        rows.append({
            "syndrome_code": syndrome,
            "new_observed": new_obs,
            "frozen_v1_0_1_observed": old_obs,
            "observed_abs_error": obs_error,
            "max_abs_null_iteration_error": max_null_error,
            "status": "PASS" if (
                obs_error <= OBSERVED_CROSSCHECK_TOL
                and max_null_error <= NULL_CROSSCHECK_TOL
            ) else "FAIL",
        })

    qa = pd.DataFrame(rows)
    if not qa["status"].eq("PASS").all():
        raise RuntimeError(
            "LOPO UBIQUITY_CORRECTED failed to reproduce frozen complete-profile v1.0.1.\n"
            + qa.to_string(index=False)
        )
    return qa


# =============================================================================
# 9. MAIN
# =============================================================================

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("[1/9] Loading and verifying frozen inputs...")
    projection, universe, cap, ref_null, ref_summary, paths, hashes = load_inputs()

    print("[2/9] Building complete 374-target weight tensor for all methods...")
    W, vector_qa = build_weight_tensor(projection, universe)
    universe_weights = build_universe_weight_table(universe, W)

    print("[3/9] Preparing frozen CAP disease vectors...")
    measured, z_primary, z_sens = disease_vectors(universe, cap)

    print("[4/9] Computing observed projection-method alignments...")
    observed = compute_observed(W, measured, z_primary)
    observed_sens = compute_observed(W, measured, z_sens) if z_sens is not None else pd.DataFrame()
    if not observed_sens.empty:
        observed_sens = observed_sens.rename(
            columns={"observed_alignment": "observed_alignment_neff_metaZ_sensitivity"}
        )

    print("[5/9] Running frozen 19-stratum paired matched null...")
    null_primary, null_sens = generate_null(universe, W, measured, z_primary, z_sens)

    print("[6/9] Hard cross-checking frozen LOPO ubiquity-corrected branch...")
    crosscheck = frozen_primary_crosscheck(
        observed, null_primary, ref_null, ref_summary
    )

    print("[7/9] Summarizing observed vs null and paired method contrasts...")
    summary = summarize_observed_vs_null(observed, null_primary)
    contrasts = paired_method_contrasts(observed, null_primary)
    trajectory = make_primary_trajectory(summary)

    sensitivity_summary = pd.DataFrame()
    if null_sens is not None and not observed_sens.empty:
        obs_for_sens = observed_sens.rename(
            columns={"observed_alignment_neff_metaZ_sensitivity": "observed_alignment"}
        )
        sensitivity_summary = summarize_observed_vs_null(obs_for_sens, null_sens)
        sensitivity_summary["CAP_metaZ_version"] = "effective_N_weighted_sensitivity"

    print("[8/9] Writing QA and reproducibility outputs...")
    invariant = invariant_qa(universe, W)

    vector_qa.to_csv(
        OUTPUT_DIR / "projection_method_ablation_weight_vector_QA_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )
    invariant.to_csv(
        OUTPUT_DIR / "projection_method_ablation_stratum_weight_invariant_QA_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )
    universe_weights.to_csv(
        OUTPUT_DIR / "projection_method_ablation_target_universe_weights_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )
    observed.to_csv(
        OUTPUT_DIR / "projection_method_ablation_observed_alignment_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )
    summary.to_csv(
        OUTPUT_DIR / "projection_method_ablation_matched_null_summary_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )
    contrasts.to_csv(
        OUTPUT_DIR / "projection_method_ablation_paired_method_contrasts_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )
    trajectory.to_csv(
        OUTPUT_DIR / "projection_method_ablation_LOPO_primary_trajectory_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )
    crosscheck.to_csv(
        OUTPUT_DIR / "projection_method_ablation_frozen_primary_crosscheck_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )

    if not sensitivity_summary.empty:
        sensitivity_summary.to_csv(
            OUTPUT_DIR / "projection_method_ablation_neff_metaZ_sensitivity_v1_0.csv",
            index=False, encoding="utf-8-sig"
        )

    if SAVE_FULL_NULL_DISTRIBUTION:
        null_long = null_to_long(null_primary, "alignment_null")
        if null_sens is not None:
            null_sens_long = null_to_long(null_sens, "alignment_neff_sensitivity_null")
            null_long = null_long.merge(
                null_sens_long,
                on=["null_iteration", "representation", "method", "syndrome_code"],
                how="left",
            )
        null_long.to_csv(
            OUTPUT_DIR / "projection_method_ablation_full_matched_null_distribution_v1_0.csv.gz",
            index=False, compression="gzip", encoding="utf-8-sig"
        )

    metadata = {
        "analysis_name": "Projection-method x matched-null ablation v1.0",
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "status": "DOWNSTREAM_PRESPECIFIED_ABLATION_DO_NOT_TUNE_TO_CAP",
        "primary_representation": PRIMARY_REPRESENTATION,
        "primary_projection_method_already_frozen_upstream": PRIMARY_METHOD,
        "projection_methods": METHOD_ORDER,
        "method_contrasts_after_minus_before": [
            {"after": a, "before": b} for a, b in METHOD_CONTRASTS
        ],
        "target_universe_n": int(len(universe)),
        "matched_strata_n": int(universe["matched_stratum"].nunique()),
        "n_null": N_NULL,
        "random_seed": RANDOM_SEED,
        "same_permutation_shared_across": [
            "all_projection_methods", "all_syndromes", "LOPO_and_FULL"
        ],
        "CAP_primary_Z": CAP_PRIMARY_Z_COLUMN,
        "CAP_sensitivity_Z": CAP_SENSITIVITY_Z_COLUMN if z_sens is not None else None,
        "multiplicity": {
            "method_upper_tail": "Holm across 4 methods within each representation x syndrome",
            "paired_method_contrasts": "Holm across 4 prespecified contrasts within each representation x syndrome",
        },
        "interpretation_guardrails": [
            "LOPO is primary; FULL is sensitivity only",
            "UBIQUITY_CORRECTED remains the frozen primary target projection regardless of ablation result",
            "Do not select a projection method based on CAP alignment",
            "Positive CAP alignment is disease-context localization, not therapeutic reversal or efficacy",
            "Matched-null context is required for interpretation",
            "No RWR is run in this module",
        ],
        "inputs": {
            key: {"path": str(paths[key]), "sha256": hashes[key]}
            for key in paths
        },
        "hard_crosscheck": {
            "LOPO_UBIQUITY_CORRECTED_reproduces_frozen_v1_0_1": True,
            "max_observed_abs_error": float(crosscheck["observed_abs_error"].max()),
            "max_null_iteration_abs_error": float(crosscheck["max_abs_null_iteration_error"].max()),
        },
    }
    with open(
        OUTPUT_DIR / "projection_method_ablation_run_metadata_v1_0.json",
        "w", encoding="utf-8"
    ) as f:
        json.dump(metadata, f, ensure_ascii=False, indent=2)

    freeze_candidate = f"""PROJECTION-METHOD x MATCHED-NULL ABLATION v1.0\nFREEZE CANDIDATE\n\nSTATUS\n------\nThis module is a downstream prespecified ablation.\nLOPO is primary; FULL is sensitivity.\nUBIQUITY_CORRECTED remains the already frozen primary target representation.\nNo CAP result may be used to select or redefine a projection method.\n\nMETHOD ORDER\n------------\n1. BINARY_UNION\n2. WEIGHTED_RAW\n3. HERB_DEGREE_NORMALIZED\n4. UBIQUITY_CORRECTED\n\nNULL\n----\nFrozen 374-target universe; frozen 19 matched strata; N={N_NULL}; seed={RANDOM_SEED}.\nThe same within-stratum target-label permutation is shared across all methods,\nall syndromes, and both LOPO/FULL representations.\n\nHARD REPRODUCIBILITY CROSS-CHECK\n--------------------------------\nLOPO UBIQUITY_CORRECTED reproduces frozen complete-profile v1.0.1.\nMaximum observed absolute error: {crosscheck['observed_abs_error'].max():.3e}\nMaximum null-iteration absolute error: {crosscheck['max_abs_null_iteration_error'].max():.3e}\n\nINTERPRETATION\n--------------\nUse this analysis to quantify sensitivity of apparent CAP localization to target\nprojection choices and matched target-background structure. It is not a method\nselection exercise, not therapeutic reversal, and not evidence of efficacy.\n"""
    with open(
        OUTPUT_DIR / "PROJECTION_METHOD_ABLATION_FREEZE_CANDIDATE_v1_0.txt",
        "w", encoding="utf-8"
    ) as f:
        f.write(freeze_candidate)

    # SHA256 manifest for all outputs except the manifest itself.
    output_files = sorted(
        p for p in OUTPUT_DIR.iterdir()
        if p.is_file() and p.name != "SHA256_MANIFEST_v1_0.csv"
    )
    manifest = pd.DataFrame([
        {
            "filename": p.name,
            "size_bytes": p.stat().st_size,
            "sha256": sha256_file(p),
        }
        for p in output_files
    ])
    manifest.to_csv(
        OUTPUT_DIR / "SHA256_MANIFEST_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )

    print("[9/9] Complete.")
    print(f"Outputs written to:\n  {OUTPUT_DIR}")
    print("\nFrozen primary cross-check:")
    print(crosscheck.to_string(index=False))
    print("\nLOPO primary trajectory:")
    print(trajectory.to_string(index=False))
    print(
        "\nIMPORTANT: Do not choose a projection method from these CAP results. "
        "UBIQUITY_CORRECTED remains the upstream-frozen primary representation."
    )


if __name__ == "__main__":
    main()
