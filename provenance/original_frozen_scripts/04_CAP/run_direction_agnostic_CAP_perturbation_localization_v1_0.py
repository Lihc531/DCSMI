# -*- coding: utf-8 -*-
"""
run_direction_agnostic_CAP_perturbation_localization_v1_0.py

Phase 12B — Direction-agnostic CAP perturbation localization v1.0
=================================================================

Scientific question
-------------------
Does the already-frozen prescription-derived molecular representation localize
preferentially to genes that are strongly perturbed in CAP, irrespective of
whether those genes are up- or down-regulated?

This analysis is downstream only. It MUST NOT modify any upstream syndrome
weights, herb-target mappings, target projection, shared/specific decomposition,
CAP meta-Z construction, or matched target strata.

Why |Z| rather than signed Z?
----------------------------
The frozen herb-target database encodes target association but not a validated
herb-specific activation/inhibition direction. Therefore the primary statistic
for this phase is direction-agnostic disease perturbation magnitude:

    L(c) = sum_g W_cg * I_CAP(g) * |Z_CAP,g|
           -----------------------------------
           sum_g W_cg * I_CAP(g)

where W_cg is the frozen raw ubiquity-corrected target weight for component c.
A larger L means that the component places more target weight on genes with
stronger CAP perturbation. It does NOT imply therapeutic reversal, beneficial
regulation, efficacy, or causal mechanism.

Prespecified primary test
-------------------------
PRIMARY PROFILE:
    SHARED

PRIMARY DISEASE SCORE:
    |Z_meta_equal| from the frozen two-cohort CAP meta-Z analysis.

PRIMARY NULL:
    Reuse verbatim the frozen 374-target universe and 19 matched strata from
    syndrome-disease alignment v1.0.1. Within each stratum, permute the complete
    target-weight vector INCLUDING zeros.

    The SAME donor permutation is applied simultaneously to all seven profiles:
        SHARED
        SPECIFIC_PDOL
        SPECIFIC_PHOL
        SPECIFIC_WHIL
        COMPLETE_PDOL
        COMPLETE_PHOL
        COMPLETE_WHIL

    This preserves each profile's total weight, global positive/zero count,
    within-stratum positive/zero count, and within-stratum weight mass.

PRIMARY ALTERNATIVE:
    L_shared,observed > matched-null expectation.

PRIMARY P VALUE:
    Upper-tail empirical P with +1 correction over 10,000 permutations.

The primary hypothesis is a SINGLE prespecified test and is not replaced by a
secondary profile if it fails.

Secondary analyses
------------------
1) The three COMPLETE LOPO profiles and three syndrome-SPECIFIC residuals are
   evaluated with the same |Z_meta_equal| statistic and same matched null.
   Their six upper-tail P values are Holm-adjusted as one secondary family.

2) Cohort-specific consistency for SHARED is evaluated using |Z_GSE103119| and
   |Z_GSE196399|. These cohorts are the same cohorts used to build the meta-Z,
   so this is NOT independent external replication. The two upper-tail P values
   are Holm-adjusted as a two-test consistency family.

3) |Z_meta_neff_weighted| is retained as a SHARED sensitivity analysis only.

4) Paired secondary contrasts use the same permutation replicate:
       SPECIFIC_s - SHARED  (3 contrasts)
       pairwise SPECIFIC differences (3 contrasts)
   and are Holm-adjusted as one six-contrast family. These are differential
   localization diagnostics, not proof of discrete syndrome mechanisms.

Hard legacy cross-check
-----------------------
To guarantee that Phase 12B uses the exact frozen matched-null machinery rather
than a subtly different permutation branch, the script simultaneously computes
SIGNED Z_meta_equal alignment under the same permutations and compares its
summary to the formally frozen Phase 8 (COMPLETE) and Phase 10
(SHARED/SPECIFIC) summaries. The |Z| result is accepted only if this legacy
cross-check passes within numerical tolerance.

Interpretation guardrails
-------------------------
- A significant SHARED upper-tail result supports preferential localization of
  the frozen shared target component to strongly perturbed CAP genes relative
  to degree- and target-ubiquity-matched TCM-ID target background.
- It does NOT establish that the shared component is the causal CAP mechanism.
- It does NOT establish therapeutic reversal or treatment efficacy.
- A non-significant primary result must not be rescued by selecting another
  profile, cohort, threshold, pathway, projection method, or network setting.
- Cohort-specific analyses are internal consistency checks, not independent
  replication, because both cohorts contribute to the frozen meta-Z.
"""

from __future__ import annotations

from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json

import numpy as np
import pandas as pd


# =============================================================================
# 1. USER SETTINGS — EDIT PATHS HERE ONLY IF NEEDED
# =============================================================================

PROJECT_ROOT = Path(r"E:\00project\Tanre Yongfei\重构")

SHARED_SPECIFIC_TARGET_FILE = (
    PROJECT_ROOT
    / "09_LOPO_shared_specific_projection_v1_0"
    / "LOPO_SHARED_SPECIFIC_primary_target_profiles_v1_0.csv"
)

COMPLETE_TARGET_FILE = (
    PROJECT_ROOT
    / "05_weighted_herb_target_projection_v1_1_FROZEN"
    / "LOPO_FROZEN_primary_target_weights_for_RWR_v1_1.csv"
)

FROZEN_TARGET_UNIVERSE_STRATA_FILE = (
    PROJECT_ROOT
    / "08_syndrome_disease_alignment_matched_null_v1_0_1"
    / "matched_null_target_universe_and_strata_v1_0_1.csv"
)

CAP_META_Z_FILE = (
    PROJECT_ROOT
    / "07_CAP_transcriptomic_metaZ_v1_0_1"
    / "CAP_CONTINUOUS_PRIMARY_equal_weight_metaZ_v1_0_1.csv"
)

# These two summaries are used ONLY as hard legacy cross-checks.
LEGACY_COMPLETE_SIGNED_SUMMARY_FILE = (
    PROJECT_ROOT
    / "08_syndrome_disease_alignment_matched_null_v1_0_1"
    / "CAP_alignment_matched_null_summary_v1_0_1.csv"
)

LEGACY_SHARED_SPECIFIC_SIGNED_SUMMARY_FILE = (
    PROJECT_ROOT
    / "10_shared_specific_CAP_alignment_matched_null_v1_0"
    / "CAP_shared_specific_alignment_matched_null_summary_v1_0.csv"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "12B_direction_agnostic_CAP_perturbation_localization_v1_0"
)

ENFORCE_EXPECTED_SHA256 = True

EXPECTED_SHA256 = {
    "shared_specific_target_profiles":
        "af66b3f49961b7a0e578e35b6c9bb86b7c84007f4c68d7efa4c286e6b9b85669",
    "complete_target_profiles":
        "9256cd3fe7ae142a91c7e081d97ba8a9225aa5aea3434c5f5673d53e5a10c9ef",
    "frozen_target_universe_strata":
        "a3871d4da3b3f4005833875f4e8d7c11e529214f6f93a95255f5acc5fffef297",
    "cap_meta_z":
        "168d9a23c06df388b69f73d92b2e7811b1f5332dfb507893511367fd752be58e",
    "legacy_complete_signed_summary":
        "f2e67bc857e7d7d5055ac3b1a9062571b15e5cc07e028a6d0bd43644578a040a",
    "legacy_shared_specific_signed_summary":
        "f60dbcf0e346922049f795b9d679dce2c3323f3b921b1a444acdfdf5cbc9765c",
}

# Frozen inference settings. Do not change after formal run.
N_NULL = 10_000
RANDOM_SEED = 20260921

RAW_WEIGHT_COLUMN = "ubiquity_corrected"
PRIMARY_Z_COLUMN = "Z_meta_equal"
COHORT_Z_COLUMNS = ["Z_GSE103119", "Z_GSE196399"]
SENSITIVITY_Z_COLUMN = "Z_meta_neff_weighted"

PROFILE_IDS = [
    "SHARED",
    "SPECIFIC_PDOL",
    "SPECIFIC_PHOL",
    "SPECIFIC_WHIL",
    "COMPLETE_PDOL",
    "COMPLETE_PHOL",
    "COMPLETE_WHIL",
]

PRIMARY_PROFILE = "SHARED"
SECONDARY_META_PROFILES = [
    "SPECIFIC_PDOL", "SPECIFIC_PHOL", "SPECIFIC_WHIL",
    "COMPLETE_PDOL", "COMPLETE_PHOL", "COMPLETE_WHIL",
]

LEGACY_TOL = 1e-10
WEIGHT_TOL = 1e-12


# =============================================================================
# 2. HELPERS
# =============================================================================

def norm_gene(x) -> str:
    if pd.isna(x):
        return ""
    return str(x).strip().upper()


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk_size)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def require_file(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing {label}: {path}")


def require_columns(df: pd.DataFrame, cols: list[str], label: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"{label} is missing columns: {missing}")


def check_sha(path: Path, key: str) -> str:
    observed = sha256_file(path)
    expected = EXPECTED_SHA256[key]
    if ENFORCE_EXPECTED_SHA256 and observed != expected:
        raise RuntimeError(
            f"Frozen input SHA256 mismatch for {key}.\n"
            f"File: {path}\nExpected: {expected}\nObserved: {observed}\n"
            "Do not silently mix analysis branches. Create a new version if the "
            "input was intentionally changed."
        )
    return observed


def empirical_summary(null_values: np.ndarray, observed: float) -> dict:
    x = np.asarray(null_values, dtype=float)
    if len(x) != N_NULL:
        raise RuntimeError(f"Expected {N_NULL} null values, observed {len(x)}.")
    if not np.isfinite(x).all():
        raise RuntimeError("Null distribution contains non-finite values.")

    mu = float(np.mean(x))
    sd = float(np.std(x, ddof=1))
    z = (float(observed) - mu) / sd if sd > 0 else np.nan
    p_upper = (1 + int(np.sum(x >= observed))) / (len(x) + 1)
    p_lower = (1 + int(np.sum(x <= observed))) / (len(x) + 1)
    obs_dev = abs(float(observed) - mu)
    p_two = (1 + int(np.sum(np.abs(x - mu) >= obs_dev))) / (len(x) + 1)

    return {
        "null_mean": mu,
        "null_sd": sd,
        "observed_minus_null_mean": float(observed) - mu,
        "Z_localization": z,
        "empirical_p_upper": p_upper,
        "empirical_p_lower": p_lower,
        "empirical_p_two_sided_centered": p_two,
        "null_q025": float(np.quantile(x, 0.025)),
        "null_q25": float(np.quantile(x, 0.25)),
        "null_median": float(np.quantile(x, 0.50)),
        "null_q75": float(np.quantile(x, 0.75)),
        "null_q975": float(np.quantile(x, 0.975)),
    }


def holm_adjust(values: pd.Series) -> pd.Series:
    p = pd.to_numeric(values, errors="coerce")
    out = pd.Series(np.nan, index=p.index, dtype=float)
    valid = p.dropna()
    if valid.empty:
        return out
    order = valid.sort_values().index.tolist()
    m = len(order)
    running = 0.0
    for rank, idx in enumerate(order, start=1):
        adj = (m - rank + 1) * float(valid.loc[idx])
        running = max(running, adj)
        out.loc[idx] = min(running, 1.0)
    return out


def profile_meta(profile_id: str) -> tuple[str, str]:
    if profile_id == "SHARED":
        return "SHARED", "ALL"
    if profile_id.startswith("SPECIFIC_"):
        return "SPECIFIC", profile_id.replace("SPECIFIC_", "")
    if profile_id.startswith("COMPLETE_"):
        return "COMPLETE", profile_id.replace("COMPLETE_", "")
    raise ValueError(profile_id)


def stratum_indices(universe: pd.DataFrame) -> dict[str, np.ndarray]:
    labels = universe["matched_stratum"].astype(str).to_numpy()
    return {
        s: np.flatnonzero(labels == s)
        for s in sorted(set(labels))
    }


# =============================================================================
# 3. LOAD + HARD-CHECK FROZEN INPUTS
# =============================================================================

def load_inputs():
    paths = {
        "shared_specific_target_profiles": SHARED_SPECIFIC_TARGET_FILE,
        "complete_target_profiles": COMPLETE_TARGET_FILE,
        "frozen_target_universe_strata": FROZEN_TARGET_UNIVERSE_STRATA_FILE,
        "cap_meta_z": CAP_META_Z_FILE,
        "legacy_complete_signed_summary": LEGACY_COMPLETE_SIGNED_SUMMARY_FILE,
        "legacy_shared_specific_signed_summary": LEGACY_SHARED_SPECIFIC_SIGNED_SUMMARY_FILE,
    }
    hashes = {}
    for key, path in paths.items():
        require_file(path, key)
        hashes[key] = check_sha(path, key)

    ss = pd.read_csv(SHARED_SPECIFIC_TARGET_FILE)
    require_columns(
        ss,
        ["component_id", "component_type", "syndrome_code", "target_id", RAW_WEIGHT_COLUMN],
        SHARED_SPECIFIC_TARGET_FILE.name,
    )
    ss = ss.copy()
    ss["target_id"] = ss["target_id"].map(norm_gene)
    ss[RAW_WEIGHT_COLUMN] = pd.to_numeric(ss[RAW_WEIGHT_COLUMN], errors="coerce")
    if ss[["component_id", "target_id"]].duplicated().any():
        raise RuntimeError("Duplicate component_id × target_id rows in shared/specific profile.")

    complete = pd.read_csv(COMPLETE_TARGET_FILE)
    require_columns(
        complete,
        ["syndrome_code", "target_id", RAW_WEIGHT_COLUMN],
        COMPLETE_TARGET_FILE.name,
    )
    complete = complete.copy()
    complete["target_id"] = complete["target_id"].map(norm_gene)
    complete[RAW_WEIGHT_COLUMN] = pd.to_numeric(complete[RAW_WEIGHT_COLUMN], errors="coerce")
    if complete[["syndrome_code", "target_id"]].duplicated().any():
        raise RuntimeError("Duplicate syndrome_code × target_id rows in complete profile.")

    universe = pd.read_csv(FROZEN_TARGET_UNIVERSE_STRATA_FILE)
    require_columns(
        universe,
        [
            "gene", "target_herb_degree", "network_degree_unweighted", "in_STRING",
            "measured_in_CAP_primary", "degree_bin", "ubiquity_bin", "matched_stratum",
            PRIMARY_Z_COLUMN,
        ],
        FROZEN_TARGET_UNIVERSE_STRATA_FILE.name,
    )
    universe = universe.copy()
    universe["gene"] = universe["gene"].map(norm_gene)
    if universe["gene"].duplicated().any():
        raise RuntimeError("Duplicate genes in frozen target universe.")
    if len(universe) != 374:
        raise RuntimeError(f"Frozen target universe must contain 374 genes; observed {len(universe)}.")
    if universe["matched_stratum"].nunique() != 19:
        raise RuntimeError(
            f"Frozen target universe must contain 19 strata; observed "
            f"{universe['matched_stratum'].nunique()}."
        )
    # Critical for exact reproduction of Phase 8/10 permutation order.
    universe = universe.sort_values("gene").reset_index(drop=True)

    cap = pd.read_csv(CAP_META_Z_FILE)
    require_columns(
        cap,
        ["gene", PRIMARY_Z_COLUMN, *COHORT_Z_COLUMNS, SENSITIVITY_Z_COLUMN],
        CAP_META_Z_FILE.name,
    )
    cap = cap.copy()
    cap["gene"] = cap["gene"].map(norm_gene)
    for c in [PRIMARY_Z_COLUMN, *COHORT_Z_COLUMNS, SENSITIVITY_Z_COLUMN]:
        cap[c] = pd.to_numeric(cap[c], errors="coerce")
    if cap["gene"].duplicated().any():
        raise RuntimeError("Duplicate genes in CAP meta-Z file.")

    # Merge cohort-specific Z values while hard-checking frozen primary Z/mask.
    merged = universe.merge(
        cap[["gene", PRIMARY_Z_COLUMN, *COHORT_Z_COLUMNS, SENSITIVITY_Z_COLUMN]],
        on="gene", how="left", suffixes=("_strata", "_cap")
    )
    a = merged[f"{PRIMARY_Z_COLUMN}_strata"].to_numpy(dtype=float)
    b = merged[f"{PRIMARY_Z_COLUMN}_cap"].to_numpy(dtype=float)
    if not np.array_equal(np.isfinite(a), np.isfinite(b)):
        raise RuntimeError("CAP measurement mask differs between frozen strata and CAP file.")
    mask = np.isfinite(a)
    if mask.any() and np.max(np.abs(a[mask] - b[mask])) > 1e-12:
        raise RuntimeError("Frozen strata Z_meta_equal differs from frozen CAP meta-Z file.")

    universe[PRIMARY_Z_COLUMN] = b
    # Cohort Z columns exist only in the CAP file; sensitivity Z exists in both
    # the frozen strata and CAP file and therefore receives merge suffixes.
    for c in COHORT_Z_COLUMNS:
        universe[c] = merged[c].to_numpy(dtype=float)
    sens_cap_col = (
        f"{SENSITIVITY_Z_COLUMN}_cap"
        if f"{SENSITIVITY_Z_COLUMN}_cap" in merged.columns
        else SENSITIVITY_Z_COLUMN
    )
    universe[SENSITIVITY_Z_COLUMN] = merged[sens_cap_col].to_numpy(dtype=float)
    if f"{SENSITIVITY_Z_COLUMN}_strata" in merged.columns:
        sa = merged[f"{SENSITIVITY_Z_COLUMN}_strata"].to_numpy(dtype=float)
        sb = universe[SENSITIVITY_Z_COLUMN].to_numpy(dtype=float)
        if not np.array_equal(np.isfinite(sa), np.isfinite(sb)):
            raise RuntimeError("CAP sensitivity-Z mask differs between frozen strata and CAP file.")
        sm = np.isfinite(sa)
        if sm.any() and np.max(np.abs(sa[sm] - sb[sm])) > 1e-12:
            raise RuntimeError("Frozen strata Z_meta_neff_weighted differs from CAP file.")

    measured = universe["measured_in_CAP_primary"].astype(bool).to_numpy()
    if not np.array_equal(measured, np.isfinite(universe[PRIMARY_Z_COLUMN].to_numpy(dtype=float))):
        raise RuntimeError("Frozen measured_in_CAP_primary mask is inconsistent with Z_meta_equal.")
    for c in COHORT_Z_COLUMNS + [SENSITIVITY_Z_COLUMN]:
        if not np.isfinite(universe.loc[measured, c].to_numpy(dtype=float)).all():
            raise RuntimeError(f"{c} is missing for genes in the frozen CAP shared-gene mask.")

    legacy_complete = pd.read_csv(LEGACY_COMPLETE_SIGNED_SUMMARY_FILE)
    legacy_ss = pd.read_csv(LEGACY_SHARED_SPECIFIC_SIGNED_SUMMARY_FILE)

    return ss, complete, universe, legacy_complete, legacy_ss, hashes


# =============================================================================
# 4. BUILD SEVEN FROZEN WEIGHT VECTORS ON THE 374-TARGET UNIVERSE
# =============================================================================

def build_weight_matrix(
    ss: pd.DataFrame,
    complete: pd.DataFrame,
    universe: pd.DataFrame,
) -> tuple[np.ndarray, pd.DataFrame]:
    genes = universe["gene"].tolist()
    gene_to_i = {g: i for i, g in enumerate(genes)}
    W = np.zeros((len(PROFILE_IDS), len(genes)), dtype=float)
    rows = []

    for p_idx, pid in enumerate(PROFILE_IDS):
        if pid == "SHARED" or pid.startswith("SPECIFIC_"):
            x = ss[ss["component_id"].eq(pid)].copy()
            if len(x) != 374:
                raise RuntimeError(f"{pid}: expected 374 zero-filled target rows; observed {len(x)}.")
        else:
            syndrome = pid.replace("COMPLETE_", "")
            x = complete[complete["syndrome_code"].eq(syndrome)].copy()

        unknown = sorted(set(x["target_id"]) - set(genes))
        if unknown:
            raise RuntimeError(f"{pid}: targets outside frozen 374-target universe: {unknown[:10]}")

        for r in x.itertuples(index=False):
            g = getattr(r, "target_id")
            w = float(getattr(r, RAW_WEIGHT_COLUMN))
            if not np.isfinite(w) or w < 0:
                raise RuntimeError(f"{pid}: invalid target weight for {g}: {w}")
            W[p_idx, gene_to_i[g]] = w

        if W[p_idx].sum() <= 0:
            raise RuntimeError(f"{pid}: non-positive total target weight.")

        family, syndrome = profile_meta(pid)
        rows.append({
            "profile_id": pid,
            "profile_family": family,
            "syndrome_code": syndrome,
            "target_universe_n": len(genes),
            "positive_weight_targets": int(np.sum(W[p_idx] > 0)),
            "zero_weight_targets": int(np.sum(W[p_idx] == 0)),
            "total_raw_target_weight": float(W[p_idx].sum()),
        })

    # Exact linear identity: COMPLETE = SHARED + matching SPECIFIC.
    for syndrome in ["PDOL", "PHOL", "WHIL"]:
        i_complete = PROFILE_IDS.index(f"COMPLETE_{syndrome}")
        i_shared = PROFILE_IDS.index("SHARED")
        i_specific = PROFILE_IDS.index(f"SPECIFIC_{syndrome}")
        err = float(np.max(np.abs(W[i_complete] - W[i_shared] - W[i_specific])))
        if err > 1e-12:
            raise RuntimeError(
                f"Frozen target linear identity failed for {syndrome}: max error={err:.3e}"
            )

    return W, pd.DataFrame(rows)


# =============================================================================
# 5. OBSERVED DIRECTION-AGNOSTIC LOCALIZATION
# =============================================================================

def score_vectors(universe: pd.DataFrame) -> dict[str, np.ndarray]:
    measured = universe["measured_in_CAP_primary"].astype(float).to_numpy()
    out = {
        "META_EQUAL_ABS_PRIMARY": measured * np.abs(
            universe[PRIMARY_Z_COLUMN].fillna(0.0).to_numpy(dtype=float)
        ),
        "GSE103119_ABS_COHORT_CONSISTENCY": measured * np.abs(
            universe["Z_GSE103119"].fillna(0.0).to_numpy(dtype=float)
        ),
        "GSE196399_ABS_COHORT_CONSISTENCY": measured * np.abs(
            universe["Z_GSE196399"].fillna(0.0).to_numpy(dtype=float)
        ),
        "META_NEFF_ABS_SENSITIVITY": measured * np.abs(
            universe[SENSITIVITY_Z_COLUMN].fillna(0.0).to_numpy(dtype=float)
        ),
        # Legacy signed score only for hard cross-check; never used as Phase 12B primary.
        "SIGNED_META_EQUAL_LEGACY_CROSSCHECK": measured * universe[
            PRIMARY_Z_COLUMN
        ].fillna(0.0).to_numpy(dtype=float),
    }
    return out


def observed_localization(
    universe: pd.DataFrame,
    W: np.ndarray,
    vectors: dict[str, np.ndarray],
) -> pd.DataFrame:
    measured = universe["measured_in_CAP_primary"].astype(float).to_numpy()
    rows = []
    for i, pid in enumerate(PROFILE_IDS):
        family, syndrome = profile_meta(pid)
        den = float(W[i] @ measured)
        total = float(W[i].sum())
        if den <= 0:
            raise RuntimeError(f"{pid}: zero CAP-measured target mass.")
        for score_name, v in vectors.items():
            num = float(W[i] @ v)
            rows.append({
                "profile_id": pid,
                "profile_family": family,
                "syndrome_code": syndrome,
                "score_name": score_name,
                "observed_localization": num / den,
                "localization_numerator": num,
                "CAP_measured_raw_target_mass": den,
                "total_raw_target_mass": total,
                "CAP_measured_fraction_of_raw_target_mass": den / total,
                "positive_weight_targets_total": int(np.sum(W[i] > 0)),
                "positive_weight_targets_measured_in_CAP": int(
                    np.sum((W[i] > 0) & (measured > 0))
                ),
            })
    return pd.DataFrame(rows)


# =============================================================================
# 6. MATCHED NULL — EXACT FROZEN 19-STRATA PERMUTATION BRANCH
# =============================================================================

def null_invariant_qa(
    universe: pd.DataFrame,
    W: np.ndarray,
) -> pd.DataFrame:
    strata = stratum_indices(universe)
    rows = []
    for i, pid in enumerate(PROFILE_IDS):
        for s, idx in strata.items():
            w = W[i, idx]
            rows.append({
                "profile_id": pid,
                "matched_stratum": s,
                "stratum_n": int(len(idx)),
                "positive_weight_n": int(np.sum(w > 0)),
                "zero_weight_n": int(np.sum(w == 0)),
                "stratum_raw_weight_sum": float(w.sum()),
                "mean_target_herb_degree": float(
                    universe.loc[idx, "target_herb_degree"].mean()
                ),
                "mean_network_degree_unweighted_zero_filled": float(
                    universe.loc[idx, "network_degree_unweighted"].fillna(0).mean()
                ),
            })
    return pd.DataFrame(rows)


def generate_null(
    universe: pd.DataFrame,
    W: np.ndarray,
    vectors: dict[str, np.ndarray],
) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    rng = np.random.default_rng(RANDOM_SEED)
    strata = stratum_indices(universe)
    measured = universe["measured_in_CAP_primary"].astype(float).to_numpy()
    base = np.arange(len(universe), dtype=int)

    score_names = list(vectors.keys())
    null_arrays = {
        name: np.empty((N_NULL, len(PROFILE_IDS)), dtype=float)
        for name in score_names
    }
    den_array = np.empty((N_NULL, len(PROFILE_IDS)), dtype=float)

    for k in range(N_NULL):
        perm = base.copy()
        for idx in strata.values():
            perm[idx] = rng.permutation(idx)

        # Same donor permutation simultaneously across all seven profiles.
        WP = W[:, perm]
        den = WP @ measured
        if np.any(den <= 0):
            raise RuntimeError(f"Null iteration {k+1}: zero CAP-measured target mass.")
        den_array[k, :] = den
        for name, v in vectors.items():
            null_arrays[name][k, :] = (WP @ v) / den

    parts = []
    for i, pid in enumerate(PROFILE_IDS):
        family, syndrome = profile_meta(pid)
        d = pd.DataFrame({
            "null_iteration": np.arange(1, N_NULL + 1),
            "profile_id": pid,
            "profile_family": family,
            "syndrome_code": syndrome,
            "CAP_measured_raw_target_mass_null": den_array[:, i],
        })
        for name in score_names:
            d[name] = null_arrays[name][:, i]
        parts.append(d)

    return pd.concat(parts, ignore_index=True), null_arrays


# =============================================================================
# 7. SUMMARIES, MULTIPLE TESTING, PAIRED CONTRASTS
# =============================================================================

def summarize_abs_results(
    observed: pd.DataFrame,
    null_arrays: dict[str, np.ndarray],
) -> pd.DataFrame:
    score_roles = {
        "META_EQUAL_ABS_PRIMARY": "PRIMARY_OR_SECONDARY_META",
        "GSE103119_ABS_COHORT_CONSISTENCY": "COHORT_CONSISTENCY",
        "GSE196399_ABS_COHORT_CONSISTENCY": "COHORT_CONSISTENCY",
        "META_NEFF_ABS_SENSITIVITY": "SENSITIVITY",
    }
    rows = []
    for score_name, role in score_roles.items():
        for i, pid in enumerate(PROFILE_IDS):
            obs = observed.loc[
                (observed["profile_id"].eq(pid))
                & (observed["score_name"].eq(score_name))
            ].iloc[0]
            inferential_role = role
            if score_name == "META_EQUAL_ABS_PRIMARY":
                inferential_role = "PRIMARY" if pid == PRIMARY_PROFILE else "SECONDARY_META"
            elif score_name in COHORT_SCORE_NAMES:
                inferential_role = (
                    "COHORT_CONSISTENCY_SHARED" if pid == PRIMARY_PROFILE
                    else "EXPLORATORY_COHORT_PROFILE"
                )
            elif score_name == "META_NEFF_ABS_SENSITIVITY":
                inferential_role = (
                    "SENSITIVITY_SHARED" if pid == PRIMARY_PROFILE
                    else "EXPLORATORY_SENSITIVITY_PROFILE"
                )

            rows.append({
                "profile_id": pid,
                "profile_family": obs["profile_family"],
                "syndrome_code": obs["syndrome_code"],
                "score_name": score_name,
                "inferential_role": inferential_role,
                "observed_localization": float(obs["observed_localization"]),
                "observed_CAP_measured_mass_fraction": float(
                    obs["CAP_measured_fraction_of_raw_target_mass"]
                ),
                "n_null": N_NULL,
                "null_model": "frozen_19_strata_same_donor_permutation_all_7_profiles",
                **empirical_summary(null_arrays[score_name][:, i], float(obs["observed_localization"])),
            })

    out = pd.DataFrame(rows)
    out["holm_p_upper_prespecified_family"] = np.nan

    # Six secondary meta profiles are one family. SHARED primary is excluded.
    idx = out.index[
        out["score_name"].eq("META_EQUAL_ABS_PRIMARY")
        & out["profile_id"].isin(SECONDARY_META_PROFILES)
    ]
    out.loc[idx, "holm_p_upper_prespecified_family"] = holm_adjust(
        out.loc[idx, "empirical_p_upper"]
    )

    # Two SHARED cohort-specific consistency checks are one family.
    idx = out.index[
        out["profile_id"].eq(PRIMARY_PROFILE)
        & out["score_name"].isin(COHORT_SCORE_NAMES)
    ]
    out.loc[idx, "holm_p_upper_prespecified_family"] = holm_adjust(
        out.loc[idx, "empirical_p_upper"]
    )
    return out


def paired_meta_contrasts(
    observed: pd.DataFrame,
    null_arrays: dict[str, np.ndarray],
) -> pd.DataFrame:
    score_name = "META_EQUAL_ABS_PRIMARY"
    arr = null_arrays[score_name]
    obs_map = observed[
        observed["score_name"].eq(score_name)
    ].set_index("profile_id")["observed_localization"].to_dict()

    contrasts = [
        ("SPECIFIC_PDOL", "SHARED"),
        ("SPECIFIC_PHOL", "SHARED"),
        ("SPECIFIC_WHIL", "SHARED"),
        ("SPECIFIC_PDOL", "SPECIFIC_PHOL"),
        ("SPECIFIC_PDOL", "SPECIFIC_WHIL"),
        ("SPECIFIC_PHOL", "SPECIFIC_WHIL"),
    ]
    rows = []
    for a, b in contrasts:
        ia, ib = PROFILE_IDS.index(a), PROFILE_IDS.index(b)
        obs_diff = float(obs_map[a] - obs_map[b])
        null_diff = arr[:, ia] - arr[:, ib]
        rows.append({
            "contrast": f"{a}_minus_{b}",
            "profile_A": a,
            "profile_B": b,
            "observed_difference": obs_diff,
            "n_null": N_NULL,
            "null_model": "paired_same_permutation_all_profiles",
            **empirical_summary(null_diff, obs_diff),
        })
    out = pd.DataFrame(rows)
    out["holm_p_two_sided_6_contrasts"] = holm_adjust(
        out["empirical_p_two_sided_centered"]
    )
    return out


# =============================================================================
# 8. HARD LEGACY SIGNED-Z CROSS-CHECK
# =============================================================================

def legacy_crosscheck(
    observed: pd.DataFrame,
    null_arrays: dict[str, np.ndarray],
    legacy_complete: pd.DataFrame,
    legacy_ss: pd.DataFrame,
) -> pd.DataFrame:
    signed_name = "SIGNED_META_EQUAL_LEGACY_CROSSCHECK"
    rows = []

    # Build reference lookup from frozen Phase 8 and Phase 10 summaries.
    refs = {}
    lc = legacy_complete[legacy_complete["method"].eq("RAW_FROZEN_TARGET_PRIMARY")].copy()
    for r in lc.itertuples(index=False):
        refs[f"COMPLETE_{r.syndrome_code}"] = r

    ls = legacy_ss[legacy_ss["method"].eq("RAW_SHARED_SPECIFIC_PRIMARY")].copy()
    for r in ls.itertuples(index=False):
        refs[str(r.component_id)] = r

    for i, pid in enumerate(PROFILE_IDS):
        if pid not in refs:
            raise RuntimeError(f"Legacy signed reference missing for {pid}.")
        ref = refs[pid]
        obs = float(observed.loc[
            observed["profile_id"].eq(pid)
            & observed["score_name"].eq(signed_name),
            "observed_localization",
        ].iloc[0])
        calc = empirical_summary(null_arrays[signed_name][:, i], obs)

        errors = {
            "observed_error": abs(obs - float(ref.observed_alignment)),
            "null_mean_error": abs(calc["null_mean"] - float(ref.null_mean)),
            "null_sd_error": abs(calc["null_sd"] - float(ref.null_sd)),
            "Z_error": abs(calc["Z_localization"] - float(ref.Z_align)),
            "p_upper_error": abs(calc["empirical_p_upper"] - float(ref.empirical_p_upper)),
            "p_lower_error": abs(calc["empirical_p_lower"] - float(ref.empirical_p_lower)),
            "p_two_error": abs(
                calc["empirical_p_two_sided_centered"] - float(ref.empirical_p_two_sided)
            ),
        }
        max_err = max(errors.values())
        rows.append({
            "profile_id": pid,
            **errors,
            "max_abs_crosscheck_error": max_err,
            "status": "PASS" if max_err <= LEGACY_TOL else "FAIL",
        })

    out = pd.DataFrame(rows)
    if not out["status"].eq("PASS").all():
        raise RuntimeError(
            "Legacy signed-Z cross-check FAILED. Phase 12B null branch does not "
            "exactly reproduce frozen Phase 8/10 matched-null inference."
        )
    return out


# =============================================================================
# 9. OUTPUT QA TABLES
# =============================================================================

def target_universe_output(
    universe: pd.DataFrame,
    W: np.ndarray,
) -> pd.DataFrame:
    out = universe[[
        "gene", "target_herb_degree", "network_degree_unweighted", "in_STRING",
        "degree_bin", "ubiquity_bin", "matched_stratum", "measured_in_CAP_primary",
        PRIMARY_Z_COLUMN, *COHORT_Z_COLUMNS, SENSITIVITY_Z_COLUMN,
    ]].copy()
    out["abs_Z_meta_equal"] = out[PRIMARY_Z_COLUMN].abs()
    out["abs_Z_GSE103119"] = out["Z_GSE103119"].abs()
    out["abs_Z_GSE196399"] = out["Z_GSE196399"].abs()
    out["abs_Z_meta_neff_weighted"] = out[SENSITIVITY_Z_COLUMN].abs()
    for i, pid in enumerate(PROFILE_IDS):
        out[f"W_{pid}"] = W[i]
    return out


def write_manifest(output_dir: Path) -> pd.DataFrame:
    rows = []
    for p in sorted(output_dir.iterdir()):
        if p.is_file() and p.name != "SHA256_MANIFEST_v1_0.csv":
            rows.append({
                "filename": p.name,
                "bytes": p.stat().st_size,
                "sha256": sha256_file(p),
            })
    m = pd.DataFrame(rows)
    m.to_csv(output_dir / "SHA256_MANIFEST_v1_0.csv", index=False, encoding="utf-8-sig")
    return m


# Global names used by summary helper.
COHORT_SCORE_NAMES = [
    "GSE103119_ABS_COHORT_CONSISTENCY",
    "GSE196399_ABS_COHORT_CONSISTENCY",
]


# =============================================================================
# 10. MAIN
# =============================================================================

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("[1/9] Loading and SHA-checking frozen inputs...")
    ss, complete, universe, legacy_complete, legacy_ss, input_hashes = load_inputs()

    print("[2/9] Building seven frozen target-weight vectors...")
    W, weight_qa = build_weight_matrix(ss, complete, universe)
    weight_qa.to_csv(
        OUTPUT_DIR / "phase12B_profile_weight_QA_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )

    print("[3/9] Building |Z| disease perturbation vectors...")
    vectors = score_vectors(universe)
    universe_out = target_universe_output(universe, W)
    universe_out.to_csv(
        OUTPUT_DIR / "phase12B_target_universe_profiles_CAP_absZ_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )

    print("[4/9] Computing observed localization...")
    observed = observed_localization(universe, W, vectors)
    observed.to_csv(
        OUTPUT_DIR / "phase12B_observed_localization_all_scores_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )

    print(f"[5/9] Running {N_NULL:,} frozen-strata matched-null permutations...")
    null_long, null_arrays = generate_null(universe, W, vectors)
    null_long.to_csv(
        OUTPUT_DIR / "phase12B_full_matched_null_distribution_v1_0.csv.gz",
        index=False, compression="gzip", encoding="utf-8-sig"
    )

    invariant_qa = null_invariant_qa(universe, W)
    invariant_qa.to_csv(
        OUTPUT_DIR / "phase12B_stratum_weight_invariant_QA_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )

    print("[6/9] Hard-crosschecking frozen signed-Z Phase 8/10 branch...")
    legacy_qa = legacy_crosscheck(
        observed, null_arrays, legacy_complete, legacy_ss
    )
    legacy_qa.to_csv(
        OUTPUT_DIR / "phase12B_signedZ_legacy_crosscheck_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )

    print("[7/9] Summarizing direction-agnostic primary and secondary results...")
    summary = summarize_abs_results(observed, null_arrays)
    summary.to_csv(
        OUTPUT_DIR / "phase12B_direction_agnostic_matched_null_summary_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )

    primary = summary[
        summary["profile_id"].eq(PRIMARY_PROFILE)
        & summary["score_name"].eq("META_EQUAL_ABS_PRIMARY")
    ].copy()
    if len(primary) != 1:
        raise RuntimeError("Primary SHARED × META_EQUAL_ABS result is not unique.")
    primary.to_csv(
        OUTPUT_DIR / "phase12B_PRIMARY_SHARED_meta_absZ_result_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )

    cohort_consistency = summary[
        summary["profile_id"].eq(PRIMARY_PROFILE)
        & summary["score_name"].isin(COHORT_SCORE_NAMES)
    ].copy()
    cohort_consistency.to_csv(
        OUTPUT_DIR / "phase12B_SHARED_cohort_consistency_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )

    paired = paired_meta_contrasts(observed, null_arrays)
    paired.to_csv(
        OUTPUT_DIR / "phase12B_paired_meta_absZ_profile_contrasts_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )

    print("[8/9] Writing metadata and result-freeze candidate...")
    primary_row = primary.iloc[0]
    metadata = {
        "analysis_name": "Phase 12B direction-agnostic CAP perturbation localization v1.0",
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "status": "FORMAL_RUN_RESULT_FREEZE_CANDIDATE",
        "primary_hypothesis": (
            "The frozen SHARED target component preferentially localizes to genes "
            "with larger absolute CAP equal-weight meta-Z than expected under the "
            "frozen degree × target-ubiquity matched TCM-ID target background."
        ),
        "primary_profile": PRIMARY_PROFILE,
        "primary_score": "weighted mean absolute Z_meta_equal over CAP-measured targets",
        "primary_alternative": "observed > matched-null expectation",
        "primary_p": "upper-tail empirical P with +1 correction",
        "n_null": N_NULL,
        "random_seed": RANDOM_SEED,
        "target_universe_n": int(len(universe)),
        "matched_strata_n": int(universe["matched_stratum"].nunique()),
        "CAP_measured_targets_in_universe": int(universe["measured_in_CAP_primary"].sum()),
        "null_permutation": (
            "same donor permutation within each frozen stratum applied simultaneously "
            "to all seven complete/shared/specific target-weight vectors including zeros"
        ),
        "cohort_note": (
            "GSE103119 and GSE196399 cohort-specific analyses are internal consistency "
            "checks, not independent external replication, because both cohorts form the meta-Z."
        ),
        "interpretation_guardrail": (
            "Absolute-Z localization measures disease perturbation magnitude only. "
            "It does not establish therapeutic reversal, efficacy, or causal mechanism."
        ),
        "input_sha256": input_hashes,
        "legacy_crosscheck_max_error": float(legacy_qa["max_abs_crosscheck_error"].max()),
        "primary_observed": float(primary_row["observed_localization"]),
        "primary_null_mean": float(primary_row["null_mean"]),
        "primary_Z_localization": float(primary_row["Z_localization"]),
        "primary_empirical_p_upper": float(primary_row["empirical_p_upper"]),
    }
    with open(
        OUTPUT_DIR / "phase12B_run_metadata_v1_0.json", "w", encoding="utf-8"
    ) as f:
        json.dump(metadata, f, ensure_ascii=False, indent=2)

    freeze_text = f"""PHASE 12B — DIRECTION-AGNOSTIC CAP PERTURBATION LOCALIZATION v1.0
RESULT FREEZE CANDIDATE

Primary question
----------------
Does the frozen SHARED target component preferentially localize to strongly
perturbed CAP genes, irrespective of direction, beyond the already frozen
STRING-degree × target-database-ubiquity matched TCM-ID target background?

Frozen primary statistic
------------------------
Weighted mean |Z_meta_equal| over CAP-measured targets using the frozen raw
ubiquity-corrected SHARED target weights.

Frozen primary null
-------------------
374-target universe; 19 frozen matched strata; N={N_NULL}; seed={RANDOM_SEED}.
Within each stratum the complete weight vector including zeros is permuted.
The same donor permutation is applied to all seven profiles in every replicate.

Frozen primary alternative
--------------------------
Observed SHARED localization > matched-null expectation.
Upper-tail empirical P with +1 correction.

Primary result
--------------
Observed = {float(primary_row['observed_localization']):.12g}
Null mean = {float(primary_row['null_mean']):.12g}
Null SD = {float(primary_row['null_sd']):.12g}
Z_localization = {float(primary_row['Z_localization']):.12g}
Upper-tail empirical P = {float(primary_row['empirical_p_upper']):.12g}

Legacy branch QA
----------------
Maximum absolute error versus frozen signed-Z Phase 8/10 summaries:
{float(legacy_qa['max_abs_crosscheck_error'].max()):.6e}

Interpretation guardrail
------------------------
This analysis is direction-agnostic. A positive result supports preferential
localization to strongly perturbed CAP genes relative to the matched target
background. It does NOT establish therapeutic reversal, treatment efficacy,
or a causal CAP mechanism. A non-significant primary result must not be rescued
by selecting a secondary profile, cohort, threshold, pathway, projection method,
or network setting.

Cohort-specific note
--------------------
GSE103119 and GSE196399 are internal cohort-consistency analyses, not independent
external replication, because both cohorts contribute to the frozen meta-Z.
"""
    with open(
        OUTPUT_DIR / "PHASE12B_RESULT_FREEZE_CANDIDATE_v1_0.txt",
        "w", encoding="utf-8"
    ) as f:
        f.write(freeze_text)

    print("[9/9] Writing SHA256 manifest...")
    manifest = write_manifest(OUTPUT_DIR)

    print("\n=== Phase 12B primary result ===")
    print(primary.to_string(index=False))
    print("\n=== SHARED cohort-specific consistency ===")
    print(cohort_consistency.to_string(index=False))
    print("\n=== Legacy signed-Z cross-check ===")
    print(legacy_qa.to_string(index=False))
    print(f"\nOutputs written to: {OUTPUT_DIR}")
    print(f"Manifest entries: {len(manifest)}")
    print(
        "\nIMPORTANT: Freeze the formal result before any pathway/module follow-up. "
        "Do not tune upstream weights, null strata, or disease scores based on this result."
    )


if __name__ == "__main__":
    main()
