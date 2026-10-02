# -*- coding: utf-8 -*-
"""
build_LOPO_shared_specific_projection_v1_0.py

Purpose
-------
Construct and freeze an observed-only shared/syndrome-specific decomposition
of the already frozen LOPO syndrome-herb representation, then project those
components through the already frozen herb-target transform.

IMPORTANT
---------
This script does NOT read CAP transcriptomics or any syndrome-disease alignment
result. It is an upstream decomposition script.

The decomposition rule was prespecified before downstream interpretation:

For every herb h, after forming the union of herbs observed in any LOPO
syndrome and setting unobserved syndrome-herb relations to ZERO:

    W_shared(h) = min_s W_hs

    W_specific(h, s) = W_hs - W_shared(h)

Therefore:
    W_hs = W_shared(h) + W_specific(h, s)

This is intentionally different from the old complete-matrix decomposition
that applied a Jeffreys posterior floor to unobserved syndrome-herb relations.
Unobserved relations are zero here. No unobserved treatment is inferred.

Frozen target projection
------------------------
The script reuses the v1.1 frozen target-layer mapping and database quantities:

    degree-normalized contribution:
        W_component(h) / |T_h|

    target IDF:
        IDF_t = log((N_H + 1)/(n_t + 1)) + 1

    ubiquity-corrected target score:
        U_t = sum_h [ W_component(h) / |T_h| * IDF_t ]

For profile comparison/alignment, each nonzero component is also normalized:

        U_sum1(t) = U_t / sum_t U_t

The UNNORMALIZED component masses are retained for magnitude diagnostics.

Critical QA
-----------
Before any shared/specific result is accepted, this script reconstructs the
original frozen LOPO target profiles from the same frozen inputs and requires
them to match:
    LOPO_FROZEN_primary_target_weights_for_RWR_v1_1.csv

It also requires:
1. exact herb-layer identity:
       original = shared + specific
2. target-layer linear identity BEFORE sum1 normalization:
       original projected score = shared projected score + specific projected score
3. degree-normalization mass conservation:
       sum target degree-normalized score = mapped herb-component mass

No CAP-driven tuning is permitted.
"""

from __future__ import annotations

from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import math
import re

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr


# =============================================================================
# 1. USER SETTINGS — EDIT PATHS HERE IF NEEDED
# =============================================================================

PROJECT_ROOT = Path(r"E:\00project\Tanre Yongfei\重构")

# The resolver searches recursively under PROJECT_ROOT by exact filename if
# these preferred paths do not exist.
LOPO_HERB_WEIGHTS_FILE = (
    PROJECT_ROOT
    / "LOPO_syndrome_herb_weights_observed_v1_0.csv"
)

FROZEN_MAPPING_DICTIONARY_FILE = (
    PROJECT_ROOT
    / "FROZEN_herb_to_target_mapping_dictionary_v1_1.csv"
)

CANONICAL_HERB_TARGET_PAIRS_FILE = (
    PROJECT_ROOT
    / "canonical_herb_target_pairs_v1_1.csv"
)

TARGET_UBIQUITY_FILE = (
    PROJECT_ROOT
    / "target_database_target_ubiquity_v1_1.csv"
)

HERB_TARGET_DEGREE_FILE = (
    PROJECT_ROOT
    / "target_database_herb_degree_v1_1.csv"
)

FROZEN_LOPO_TARGET_PROFILE_FILE = (
    PROJECT_ROOT
    / "LOPO_FROZEN_primary_target_weights_for_RWR_v1_1.csv"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "09_LOPO_shared_specific_projection_v1_0"
)

# Frozen input hashes. These prevent accidental branch mixing.
ENFORCE_EXPECTED_SHA256 = True

EXPECTED_SHA256 = {
    "LOPO_herb_weights":
        "4470442be38c5a7c92f191f790ee9d2d5799b707e8526b2f0d7239f7eec84d11",
    "frozen_mapping_dictionary":
        "709fac26761667242c5754e4bc7c67be047f9398aea010542e25d195eea95e66",
    "canonical_herb_target_pairs":
        "a6b283fe945885046b78141a94b6474d1acdef41a4ace8e941096fee95dff983",
    "target_ubiquity":
        "7a714cbfadd44033aabcc7178d4ff045bbc826f6cc5606705acc76fcad002bba",
    "herb_target_degree":
        "8f3e2031674a959e5f6af57b5506c54d4842e7e7893db23f0659804a07964bce",
    "frozen_LOPO_target_profile":
        "9256cd3fe7ae142a91c7e081d97ba8a9225aa5aea3434c5f5673d53e5a10c9ef",
}

SYNDROME_CODES = ["PDOL", "PHOL", "WHIL"]

# Frozen target projection field.
PRIMARY_TARGET_SCORE = "ubiquity_corrected"
PRIMARY_TARGET_SCORE_SUM1 = "ubiquity_corrected_sum1"

# Numerical QA tolerances.
HERB_IDENTITY_TOL = 1e-14
TARGET_RECONSTRUCTION_TOL = 1e-12
TARGET_LINEAR_IDENTITY_TOL = 1e-12
MASS_CONSERVATION_TOL = 1e-12
NORMALIZED_SUM_TOL = 1e-12

TOP_K_VALUES = [20, 50, 100]


# =============================================================================
# 2. HELPERS
# =============================================================================

def norm_text(x) -> str:
    if pd.isna(x):
        return ""
    return re.sub(r"\s+", "", str(x).strip())


def norm_gene(x) -> str:
    return norm_text(x).upper()


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
    Use preferred path if present; otherwise recursively search PROJECT_ROOT
    for the exact filename.

    Multiple identical copies are allowed. Multiple non-identical copies stop
    the run rather than silently selecting an analysis branch.
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
    if len(set(hashes.values())) == 1:
        chosen = matches[0]
        print(
            f"  Multiple identical copies found for {label}; using:\n"
            f"    {chosen}"
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
    expected = EXPECTED_SHA256[key]
    if ENFORCE_EXPECTED_SHA256 and observed != expected:
        raise RuntimeError(
            f"Frozen input SHA256 mismatch for {key}.\n"
            f"File: {path}\n"
            f"Expected: {expected}\n"
            f"Observed: {observed}\n"
            "Do not mix analysis branches. Create a new version if the "
            "upstream frozen input is intentionally changed."
        )
    return observed


def safe_corr(x: np.ndarray, y: np.ndarray, method: str) -> float:
    if np.allclose(np.std(x), 0) or np.allclose(np.std(y), 0):
        return np.nan
    if method == "spearman":
        return float(spearmanr(x, y).statistic)
    if method == "pearson":
        return float(pearsonr(x, y).statistic)
    raise ValueError(method)


# =============================================================================
# 3. LOAD AND VERIFY FROZEN INPUTS
# =============================================================================

def load_inputs():
    paths = {
        "LOPO_herb_weights":
            resolve_input(LOPO_HERB_WEIGHTS_FILE, "LOPO observed herb weights"),
        "frozen_mapping_dictionary":
            resolve_input(FROZEN_MAPPING_DICTIONARY_FILE, "frozen herb mapping dictionary"),
        "canonical_herb_target_pairs":
            resolve_input(CANONICAL_HERB_TARGET_PAIRS_FILE, "canonical herb-target pairs"),
        "target_ubiquity":
            resolve_input(TARGET_UBIQUITY_FILE, "target ubiquity"),
        "herb_target_degree":
            resolve_input(HERB_TARGET_DEGREE_FILE, "herb target degree"),
        "frozen_LOPO_target_profile":
            resolve_input(FROZEN_LOPO_TARGET_PROFILE_FILE, "frozen LOPO target profile"),
    }

    hashes = {k: verify_hash(v, k) for k, v in paths.items()}

    w = pd.read_csv(paths["LOPO_herb_weights"])
    require_columns(
        w,
        [
            "representation",
            "canonical_name",
            "syndrome_code",
            "primary_weight_W",
            "observed_relation",
        ],
        paths["LOPO_herb_weights"].name,
    )
    w = w.copy()
    w["canonical_name"] = w["canonical_name"].map(norm_text)
    w["syndrome_code"] = w["syndrome_code"].astype(str).str.strip()
    w["primary_weight_W"] = pd.to_numeric(
        w["primary_weight_W"], errors="coerce"
    )
    w["observed_relation"] = pd.to_numeric(
        w["observed_relation"], errors="coerce"
    )
    w = w[
        w["representation"].eq("LOPO")
        & w["syndrome_code"].isin(SYNDROME_CODES)
        & w["canonical_name"].ne("")
        & w["primary_weight_W"].notna()
        & w["primary_weight_W"].gt(0)
        & w["observed_relation"].eq(1)
    ].copy()

    if w.duplicated(["canonical_name", "syndrome_code"]).any():
        dup = w[w.duplicated(
            ["canonical_name", "syndrome_code"], keep=False
        )]
        raise RuntimeError(
            "Duplicate observed LOPO herb-syndrome rows:\n"
            f"{dup.head(20).to_string(index=False)}"
        )

    mapping = pd.read_csv(paths["frozen_mapping_dictionary"])
    require_columns(
        mapping,
        [
            "representation",
            "syndrome_code",
            "canonical_name",
            "target_mapping_name",
            "mapping_method",
            "mapped_to_target_db",
        ],
        paths["frozen_mapping_dictionary"].name,
    )
    mapping = mapping[mapping["representation"].eq("LOPO")].copy()
    mapping["canonical_name"] = mapping["canonical_name"].map(norm_text)
    mapping["target_mapping_name"] = mapping["target_mapping_name"].map(norm_text)
    mapping["mapped_to_target_db"] = pd.to_numeric(
        mapping["mapped_to_target_db"], errors="coerce"
    ).fillna(0).astype(int)

    # The frozen LOPO mapping must be syndrome-independent for the same herb.
    mapping_unique = (
        mapping[
            [
                "canonical_name",
                "target_mapping_name",
                "mapping_method",
                "mapped_to_target_db",
            ]
        ]
        .drop_duplicates()
        .copy()
    )
    if mapping_unique["canonical_name"].duplicated().any():
        bad = mapping_unique[
            mapping_unique["canonical_name"].duplicated(keep=False)
        ]
        raise RuntimeError(
            "Conflicting frozen target mappings for the same LOPO herb:\n"
            f"{bad.to_string(index=False)}"
        )

    pairs = pd.read_csv(paths["canonical_herb_target_pairs"])
    require_columns(
        pairs,
        ["canonical_name", "target_id"],
        paths["canonical_herb_target_pairs"].name,
    )
    pairs = pairs[["canonical_name", "target_id"]].copy()
    pairs["canonical_name"] = pairs["canonical_name"].map(norm_text)
    pairs["target_id"] = pairs["target_id"].map(norm_gene)
    pairs = pairs[
        pairs["canonical_name"].ne("")
        & pairs["target_id"].ne("")
    ].drop_duplicates()

    target_stats = pd.read_csv(paths["target_ubiquity"])
    require_columns(
        target_stats,
        ["target_id", "target_herb_degree", "N_database_herbs", "target_idf"],
        paths["target_ubiquity"].name,
    )
    target_stats = target_stats.copy()
    target_stats["target_id"] = target_stats["target_id"].map(norm_gene)
    target_stats["target_herb_degree"] = pd.to_numeric(
        target_stats["target_herb_degree"], errors="coerce"
    )
    target_stats["target_idf"] = pd.to_numeric(
        target_stats["target_idf"], errors="coerce"
    )
    if target_stats["target_id"].duplicated().any():
        raise RuntimeError("Duplicate targets in target ubiquity table.")

    herb_degree = pd.read_csv(paths["herb_target_degree"])
    require_columns(
        herb_degree,
        ["canonical_name", "herb_target_degree"],
        paths["herb_target_degree"].name,
    )
    herb_degree = herb_degree.copy()
    herb_degree["canonical_name"] = herb_degree["canonical_name"].map(norm_text)
    herb_degree["herb_target_degree"] = pd.to_numeric(
        herb_degree["herb_target_degree"], errors="coerce"
    )
    if herb_degree["canonical_name"].duplicated().any():
        raise RuntimeError("Duplicate herbs in herb target-degree table.")

    frozen = pd.read_csv(paths["frozen_LOPO_target_profile"])
    require_columns(
        frozen,
        [
            "syndrome_code",
            "target_id",
            "ubiquity_corrected",
            "ubiquity_corrected_sum1",
        ],
        paths["frozen_LOPO_target_profile"].name,
    )
    frozen = frozen.copy()
    frozen["syndrome_code"] = frozen["syndrome_code"].astype(str).str.strip()
    frozen["target_id"] = frozen["target_id"].map(norm_gene)
    frozen["ubiquity_corrected"] = pd.to_numeric(
        frozen["ubiquity_corrected"], errors="coerce"
    )
    frozen["ubiquity_corrected_sum1"] = pd.to_numeric(
        frozen["ubiquity_corrected_sum1"], errors="coerce"
    )

    return (
        w,
        mapping_unique,
        pairs,
        target_stats,
        herb_degree,
        frozen,
        paths,
        hashes,
    )


# =============================================================================
# 4. OBSERVED-ONLY HERB DECOMPOSITION
# =============================================================================

def build_observed_zero_filled_decomposition(
    w: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build the union of observed herbs and fill an unobserved syndrome-herb
    relation with ZERO before taking the across-syndrome minimum.
    """
    pivot = (
        w.pivot(
            index="canonical_name",
            columns="syndrome_code",
            values="primary_weight_W",
        )
        .reindex(columns=SYNDROME_CODES)
        .fillna(0.0)
        .sort_index()
    )

    out = pivot.reset_index().copy()
    out = out.rename(
        columns={s: f"W_original_{s}" for s in SYNDROME_CODES}
    )

    original_cols = [f"W_original_{s}" for s in SYNDROME_CODES]
    out["W_shared"] = out[original_cols].min(axis=1)

    for s in SYNDROME_CODES:
        out[f"W_specific_{s}"] = (
            out[f"W_original_{s}"] - out["W_shared"]
        )

    out["n_syndromes_observed"] = (
        out[original_cols].gt(0).sum(axis=1)
    )
    out["shared_positive"] = out["W_shared"].gt(0).astype(int)

    # Exact identity QA.
    max_err = 0.0
    for s in SYNDROME_CODES:
        err = np.max(
            np.abs(
                out[f"W_original_{s}"].to_numpy()
                - out["W_shared"].to_numpy()
                - out[f"W_specific_{s}"].to_numpy()
            )
        )
        max_err = max(max_err, float(err))

    if max_err > HERB_IDENTITY_TOL:
        raise RuntimeError(
            f"Herb-layer shared+specific identity failed: {max_err:.3e}"
        )

    # Long component table.
    parts = []

    shared = out[
        ["canonical_name", "W_shared"]
    ].rename(columns={"W_shared": "component_weight"})
    shared["component_type"] = "SHARED"
    shared["syndrome_code"] = "ALL"
    shared["component_id"] = "SHARED"
    parts.append(shared)

    for s in SYNDROME_CODES:
        x = out[
            ["canonical_name", f"W_specific_{s}"]
        ].rename(
            columns={f"W_specific_{s}": "component_weight"}
        )
        x["component_type"] = "SPECIFIC"
        x["syndrome_code"] = s
        x["component_id"] = f"SPECIFIC_{s}"
        parts.append(x)

    component_long = pd.concat(parts, ignore_index=True)
    component_long = component_long[
        component_long["component_weight"].gt(0)
    ].copy()

    return out, component_long


# =============================================================================
# 5. FROZEN TARGET PROJECTION
# =============================================================================

def attach_frozen_mapping(
    herb_weights: pd.DataFrame,
    mapping_unique: pd.DataFrame,
) -> pd.DataFrame:
    x = herb_weights.merge(
        mapping_unique,
        on="canonical_name",
        how="left",
        validate="many_to_one",
    )

    missing_mapping_record = x["mapped_to_target_db"].isna()
    if missing_mapping_record.any():
        bad = sorted(
            x.loc[missing_mapping_record, "canonical_name"].unique()
        )
        raise RuntimeError(
            "Herbs are absent from the frozen LOPO mapping dictionary:\n"
            f"{bad}"
        )

    return x


def project_components(
    component_long: pd.DataFrame,
    mapping_unique: pd.DataFrame,
    pairs: pd.DataFrame,
    target_stats: pd.DataFrame,
    herb_degree: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Project SHARED and three SPECIFIC components through the frozen v1.1
    degree-normalized + target-IDF transform.
    """
    mapped = attach_frozen_mapping(component_long, mapping_unique)

    # Component-level mapping QA before target expansion.
    qa_rows = []
    for component_id, g in mapped.groupby("component_id", sort=False):
        raw_mass = float(g["component_weight"].sum())
        mapped_mask = g["mapped_to_target_db"].eq(1)
        mapped_mass = float(
            g.loc[mapped_mask, "component_weight"].sum()
        )

        qa_rows.append({
            "component_id": component_id,
            "component_type": g["component_type"].iloc[0],
            "syndrome_code": g["syndrome_code"].iloc[0],
            "positive_herbs": int(len(g)),
            "mapped_positive_herbs": int(mapped_mask.sum()),
            "unmapped_positive_herbs": int((~mapped_mask).sum()),
            "raw_component_herb_mass": raw_mass,
            "mapped_component_herb_mass": mapped_mass,
            "mapped_weight_fraction":
                mapped_mass / raw_mass if raw_mass > 0 else np.nan,
        })

    component_qa = pd.DataFrame(qa_rows)

    mapped_pos = mapped[mapped["mapped_to_target_db"].eq(1)].copy()

    pairs_join = pairs.rename(
        columns={"canonical_name": "target_mapping_name"}
    )
    degree_join = herb_degree.rename(
        columns={"canonical_name": "target_mapping_name"}
    )

    x = mapped_pos.merge(
        pairs_join,
        on="target_mapping_name",
        how="inner",
        validate="many_to_many",
    )
    x = x.merge(
        degree_join,
        on="target_mapping_name",
        how="left",
        validate="many_to_one",
    )
    x = x.merge(
        target_stats[
            ["target_id", "target_herb_degree", "target_idf"]
        ],
        on="target_id",
        how="left",
        validate="many_to_one",
    )

    if x[
        ["herb_target_degree", "target_herb_degree", "target_idf"]
    ].isna().any().any():
        raise RuntimeError(
            "Missing frozen degree/IDF metadata after component projection join."
        )

    x["degree_norm_contribution"] = (
        x["component_weight"] / x["herb_target_degree"]
    )
    x["ubiquity_corrected_contribution"] = (
        x["degree_norm_contribution"] * x["target_idf"]
    )

    contribution_cols = [
        "component_id",
        "component_type",
        "syndrome_code",
        "canonical_name",
        "target_mapping_name",
        "mapping_method",
        "component_weight",
        "target_id",
        "herb_target_degree",
        "target_herb_degree",
        "target_idf",
        "degree_norm_contribution",
        "ubiquity_corrected_contribution",
    ]
    contributions = x[contribution_cols].copy()

    grouped = (
        x.groupby(
            ["component_id", "component_type", "syndrome_code", "target_id"],
            as_index=False,
        )
        .agg(
            contributing_herbs=("canonical_name", "nunique"),
            herb_degree_normalized=("degree_norm_contribution", "sum"),
            ubiquity_corrected=("ubiquity_corrected_contribution", "sum"),
        )
    )

    component_meta = (
        component_long[
            ["component_id", "component_type", "syndrome_code"]
        ]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    all_targets = sorted(target_stats["target_id"].unique())

    grid = (
        component_meta.assign(_key=1)
        .merge(
            pd.DataFrame({"target_id": all_targets, "_key": 1}),
            on="_key",
            how="inner",
        )
        .drop(columns="_key")
    )

    scores = grid.merge(
        grouped,
        on=["component_id", "component_type", "syndrome_code", "target_id"],
        how="left",
    )

    for c in [
        "contributing_herbs",
        "herb_degree_normalized",
        "ubiquity_corrected",
    ]:
        scores[c] = scores[c].fillna(0.0)

    scores = scores.merge(
        target_stats[
            ["target_id", "target_herb_degree", "N_database_herbs", "target_idf"]
        ],
        on="target_id",
        how="left",
        validate="many_to_one",
    )

    for score in ["herb_degree_normalized", "ubiquity_corrected"]:
        denom = scores.groupby("component_id")[score].transform("sum")
        scores[score + "_sum1"] = np.where(
            denom > 0,
            scores[score] / denom,
            0.0,
        )

    # Degree-normalization mass conservation.
    mass = (
        scores.groupby("component_id", as_index=False)
        .agg(
            projected_degree_norm_mass=("herb_degree_normalized", "sum"),
            projected_ubiquity_mass=("ubiquity_corrected", "sum"),
            normalized_target_sum=("ubiquity_corrected_sum1", "sum"),
            nonzero_targets=("ubiquity_corrected", lambda z: int((z > 0).sum())),
        )
    )
    component_qa = component_qa.merge(
        mass,
        on="component_id",
        how="left",
        validate="one_to_one",
    )
    component_qa["degree_norm_mass_error"] = (
        component_qa["projected_degree_norm_mass"]
        - component_qa["mapped_component_herb_mass"]
    ).abs()

    if component_qa["degree_norm_mass_error"].max() > MASS_CONSERVATION_TOL:
        raise RuntimeError(
            "Degree-normalization mass conservation failed."
        )

    if (
        component_qa["normalized_target_sum"].sub(1.0).abs().max()
        > NORMALIZED_SUM_TOL
    ):
        raise RuntimeError(
            "At least one shared/specific target profile does not sum to 1."
        )

    return scores, contributions, component_qa


# =============================================================================
# 6. RECONSTRUCT ORIGINAL FROZEN LOPO PROFILES AS A HARD CROSS-CHECK
# =============================================================================

def build_original_component_long(
    zero_filled: pd.DataFrame,
) -> pd.DataFrame:
    parts = []
    for s in SYNDROME_CODES:
        x = zero_filled[
            ["canonical_name", f"W_original_{s}"]
        ].rename(
            columns={f"W_original_{s}": "component_weight"}
        )
        x = x[x["component_weight"].gt(0)].copy()
        x["component_id"] = f"ORIGINAL_{s}"
        x["component_type"] = "ORIGINAL"
        x["syndrome_code"] = s
        parts.append(x)
    return pd.concat(parts, ignore_index=True)


def frozen_profile_crosscheck(
    original_scores: pd.DataFrame,
    frozen: pd.DataFrame,
    target_stats: pd.DataFrame,
) -> pd.DataFrame:
    all_targets = sorted(target_stats["target_id"].unique())

    frozen_grid = pd.MultiIndex.from_product(
        [SYNDROME_CODES, all_targets],
        names=["syndrome_code", "target_id"],
    ).to_frame(index=False)

    frozen_full = frozen_grid.merge(
        frozen[
            [
                "syndrome_code",
                "target_id",
                "ubiquity_corrected",
                "ubiquity_corrected_sum1",
            ]
        ],
        on=["syndrome_code", "target_id"],
        how="left",
    ).fillna({
        "ubiquity_corrected": 0.0,
        "ubiquity_corrected_sum1": 0.0,
    })

    reconstructed = original_scores[
        [
            "syndrome_code",
            "target_id",
            "ubiquity_corrected",
            "ubiquity_corrected_sum1",
        ]
    ].rename(
        columns={
            "ubiquity_corrected": "reconstructed_ubiquity_corrected",
            "ubiquity_corrected_sum1":
                "reconstructed_ubiquity_corrected_sum1",
        }
    )

    cmp = frozen_full.merge(
        reconstructed,
        on=["syndrome_code", "target_id"],
        how="left",
        validate="one_to_one",
    )

    cmp["abs_error_raw"] = (
        cmp["ubiquity_corrected"]
        - cmp["reconstructed_ubiquity_corrected"]
    ).abs()
    cmp["abs_error_sum1"] = (
        cmp["ubiquity_corrected_sum1"]
        - cmp["reconstructed_ubiquity_corrected_sum1"]
    ).abs()

    qa = (
        cmp.groupby("syndrome_code", as_index=False)
        .agg(
            max_abs_error_raw=("abs_error_raw", "max"),
            max_abs_error_sum1=("abs_error_sum1", "max"),
            frozen_positive_targets=(
                "ubiquity_corrected",
                lambda z: int((z > 0).sum()),
            ),
            reconstructed_positive_targets=(
                "reconstructed_ubiquity_corrected",
                lambda z: int((z > 0).sum()),
            ),
        )
    )

    if max(
        qa["max_abs_error_raw"].max(),
        qa["max_abs_error_sum1"].max(),
    ) > TARGET_RECONSTRUCTION_TOL:
        raise RuntimeError(
            "Reconstructed original LOPO target profiles do not match the "
            "frozen v1.1 target profiles."
        )

    return qa


# =============================================================================
# 7. LINEAR IDENTITY: ORIGINAL = SHARED + SPECIFIC BEFORE NORMALIZATION
# =============================================================================

def target_linear_identity_qa(
    original_scores: pd.DataFrame,
    component_scores: pd.DataFrame,
) -> pd.DataFrame:
    shared = (
        component_scores[
            component_scores["component_id"].eq("SHARED")
        ][["target_id", "herb_degree_normalized", "ubiquity_corrected"]]
        .rename(
            columns={
                "herb_degree_normalized": "shared_degree",
                "ubiquity_corrected": "shared_ubiquity",
            }
        )
    )

    rows = []
    for s in SYNDROME_CODES:
        orig = original_scores[
            original_scores["syndrome_code"].eq(s)
        ][
            ["target_id", "herb_degree_normalized", "ubiquity_corrected"]
        ].rename(
            columns={
                "herb_degree_normalized": "original_degree",
                "ubiquity_corrected": "original_ubiquity",
            }
        )

        spec = component_scores[
            component_scores["component_id"].eq(f"SPECIFIC_{s}")
        ][
            ["target_id", "herb_degree_normalized", "ubiquity_corrected"]
        ].rename(
            columns={
                "herb_degree_normalized": "specific_degree",
                "ubiquity_corrected": "specific_ubiquity",
            }
        )

        x = orig.merge(shared, on="target_id").merge(spec, on="target_id")
        degree_err = np.max(
            np.abs(
                x["original_degree"]
                - x["shared_degree"]
                - x["specific_degree"]
            )
        )
        ubiq_err = np.max(
            np.abs(
                x["original_ubiquity"]
                - x["shared_ubiquity"]
                - x["specific_ubiquity"]
            )
        )

        rows.append({
            "syndrome_code": s,
            "max_abs_degree_normalized_identity_error": float(degree_err),
            "max_abs_ubiquity_corrected_identity_error": float(ubiq_err),
        })

    qa = pd.DataFrame(rows)

    if (
        qa[
            [
                "max_abs_degree_normalized_identity_error",
                "max_abs_ubiquity_corrected_identity_error",
            ]
        ].to_numpy().max()
        > TARGET_LINEAR_IDENTITY_TOL
    ):
        raise RuntimeError(
            "Target-layer original = shared + specific linear identity failed."
        )

    return qa


# =============================================================================
# 8. TARGET-PROFILE DESCRIPTORS — NO CAP USED
# =============================================================================

def component_pairwise_similarity(
    scores: pd.DataFrame,
) -> pd.DataFrame:
    specific = scores[
        scores["component_type"].eq("SPECIFIC")
    ].copy()

    pivot = specific.pivot(
        index="target_id",
        columns="syndrome_code",
        values="ubiquity_corrected_sum1",
    ).fillna(0.0)

    pairs = [
        ("PDOL", "PHOL"),
        ("PDOL", "WHIL"),
        ("PHOL", "WHIL"),
    ]
    rows = []

    for a, b in pairs:
        xa = pivot[a].to_numpy()
        xb = pivot[b].to_numpy()

        row = {
            "component_type": "SPECIFIC",
            "syndrome_A": a,
            "syndrome_B": b,
            "spearman_all_targets":
                safe_corr(xa, xb, "spearman"),
            "pearson_all_targets":
                safe_corr(xa, xb, "pearson"),
        }

        for k in TOP_K_VALUES:
            top_a = set(
                pivot[a].sort_values(ascending=False).head(k).index
            )
            top_b = set(
                pivot[b].sort_values(ascending=False).head(k).index
            )
            union = top_a | top_b
            row[f"top{k}_jaccard"] = (
                len(top_a & top_b) / len(union)
                if union else np.nan
            )

        rows.append(row)

    return pd.DataFrame(rows)


def top_targets(scores: pd.DataFrame, top_n: int = 100) -> pd.DataFrame:
    x = scores[scores["ubiquity_corrected"].gt(0)].copy()
    x["rank_within_component"] = (
        x.groupby("component_id")["ubiquity_corrected_sum1"]
        .rank(method="first", ascending=False)
    )
    return (
        x[x["rank_within_component"].le(top_n)]
        .sort_values(["component_id", "rank_within_component"])
        .reset_index(drop=True)
    )


# =============================================================================
# 9. MAIN
# =============================================================================

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 84)
    print("LOPO OBSERVED-ONLY SHARED / SYNDROME-SPECIFIC PROJECTION v1.0")
    print("NO CAP INPUT IS READ")
    print("=" * 84)

    print("\n[1/8] Resolving and verifying frozen upstream inputs...")
    (
        w,
        mapping_unique,
        pairs,
        target_stats,
        herb_degree,
        frozen,
        input_paths,
        input_hashes,
    ) = load_inputs()

    print(f"  observed LOPO herb-syndrome rows: {len(w):,}")
    print(f"  unique observed LOPO herbs: {w['canonical_name'].nunique():,}")
    print(f"  frozen target universe: {len(target_stats):,}")

    print("\n[2/8] Building observed-only zero-filled herb decomposition...")
    zero_filled, component_long = build_observed_zero_filled_decomposition(w)

    shared_n = int((zero_filled["W_shared"] > 0).sum())
    shared_mass = float(zero_filled["W_shared"].sum())

    print(f"  union herbs: {len(zero_filled):,}")
    print(f"  herbs observed in all 3 syndromes / shared-positive: {shared_n:,}")
    print(f"  raw shared herb mass: {shared_mass:.6f}")

    for s in SYNDROME_CODES:
        total = float(zero_filled[f"W_original_{s}"].sum())
        spec = float(zero_filled[f"W_specific_{s}"].sum())
        print(
            f"  {s}: original mass={total:.6f}; "
            f"specific mass={spec:.6f}; "
            f"shared fraction={shared_mass / total:.3%}"
        )

    print("\n[3/8] Reconstructing original frozen LOPO target profiles...")
    original_long = build_original_component_long(zero_filled)
    original_scores, original_contrib, original_component_qa = project_components(
        original_long,
        mapping_unique,
        pairs,
        target_stats,
        herb_degree,
    )

    frozen_crosscheck = frozen_profile_crosscheck(
        original_scores,
        frozen,
        target_stats,
    )
    print(frozen_crosscheck.to_string(index=False))

    print("\n[4/8] Projecting SHARED and SPECIFIC components...")
    component_scores, component_contrib, component_qa = project_components(
        component_long,
        mapping_unique,
        pairs,
        target_stats,
        herb_degree,
    )
    print(
        component_qa[
            [
                "component_id",
                "positive_herbs",
                "mapped_positive_herbs",
                "raw_component_herb_mass",
                "mapped_component_herb_mass",
                "mapped_weight_fraction",
                "nonzero_targets",
            ]
        ].to_string(index=False)
    )

    print("\n[5/8] Verifying target-layer linear identity...")
    linear_qa = target_linear_identity_qa(
        original_scores,
        component_scores,
    )
    print(linear_qa.to_string(index=False))

    print("\n[6/8] Computing CAP-independent component descriptors...")
    pairwise = component_pairwise_similarity(component_scores)
    top100 = top_targets(component_scores, top_n=100)
    print(pairwise.to_string(index=False))

    print("\n[7/8] Writing outputs...")

    zero_filled.to_csv(
        OUTPUT_DIR / "LOPO_observed_zero_filled_herb_matrix_v1_0.csv",
        index=False,
        encoding="utf-8-sig",
    )

    component_long.to_csv(
        OUTPUT_DIR / "LOPO_shared_specific_herb_components_v1_0.csv",
        index=False,
        encoding="utf-8-sig",
    )

    component_scores.to_csv(
        OUTPUT_DIR / "LOPO_shared_specific_target_projection_all_targets_v1_0.csv",
        index=False,
        encoding="utf-8-sig",
    )

    # Downstream-ready primary target profiles, including zero targets so that
    # the frozen 374-target universe remains explicit.
    component_scores[
        [
            "component_id",
            "component_type",
            "syndrome_code",
            "target_id",
            "ubiquity_corrected",
            "ubiquity_corrected_sum1",
            "contributing_herbs",
            "target_herb_degree",
            "target_idf",
        ]
    ].to_csv(
        OUTPUT_DIR / "LOPO_SHARED_SPECIFIC_primary_target_profiles_v1_0.csv",
        index=False,
        encoding="utf-8-sig",
    )

    component_contrib.to_csv(
        OUTPUT_DIR / "LOPO_shared_specific_herb_target_contributions_v1_0.csv",
        index=False,
        encoding="utf-8-sig",
    )

    component_qa.to_csv(
        OUTPUT_DIR / "LOPO_shared_specific_component_mapping_mass_QA_v1_0.csv",
        index=False,
        encoding="utf-8-sig",
    )

    frozen_crosscheck.to_csv(
        OUTPUT_DIR / "LOPO_original_frozen_target_reconstruction_QA_v1_0.csv",
        index=False,
        encoding="utf-8-sig",
    )

    linear_qa.to_csv(
        OUTPUT_DIR / "LOPO_shared_specific_target_linear_identity_QA_v1_0.csv",
        index=False,
        encoding="utf-8-sig",
    )

    pairwise.to_csv(
        OUTPUT_DIR / "LOPO_specific_target_pairwise_similarity_v1_0.csv",
        index=False,
        encoding="utf-8-sig",
    )

    top100.to_csv(
        OUTPUT_DIR / "LOPO_shared_specific_top100_targets_v1_0.csv",
        index=False,
        encoding="utf-8-sig",
    )

    print("\n[8/8] Writing metadata and freeze-candidate statement...")

    herb_summary_rows = []
    for s in SYNDROME_CODES:
        total = float(zero_filled[f"W_original_{s}"].sum())
        specific = float(zero_filled[f"W_specific_{s}"].sum())
        herb_summary_rows.append({
            "syndrome_code": s,
            "original_positive_herbs":
                int((zero_filled[f"W_original_{s}"] > 0).sum()),
            "original_raw_herb_mass": total,
            "shared_positive_herbs": shared_n,
            "shared_raw_herb_mass": shared_mass,
            "specific_positive_herbs":
                int((zero_filled[f"W_specific_{s}"] > 0).sum()),
            "specific_raw_herb_mass": specific,
            "shared_fraction_of_original_raw_herb_mass":
                shared_mass / total,
            "specific_fraction_of_original_raw_herb_mass":
                specific / total,
        })

    herb_summary = pd.DataFrame(herb_summary_rows)
    herb_summary.to_csv(
        OUTPUT_DIR / "LOPO_shared_specific_herb_mass_summary_v1_0.csv",
        index=False,
        encoding="utf-8-sig",
    )

    metadata = {
        "analysis_name":
            "LOPO observed-only shared/syndrome-specific projection v1.0",
        "created_utc":
            datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "status":
            "UPSTREAM_FREEZE_CANDIDATE_BEFORE_SHARED_SPECIFIC_CAP_ALIGNMENT",
        "CAP_transcriptomics_read_by_this_script": False,
        "CAP_alignment_results_read_by_this_script": False,
        "decomposition_rule": {
            "unobserved_syndrome_herb_relation": 0,
            "shared":
                "W_shared(h) = min_s W_hs after zero-filling unobserved relations",
            "specific":
                "W_specific(h,s) = W_hs - W_shared(h)",
            "identity":
                "W_hs = W_shared(h) + W_specific(h,s)",
            "important":
                "No Jeffreys/prior floor is assigned to unobserved syndrome-herb relations.",
        },
        "projection_rule": {
            "mapping":
                "frozen v1.1 target mapping dictionary; no remapping or fuzzy matching",
            "degree_normalization":
                "component_weight / frozen herb_target_degree",
            "target_ubiquity_correction":
                "multiply by frozen target_idf",
            "target_idf":
                "log((N_H+1)/(n_t+1))+1 from frozen v1.1 target universe",
            "normalized_profile":
                "ubiquity_corrected_sum1",
            "raw_mass_preserved_for_diagnostics": True,
        },
        "shared_positive_herbs": shared_n,
        "shared_raw_herb_mass": shared_mass,
        "target_universe_n": int(len(target_stats)),
        "input_paths": {k: str(v) for k, v in input_paths.items()},
        "input_sha256": input_hashes,
        "hard_QA": {
            "max_original_frozen_profile_reconstruction_error":
                float(
                    frozen_crosscheck[
                        ["max_abs_error_raw", "max_abs_error_sum1"]
                    ].to_numpy().max()
                ),
            "max_target_linear_identity_error":
                float(
                    linear_qa[
                        [
                            "max_abs_degree_normalized_identity_error",
                            "max_abs_ubiquity_corrected_identity_error",
                        ]
                    ].to_numpy().max()
                ),
            "max_degree_normalization_mass_error":
                float(component_qa["degree_norm_mass_error"].max()),
        },
        "next_step_policy":
            "Audit and freeze these profiles before using CAP. Shared/specific CAP alignment must not modify this decomposition.",
    }

    with open(
        OUTPUT_DIR / "LOPO_shared_specific_projection_run_metadata_v1_0.json",
        "w",
        encoding="utf-8",
    ) as f:
        json.dump(metadata, f, ensure_ascii=False, indent=2)

    freeze_statement = f"""LOPO OBSERVED-ONLY SHARED / SYNDROME-SPECIFIC PROJECTION v1.0
FREEZE CANDIDATE — AUDIT BEFORE CAP ALIGNMENT

STATUS
------
UPSTREAM FREEZE CANDIDATE.
This script reads NO CAP transcriptomic data and NO syndrome-disease alignment
result.

HERB-LAYER RULE
---------------
The union of observed LOPO herbs is formed first.
An unobserved syndrome-herb relation is assigned exactly ZERO.

For every herb h:
    W_shared(h) = min_s W_hs
    W_specific(h,s) = W_hs - W_shared(h)

Therefore:
    W_hs = W_shared(h) + W_specific(h,s)

No unobserved treatment is inferred.
No Jeffreys/prior floor is assigned to an unobserved syndrome-herb relation.

CURRENT DECOMPOSITION
---------------------
Union herbs: {len(zero_filled)}
Shared-positive herbs: {shared_n}
Raw shared herb mass: {shared_mass:.12f}

{herb_summary.to_string(index=False)}

TARGET PROJECTION
-----------------
Frozen v1.1 mapping and target database are reused without modification.

For each mapped herb:
    degree contribution = component_weight / |T_h|

For each target:
    IDF_t = frozen v1.1 target IDF
    ubiquity-corrected score =
        sum_h [component_weight / |T_h| * IDF_t]

Each nonzero SHARED/SPECIFIC target profile is also normalized to sum 1 for
profile comparison and later alignment. Raw component masses and raw projected
scores remain saved for magnitude diagnostics.

HARD QA
-------
Original frozen LOPO target-profile reconstruction max error:
{metadata["hard_QA"]["max_original_frozen_profile_reconstruction_error"]:.3e}

Target-layer original = shared + specific max linear identity error:
{metadata["hard_QA"]["max_target_linear_identity_error"]:.3e}

Degree-normalization mass-conservation max error:
{metadata["hard_QA"]["max_degree_normalization_mass_error"]:.3e}

FREEZE POLICY
-------------
1. Do not alter zero-filling, min-shared decomposition, frozen herb-target
   mapping, herb target degree, target IDF, or sum1 normalization after CAP is
   inspected.
2. The shared component is common to all syndromes by construction.
3. Specific components are residuals after removing the shared herb layer.
4. Magnitude comparisons must retain raw component mass; normalized profiles
   describe composition and must not be interpreted as equal raw magnitude.
5. Do not construct shared/specific components by taking a minimum directly
   in target space; decomposition is defined at the herb layer.
6. After audit, freeze this upstream layer before shared/specific CAP alignment.
"""

    (
        OUTPUT_DIR
        / "LOPO_SHARED_SPECIFIC_PROJECTION_FREEZE_CANDIDATE_v1_0.txt"
    ).write_text(freeze_statement, encoding="utf-8")

    print("\n" + "=" * 84)
    print("Shared/specific projection completed.")
    print("STATUS: UPSTREAM FREEZE CANDIDATE — audit before CAP alignment.")
    print("=" * 84)
    print(f"\nOutputs:\n{OUTPUT_DIR}")


if __name__ == "__main__":
    main()
