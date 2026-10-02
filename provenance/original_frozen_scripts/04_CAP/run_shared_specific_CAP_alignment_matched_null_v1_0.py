# -*- coding: utf-8 -*-
"""
run_shared_specific_CAP_alignment_matched_null_v1_0.py

Purpose
-------
Evaluate the already FROZEN observed-only SHARED and syndrome-SPECIFIC target
profiles against the already FROZEN CAP continuous meta-Z program using the
same matched-null principles as the frozen complete-profile analysis v1.0.1.

This script is DOWNSTREAM ONLY. It must not modify the upstream decomposition,
target mapping, target IDF, CAP meta-Z, STRING network, or null strata.

Frozen components
-----------------
The input target profiles contain four components over the same 374-target
TCM-ID universe:

    SHARED
    SPECIFIC_PDOL
    SPECIFIC_PHOL
    SPECIFIC_WHIL

They were constructed upstream, without CAP, from:

    W_shared(h) = min_s W_hs
    W_specific(h,s) = W_hs - W_shared(h)

with unobserved syndrome-herb relations set to zero.

Primary component alignment
---------------------------
Use the UNNORMALIZED frozen ubiquity-corrected target score U_cg. Scaling does
not affect the weighted-mean alignment, but retaining raw U preserves component
magnitude and allows exact decomposition of each original syndrome profile.

For component c:

    A_raw(c) = sum_g U_cg I_CAP(g) Z_CAP,g
               --------------------------------
               sum_g U_cg I_CAP(g)

The denominator and its fraction of the component's total projected target
mass are reported explicitly.

Exact original-profile decomposition
------------------------------------
Because the target projection is linear BEFORE sum1 normalization:

    U_original,s = U_shared + U_specific,s

therefore, on CAP-measured genes:

    N_original,s = N_shared + N_specific,s
    D_original,s = D_shared + D_specific,s
    A_original,s = (N_shared + N_specific,s) / (D_shared + D_specific,s)

The script reports the shared and specific contributions to the original
alignment while preserving their measured raw target masses.

Matched null
------------
The target universe and 19 matched strata are NOT rebuilt. They are read
verbatim from the formally frozen complete-profile alignment v1.0.1:

    matched_null_target_universe_and_strata_v1_0_1.csv

Within every matched stratum, the entire RAW component target-weight vectors,
INCLUDING zeros, are permuted.

CRITICAL: the SAME donor permutation is applied simultaneously to SHARED and
all three SPECIFIC components in every replicate. Consequently:

    perm(U_original,s)
      = perm(U_shared) + perm(U_specific,s)

exactly. This preserves the upstream linear decomposition under the null and
allows the reconstructed complete-profile null to reproduce the already frozen
complete-profile null distribution.

Each component null preserves exactly:
- total raw projected target mass;
- global positive/zero target count;
- positive/zero count within each frozen matched stratum;
- raw target-weight mass within each frozen matched stratum;
- the 374-target database universe.

The frozen strata already match:
- value-based STRING unweighted-degree quartiles, with ties kept together;
- value-based target-herb-ubiquity quartiles over the complete target universe;
- NO_STRING as a separate degree category while retaining ubiquity strata.

Primary inference
-----------------
For SHARED and each SPECIFIC component, report:

    Z_align = (A_obs - mean(A_null)) / sd(A_null)

plus empirical upper-tail, lower-tail, and two-sided P values with +1
correction.

Paired contrasts
----------------
Because the same stratum permutation is shared across all components, the
script also evaluates paired null distributions for:

1. each SPECIFIC component minus SHARED;
2. pairwise differences among the three SPECIFIC components.

These contrasts describe differential disease-context localization. They do
NOT imply discrete molecular identities, therapeutic reversal, or efficacy.

Secondary STRING RWR
--------------------
The frozen STRING v12 high-confidence network and restart r=0.5 are reused only
as a SECONDARY network-localization layer. Component and null RWR alignments are
computed exactly via the same adjoint stationary-RWR identity used in the
frozen complete-profile analysis.

Hard cross-checks
-----------------
The script requires all of the following before accepting output:
1. shared + specific raw target vectors reconstruct each frozen original LOPO
   target vector;
2. reconstructed original observed RAW alignment reproduces the frozen
   complete-profile observed RAW result;
3. reconstructed original observed RWR alignment reproduces the frozen
   complete-profile observed RWR result;
4. reconstructed complete-profile RAW and RWR null distributions reproduce the
   already frozen complete-profile v1.0.1 null distribution.

Interpretation
--------------
Positive CAP meta-Z means higher expression in CAP than healthy control.
Negative CAP meta-Z means lower expression in CAP than healthy control.
This is disease-context localization, NOT therapeutic reversal and NOT
evidence of treatment efficacy.
"""

from __future__ import annotations

from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import re

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix, csr_matrix, diags


# =============================================================================
# 1. USER SETTINGS — EDIT PATHS HERE IF NEEDED
# =============================================================================

PROJECT_ROOT = Path(r"E:\00project\Tanre Yongfei\重构")

SHARED_SPECIFIC_TARGET_FILE = (
    PROJECT_ROOT
    / "09_LOPO_shared_specific_projection_v1_0"
    / "LOPO_SHARED_SPECIFIC_primary_target_profiles_v1_0.csv"
)

COMPONENT_MASS_QA_FILE = (
    PROJECT_ROOT
    / "09_LOPO_shared_specific_projection_v1_0"
    / "LOPO_shared_specific_component_mapping_mass_QA_v1_0.csv"
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

STRING_CLEAN_EDGE_FILE = (
    PROJECT_ROOT
    / "06_STRING_RWR_frozen_v1_0"
    / "STRING_clean_high_confidence_edges_v1_0.csv"
)

COMPLETE_FROZEN_TARGET_PROFILE_FILE = (
    PROJECT_ROOT
    / "05_weighted_herb_target_projection_v1_1_FROZEN"
    / "LOPO_FROZEN_primary_target_weights_for_RWR_v1_1.csv"
)

COMPLETE_OBSERVED_RAW_FILE = (
    PROJECT_ROOT
    / "08_syndrome_disease_alignment_matched_null_v1_0_1"
    / "CAP_alignment_observed_RAW_primary_v1_0_1.csv"
)

COMPLETE_OBSERVED_RWR_FILE = (
    PROJECT_ROOT
    / "08_syndrome_disease_alignment_matched_null_v1_0_1"
    / "CAP_alignment_observed_RWR_secondary_v1_0_1.csv"
)

COMPLETE_NULL_DISTRIBUTION_FILE = (
    PROJECT_ROOT
    / "08_syndrome_disease_alignment_matched_null_v1_0_1"
    / "CAP_alignment_full_matched_null_distribution_v1_0_1.csv.gz"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "10_shared_specific_CAP_alignment_matched_null_v1_0"
)

ENFORCE_EXPECTED_SHA256 = True

EXPECTED_SHA256 = {
    "shared_specific_target_profiles":
        "af66b3f49961b7a0e578e35b6c9bb86b7c84007f4c68d7efa4c286e6b9b85669",
    "component_mass_QA":
        "c1307f444c4f35eeeff07eee280a73f32287ee4d08fefbf440bb6f60d45d71b7",
    "frozen_target_universe_strata":
        "a3871d4da3b3f4005833875f4e8d7c11e529214f6f93a95255f5acc5fffef297",
    "cap_meta_z":
        "168d9a23c06df388b69f73d92b2e7811b1f5332dfb507893511367fd752be58e",
    "string_clean_edges":
        "deb2d74bdf0bcb7f249883c39d18bfecab9c9e9f5567612e0fc1fdc6db5ac6ba",
    "complete_frozen_target_profile":
        "9256cd3fe7ae142a91c7e081d97ba8a9225aa5aea3434c5f5673d53e5a10c9ef",
    "complete_observed_raw":
        "4e7689a56fdff4512f1bcdde1d884bbef90a4aa53bada6b0eafa6fe085d3df6e",
    "complete_observed_rwr":
        "add68d920ff8bb9d7e80b289546ec2137f3d873b2abb5afed71c20117cbfc24e",
    "complete_null_distribution":
        "7062555f8443581a46c3d328d77afa8d1a08095db411ea60e4e68f6cdf744f09",
}

COMPONENT_IDS = [
    "SHARED",
    "SPECIFIC_PDOL",
    "SPECIFIC_PHOL",
    "SPECIFIC_WHIL",
]

SPECIFIC_COMPONENT_BY_SYNDROME = {
    "PDOL": "SPECIFIC_PDOL",
    "PHOL": "SPECIFIC_PHOL",
    "WHIL": "SPECIFIC_WHIL",
}
SYNDROME_CODES = ["PDOL", "PHOL", "WHIL"]

RAW_TARGET_WEIGHT_COLUMN = "ubiquity_corrected"
NORMALIZED_TARGET_WEIGHT_COLUMN = "ubiquity_corrected_sum1"
CAP_PRIMARY_Z_COLUMN = "Z_meta_equal"
CAP_SENSITIVITY_Z_COLUMN = "Z_meta_neff_weighted"
PRIMARY_RWR_RESTART = 0.5

N_NULL = 10_000
RANDOM_SEED = 20260921

ADJOINT_TOL_MAXABS = 1e-12
ADJOINT_MAX_ITER = 10_000

WEIGHT_SUM_TOL = 1e-10
RECONSTRUCTION_TOL = 1e-12
OBSERVED_CROSSCHECK_TOL = 1e-10
NULL_CROSSCHECK_TOL = 1e-10
CAP_CROSSCHECK_TOL = 1e-12

SAVE_FULL_NULL_DISTRIBUTION = True


# =============================================================================
# 2. HELPERS
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
            f"{label}: could not find '{preferred.name}' under:\n{PROJECT_ROOT}"
        )

    if len(matches) == 1:
        print(f"  Auto-resolved {label}:\n    {matches[0]}")
        return matches[0]

    hashes = {str(p): sha256_file(p) for p in matches}
    if len(set(hashes.values())) == 1:
        chosen = matches[0]
        print(f"  Multiple identical copies found for {label}; using:\n    {chosen}")
        return chosen

    expected = EXPECTED_SHA256.get(label)
    if expected is not None:
        matching = [Path(p) for p, h in hashes.items() if h == expected]
        if len(matching) == 1:
            print(
                f"  Multiple non-identical copies found for {label}; "
                f"selected the unique frozen-hash match:\n    {matching[0]}"
            )
            return matching[0]

    details = "\n".join(f"  {p}\n    SHA256={h}" for p, h in hashes.items())
    raise RuntimeError(
        f"{label}: multiple non-identical files named '{preferred.name}' were found. "
        f"Resolve manually.\n{details}"
    )


def verify_hash(path: Path, key: str) -> str:
    observed = sha256_file(path)
    expected = EXPECTED_SHA256[key]
    if ENFORCE_EXPECTED_SHA256 and observed != expected:
        raise RuntimeError(
            f"Frozen input SHA256 mismatch for {key}.\n"
            f"File: {path}\nExpected: {expected}\nObserved: {observed}\n"
            "Do not mix analysis branches. If the upstream frozen input was "
            "intentionally changed, create a new downstream version."
        )
    return observed


def empirical_p_values(null_values: np.ndarray, observed: float) -> dict:
    x = np.asarray(null_values, dtype=float)
    if not np.isfinite(x).all():
        raise RuntimeError("Null distribution contains non-finite values.")
    n = len(x)
    mu = float(np.mean(x))
    sd = float(np.std(x, ddof=1))
    z_align = (float(observed) - mu) / sd if sd > 0 else np.nan
    p_upper = (1 + int(np.sum(x >= observed))) / (n + 1)
    p_lower = (1 + int(np.sum(x <= observed))) / (n + 1)
    obs_dev = abs(float(observed) - mu)
    p_two = (1 + int(np.sum(np.abs(x - mu) >= obs_dev))) / (n + 1)
    return {
        "null_mean": mu,
        "null_sd": sd,
        "Z_align": z_align,
        "empirical_p_upper": p_upper,
        "empirical_p_lower": p_lower,
        "empirical_p_two_sided": p_two,
        "null_q025": float(np.quantile(x, 0.025)),
        "null_q25": float(np.quantile(x, 0.25)),
        "null_median": float(np.quantile(x, 0.50)),
        "null_q75": float(np.quantile(x, 0.75)),
        "null_q975": float(np.quantile(x, 0.975)),
    }


# =============================================================================
# 3. LOAD FROZEN INPUTS
# =============================================================================

def load_inputs():
    preferred = {
        "shared_specific_target_profiles": SHARED_SPECIFIC_TARGET_FILE,
        "component_mass_QA": COMPONENT_MASS_QA_FILE,
        "frozen_target_universe_strata": FROZEN_TARGET_UNIVERSE_STRATA_FILE,
        "cap_meta_z": CAP_META_Z_FILE,
        "string_clean_edges": STRING_CLEAN_EDGE_FILE,
        "complete_frozen_target_profile": COMPLETE_FROZEN_TARGET_PROFILE_FILE,
        "complete_observed_raw": COMPLETE_OBSERVED_RAW_FILE,
        "complete_observed_rwr": COMPLETE_OBSERVED_RWR_FILE,
        "complete_null_distribution": COMPLETE_NULL_DISTRIBUTION_FILE,
    }

    paths = {k: resolve_input(v, k) for k, v in preferred.items()}
    hashes = {k: verify_hash(v, k) for k, v in paths.items()}

    profiles = pd.read_csv(paths["shared_specific_target_profiles"])
    require_columns(
        profiles,
        [
            "component_id", "component_type", "syndrome_code", "target_id",
            RAW_TARGET_WEIGHT_COLUMN, NORMALIZED_TARGET_WEIGHT_COLUMN,
        ],
        paths["shared_specific_target_profiles"].name,
    )
    profiles = profiles.copy()
    profiles["component_id"] = profiles["component_id"].astype(str).str.strip()
    profiles["target_id"] = profiles["target_id"].map(norm_gene)
    profiles[RAW_TARGET_WEIGHT_COLUMN] = pd.to_numeric(
        profiles[RAW_TARGET_WEIGHT_COLUMN], errors="coerce"
    )
    profiles[NORMALIZED_TARGET_WEIGHT_COLUMN] = pd.to_numeric(
        profiles[NORMALIZED_TARGET_WEIGHT_COLUMN], errors="coerce"
    )
    profiles = profiles[profiles["component_id"].isin(COMPONENT_IDS)].copy()

    mass_qa = pd.read_csv(paths["component_mass_QA"])
    require_columns(
        mass_qa,
        [
            "component_id", "raw_component_herb_mass", "mapped_component_herb_mass",
            "mapped_weight_fraction", "projected_ubiquity_mass", "nonzero_targets",
        ],
        paths["component_mass_QA"].name,
    )

    universe = pd.read_csv(paths["frozen_target_universe_strata"])
    require_columns(
        universe,
        [
            "gene", "target_herb_degree", "target_idf",
            "network_degree_unweighted", "network_degree_weighted", "in_STRING",
            CAP_PRIMARY_Z_COLUMN, "measured_in_CAP_primary", "matched_stratum",
        ],
        paths["frozen_target_universe_strata"].name,
    )
    universe = universe.copy()
    universe["gene"] = universe["gene"].map(norm_gene)
    universe = universe.sort_values("gene").reset_index(drop=True)
    if universe["gene"].duplicated().any():
        raise RuntimeError("Duplicate genes in frozen target-universe strata file.")

    cap = pd.read_csv(paths["cap_meta_z"])
    require_columns(cap, ["gene", CAP_PRIMARY_Z_COLUMN], paths["cap_meta_z"].name)
    cap = cap.copy()
    cap["gene"] = cap["gene"].map(norm_gene)
    cap[CAP_PRIMARY_Z_COLUMN] = pd.to_numeric(cap[CAP_PRIMARY_Z_COLUMN], errors="coerce")
    if CAP_SENSITIVITY_Z_COLUMN in cap.columns:
        cap[CAP_SENSITIVITY_Z_COLUMN] = pd.to_numeric(
            cap[CAP_SENSITIVITY_Z_COLUMN], errors="coerce"
        )
    cap = cap[cap["gene"].ne("") & cap[CAP_PRIMARY_Z_COLUMN].notna()].copy()
    if cap["gene"].duplicated().any():
        raise RuntimeError("Duplicate genes in frozen CAP meta-Z file.")

    edges = pd.read_csv(paths["string_clean_edges"])
    require_columns(edges, ["u", "v", "edge_weight"], paths["string_clean_edges"].name)
    edges = edges.copy()
    edges["u"] = edges["u"].map(norm_gene)
    edges["v"] = edges["v"].map(norm_gene)
    edges["edge_weight"] = pd.to_numeric(edges["edge_weight"], errors="coerce")
    edges = edges[
        edges["u"].ne("") & edges["v"].ne("")
        & edges["edge_weight"].notna() & edges["edge_weight"].gt(0)
        & edges["u"].ne(edges["v"])
    ].copy()

    complete = pd.read_csv(paths["complete_frozen_target_profile"])
    require_columns(
        complete,
        ["syndrome_code", "target_id", "ubiquity_corrected", "ubiquity_corrected_sum1"],
        paths["complete_frozen_target_profile"].name,
    )
    complete = complete.copy()
    complete["syndrome_code"] = complete["syndrome_code"].astype(str).str.strip()
    complete["target_id"] = complete["target_id"].map(norm_gene)

    complete_raw_obs = pd.read_csv(paths["complete_observed_raw"])
    require_columns(
        complete_raw_obs,
        ["syndrome_code", "observed_alignment"],
        paths["complete_observed_raw"].name,
    )

    complete_rwr_obs = pd.read_csv(paths["complete_observed_rwr"])
    require_columns(
        complete_rwr_obs,
        ["syndrome_code", "observed_alignment"],
        paths["complete_observed_rwr"].name,
    )

    complete_null = pd.read_csv(paths["complete_null_distribution"])
    require_columns(
        complete_null,
        ["null_iteration", "syndrome_code", "raw_alignment_null", "rwr_alignment_null"],
        paths["complete_null_distribution"].name,
    )

    return (
        profiles, mass_qa, universe, cap, edges, complete,
        complete_raw_obs, complete_rwr_obs, complete_null, paths, hashes,
    )


# =============================================================================
# 4. HARD INPUT CROSS-CHECKS + COMPONENT MATRICES
# =============================================================================

def validate_cap_against_frozen_universe(universe: pd.DataFrame, cap: pd.DataFrame) -> float:
    c = cap[["gene", CAP_PRIMARY_Z_COLUMN]].rename(
        columns={CAP_PRIMARY_Z_COLUMN: "cap_file_z"}
    )
    x = universe[["gene", CAP_PRIMARY_Z_COLUMN]].merge(c, on="gene", how="left")
    measured = x[CAP_PRIMARY_Z_COLUMN].notna()
    if x.loc[measured, "cap_file_z"].isna().any():
        raise RuntimeError("A CAP-measured target in the frozen universe is absent from CAP meta-Z.")
    err = float(np.max(np.abs(
        x.loc[measured, CAP_PRIMARY_Z_COLUMN].to_numpy(dtype=float)
        - x.loc[measured, "cap_file_z"].to_numpy(dtype=float)
    ))) if measured.any() else 0.0
    if err > CAP_CROSSCHECK_TOL:
        raise RuntimeError(f"CAP meta-Z mismatch versus frozen target universe: {err:.3e}")
    return err


def build_component_matrices(
    profiles: pd.DataFrame,
    mass_qa: pd.DataFrame,
    universe: pd.DataFrame,
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    genes = universe["gene"].tolist()
    gene_set = set(genes)
    W_raw = np.zeros((len(COMPONENT_IDS), len(genes)), dtype=float)
    W_sum1 = np.zeros_like(W_raw)
    qa_rows = []

    for i, cid in enumerate(COMPONENT_IDS):
        x = profiles[profiles["component_id"].eq(cid)].copy()
        if len(x) != len(universe):
            raise RuntimeError(f"{cid}: expected {len(universe)} target rows, found {len(x)}.")
        if x["target_id"].duplicated().any():
            raise RuntimeError(f"{cid}: duplicate target IDs.")
        if set(x["target_id"]) != gene_set:
            raise RuntimeError(f"{cid}: target universe differs from frozen 374-target universe.")

        raw_map = dict(zip(x["target_id"], x[RAW_TARGET_WEIGHT_COLUMN]))
        sum1_map = dict(zip(x["target_id"], x[NORMALIZED_TARGET_WEIGHT_COLUMN]))
        W_raw[i, :] = np.array([float(raw_map[g]) for g in genes])
        W_sum1[i, :] = np.array([float(sum1_map[g]) for g in genes])

        if np.any(W_raw[i] < -WEIGHT_SUM_TOL) or np.any(W_sum1[i] < -WEIGHT_SUM_TOL):
            raise RuntimeError(f"{cid}: negative target weights detected.")
        if abs(float(W_sum1[i].sum()) - 1.0) > WEIGHT_SUM_TOL:
            raise RuntimeError(f"{cid}: normalized target profile does not sum to 1.")

        q = mass_qa[mass_qa["component_id"].eq(cid)]
        if len(q) != 1:
            raise RuntimeError(f"{cid}: expected exactly one component mass QA row.")
        q = q.iloc[0]
        raw_mass_error = abs(float(W_raw[i].sum()) - float(q["projected_ubiquity_mass"]))
        pos = int(np.sum(W_raw[i] > 0))
        if raw_mass_error > RECONSTRUCTION_TOL:
            raise RuntimeError(f"{cid}: projected target mass mismatch: {raw_mass_error:.3e}")
        if pos != int(q["nonzero_targets"]):
            raise RuntimeError(f"{cid}: nonzero target count mismatch.")

        qa_rows.append({
            "component_id": cid,
            "raw_projected_target_mass": float(W_raw[i].sum()),
            "normalized_target_mass": float(W_sum1[i].sum()),
            "nonzero_targets": pos,
            "raw_component_herb_mass": float(q["raw_component_herb_mass"]),
            "mapped_component_herb_mass": float(q["mapped_component_herb_mass"]),
            "mapped_weight_fraction": float(q["mapped_weight_fraction"]),
            "projected_target_mass_error_vs_upstream_QA": raw_mass_error,
        })

    return W_raw, W_sum1, pd.DataFrame(qa_rows)


def build_complete_raw_matrix(
    complete: pd.DataFrame,
    universe: pd.DataFrame,
) -> np.ndarray:
    genes = universe["gene"].tolist()
    W = np.zeros((len(SYNDROME_CODES), len(genes)), dtype=float)
    for i, s in enumerate(SYNDROME_CODES):
        x = complete[complete["syndrome_code"].eq(s)].copy()
        d = dict(zip(x["target_id"], x["ubiquity_corrected"]))
        outside = set(d) - set(genes)
        if outside:
            raise RuntimeError(f"{s}: complete frozen targets outside frozen universe: {sorted(outside)[:10]}")
        W[i, :] = np.array([float(d.get(g, 0.0)) for g in genes])
    return W


def component_reconstruction_qa(
    W_raw: np.ndarray,
    W_complete: np.ndarray,
) -> pd.DataFrame:
    shared = W_raw[COMPONENT_IDS.index("SHARED")]
    rows = []
    for i, s in enumerate(SYNDROME_CODES):
        spec = W_raw[COMPONENT_IDS.index(SPECIFIC_COMPONENT_BY_SYNDROME[s])]
        recon = shared + spec
        err = float(np.max(np.abs(recon - W_complete[i])))
        mass_err = abs(float(recon.sum()) - float(W_complete[i].sum()))
        rows.append({
            "syndrome_code": s,
            "max_abs_target_weight_reconstruction_error": err,
            "raw_target_mass_reconstruction_error": mass_err,
            "complete_raw_target_mass": float(W_complete[i].sum()),
            "shared_raw_target_mass": float(shared.sum()),
            "specific_raw_target_mass": float(spec.sum()),
        })
    qa = pd.DataFrame(rows)
    if qa[[
        "max_abs_target_weight_reconstruction_error",
        "raw_target_mass_reconstruction_error",
    ]].to_numpy().max() > RECONSTRUCTION_TOL:
        raise RuntimeError("Shared + specific does not reconstruct frozen complete target profiles.")
    return qa


# =============================================================================
# 5. OBSERVED RAW COMPONENT ALIGNMENT + EXACT MIXTURE DECOMPOSITION
# =============================================================================

def observed_raw_components(
    universe: pd.DataFrame,
    W_raw: np.ndarray,
    mass_qa: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    measured = universe["measured_in_CAP_primary"].astype(float).to_numpy()
    z = universe[CAP_PRIMARY_Z_COLUMN].fillna(0.0).to_numpy(dtype=float)
    z_sens = (
        universe[CAP_SENSITIVITY_Z_COLUMN].fillna(0.0).to_numpy(dtype=float)
        if CAP_SENSITIVITY_Z_COLUMN in universe.columns else None
    )

    rows = []
    contributions = []
    raw_cache = {}

    for i, cid in enumerate(COMPONENT_IDS):
        w = W_raw[i]
        total_mass = float(w.sum())
        den = float(w @ measured)
        num = float(w @ z)
        if den <= 0:
            raise RuntimeError(f"{cid}: zero CAP-measured target mass.")
        alignment = num / den
        measured_fraction = den / total_mass

        sens = np.nan
        if z_sens is not None:
            sens = float(w @ z_sens) / den

        q = mass_qa[mass_qa["component_id"].eq(cid)].iloc[0]
        rows.append({
            "method": "RAW_SHARED_SPECIFIC_PRIMARY",
            "component_id": cid,
            "component_type": str(q["component_type"]),
            "syndrome_code": str(q["syndrome_code"]),
            "observed_alignment": alignment,
            "alignment_raw_numerator": num,
            "CAP_measured_raw_target_mass": den,
            "total_raw_projected_target_mass": total_mass,
            "CAP_measured_fraction_of_raw_target_mass": measured_fraction,
            "CAP_unmeasured_fraction_of_raw_target_mass": 1.0 - measured_fraction,
            "raw_component_herb_mass": float(q["raw_component_herb_mass"]),
            "mapped_component_herb_mass": float(q["mapped_component_herb_mass"]),
            "mapped_weight_fraction": float(q["mapped_weight_fraction"]),
            "positive_weight_targets_total": int(np.sum(w > 0)),
            "positive_weight_targets_measured_in_CAP": int(np.sum((w > 0) & (measured > 0))),
            "observed_alignment_neff_weighted_sensitivity": sens,
        })

        c = universe[[
            "gene", "target_herb_degree", "target_idf",
            "network_degree_unweighted", "network_degree_weighted", "in_STRING",
            "matched_stratum", CAP_PRIMARY_Z_COLUMN, "measured_in_CAP_primary",
        ]].copy()
        c["component_id"] = cid
        c["raw_component_target_weight"] = w
        c["within_component_sum1_weight"] = w / total_mass
        c["CAP_measured_component_weight"] = np.where(measured > 0, w / den, 0.0)
        c["alignment_contribution"] = c["CAP_measured_component_weight"] * c[CAP_PRIMARY_Z_COLUMN].fillna(0.0)
        contributions.append(c)

        raw_cache[cid] = {"num": num, "den": den, "alignment": alignment, "total": total_mass}

    return pd.DataFrame(rows), pd.concat(contributions, ignore_index=True), raw_cache


def mixture_decomposition_observed(raw_cache: dict) -> pd.DataFrame:
    sh = raw_cache["SHARED"]
    rows = []
    for s in SYNDROME_CODES:
        sp = raw_cache[SPECIFIC_COMPONENT_BY_SYNDROME[s]]
        total_den = sh["den"] + sp["den"]
        total_num = sh["num"] + sp["num"]
        total_target_mass = sh["total"] + sp["total"]
        alignment = total_num / total_den
        rows.append({
            "syndrome_code": s,
            "reconstructed_complete_alignment": alignment,
            "complete_CAP_measured_raw_target_mass": total_den,
            "complete_raw_projected_target_mass": total_target_mass,
            "shared_fraction_of_total_projected_target_mass": sh["total"] / total_target_mass,
            "specific_fraction_of_total_projected_target_mass": sp["total"] / total_target_mass,
            "shared_fraction_of_CAP_measured_raw_target_mass": sh["den"] / total_den,
            "specific_fraction_of_CAP_measured_raw_target_mass": sp["den"] / total_den,
            "shared_alignment": sh["alignment"],
            "specific_alignment": sp["alignment"],
            "shared_contribution_to_complete_alignment": sh["num"] / total_den,
            "specific_contribution_to_complete_alignment": sp["num"] / total_den,
            "contribution_sum": total_num / total_den,
        })
    return pd.DataFrame(rows)


# =============================================================================
# 6. STRING RWR ADJOINT
# =============================================================================

def build_string_network(edges: pd.DataFrame) -> dict:
    nodes = sorted(set(edges["u"]) | set(edges["v"]))
    index = {g: i for i, g in enumerate(nodes)}
    i = edges["u"].map(index).to_numpy(dtype=int)
    j = edges["v"].map(index).to_numpy(dtype=int)
    w = edges["edge_weight"].to_numpy(dtype=float)
    A = coo_matrix(
        (np.concatenate([w, w]), (np.concatenate([i, j]), np.concatenate([j, i]))),
        shape=(len(nodes), len(nodes)), dtype=float,
    ).tocsr()
    degree_weighted = np.asarray(A.sum(axis=1)).ravel()
    degree_unweighted = np.diff(A.indptr).astype(int)
    if np.any(degree_weighted <= 0):
        raise RuntimeError("Zero weighted-degree node in STRING network.")
    P = (diags(1.0 / degree_weighted, offsets=0, format="csr") @ A).tocsr()
    row_sums = np.asarray(P.sum(axis=1)).ravel()
    max_row_error = float(np.max(np.abs(row_sums - 1.0)))
    if max_row_error > 1e-12:
        raise RuntimeError("STRING transition matrix is not row-stochastic.")
    return {
        "nodes": nodes, "index": index, "P": P,
        "degree_unweighted": degree_unweighted,
        "degree_weighted": degree_weighted,
        "max_row_sum_error": max_row_error,
    }


def solve_adjoint_rwr(P: csr_matrix, vector: np.ndarray, restart: float) -> tuple[np.ndarray, dict]:
    x = np.asarray(vector, dtype=float).reshape(-1)
    q = np.zeros_like(x)
    for iteration in range(1, ADJOINT_MAX_ITER + 1):
        q_new = restart * x + (1.0 - restart) * (P @ q)
        delta = float(np.max(np.abs(q_new - q)))
        q = q_new
        if delta < ADJOINT_TOL_MAXABS:
            return q, {
                "iterations": int(iteration),
                "final_max_abs_delta": delta,
                "converged": True,
            }
    raise RuntimeError(f"Adjoint RWR failed to converge; last delta={delta:.3e}")


def build_adjoint_vectors(network: dict, cap: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, dict]:
    cap_z = dict(zip(cap["gene"], cap[CAP_PRIMARY_Z_COLUMN]))
    z = np.zeros(len(network["nodes"]), dtype=float)
    measured = np.zeros(len(network["nodes"]), dtype=float)
    for i, gene in enumerate(network["nodes"]):
        if gene in cap_z:
            z[i] = float(cap_z[gene])
            measured[i] = 1.0
    qz, qaz = solve_adjoint_rwr(network["P"], z, PRIMARY_RWR_RESTART)
    qm, qam = solve_adjoint_rwr(network["P"], measured, PRIMARY_RWR_RESTART)
    qa = {
        "CAP_genes_in_STRING": int(measured.sum()),
        "adjoint_Z_iterations": qaz["iterations"],
        "adjoint_Z_final_max_abs_delta": qaz["final_max_abs_delta"],
        "adjoint_measurement_iterations": qam["iterations"],
        "adjoint_measurement_final_max_abs_delta": qam["final_max_abs_delta"],
    }
    return qz, qm, qa


def map_adjoint_to_target_universe(
    universe: pd.DataFrame,
    network: dict,
    qz: np.ndarray,
    qm: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    qzu = np.zeros(len(universe), dtype=float)
    qmu = np.zeros(len(universe), dtype=float)
    for i, gene in enumerate(universe["gene"]):
        j = network["index"].get(gene)
        if j is not None:
            qzu[i] = qz[j]
            qmu[i] = qm[j]
    return qzu, qmu


def observed_rwr_components(
    W_raw: np.ndarray,
    universe: pd.DataFrame,
    qzu: np.ndarray,
    qmu: np.ndarray,
    mass_qa: pd.DataFrame,
) -> tuple[pd.DataFrame, dict]:
    in_string = universe["in_STRING"].astype(float).to_numpy()
    rows = []
    cache = {}
    for i, cid in enumerate(COMPONENT_IDS):
        w = W_raw[i]
        total = float(w.sum())
        num = float(w @ qzu)
        den = float(w @ qmu)
        if den <= 0:
            raise RuntimeError(f"{cid}: zero RWR CAP-measured mass.")
        align = num / den
        q = mass_qa[mass_qa["component_id"].eq(cid)].iloc[0]
        rows.append({
            "method": "STRING_RWR_SHARED_SPECIFIC_SECONDARY",
            "component_id": cid,
            "component_type": str(q["component_type"]),
            "syndrome_code": str(q["syndrome_code"]),
            "observed_alignment": align,
            "adjoint_alignment_numerator_scaled": num,
            "adjoint_CAP_measured_RWR_mass_scaled": den,
            "CAP_measured_RWR_mass_fraction_if_seed_sum1": den / total,
            "STRING_seed_raw_target_mass_fraction": float(w @ in_string) / total,
            "total_raw_projected_target_mass": total,
        })
        cache[cid] = {"num": num, "den": den, "alignment": align, "total": total}
    return pd.DataFrame(rows), cache


# =============================================================================
# 7. MATCHED NULL — SAME PERMUTATION FOR ALL FOUR COMPONENTS
# =============================================================================

def stratum_indices(universe: pd.DataFrame) -> dict[str, np.ndarray]:
    strata = universe["matched_stratum"].astype(str).to_numpy()
    return {
        s: np.flatnonzero(strata == s)
        for s in sorted(set(strata))
    }


def null_invariant_qa(
    universe: pd.DataFrame,
    W_raw: np.ndarray,
    strata: dict[str, np.ndarray],
) -> pd.DataFrame:
    rows = []
    for i, cid in enumerate(COMPONENT_IDS):
        w = W_raw[i]
        for name, idx in strata.items():
            rows.append({
                "component_id": cid,
                "matched_stratum": name,
                "stratum_n": int(len(idx)),
                "positive_weight_n": int(np.sum(w[idx] > 0)),
                "zero_weight_n": int(np.sum(w[idx] == 0)),
                "stratum_raw_target_weight_sum": float(w[idx].sum()),
                "mean_network_degree_unweighted": float(
                    universe.loc[idx, "network_degree_unweighted"].fillna(0).mean()
                ),
                "mean_target_herb_degree": float(
                    universe.loc[idx, "target_herb_degree"].mean()
                ),
            })
    return pd.DataFrame(rows)


def generate_matched_null(
    universe: pd.DataFrame,
    W_raw: np.ndarray,
    qzu: np.ndarray,
    qmu: np.ndarray,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    rng = np.random.default_rng(RANDOM_SEED)
    strata = stratum_indices(universe)
    measured = universe["measured_in_CAP_primary"].astype(float).to_numpy()
    z = universe[CAP_PRIMARY_Z_COLUMN].fillna(0.0).to_numpy(dtype=float)
    z_sens = (
        universe[CAP_SENSITIVITY_Z_COLUMN].fillna(0.0).to_numpy(dtype=float)
        if CAP_SENSITIVITY_Z_COLUMN in universe.columns else None
    )

    n_c = len(COMPONENT_IDS)
    raw_null = np.empty((N_NULL, n_c), dtype=float)
    rwr_null = np.empty((N_NULL, n_c), dtype=float)
    raw_den_null = np.empty((N_NULL, n_c), dtype=float)
    rwr_den_null = np.empty((N_NULL, n_c), dtype=float)
    raw_num_null = np.empty((N_NULL, n_c), dtype=float)
    rwr_num_null = np.empty((N_NULL, n_c), dtype=float)
    sens_null = np.empty((N_NULL, n_c), dtype=float) if z_sens is not None else None

    # Exact reconstructed complete-profile nulls for hard cross-checking.
    recon_raw = {s: np.empty(N_NULL, dtype=float) for s in SYNDROME_CODES}
    recon_rwr = {s: np.empty(N_NULL, dtype=float) for s in SYNDROME_CODES}

    base = np.arange(len(universe), dtype=int)
    shared_idx = COMPONENT_IDS.index("SHARED")

    for k in range(N_NULL):
        perm = base.copy()
        for idx in strata.values():
            perm[idx] = rng.permutation(idx)

        # Same donor permutation for SHARED and every SPECIFIC component.
        WP = W_raw[:, perm]

        raw_den = WP @ measured
        raw_num = WP @ z
        if np.any(raw_den <= 0):
            raise RuntimeError("A null replicate has zero CAP-measured raw target mass.")
        raw_null[k, :] = raw_num / raw_den
        raw_den_null[k, :] = raw_den
        raw_num_null[k, :] = raw_num

        if z_sens is not None:
            sens_null[k, :] = (WP @ z_sens) / raw_den

        rwr_den = WP @ qmu
        rwr_num = WP @ qzu
        if np.any(rwr_den <= 0):
            raise RuntimeError("A null replicate has zero CAP-measured RWR mass.")
        rwr_null[k, :] = rwr_num / rwr_den
        rwr_den_null[k, :] = rwr_den
        rwr_num_null[k, :] = rwr_num

        for s in SYNDROME_CODES:
            sp_idx = COMPONENT_IDS.index(SPECIFIC_COMPONENT_BY_SYNDROME[s])
            recon_raw[s][k] = (
                raw_num[shared_idx] + raw_num[sp_idx]
            ) / (
                raw_den[shared_idx] + raw_den[sp_idx]
            )
            recon_rwr[s][k] = (
                rwr_num[shared_idx] + rwr_num[sp_idx]
            ) / (
                rwr_den[shared_idx] + rwr_den[sp_idx]
            )

    parts = []
    for i, cid in enumerate(COMPONENT_IDS):
        d = pd.DataFrame({
            "null_iteration": np.arange(1, N_NULL + 1),
            "component_id": cid,
            "raw_alignment_null": raw_null[:, i],
            "rwr_alignment_null": rwr_null[:, i],
            "raw_CAP_measured_target_mass_null": raw_den_null[:, i],
            "raw_alignment_numerator_null": raw_num_null[:, i],
            "rwr_CAP_measured_mass_scaled_null": rwr_den_null[:, i],
            "rwr_alignment_numerator_scaled_null": rwr_num_null[:, i],
        })
        if sens_null is not None:
            d["raw_alignment_neff_sensitivity_null"] = sens_null[:, i]
        parts.append(d)

    null_long = pd.concat(parts, ignore_index=True)
    recon = {"raw": recon_raw, "rwr": recon_rwr}
    return null_long, null_invariant_qa(universe, W_raw, strata), recon


# =============================================================================
# 8. SUMMARIES + PAIRED CONTRASTS
# =============================================================================

def summarize_components(
    raw_obs: pd.DataFrame,
    rwr_obs: pd.DataFrame,
    null_long: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for cid in COMPONENT_IDS:
        ro = raw_obs[raw_obs["component_id"].eq(cid)].iloc[0]
        rn = null_long[null_long["component_id"].eq(cid)]["raw_alignment_null"].to_numpy()
        rows.append({
            "method": "RAW_SHARED_SPECIFIC_PRIMARY",
            "component_id": cid,
            "component_type": ro["component_type"],
            "syndrome_code": ro["syndrome_code"],
            "observed_alignment": float(ro["observed_alignment"]),
            "observed_CAP_measured_mass_fraction": float(ro["CAP_measured_fraction_of_raw_target_mass"]),
            "n_null": N_NULL,
            "null_model": "frozen_19_strata_same_permutation_all_components_raw_weight_permutation",
            **empirical_p_values(rn, float(ro["observed_alignment"])),
        })

        rw = rwr_obs[rwr_obs["component_id"].eq(cid)].iloc[0]
        wn = null_long[null_long["component_id"].eq(cid)]["rwr_alignment_null"].to_numpy()
        rows.append({
            "method": "STRING_RWR_SHARED_SPECIFIC_SECONDARY",
            "component_id": cid,
            "component_type": rw["component_type"],
            "syndrome_code": rw["syndrome_code"],
            "observed_alignment": float(rw["observed_alignment"]),
            "observed_CAP_measured_mass_fraction": float(rw["CAP_measured_RWR_mass_fraction_if_seed_sum1"]),
            "n_null": N_NULL,
            "null_model": "same_frozen_matched_permutation_then_exact_adjoint_RWR",
            **empirical_p_values(wn, float(rw["observed_alignment"])),
        })
    return pd.DataFrame(rows)


def summarize_neff_sensitivity(raw_obs: pd.DataFrame, null_long: pd.DataFrame) -> pd.DataFrame:
    if "raw_alignment_neff_sensitivity_null" not in null_long.columns:
        return pd.DataFrame()
    rows = []
    for cid in COMPONENT_IDS:
        obs = float(raw_obs.loc[
            raw_obs["component_id"].eq(cid),
            "observed_alignment_neff_weighted_sensitivity",
        ].iloc[0])
        null = null_long.loc[
            null_long["component_id"].eq(cid),
            "raw_alignment_neff_sensitivity_null",
        ].to_numpy()
        rows.append({
            "method": "RAW_SHARED_SPECIFIC_CAP_NEFF_META_Z_SENSITIVITY",
            "component_id": cid,
            "observed_alignment": obs,
            "n_null": N_NULL,
            **empirical_p_values(null, obs),
        })
    return pd.DataFrame(rows)


def paired_component_contrasts(
    raw_obs: pd.DataFrame,
    rwr_obs: pd.DataFrame,
    null_long: pd.DataFrame,
) -> pd.DataFrame:
    contrast_pairs = [
        ("SPECIFIC_PDOL", "SHARED", "SPECIFIC_MINUS_SHARED"),
        ("SPECIFIC_PHOL", "SHARED", "SPECIFIC_MINUS_SHARED"),
        ("SPECIFIC_WHIL", "SHARED", "SPECIFIC_MINUS_SHARED"),
        ("SPECIFIC_PDOL", "SPECIFIC_PHOL", "SPECIFIC_PAIRWISE"),
        ("SPECIFIC_PDOL", "SPECIFIC_WHIL", "SPECIFIC_PAIRWISE"),
        ("SPECIFIC_PHOL", "SPECIFIC_WHIL", "SPECIFIC_PAIRWISE"),
    ]

    raw_map = dict(zip(raw_obs["component_id"], raw_obs["observed_alignment"]))
    rwr_map = dict(zip(rwr_obs["component_id"], rwr_obs["observed_alignment"]))
    wide_raw = null_long.pivot(index="null_iteration", columns="component_id", values="raw_alignment_null")
    wide_rwr = null_long.pivot(index="null_iteration", columns="component_id", values="rwr_alignment_null")

    rows = []
    for method, obs_map, wide in [
        ("RAW_SHARED_SPECIFIC_PRIMARY", raw_map, wide_raw),
        ("STRING_RWR_SHARED_SPECIFIC_SECONDARY", rwr_map, wide_rwr),
    ]:
        for a, b, ctype in contrast_pairs:
            obs = float(obs_map[a] - obs_map[b])
            null = wide[a].to_numpy() - wide[b].to_numpy()
            rows.append({
                "method": method,
                "contrast_type": ctype,
                "component_A": a,
                "component_B": b,
                "observed_A_minus_B": obs,
                "n_null": N_NULL,
                **empirical_p_values(null, obs),
            })
    return pd.DataFrame(rows)


# =============================================================================
# 9. HARD RECONSTRUCTION OF PREVIOUS COMPLETE-PROFILE RESULTS
# =============================================================================

def observed_complete_crosscheck(
    mixture: pd.DataFrame,
    rwr_cache: dict,
    complete_raw_obs: pd.DataFrame,
    complete_rwr_obs: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    sh = rwr_cache["SHARED"]
    for s in SYNDROME_CODES:
        sp = rwr_cache[SPECIFIC_COMPONENT_BY_SYNDROME[s]]
        reconstructed_rwr = (sh["num"] + sp["num"]) / (sh["den"] + sp["den"])
        reconstructed_raw = float(mixture.loc[
            mixture["syndrome_code"].eq(s), "reconstructed_complete_alignment"
        ].iloc[0])
        frozen_raw = float(complete_raw_obs.loc[
            complete_raw_obs["syndrome_code"].eq(s), "observed_alignment"
        ].iloc[0])
        frozen_rwr = float(complete_rwr_obs.loc[
            complete_rwr_obs["syndrome_code"].eq(s), "observed_alignment"
        ].iloc[0])
        rows.append({
            "syndrome_code": s,
            "reconstructed_complete_RAW_alignment": reconstructed_raw,
            "frozen_complete_RAW_alignment": frozen_raw,
            "abs_error_RAW": abs(reconstructed_raw - frozen_raw),
            "reconstructed_complete_RWR_alignment": reconstructed_rwr,
            "frozen_complete_RWR_alignment": frozen_rwr,
            "abs_error_RWR": abs(reconstructed_rwr - frozen_rwr),
        })
    qa = pd.DataFrame(rows)
    if qa[["abs_error_RAW", "abs_error_RWR"]].to_numpy().max() > OBSERVED_CROSSCHECK_TOL:
        raise RuntimeError("Reconstructed complete observed alignment does not match frozen v1.0.1.")
    return qa


def previous_null_crosscheck(recon: dict, complete_null: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for s in SYNDROME_CODES:
        old = complete_null[complete_null["syndrome_code"].eq(s)].sort_values("null_iteration")
        if len(old) != N_NULL:
            raise RuntimeError(f"{s}: previous complete null has {len(old)} rows, expected {N_NULL}.")
        raw_old = old["raw_alignment_null"].to_numpy(dtype=float)
        rwr_old = old["rwr_alignment_null"].to_numpy(dtype=float)
        raw_err = np.abs(recon["raw"][s] - raw_old)
        rwr_err = np.abs(recon["rwr"][s] - rwr_old)
        rows.append({
            "syndrome_code": s,
            "max_abs_error_reconstructed_RAW_null": float(raw_err.max()),
            "mean_abs_error_reconstructed_RAW_null": float(raw_err.mean()),
            "max_abs_error_reconstructed_RWR_null": float(rwr_err.max()),
            "mean_abs_error_reconstructed_RWR_null": float(rwr_err.mean()),
        })
    qa = pd.DataFrame(rows)
    if qa[[
        "max_abs_error_reconstructed_RAW_null",
        "max_abs_error_reconstructed_RWR_null",
    ]].to_numpy().max() > NULL_CROSSCHECK_TOL:
        raise RuntimeError(
            "Reconstructed complete-profile null does not reproduce frozen v1.0.1 null. "
            "Do not proceed: null method or randomization order has drifted."
        )
    return qa


# =============================================================================
# 10. MAIN
# =============================================================================

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 88)
    print("SHARED / SPECIFIC CAP ALIGNMENT + MATCHED NULL v1.0")
    print("RAW COMPONENT ALIGNMENT = PRIMARY; STRING RWR = SECONDARY")
    print("=" * 88)

    print("\n[1/10] Resolving and verifying all frozen inputs...")
    (
        profiles, mass_qa, universe, cap, edges, complete,
        complete_raw_obs, complete_rwr_obs, complete_null,
        input_paths, input_hashes,
    ) = load_inputs()
    print(f"  frozen target universe: {len(universe):,}")
    print(f"  frozen CAP continuous universe: {len(cap):,}")
    print(f"  frozen matched strata: {universe['matched_stratum'].nunique():,}")
    print(f"  frozen STRING edges: {len(edges):,}")

    print("\n[2/10] Cross-checking CAP values and component target vectors...")
    cap_err = validate_cap_against_frozen_universe(universe, cap)
    W_raw, W_sum1, component_input_qa = build_component_matrices(profiles, mass_qa, universe)
    W_complete = build_complete_raw_matrix(complete, universe)
    reconstruction_qa = component_reconstruction_qa(W_raw, W_complete)
    print(f"  CAP target-universe cross-check max error: {cap_err:.3e}")
    print(reconstruction_qa.to_string(index=False))

    print("\n[3/10] Computing observed RAW shared/specific alignment...")
    raw_obs, raw_contrib, raw_cache = observed_raw_components(universe, W_raw, mass_qa)
    mixture = mixture_decomposition_observed(raw_cache)
    print(raw_obs[[
        "component_id", "observed_alignment", "CAP_measured_fraction_of_raw_target_mass",
        "total_raw_projected_target_mass",
    ]].to_string(index=False))

    print("\n[4/10] Building frozen STRING network and adjoint disease vectors...")
    network = build_string_network(edges)
    qz, qm, adjoint_qa = build_adjoint_vectors(network, cap)
    qzu, qmu = map_adjoint_to_target_universe(universe, network, qz, qm)
    print(f"  STRING nodes: {len(network['nodes']):,}")
    print(f"  transition max row-sum error: {network['max_row_sum_error']:.3e}")

    print("\n[5/10] Computing observed secondary RWR component alignment...")
    rwr_obs, rwr_cache = observed_rwr_components(W_raw, universe, qzu, qmu, mass_qa)
    print(rwr_obs[[
        "component_id", "observed_alignment", "CAP_measured_RWR_mass_fraction_if_seed_sum1",
    ]].to_string(index=False))

    print("\n[6/10] Cross-checking reconstructed complete observed results...")
    observed_crosscheck = observed_complete_crosscheck(
        mixture, rwr_cache, complete_raw_obs, complete_rwr_obs
    )
    print(observed_crosscheck.to_string(index=False))

    print(f"\n[7/10] Running {N_NULL:,} frozen-strata matched-null permutations...")
    null_long, invariant_qa, recon_null = generate_matched_null(
        universe, W_raw, qzu, qmu
    )

    print("\n[8/10] Reconstructing the previously frozen complete-profile null...")
    previous_null_qa = previous_null_crosscheck(recon_null, complete_null)
    print(previous_null_qa.to_string(index=False))

    print("\n[9/10] Summarizing component inference and paired contrasts...")
    summary = summarize_components(raw_obs, rwr_obs, null_long)
    contrasts = paired_component_contrasts(raw_obs, rwr_obs, null_long)
    neff = summarize_neff_sensitivity(raw_obs, null_long)
    print(summary[[
        "method", "component_id", "observed_alignment", "null_mean", "Z_align",
        "empirical_p_upper", "empirical_p_two_sided",
    ]].to_string(index=False))

    print("\n[10/10] Writing reproducible outputs...")

    raw_obs.to_csv(
        OUTPUT_DIR / "CAP_shared_specific_alignment_observed_RAW_primary_v1_0.csv",
        index=False, encoding="utf-8-sig",
    )
    rwr_obs.to_csv(
        OUTPUT_DIR / "CAP_shared_specific_alignment_observed_RWR_secondary_v1_0.csv",
        index=False, encoding="utf-8-sig",
    )
    summary.to_csv(
        OUTPUT_DIR / "CAP_shared_specific_alignment_matched_null_summary_v1_0.csv",
        index=False, encoding="utf-8-sig",
    )
    contrasts.to_csv(
        OUTPUT_DIR / "CAP_shared_specific_alignment_component_contrasts_v1_0.csv",
        index=False, encoding="utf-8-sig",
    )
    mixture.to_csv(
        OUTPUT_DIR / "CAP_shared_specific_alignment_mixture_decomposition_v1_0.csv",
        index=False, encoding="utf-8-sig",
    )
    neff.to_csv(
        OUTPUT_DIR / "CAP_shared_specific_alignment_neff_metaZ_sensitivity_v1_0.csv",
        index=False, encoding="utf-8-sig",
    )
    raw_contrib.to_csv(
        OUTPUT_DIR / "CAP_shared_specific_alignment_RAW_gene_contributions_v1_0.csv",
        index=False, encoding="utf-8-sig",
    )
    component_input_qa.to_csv(
        OUTPUT_DIR / "shared_specific_component_input_QA_v1_0.csv",
        index=False, encoding="utf-8-sig",
    )
    reconstruction_qa.to_csv(
        OUTPUT_DIR / "shared_specific_complete_target_reconstruction_QA_v1_0.csv",
        index=False, encoding="utf-8-sig",
    )
    observed_crosscheck.to_csv(
        OUTPUT_DIR / "shared_specific_complete_observed_alignment_crosscheck_QA_v1_0.csv",
        index=False, encoding="utf-8-sig",
    )
    previous_null_qa.to_csv(
        OUTPUT_DIR / "shared_specific_previous_complete_null_reconstruction_QA_v1_0.csv",
        index=False, encoding="utf-8-sig",
    )
    invariant_qa.to_csv(
        OUTPUT_DIR / "shared_specific_matched_null_weight_invariant_QA_v1_0.csv",
        index=False, encoding="utf-8-sig",
    )

    # Save the exact frozen universe together with component raw weights for audit.
    audit_universe = universe.copy()
    for i, cid in enumerate(COMPONENT_IDS):
        audit_universe[f"raw_weight_{cid}"] = W_raw[i]
        audit_universe[f"sum1_weight_{cid}"] = W_sum1[i]
    audit_universe.to_csv(
        OUTPUT_DIR / "shared_specific_target_universe_weights_and_frozen_strata_v1_0.csv",
        index=False, encoding="utf-8-sig",
    )

    if SAVE_FULL_NULL_DISTRIBUTION:
        null_long.to_csv(
            OUTPUT_DIR / "CAP_shared_specific_alignment_full_matched_null_distribution_v1_0.csv.gz",
            index=False, compression="gzip",
        )

    metadata = {
        "analysis_name": "shared/specific CAP alignment + matched null v1.0",
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "status": "FREEZE_CANDIDATE_AFTER_AUDIT",
        "primary_representation": "frozen raw ubiquity-corrected shared/specific target components",
        "primary_disease_vector": CAP_PRIMARY_Z_COLUMN,
        "secondary_network_layer": {
            "method": "STRING stationary RWR via exact adjoint",
            "restart_probability": PRIMARY_RWR_RESTART,
            "role": "secondary network-localization/sensitivity only",
        },
        "null": {
            "n_replicates": N_NULL,
            "random_seed": RANDOM_SEED,
            "strata_source": str(input_paths["frozen_target_universe_strata"]),
            "n_frozen_strata": int(universe["matched_stratum"].nunique()),
            "permutation": "same within-stratum donor permutation applied to SHARED and all SPECIFIC components",
            "preserves": [
                "component total raw projected target mass",
                "component global positive/zero count",
                "component positive/zero count within each frozen stratum",
                "component raw target-weight mass within each frozen stratum",
                "shared+specific linear reconstruction under every null permutation",
            ],
            "empirical_p_correction": "(count+1)/(N+1)",
        },
        "hard_QA": {
            "CAP_target_universe_max_abs_error": cap_err,
            "max_shared_plus_specific_target_reconstruction_error": float(
                reconstruction_qa["max_abs_target_weight_reconstruction_error"].max()
            ),
            "max_complete_observed_RAW_crosscheck_error": float(
                observed_crosscheck["abs_error_RAW"].max()
            ),
            "max_complete_observed_RWR_crosscheck_error": float(
                observed_crosscheck["abs_error_RWR"].max()
            ),
            "max_previous_complete_RAW_null_reconstruction_error": float(
                previous_null_qa["max_abs_error_reconstructed_RAW_null"].max()
            ),
            "max_previous_complete_RWR_null_reconstruction_error": float(
                previous_null_qa["max_abs_error_reconstructed_RWR_null"].max()
            ),
            "STRING_transition_max_row_sum_error": network["max_row_sum_error"],
            **adjoint_qa,
        },
        "input_paths": {k: str(v) for k, v in input_paths.items()},
        "input_sha256": input_hashes,
        "interpretation_guardrails": [
            "positive CAP alignment is disease-context localization, not therapeutic reversal",
            "component sum1 normalization does not imply equal component magnitude",
            "specific profiles are residual weighted representations, not discrete molecular identities",
            "no upstream parameter may be retuned after inspecting these CAP results",
        ],
    }

    with open(
        OUTPUT_DIR / "shared_specific_CAP_alignment_matched_null_run_metadata_v1_0.json",
        "w", encoding="utf-8",
    ) as f:
        json.dump(metadata, f, ensure_ascii=False, indent=2)

    freeze_text = f"""SHARED / SPECIFIC CAP ALIGNMENT + MATCHED NULL v1.0
FREEZE CANDIDATE — AUDIT BEFORE FINAL FREEZE

STATUS
------
DOWNSTREAM FREEZE CANDIDATE.
Do not modify upstream shared/specific profiles, CAP meta-Z, frozen 19 null
strata, STRING network, or restart probability after inspecting these results.

PRIMARY COMPONENT STATISTIC
---------------------------
For frozen raw ubiquity-corrected component target weights U_cg:

A_raw(c) = sum_g U_cg I_CAP(g) Z_CAP,g / sum_g U_cg I_CAP(g)

Raw component target mass is retained. The weighted-mean alignment is scale
invariant, but raw mass is required for exact shared/specific decomposition.

NULL
----
N = {N_NULL}
Random seed = {RANDOM_SEED}
Frozen matched strata = {universe['matched_stratum'].nunique()}

The SAME within-stratum donor permutation is applied to SHARED and all three
SPECIFIC components. Therefore every null replicate exactly preserves:

perm(U_original,s) = perm(U_shared) + perm(U_specific,s)

and the reconstructed complete-profile null can be checked against the already
frozen v1.0.1 complete-profile null distribution.

HARD QA
-------
CAP target-universe max error:
{cap_err:.3e}

Shared + specific target reconstruction max error:
{reconstruction_qa['max_abs_target_weight_reconstruction_error'].max():.3e}

Complete observed RAW cross-check max error:
{observed_crosscheck['abs_error_RAW'].max():.3e}

Complete observed RWR cross-check max error:
{observed_crosscheck['abs_error_RWR'].max():.3e}

Previous complete RAW null reconstruction max error:
{previous_null_qa['max_abs_error_reconstructed_RAW_null'].max():.3e}

Previous complete RWR null reconstruction max error:
{previous_null_qa['max_abs_error_reconstructed_RWR_null'].max():.3e}

INTERPRETATION POLICY
---------------------
1. RAW shared/specific target alignment is primary.
2. STRING RWR r=0.5 is secondary network-localization only.
3. Effective-N CAP meta-Z is sensitivity only.
4. SHARED and SPECIFIC component alignments are evaluated against their own
   frozen-strata matched null distributions.
5. SPECIFIC-minus-SHARED and pairwise SPECIFIC contrasts use paired nulls from
   the same target-label permutation.
6. Raw component mass must be reported alongside normalized/compositional
   alignment; sum1 profiles do not imply equal component magnitude.
7. Positive CAP meta-Z alignment is disease-context localization, NOT
   therapeutic reversal and NOT treatment efficacy.
8. Do not tune upstream parameters, null strata, STRING settings, or component
   definitions based on significance.
"""
    (
        OUTPUT_DIR / "SHARED_SPECIFIC_CAP_ALIGNMENT_FREEZE_CANDIDATE_v1_0.txt"
    ).write_text(freeze_text, encoding="utf-8")

    print("\n" + "=" * 88)
    print("Completed. STATUS: FREEZE CANDIDATE — audit outputs before final freeze.")
    print("=" * 88)
    print(f"\nOutputs:\n{OUTPUT_DIR}")


if __name__ == "__main__":
    main()
