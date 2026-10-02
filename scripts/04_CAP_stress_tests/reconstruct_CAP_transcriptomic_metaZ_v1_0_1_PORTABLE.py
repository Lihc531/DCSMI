# -*- coding: utf-8 -*-
"""
reconstruct_CAP_transcriptomic_metaZ_v1_0_1.py

Disease-layer reconstruction for community-acquired pneumonia (CAP) using two
independently processed blood transcriptomic cohorts:

    GSE103119  (microarray; CAP vs healthy control)
    GSE196399  (RNA-seq; severe CAP vs healthy control)

PRIMARY DISEASE REPRESENTATION
------------------------------
For every gene tested in BOTH cohorts:

1) Convert each cohort's two-sided differential-expression P value to a signed
   standard-normal statistic using the sign of logFC:

       Z_i = sign(logFC_i) * Phi^{-1}(1 - P_i/2)

2) Combine the two cohort Z statistics with EQUAL-WEIGHT Stouffer meta-analysis:

       Z_meta = (Z_103119 + Z_196399) / sqrt(2)

3) Convert Z_meta to a two-sided meta P value and apply BH FDR across the shared
   analysis universe.

The CONTINUOUS equal-weight Z_meta vector over all shared genes is the primary
CAP disease program for downstream syndrome-disease alignment. It is NOT
filtered by FDR before continuous alignment.

DIRECTION-CONSISTENT MODULE
---------------------------
A descriptive/set-based CAP module is additionally defined as:

    meta_FDR < 0.05 AND sign(logFC_103119) == sign(logFC_196399)

This module is useful for pathway/module summaries, but it does not replace the
continuous Z_meta vector in the primary alignment statistic.

PREDECLARED SENSITIVITY
-----------------------
A sample-size-weighted Stouffer statistic is calculated using the two-group
EFFECTIVE sample size

    n_eff = n_case * n_control / (n_case + n_control)
    w_i   = sqrt(n_eff_i)

and

    Z_weighted = sum_i(w_i Z_i) / sqrt(sum_i w_i^2)

This is a sensitivity analysis only. Equal-weight Stouffer remains primary.

LEGACY / STRICT REPLICATION MODULE
----------------------------------
For continuity with the earlier CAP analysis, the script also reports a strict
replicated-DEG set:

    cohort FDR < 0.05 in BOTH cohorts
    |logFC| >= 1 in BOTH cohorts
    same direction in BOTH cohorts

This set is descriptive/sensitivity only and is NOT used to tune the meta-Z
method.

IMPORTANT METHOD FREEZE RULES
-----------------------------
- Do not use syndrome scores, STRING/RWR results, or CAP alignment outcomes to
  choose the meta-Z method or thresholds in this script.
- Do not force the shared-gene universe to reproduce an older row count. The
  universe is reconstructed from the CURRENT gene-level all-results files and
  is reported explicitly.
- Do not meta-analyze raw logFC across these platforms; the primary cross-cohort
  combination is performed in signed-Z space.
- Do not interpret positive Z as therapeutic benefit. Positive Z means higher
  expression in CAP relative to healthy control; negative Z means lower.
"""

from __future__ import annotations

from pathlib import Path
import os
from datetime import datetime, timezone
import hashlib
import json
import math
import re
import warnings

import numpy as np
import pandas as pd
from scipy.stats import norm, spearmanr, pearsonr


# =============================================================================
# 1. USER SETTINGS — EDIT PATHS HERE
# =============================================================================

PROJECT_ROOT = Path(os.environ.get("DCSMI_PROJECT_ROOT", Path(__file__).resolve().parents[2]))

# Prefer the gene-level ALL-results tables, not prefiltered DEG tables.
GSE103119_FILE = PROJECT_ROOT / "R" / "GSE103119_gene_level_all.csv"
GSE196399_FILE = PROJECT_ROOT / "R" / "GSE196399_gene_level_all.csv"

# Optional fallback locations if your files are directly under PROJECT_ROOT.
GSE103119_FALLBACKS = [
    PROJECT_ROOT / "GSE103119_gene_level_all.csv",
    PROJECT_ROOT / "CAP" / "GSE103119_gene_level_all.csv",
]
GSE196399_FALLBACKS = [
    PROJECT_ROOT / "GSE196399_gene_level_all.csv",
    PROJECT_ROOT / "CAP" / "GSE196399_gene_level_all.csv",
]

OUTPUT_DIR = PROJECT_ROOT / "07_CAP_transcriptomic_metaZ_v1_0_1"

# Column names in the current gene-level files.
# Gene column is auto-detected because the current exported tables use SYMBOL,
# while some earlier intermediate tables used gene.
GENE_COLUMN = None
GENE_COLUMN_CANDIDATES = [
    "gene", "SYMBOL", "symbol", "Gene", "gene_symbol", "GeneSymbol"
]
LOGFC_COLUMN = "logFC"
P_COLUMN = "P.Value"
FDR_COLUMN = "adj.P.Val"
T_COLUMN = "t"  # retained for audit; not required for meta-Z if absent

# Current preprocessing contrast is CAP - Control in both cohorts.
# Keep +1 only if positive logFC means higher expression in CAP.
GSE103119_DIRECTION_MULTIPLIER = +1
GSE196399_DIRECTION_MULTIPLIER = +1

# Sample counts used ONLY for the predeclared weighted-Stouffer sensitivity.
# These should match the samples used to generate the gene-level tables.
SAMPLE_COUNTS = {
    "GSE103119": {"case": 152, "control": 20},
    "GSE196399": {"case": 56, "control": 21},
}

# Primary descriptive module threshold.
META_FDR_THRESHOLD = 0.05

# Strict replicated-DEG sensitivity threshold (legacy continuity only).
STRICT_COHORT_FDR_THRESHOLD = 0.05
STRICT_ABS_LOGFC_THRESHOLD = 1.0

# Numerical handling for two-sided P -> signed Z.
# Values smaller than this are clipped only to avoid +/-inf.
MIN_P_FOR_Z = 1e-300

# Gene symbols are normalized to uppercase to match the target/PPI layers.
UPPERCASE_GENE_SYMBOLS = True

# Canonical outputs are CSV/JSON/TXT; no Excel dependency is required.


# =============================================================================
# 2. HELPERS
# =============================================================================

def norm_gene(x) -> str:
    if pd.isna(x):
        return ""
    s = re.sub(r"\s+", "", str(x).strip())
    return s.upper() if UPPERCASE_GENE_SYMBOLS else s


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk_size)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def resolve_input(primary: Path, fallbacks: list[Path], label: str) -> Path:
    candidates = [primary] + list(fallbacks)
    existing = []
    seen = set()
    for p in candidates:
        key = str(p.resolve()) if p.exists() else str(p)
        if key in seen:
            continue
        seen.add(key)
        if p.exists():
            existing.append(p)

    if not existing:
        searched = "\n".join(f"  - {p}" for p in candidates)
        raise FileNotFoundError(
            f"Could not find {label}. Searched:\n{searched}"
        )

    # If multiple copies exist, protect against silent version drift.
    if len(existing) > 1:
        hashes = {p: sha256_file(p) for p in existing}
        unique_hashes = set(hashes.values())
        if len(unique_hashes) > 1:
            detail = "\n".join(f"  {p}: {h}" for p, h in hashes.items())
            raise RuntimeError(
                f"Multiple non-identical copies of {label} were found.\n"
                f"Please keep/point to one intended version:\n{detail}"
            )

    return existing[0]


def require_columns(df: pd.DataFrame, cols: list[str], name: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"{name} missing required columns: {missing}\n"
            f"Available columns: {list(df.columns)}"
        )


def detect_gene_column(df: pd.DataFrame, name: str) -> str:
    """Deterministically resolve the gene-symbol column."""
    if GENE_COLUMN is not None:
        if GENE_COLUMN not in df.columns:
            raise ValueError(
                f"{name} configured gene column '{GENE_COLUMN}' not found.\n"
                f"Available columns: {list(df.columns)}"
            )
        return GENE_COLUMN

    # Exact matching first, in declared priority order.
    for c in GENE_COLUMN_CANDIDATES:
        if c in df.columns:
            return c

    # Then case-insensitive exact matching.
    lower_map = {str(c).strip().lower(): c for c in df.columns}
    for c in GENE_COLUMN_CANDIDATES:
        if c.lower() in lower_map:
            return lower_map[c.lower()]

    raise ValueError(
        f"{name}: could not identify a gene-symbol column.\n"
        f"Tried: {GENE_COLUMN_CANDIDATES}\n"
        f"Available columns: {list(df.columns)}"
    )


def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR with NA-safe handling."""
    p = np.asarray(pvalues, dtype=float)
    out = np.full(len(p), np.nan, dtype=float)
    valid = np.isfinite(p)
    pv = p[valid]
    if len(pv) == 0:
        return out

    order = np.argsort(pv)
    ranked = pv[order]
    m = len(ranked)
    q = ranked * m / np.arange(1, m + 1, dtype=float)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0.0, 1.0)

    restored = np.empty(m, dtype=float)
    restored[order] = q
    out[valid] = restored
    return out


def signed_z_from_two_sided_p(logfc: np.ndarray, pvalue: np.ndarray) -> np.ndarray:
    logfc = np.asarray(logfc, dtype=float)
    p = np.asarray(pvalue, dtype=float)

    if np.any(~np.isfinite(logfc)) or np.any(~np.isfinite(p)):
        raise ValueError("signed_z_from_two_sided_p received non-finite values.")
    if np.any((p < 0) | (p > 1)):
        raise ValueError("P values must be within [0, 1].")

    p_clip = np.clip(p, MIN_P_FOR_Z, 1.0)
    magnitude = norm.isf(p_clip / 2.0)
    sign = np.sign(logfc)
    # A true logFC==0 gets Z=0 regardless of P.
    return sign * magnitude


def two_sided_p_from_z(z: np.ndarray) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    return 2.0 * norm.sf(np.abs(z))


def effective_two_group_n(n_case: int, n_control: int) -> float:
    if n_case <= 0 or n_control <= 0:
        raise ValueError("Sample counts must be positive.")
    return (n_case * n_control) / (n_case + n_control)


def safe_spearman(x, y) -> tuple[float, float]:
    r = spearmanr(x, y)
    return float(r.statistic), float(r.pvalue)


def safe_pearson(x, y) -> tuple[float, float]:
    r = pearsonr(x, y)
    return float(r.statistic), float(r.pvalue)


# =============================================================================
# 3. LOAD AND STANDARDIZE EACH COHORT
# =============================================================================

def load_cohort(
    path: Path,
    cohort: str,
    direction_multiplier: int,
) -> tuple[pd.DataFrame, dict]:
    raw = pd.read_csv(path)

    detected_gene_column = detect_gene_column(raw, path.name)
    required_non_gene = [LOGFC_COLUMN, P_COLUMN, FDR_COLUMN]
    require_columns(raw, required_non_gene, path.name)

    keep_cols = [detected_gene_column] + required_non_gene
    if T_COLUMN in raw.columns:
        keep_cols.append(T_COLUMN)

    x = raw[keep_cols].copy()
    rename = {
        detected_gene_column: "gene",
        LOGFC_COLUMN: "logFC",
        P_COLUMN: "p_value",
        FDR_COLUMN: "cohort_fdr",
    }
    if T_COLUMN in x.columns:
        rename[T_COLUMN] = "t_statistic"
    x = x.rename(columns=rename)

    x["gene"] = x["gene"].map(norm_gene)
    for c in ["logFC", "p_value", "cohort_fdr"]:
        x[c] = pd.to_numeric(x[c], errors="coerce")
    if "t_statistic" in x.columns:
        x["t_statistic"] = pd.to_numeric(x["t_statistic"], errors="coerce")

    n_input = len(x)
    x = x[
        x["gene"].ne("")
        & x["logFC"].notna()
        & x["p_value"].notna()
        & x["cohort_fdr"].notna()
    ].copy()

    x = x[
        np.isfinite(x["logFC"])
        & np.isfinite(x["p_value"])
        & np.isfinite(x["cohort_fdr"])
        & x["p_value"].between(0, 1, inclusive="both")
        & x["cohort_fdr"].between(0, 1, inclusive="both")
    ].copy()

    # Orient both cohorts to CAP - Control before signed-Z construction.
    if direction_multiplier not in {-1, +1}:
        raise ValueError("Direction multiplier must be +1 or -1.")
    x["logFC"] = x["logFC"] * direction_multiplier
    if "t_statistic" in x.columns:
        x["t_statistic"] = x["t_statistic"] * direction_multiplier

    # The current files should already be gene-level unique. Do not silently
    # collapse duplicates because that would change the disease universe.
    dup = x["gene"].duplicated(keep=False)
    if dup.any():
        examples = x.loc[dup, "gene"].drop_duplicates().head(20).tolist()
        raise RuntimeError(
            f"{cohort} contains duplicated gene symbols after normalization. "
            f"Examples: {examples}. Resolve upstream rather than silently "
            "choosing a row here."
        )

    x["signed_z"] = signed_z_from_two_sided_p(
        x["logFC"].to_numpy(), x["p_value"].to_numpy()
    )

    # Audit whether supplied P values and t statistics have consistent signs.
    # We do not reconstruct P from t because moderated t degrees of freedom are
    # not carried in the current gene-level tables.
    if "t_statistic" in x.columns:
        nonzero = (x["logFC"] != 0) & (x["t_statistic"] != 0)
        sign_match = (
            np.sign(x.loc[nonzero, "logFC"])
            == np.sign(x.loc[nonzero, "t_statistic"])
        )
        t_logfc_sign_agreement = float(sign_match.mean()) if len(sign_match) else np.nan
    else:
        t_logfc_sign_agreement = np.nan

    stats = {
        "cohort": cohort,
        "input_file": str(path),
        "input_sha256": sha256_file(path),
        "detected_gene_column": detected_gene_column,
        "rows_input": int(n_input),
        "valid_unique_genes": int(len(x)),
        "direction_multiplier": int(direction_multiplier),
        "min_p_value": float(x["p_value"].min()),
        "max_abs_signed_z": float(np.abs(x["signed_z"]).max()),
        "t_logFC_sign_agreement": t_logfc_sign_agreement,
    }

    suffix_map = {
        "logFC": f"logFC_{cohort}",
        "p_value": f"p_{cohort}",
        "cohort_fdr": f"fdr_{cohort}",
        "signed_z": f"Z_{cohort}",
        "t_statistic": f"t_{cohort}",
    }
    x = x.rename(columns=suffix_map)
    return x, stats


# =============================================================================
# 4. META-Z RECONSTRUCTION
# =============================================================================

def reconstruct_meta(
    g103: pd.DataFrame,
    g196: pd.DataFrame,
) -> tuple[pd.DataFrame, dict]:
    m = g103.merge(g196, on="gene", how="inner", validate="one_to_one")
    if m.empty:
        raise RuntimeError("No shared genes between the two CAP cohorts.")

    z1 = m["Z_GSE103119"].to_numpy(dtype=float)
    z2 = m["Z_GSE196399"].to_numpy(dtype=float)

    # PRIMARY: equal-weight Stouffer.
    m["Z_meta_equal"] = (z1 + z2) / math.sqrt(2.0)
    m["P_meta_equal"] = two_sided_p_from_z(m["Z_meta_equal"].to_numpy())
    m["FDR_meta_equal"] = bh_fdr(m["P_meta_equal"].to_numpy())

    # Direction consistency is evaluated on oriented CAP-vs-Control logFC.
    sign1 = np.sign(m["logFC_GSE103119"].to_numpy(dtype=float))
    sign2 = np.sign(m["logFC_GSE196399"].to_numpy(dtype=float))
    m["direction_consistent"] = (sign1 == sign2) & (sign1 != 0)
    m["cohort_direction_GSE103119"] = np.where(
        sign1 > 0, "UP", np.where(sign1 < 0, "DOWN", "ZERO")
    )
    m["cohort_direction_GSE196399"] = np.where(
        sign2 > 0, "UP", np.where(sign2 < 0, "DOWN", "ZERO")
    )
    zm = m["Z_meta_equal"].to_numpy(dtype=float)
    m["meta_direction_equal"] = np.where(
        zm > 0, "UP", np.where(zm < 0, "DOWN", "ZERO")
    )

    # Useful descriptive discordance magnitude; not a formal heterogeneity test.
    m["abs_Z_difference_between_cohorts"] = np.abs(z1 - z2)
    m["min_abs_cohort_Z"] = np.minimum(np.abs(z1), np.abs(z2))

    # PREDECLARED SENSITIVITY: effective-sample-size weighted Stouffer.
    n103 = SAMPLE_COUNTS["GSE103119"]
    n196 = SAMPLE_COUNTS["GSE196399"]
    neff103 = effective_two_group_n(n103["case"], n103["control"])
    neff196 = effective_two_group_n(n196["case"], n196["control"])
    w103 = math.sqrt(neff103)
    w196 = math.sqrt(neff196)
    denom = math.sqrt(w103**2 + w196**2)

    m["Z_meta_neff_weighted"] = (w103 * z1 + w196 * z2) / denom
    m["P_meta_neff_weighted"] = two_sided_p_from_z(
        m["Z_meta_neff_weighted"].to_numpy()
    )
    m["FDR_meta_neff_weighted"] = bh_fdr(
        m["P_meta_neff_weighted"].to_numpy()
    )

    # Descriptive direction-consistent meta-FDR module.
    m["in_direction_consistent_metaFDR05_module"] = (
        m["direction_consistent"]
        & (m["FDR_meta_equal"] < META_FDR_THRESHOLD)
    )

    # Strict replicated-DEG sensitivity / continuity set.
    m["in_strict_replicated_DEG_module"] = (
        m["direction_consistent"]
        & (m["fdr_GSE103119"] < STRICT_COHORT_FDR_THRESHOLD)
        & (m["fdr_GSE196399"] < STRICT_COHORT_FDR_THRESHOLD)
        & (m["logFC_GSE103119"].abs() >= STRICT_ABS_LOGFC_THRESHOLD)
        & (m["logFC_GSE196399"].abs() >= STRICT_ABS_LOGFC_THRESHOLD)
    )

    # Rank by disease-program magnitude, retaining direction.
    m["rank_abs_Z_meta_equal"] = (
        m["Z_meta_equal"].abs().rank(method="min", ascending=False).astype(int)
    )
    m["rank_Z_meta_equal_desc"] = (
        m["Z_meta_equal"].rank(method="min", ascending=False).astype(int)
    )
    m["rank_Z_meta_equal_asc"] = (
        m["Z_meta_equal"].rank(method="min", ascending=True).astype(int)
    )

    # Sort for human inspection only. Downstream joins must use gene identifiers.
    m = m.sort_values(
        ["FDR_meta_equal", "P_meta_equal", "gene"],
        ascending=[True, True, True],
    ).reset_index(drop=True)

    # Cross-cohort reproducibility summaries.
    rho_logfc, p_rho_logfc = safe_spearman(
        m["logFC_GSE103119"], m["logFC_GSE196399"]
    )
    pear_logfc, p_pear_logfc = safe_pearson(
        m["logFC_GSE103119"], m["logFC_GSE196399"]
    )
    rho_z, p_rho_z = safe_spearman(m["Z_GSE103119"], m["Z_GSE196399"])
    pear_z, p_pear_z = safe_pearson(m["Z_GSE103119"], m["Z_GSE196399"])
    rho_meta_sens, p_meta_sens = safe_spearman(
        m["Z_meta_equal"], m["Z_meta_neff_weighted"]
    )

    stats = {
        "shared_gene_universe": int(len(m)),
        "logFC_spearman": rho_logfc,
        "logFC_spearman_p": p_rho_logfc,
        "logFC_pearson": pear_logfc,
        "logFC_pearson_p": p_pear_logfc,
        "signedZ_spearman": rho_z,
        "signedZ_spearman_p": p_rho_z,
        "signedZ_pearson": pear_z,
        "signedZ_pearson_p": p_pear_z,
        "direction_concordance_fraction": float(m["direction_consistent"].mean()),
        "metaFDR05_total": int((m["FDR_meta_equal"] < META_FDR_THRESHOLD).sum()),
        "direction_consistent_metaFDR05_total": int(
            m["in_direction_consistent_metaFDR05_module"].sum()
        ),
        "direction_consistent_metaFDR05_up": int(
            (
                m["in_direction_consistent_metaFDR05_module"]
                & (m["Z_meta_equal"] > 0)
            ).sum()
        ),
        "direction_consistent_metaFDR05_down": int(
            (
                m["in_direction_consistent_metaFDR05_module"]
                & (m["Z_meta_equal"] < 0)
            ).sum()
        ),
        "strict_replicated_DEG_total": int(m["in_strict_replicated_DEG_module"].sum()),
        "strict_replicated_DEG_up": int(
            (
                m["in_strict_replicated_DEG_module"]
                & (m["Z_meta_equal"] > 0)
            ).sum()
        ),
        "strict_replicated_DEG_down": int(
            (
                m["in_strict_replicated_DEG_module"]
                & (m["Z_meta_equal"] < 0)
            ).sum()
        ),
        "effective_n_GSE103119": neff103,
        "effective_n_GSE196399": neff196,
        "weighted_stouffer_weight_GSE103119": w103,
        "weighted_stouffer_weight_GSE196399": w196,
        "equal_vs_neff_weighted_metaZ_spearman": rho_meta_sens,
        "equal_vs_neff_weighted_metaZ_spearman_p": p_meta_sens,
    }

    return m, stats


# =============================================================================
# 5. OUTPUT TABLES
# =============================================================================

def build_continuous_primary(meta: pd.DataFrame) -> pd.DataFrame:
    """
    Minimal disease vector for downstream syndrome-disease alignment.

    IMPORTANT: this table contains ALL shared genes, not only significant genes.
    """
    cols = [
        "gene",
        "Z_meta_equal",
        "P_meta_equal",
        "FDR_meta_equal",
        "meta_direction_equal",
        "direction_consistent",
        "Z_GSE103119",
        "Z_GSE196399",
        "logFC_GSE103119",
        "logFC_GSE196399",
        "Z_meta_neff_weighted",
        "FDR_meta_neff_weighted",
    ]
    return meta[cols].copy()


def build_module_table(meta: pd.DataFrame, flag: str) -> pd.DataFrame:
    x = meta[meta[flag]].copy()
    return x.sort_values(
        ["Z_meta_equal"], ascending=False
    ).reset_index(drop=True)


def build_top_gene_table(meta: pd.DataFrame, n: int = 200) -> pd.DataFrame:
    up = meta.nlargest(n, "Z_meta_equal").copy()
    up["top_set"] = f"top_{n}_CAP_UP_by_metaZ"
    down = meta.nsmallest(n, "Z_meta_equal").copy()
    down["top_set"] = f"top_{n}_CAP_DOWN_by_metaZ"
    return pd.concat([up, down], ignore_index=True)


# =============================================================================
# 6. MAIN
# =============================================================================

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("CAP TRANSCRIPTOMIC META-Z RECONSTRUCTION v1.0.1")
    print("Primary disease vector: equal-weight signed-Z Stouffer meta-analysis")
    print("=" * 80)

    f103 = resolve_input(
        GSE103119_FILE, GSE103119_FALLBACKS, "GSE103119 gene-level all-results"
    )
    f196 = resolve_input(
        GSE196399_FILE, GSE196399_FALLBACKS, "GSE196399 gene-level all-results"
    )

    print(f"\n[1/7] GSE103119 input:\n  {f103}")
    print(f"[2/7] GSE196399 input:\n  {f196}")

    g103, qa103 = load_cohort(
        f103, "GSE103119", GSE103119_DIRECTION_MULTIPLIER
    )
    g196, qa196 = load_cohort(
        f196, "GSE196399", GSE196399_DIRECTION_MULTIPLIER
    )

    print("\n[3/7] Cohort-level QA")
    print(pd.DataFrame([qa103, qa196]).to_string(index=False))

    print("\n[4/7] Reconstructing shared-gene signed-Z meta-analysis...")
    meta, meta_stats = reconstruct_meta(g103, g196)

    print(f"    shared genes: {meta_stats['shared_gene_universe']:,}")
    print(f"    logFC Spearman: {meta_stats['logFC_spearman']:.6f}")
    print(f"    signed-Z Spearman: {meta_stats['signedZ_spearman']:.6f}")
    print(
        f"    direction concordance: "
        f"{meta_stats['direction_concordance_fraction']:.3%}"
    )
    print(
        "    direction-consistent meta-FDR<0.05 module: "
        f"{meta_stats['direction_consistent_metaFDR05_total']:,} "
        f"(UP {meta_stats['direction_consistent_metaFDR05_up']:,}; "
        f"DOWN {meta_stats['direction_consistent_metaFDR05_down']:,})"
    )
    print(
        "    strict replicated DEG sensitivity module: "
        f"{meta_stats['strict_replicated_DEG_total']:,} "
        f"(UP {meta_stats['strict_replicated_DEG_up']:,}; "
        f"DOWN {meta_stats['strict_replicated_DEG_down']:,})"
    )
    print(
        "    equal-vs-neff-weighted meta-Z Spearman: "
        f"{meta_stats['equal_vs_neff_weighted_metaZ_spearman']:.6f}"
    )

    print("\n[5/7] Building downstream-ready disease vectors/modules...")
    continuous = build_continuous_primary(meta)
    module_meta = build_module_table(
        meta, "in_direction_consistent_metaFDR05_module"
    )
    module_strict = build_module_table(
        meta, "in_strict_replicated_DEG_module"
    )
    top200 = build_top_gene_table(meta, n=200)

    print("\n[6/7] Writing outputs...")
    meta.to_csv(
        OUTPUT_DIR / "CAP_metaZ_all_shared_genes_v1_0_1.csv",
        index=False, encoding="utf-8-sig"
    )
    continuous.to_csv(
        OUTPUT_DIR / "CAP_CONTINUOUS_PRIMARY_equal_weight_metaZ_v1_0_1.csv",
        index=False, encoding="utf-8-sig"
    )
    module_meta.to_csv(
        OUTPUT_DIR / "CAP_direction_consistent_metaFDR05_module_v1_0_1.csv",
        index=False, encoding="utf-8-sig"
    )
    module_strict.to_csv(
        OUTPUT_DIR / "CAP_strict_replicated_DEG_sensitivity_module_v1_0_1.csv",
        index=False, encoding="utf-8-sig"
    )
    top200.to_csv(
        OUTPUT_DIR / "CAP_top200_up_down_by_equal_metaZ_v1_0_1.csv",
        index=False, encoding="utf-8-sig"
    )

    cohort_qa = pd.DataFrame([qa103, qa196])
    cohort_qa.to_csv(
        OUTPUT_DIR / "CAP_metaZ_cohort_input_QA_v1_0_1.csv",
        index=False, encoding="utf-8-sig"
    )
    pd.DataFrame([meta_stats]).to_csv(
        OUTPUT_DIR / "CAP_metaZ_reproducibility_summary_v1_0_1.csv",
        index=False, encoding="utf-8-sig"
    )

    # Compact rank-stability table for the predeclared weighted sensitivity.
    sensitivity = meta[[
        "gene",
        "Z_meta_equal",
        "Z_meta_neff_weighted",
        "FDR_meta_equal",
        "FDR_meta_neff_weighted",
        "direction_consistent",
    ]].copy()
    sensitivity["rank_equal_abs"] = (
        sensitivity["Z_meta_equal"].abs()
        .rank(method="min", ascending=False).astype(int)
    )
    sensitivity["rank_neff_weighted_abs"] = (
        sensitivity["Z_meta_neff_weighted"].abs()
        .rank(method="min", ascending=False).astype(int)
    )
    sensitivity["absolute_rank_shift"] = (
        sensitivity["rank_equal_abs"]
        - sensitivity["rank_neff_weighted_abs"]
    ).abs()
    sensitivity.to_csv(
        OUTPUT_DIR / "CAP_equal_vs_neff_weighted_metaZ_sensitivity_v1_0_1.csv",
        index=False, encoding="utf-8-sig"
    )

    print("\n[7/7] Writing metadata and freeze candidate statement...")
    metadata = {
        "analysis_name": "CAP transcriptomic signed-Z meta-analysis v1.0.1",
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "status": "DISEASE_LAYER_FREEZE_CANDIDATE_REQUIRES_OUTPUT_AUDIT",
        "disease_context": "community-acquired pneumonia",
        "cohorts": {
            "GSE103119": qa103,
            "GSE196399": qa196,
        },
        "contrast_orientation": "CAP minus healthy control in both cohorts",
        "primary_gene_universe": (
            "intersection of valid unique gene symbols in the two CURRENT "
            "gene-level all-results files"
        ),
        "primary_meta_method": {
            "cohort_statistic": (
                "signed Z = sign(logFC) * Phi^-1(1 - two-sided P/2)"
            ),
            "combination": "equal-weight Stouffer",
            "formula": "Z_meta=(Z_GSE103119+Z_GSE196399)/sqrt(2)",
            "multiple_testing": "Benjamini-Hochberg over shared gene universe",
            "continuous_vector": "ALL shared genes; no FDR filtering",
        },
        "direction_consistent_module": {
            "meta_FDR_threshold": META_FDR_THRESHOLD,
            "requires_same_logFC_sign": True,
            "role": "descriptive/set-based module; not primary continuous alignment vector",
        },
        "weighted_sensitivity": {
            "method": "effective-two-group-sample-size weighted Stouffer",
            "sample_counts": SAMPLE_COUNTS,
            "effective_n_formula": "n_case*n_control/(n_case+n_control)",
            "role": "predeclared sensitivity only",
        },
        "strict_replicated_DEG_sensitivity": {
            "both_cohort_FDR_below": STRICT_COHORT_FDR_THRESHOLD,
            "both_cohort_abs_logFC_at_least": STRICT_ABS_LOGFC_THRESHOLD,
            "requires_same_direction": True,
            "role": "legacy continuity / sensitivity only",
        },
        "summary": meta_stats,
        "important_interpretation": (
            "Positive meta-Z means higher expression in CAP than healthy control; "
            "negative meta-Z means lower expression. This is disease direction, "
            "not therapeutic reversal or efficacy."
        ),
        "downstream_freeze_policy": (
            "Do not alter the meta-Z method, cohort direction, FDR threshold, "
            "or gene universe based on syndrome-target alignment results. "
            "Audit this output first, then freeze the disease layer before "
            "constructing observed-vs-null syndrome-disease alignment tests."
        ),
    }

    with open(
        OUTPUT_DIR / "CAP_transcriptomic_metaZ_run_metadata_v1_0_1.json",
        "w", encoding="utf-8"
    ) as f:
        json.dump(metadata, f, ensure_ascii=False, indent=2)

    statement = f"""CAP TRANSCRIPTOMIC META-Z v1.0.1 — FREEZE CANDIDATE

PRIMARY DISEASE REPRESENTATION
------------------------------
Cohorts: GSE103119 and GSE196399
Contrast: CAP - healthy control in both cohorts
Shared-gene universe: {meta_stats['shared_gene_universe']}

Per-cohort statistic:
Z_i = sign(logFC_i) * Phi^(-1)(1 - P_i/2)

Primary equal-weight Stouffer statistic:
Z_meta = (Z_GSE103119 + Z_GSE196399) / sqrt(2)

Primary downstream disease vector:
CAP_CONTINUOUS_PRIMARY_equal_weight_metaZ_v1_0_1.csv
Use ALL shared genes; do not prefilter the continuous vector by FDR.

Direction-consistent descriptive module:
meta FDR < {META_FDR_THRESHOLD} AND same logFC sign in both cohorts
N = {meta_stats['direction_consistent_metaFDR05_total']}
UP = {meta_stats['direction_consistent_metaFDR05_up']}
DOWN = {meta_stats['direction_consistent_metaFDR05_down']}

Strict replicated-DEG sensitivity module:
both cohort FDR < {STRICT_COHORT_FDR_THRESHOLD}
both |logFC| >= {STRICT_ABS_LOGFC_THRESHOLD}
same direction
N = {meta_stats['strict_replicated_DEG_total']}
UP = {meta_stats['strict_replicated_DEG_up']}
DOWN = {meta_stats['strict_replicated_DEG_down']}

Cross-cohort reproducibility:
logFC Spearman = {meta_stats['logFC_spearman']:.6f}
signed-Z Spearman = {meta_stats['signedZ_spearman']:.6f}
direction concordance = {meta_stats['direction_concordance_fraction']:.6f}

effective-N weighted sensitivity vs equal-weight meta-Z Spearman =
{meta_stats['equal_vs_neff_weighted_metaZ_spearman']:.6f}

STATUS
------
FREEZE CANDIDATE ONLY. Audit current gene-universe size, cross-cohort
reproducibility, direction-consistent module size, and sensitivity stability.
Do not use syndrome alignment results to change this disease-layer method.
"""
    (OUTPUT_DIR / "CAP_META_Z_FREEZE_CANDIDATE_v1_0_1.txt").write_text(
        statement, encoding="utf-8"
    )

    print("\n" + "=" * 80)
    print("CAP meta-Z reconstruction completed.")
    print("STATUS: FREEZE CANDIDATE — audit before downstream alignment/null testing.")
    print("=" * 80)
    print(f"\nPrimary downstream file:\n{OUTPUT_DIR / 'CAP_CONTINUOUS_PRIMARY_equal_weight_metaZ_v1_0_1.csv'}")
    print(f"\nAll outputs:\n{OUTPUT_DIR}")


if __name__ == "__main__":
    main()
