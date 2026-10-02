# -*- coding: utf-8 -*-
"""
weighted_herb_target_projection_v1_1_FROZEN.py

Project frozen syndrome-herb weights to molecular targets WITHOUT using CAP outcomes.

Primary purpose
---------------
Quantify and reduce target-space collapse caused by:
1) herbs with very large target lists (herb promiscuity), and
2) targets linked to many herbs in the target database (target ubiquity).

Inputs
------
1) Frozen syndrome-herb weights:
   - FULL_syndrome_herb_weights_observed_v1_0.csv
   - LOPO_syndrome_herb_weights_observed_v1_0.csv

2) A herb-target long table (CSV/XLSX/TSV), one herb-target association per row.
   The script tries to auto-detect common herb and target column names.
   You can force the names in USER SETTINGS.

Frozen projection definitions
-----------------------------
BINARY_UNION:
    target present if any observed syndrome herb maps to target.
    Diagnostic baseline only.

WEIGHTED_RAW:
    T_ts = sum_h W_hs * A_ht
    Diagnostic baseline; highly sensitive to herb target degree.

HERB_DEGREE_NORM:
    T_ts = sum_h W_hs * A_ht / |T_h|
    Each mapped herb distributes its syndrome weight across its targets.

UBIQUITY_CORRECTED (primary candidate for downstream RWR):
    IDF_t = log((N_H + 1)/(n_t + 1)) + 1
    T_ts  = sum_h W_hs * (A_ht / |T_h|) * IDF_t

where:
    W_hs = frozen article-balanced syndrome-herb weight
    A_ht = binary herb-target association
    |T_h| = number of unique targets linked to herb h
    N_H = number of unique mapped herbs in the target database after canonical mapping
    n_t = number of unique herbs linked to target t

No CAP transcriptomic data are read or used.
"""

from pathlib import Path
import json
import hashlib
import re
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


# =============================================================================
# 1. USER SETTINGS — EDIT FILE PATHS HERE
# =============================================================================

FROZEN_WEIGHTS_DIR = Path(
    r"E:\00project\Tanre Yongfei\重构\04_frozen_syndrome_herb_weights_v1_0"
)

FULL_WEIGHTS_FILE = (
    FROZEN_WEIGHTS_DIR / "FULL_syndrome_herb_weights_observed_v1_0.csv"
)

LOPO_WEIGHTS_FILE = (
    FROZEN_WEIGHTS_DIR / "LOPO_syndrome_herb_weights_observed_v1_0.csv"
)

# IMPORTANT:
# Point this to your CURRENT herb-target long table.
# CSV, TSV/TXT, and XLSX are supported.
HERB_TARGET_FILE = Path(
    r"E:\00project\Tanre Yongfei\重构\05B_TCMID_web_completion_v1_0\HERB_TARGET_LONG_TABLE_TCMID_web_completed_v1_0.xlsx"
)

# For XLSX only. Set to None to use the first sheet.
HERB_TARGET_SHEET = "targets_long"

OUTPUT_DIR = Path(
    r"E:\00project\Tanre Yongfei\重构\05_weighted_herb_target_projection_v1_1_FROZEN"
)

# -----------------------------------------------------------------------------
# Column names
# -----------------------------------------------------------------------------
# Leave as None for auto-detection.
# If auto-detection fails, set explicitly, e.g.:
# HERB_COLUMN = "canonical_name"
# TARGET_COLUMN = "Gene"
HERB_COLUMN = "Herb_ChineseName"
TARGET_COLUMN = "Gene_Symbol"

# Optional evidence/source column. This script does NOT automatically filter
# database evidence because schemas differ. If you have already curated the
# target table, leave this None.
TARGET_EVIDENCE_COLUMN = None

# Optional value filter, used only if TARGET_EVIDENCE_COLUMN is set.
# Example: {"experiment", "literature"}
TARGET_EVIDENCE_ALLOWED = None

# -----------------------------------------------------------------------------
# Herb-name normalization
# -----------------------------------------------------------------------------
# Use only deterministic, pre-declared mappings. Never infer names from results.
# Add aliases here only when you have verified them.
HERB_ALIAS_MAP = {
    "贝母": "川贝母",
    "蛤壳": "海蛤壳",
    "茅葙子": "青葙子",
    "紫苑": "紫菀",
    "蜜紫苑": "蜜紫菀",
    "浙贝": "浙贝母",
    "川贝": "川贝母",
    "瓜萎": "瓜蒌",
    "瓜萎皮": "瓜蒌皮",
    "瓜萎子": "瓜蒌子",
    "全瓜萎": "瓜蒌",
    "全瓜蒌": "瓜蒌",
    "瓜蒌仁": "瓜蒌子",
    "瓜萎仁": "瓜蒌子",
    "蝉衣": "蝉蜕",
    "苏子": "紫苏子",
    "冬花": "款冬花",
    "旋复花": "旋覆花",
    "山萸肉": "山茱萸",
    "熟地": "熟地黄",
    "黑脂麻": "黑芝麻",
    "山栀": "栀子",
    "山栀子": "栀子",
    "天花": "天花粉",
    "杭菊花": "菊花",
    "北柴胡": "柴胡",
    "川朴": "厚朴",
    "川厚朴": "厚朴",
    "白蔻仁": "豆蔻",
    "白豆蔻": "豆蔻",
    "生石膏": "石膏",
    "生甘草": "甘草",
    "生黄芪": "黄芪",
    "炙黄芪": "黄芪",
    "生大黄": "大黄",
    "生杜仲": "杜仲",
    "生牡蛎": "牡蛎",
    "生薏苡仁": "薏苡仁",
    "生麻黄": "麻黄",
    "茅苍术": "苍术",
    "青蒿草": "青蒿",
    "干芦根": "芦根",
    "冬瓜仁": "冬瓜子",
    "麦门冬": "麦冬",
    "炙桑皮": "桑白皮",
    "炙桑白皮": "桑白皮",
    "蜜桑白皮": "桑白皮",
    "姜竹茹": "竹茹",
    "夏枯全草": "夏枯草",
}

# -----------------------------------------------------------------------------
# Frozen target-layer processed-form -> parent mappings
# -----------------------------------------------------------------------------
# IMPORTANT: These mappings are used ONLY when matching a frozen prescription
# herb to an existing molecular target profile. The prescription-layer
# canonical_name is preserved in all weights, contributions, and QA outputs.
# They were resolved upstream before CAP analysis and MUST NOT be expanded based
# on downstream disease-alignment results.
PROCESSED_FORM_PARENT_MAP = {
    "制陈皮": "陈皮",
    "炙百部": "百部",
    "焦山楂": "山楂",
}

# Freeze identity. Any change to mapping rules, target database, IDF reference
# universe, or projection formula requires a new version rather than overwriting
# v1.1.
FREEZE_VERSION = "v1.1"
FREEZE_STATUS = "FROZEN_UPSTREAM_BEFORE_CAP_ALIGNMENT"
IDF_REFERENCE_UNIVERSE = "study_target_database_herbs_with_at_least_one_valid_target"


# -----------------------------------------------------------------------------
# Analysis settings
# -----------------------------------------------------------------------------
SYNDROME_CODES = ["PDOL", "PHOL", "WHIL"]
REPRESENTATIONS = ["FULL", "LOPO"]

# Target ID handling.
# For gene symbols, True is usually appropriate.
UPPERCASE_TARGET_IDS = True

# Minimum target degree is deliberately 1: do not tune this using CAP results.
MIN_TARGETS_PER_HERB = 1

# Optional hard cap for diagnostic sensitivity only.
# Primary projection does NOT remove promiscuous herbs by degree.
MAX_TARGETS_PER_HERB = None

# Numerical normalization:
# In addition to raw projection scores, each syndrome vector is normalized to sum 1.
NORMALIZE_TARGET_VECTOR_TO_SUM1 = True

# Pairwise top-k overlap diagnostics.
TOP_K_VALUES = [50, 100, 200]

# Save a human-readable XLSX QA summary.
WRITE_EXCEL_SUMMARY = True


# =============================================================================
# 2. HELPERS
# =============================================================================

def norm_text(x):
    if pd.isna(x):
        return ""
    return str(x).strip()


def norm_herb(x):
    """Prescription/database name normalization: aliases only, no parent collapse."""
    s = norm_text(x)
    s = re.sub(r"\s+", "", s)
    return HERB_ALIAS_MAP.get(s, s)


def target_mapping_name(canonical_name, available_target_herbs):
    """
    Resolve a frozen prescription herb to a target-database herb.

    Priority:
      1) exact/verified-alias canonical name already present in target DB;
      2) frozen processed-form -> parent mapping, only if parent exists in DB;
      3) unresolved.

    Returns (mapping_name, mapping_method).
    """
    h = norm_herb(canonical_name)
    if h in available_target_herbs:
        return h, "exact_or_verified_alias"

    parent = PROCESSED_FORM_PARENT_MAP.get(h, "")
    parent = norm_herb(parent) if parent else ""
    if parent and parent in available_target_herbs:
        return parent, "processed_form_to_parent"

    return "", "unmapped"


def norm_target(x):
    s = norm_text(x)
    s = re.sub(r"\s+", "", s)
    if UPPERCASE_TARGET_IDS:
        s = s.upper()
    return s


def sha256_file(path: Path, chunk_size=1024 * 1024):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk_size)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def read_table(path: Path, sheet_name=None):
    if not path.exists():
        raise FileNotFoundError(f"Input file not found:\n{path}")

    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xlsm", ".xls"}:
        return pd.read_excel(path, sheet_name=0 if sheet_name is None else sheet_name)
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")
    raise ValueError(f"Unsupported input format: {suffix}")


def detect_column(df, explicit, candidates, kind):
    if explicit is not None:
        if explicit not in df.columns:
            raise ValueError(
                f"{kind} column '{explicit}' not found.\n"
                f"Available columns: {list(df.columns)}"
            )
        return explicit

    lower_map = {str(c).strip().lower(): c for c in df.columns}
    for c in candidates:
        if c.lower() in lower_map:
            return lower_map[c.lower()]

    # relaxed contains matching
    for original in df.columns:
        low = str(original).strip().lower()
        if any(c.lower() in low for c in candidates):
            return original

    raise ValueError(
        f"Could not auto-detect {kind} column.\n"
        f"Available columns: {list(df.columns)}\n"
        f"Please set the column explicitly in USER SETTINGS."
    )


def require_columns(df, cols, name):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"{name} missing required columns: {missing}\n"
            f"Available: {list(df.columns)}"
        )


def load_weights(path: Path, representation: str):
    df = pd.read_csv(path)
    require_columns(
        df,
        ["canonical_name", "syndrome_code", "primary_weight_W", "article_support_m"],
        path.name,
    )
    df = df.copy()
    df["canonical_name"] = df["canonical_name"].map(norm_herb)
    df["syndrome_code"] = df["syndrome_code"].map(norm_text)
    df["primary_weight_W"] = pd.to_numeric(df["primary_weight_W"], errors="coerce")
    df["article_support_m"] = pd.to_numeric(df["article_support_m"], errors="coerce")

    df = df[
        df["syndrome_code"].isin(SYNDROME_CODES)
        & df["canonical_name"].ne("")
        & df["article_support_m"].gt(0)
        & df["primary_weight_W"].notna()
    ].copy()

    # After alias normalization, merge any deterministic duplicates conservatively.
    # Weight should be identical only if the canonical source already agrees; max avoids
    # double-counting an alias as two herbs.
    agg = (
        df.groupby(["canonical_name", "syndrome_code"], as_index=False)
        .agg(
            primary_weight_W=("primary_weight_W", "max"),
            article_support_m=("article_support_m", "max"),
        )
    )
    agg["representation"] = representation
    return agg


def load_herb_target():
    raw = read_table(HERB_TARGET_FILE, HERB_TARGET_SHEET)

    herb_col = detect_column(
        raw,
        HERB_COLUMN,
        [
            "canonical_name", "herb", "herb_name", "herbname",
            "Chinese_name", "Chinese name", "中药", "中药名", "药材", "药材名"
        ],
        "herb",
    )
    target_col = detect_column(
        raw,
        TARGET_COLUMN,
        [
            "target", "target_gene", "targetgene", "gene", "gene_symbol",
            "genesymbol", "symbol", "target_name", "靶点", "靶基因"
        ],
        "target",
    )

    x = raw.copy()

    if TARGET_EVIDENCE_COLUMN is not None:
        if TARGET_EVIDENCE_COLUMN not in x.columns:
            raise ValueError(
                f"TARGET_EVIDENCE_COLUMN '{TARGET_EVIDENCE_COLUMN}' not found."
            )
        if TARGET_EVIDENCE_ALLOWED is not None:
            allowed = {str(v).strip() for v in TARGET_EVIDENCE_ALLOWED}
            x = x[
                x[TARGET_EVIDENCE_COLUMN].astype(str).str.strip().isin(allowed)
            ].copy()

    x = x[[herb_col, target_col]].copy()
    x.columns = ["herb_raw", "target_raw"]
    x["canonical_name"] = x["herb_raw"].map(norm_herb)
    x["target_id"] = x["target_raw"].map(norm_target)

    x = x[
        x["canonical_name"].ne("")
        & x["target_id"].ne("")
    ].copy()

    # Binary association: duplicate database rows must not multiply evidence.
    x = x.drop_duplicates(["canonical_name", "target_id"]).reset_index(drop=True)

    herb_degree = x.groupby("canonical_name")["target_id"].nunique()
    keep = herb_degree[herb_degree >= MIN_TARGETS_PER_HERB].index

    if MAX_TARGETS_PER_HERB is not None:
        keep = herb_degree[
            (herb_degree >= MIN_TARGETS_PER_HERB)
            & (herb_degree <= MAX_TARGETS_PER_HERB)
        ].index

    x = x[x["canonical_name"].isin(keep)].copy()

    return x, herb_col, target_col


def build_target_database_stats(ht: pd.DataFrame):
    herb_degree = (
        ht.groupby("canonical_name")["target_id"]
        .nunique()
        .rename("herb_target_degree")
    )

    target_degree = (
        ht.groupby("target_id")["canonical_name"]
        .nunique()
        .rename("target_herb_degree")
    )

    n_herbs = int(ht["canonical_name"].nunique())

    target_stats = target_degree.reset_index()
    target_stats["N_database_herbs"] = n_herbs
    target_stats["target_idf"] = np.log(
        (n_herbs + 1.0) / (target_stats["target_herb_degree"] + 1.0)
    ) + 1.0

    herb_stats = herb_degree.reset_index()

    return herb_stats, target_stats


def build_weight_target_mapping(weights: pd.DataFrame, ht: pd.DataFrame, representation: str):
    """One row per observed frozen herb-syndrome relation with target-layer mapping."""
    available = set(ht["canonical_name"].unique())
    x = weights.copy()
    resolved = x["canonical_name"].map(lambda h: target_mapping_name(h, available))
    x["target_mapping_name"] = [r[0] for r in resolved]
    x["mapping_method"] = [r[1] for r in resolved]
    x["mapped_to_target_db"] = x["target_mapping_name"].ne("").astype(int)
    x["representation"] = representation
    return x


def mapping_qa(weights: pd.DataFrame, ht: pd.DataFrame, representation: str):
    mapped = build_weight_target_mapping(weights, ht, representation)
    rows = []

    for syndrome in SYNDROME_CODES:
        w = mapped[mapped["syndrome_code"].eq(syndrome)].copy()
        total_n = int(w["canonical_name"].nunique())
        mapped_n = int(w.loc[w["mapped_to_target_db"].eq(1), "canonical_name"].nunique())
        total_weight = float(w["primary_weight_W"].sum())
        mapped_weight = float(w.loc[w["mapped_to_target_db"].eq(1), "primary_weight_W"].sum())

        rows.append({
            "representation": representation,
            "syndrome_code": syndrome,
            "observed_herbs": total_n,
            "mapped_herbs": mapped_n,
            "unmapped_herbs": total_n - mapped_n,
            "herb_mapping_fraction": mapped_n / total_n if total_n else np.nan,
            "total_frozen_weight": total_weight,
            "mapped_frozen_weight": mapped_weight,
            "weight_mapping_fraction": mapped_weight / total_weight if total_weight else np.nan,
            "n_processed_form_to_parent": int(
                w.loc[w["mapping_method"].eq("processed_form_to_parent"), "canonical_name"].nunique()
            ),
        })

    return pd.DataFrame(rows)


def unmapped_table(weights: pd.DataFrame, ht: pd.DataFrame, representation: str):
    x = build_weight_target_mapping(weights, ht, representation)
    x = x[x["mapped_to_target_db"].eq(0)].copy()
    return x.sort_values(
        ["syndrome_code", "primary_weight_W", "canonical_name"],
        ascending=[True, False, True],
    )


def project_one(weights: pd.DataFrame,
                ht: pd.DataFrame,
                herb_stats: pd.DataFrame,
                target_stats: pd.DataFrame,
                representation: str):
    """
    Produce target-level scores for all four projection variants.
    """
    # Preserve prescription-layer canonical_name and resolve a separate target-layer
    # mapping name. Processed forms therefore inherit parent targets without being
    # collapsed in the frozen syndrome-herb representation.
    mapped_w = build_weight_target_mapping(weights, ht, representation)
    mapped_w = mapped_w[mapped_w["mapped_to_target_db"].eq(1)].copy()

    ht_for_join = ht.rename(columns={"canonical_name": "target_mapping_name"})
    herb_stats_for_join = herb_stats.rename(columns={"canonical_name": "target_mapping_name"})

    x = mapped_w.merge(ht_for_join, on="target_mapping_name", how="inner")
    x = x.merge(herb_stats_for_join, on="target_mapping_name", how="left")
    x = x.merge(target_stats[["target_id", "target_herb_degree", "target_idf"]],
                on="target_id", how="left")

    x["binary_contribution"] = 1.0
    x["weighted_raw_contribution"] = x["primary_weight_W"]
    x["degree_norm_contribution"] = (
        x["primary_weight_W"] / x["herb_target_degree"]
    )
    x["ubiquity_corrected_contribution"] = (
        x["degree_norm_contribution"] * x["target_idf"]
    )

    # Save herb-target contributions for audit.
    contribution_cols = [
        "representation", "syndrome_code", "canonical_name",
        "target_mapping_name", "mapping_method",
        "primary_weight_W", "article_support_m", "target_id",
        "herb_target_degree", "target_herb_degree", "target_idf",
        "weighted_raw_contribution",
        "degree_norm_contribution",
        "ubiquity_corrected_contribution",
    ]
    x["representation"] = representation
    contributions = x[contribution_cols].copy()

    # Aggregate target scores.
    grouped = (
        x.groupby(["syndrome_code", "target_id"], as_index=False)
        .agg(
            contributing_herbs=("canonical_name", "nunique"),
            weighted_raw=("weighted_raw_contribution", "sum"),
            herb_degree_normalized=("degree_norm_contribution", "sum"),
            ubiquity_corrected=("ubiquity_corrected_contribution", "sum"),
            target_herb_degree=("target_herb_degree", "first"),
            target_idf=("target_idf", "first"),
        )
    )
    grouped["binary_union"] = 1

    # Complete target universe across syndromes, fill absent targets with zero.
    all_targets = sorted(ht["target_id"].unique())
    grid = pd.MultiIndex.from_product(
        [SYNDROME_CODES, all_targets],
        names=["syndrome_code", "target_id"]
    ).to_frame(index=False)

    out = grid.merge(
        grouped,
        on=["syndrome_code", "target_id"],
        how="left"
    )

    for c in [
        "contributing_herbs", "weighted_raw", "herb_degree_normalized",
        "ubiquity_corrected", "binary_union"
    ]:
        out[c] = out[c].fillna(0)

    out = out.merge(
        target_stats[["target_id", "target_herb_degree", "target_idf"]],
        on="target_id",
        how="left",
        suffixes=("", "_db"),
    )

    # Resolve duplicated metadata from merge.
    for c in ["target_herb_degree", "target_idf"]:
        db = c + "_db"
        if db in out.columns:
            out[c] = out[c].fillna(out[db])
            out = out.drop(columns=[db])

    if NORMALIZE_TARGET_VECTOR_TO_SUM1:
        for score in ["weighted_raw", "herb_degree_normalized", "ubiquity_corrected"]:
            denom = out.groupby("syndrome_code")[score].transform("sum")
            out[score + "_sum1"] = np.where(
                denom > 0, out[score] / denom, 0.0
            )

    out["representation"] = representation

    return out, contributions


def pairwise_similarity(target_scores: pd.DataFrame, representation: str):
    score_cols = [
        "binary_union",
        "weighted_raw_sum1",
        "herb_degree_normalized_sum1",
        "ubiquity_corrected_sum1",
    ]
    score_cols = [c for c in score_cols if c in target_scores.columns]

    rows = []
    pairs = [("PDOL", "PHOL"), ("PDOL", "WHIL"), ("PHOL", "WHIL")]

    for metric in score_cols:
        pivot = target_scores.pivot(
            index="target_id",
            columns="syndrome_code",
            values=metric,
        ).fillna(0)

        for a, b in pairs:
            va = pivot[a].values
            vb = pivot[b].values

            rho, p = spearmanr(va, vb)
            if np.std(va) == 0 or np.std(vb) == 0:
                pearson = np.nan
            else:
                pearson = float(np.corrcoef(va, vb)[0, 1])

            positive_a = set(pivot.index[pivot[a] > 0])
            positive_b = set(pivot.index[pivot[b] > 0])
            union = positive_a | positive_b
            inter = positive_a & positive_b
            jaccard = len(inter) / len(union) if union else np.nan

            rows.append({
                "representation": representation,
                "metric": metric,
                "syndrome_A": a,
                "syndrome_B": b,
                "spearman": rho,
                "spearman_p": p,
                "pearson": pearson,
                "positive_target_overlap": len(inter),
                "positive_target_union": len(union),
                "positive_target_jaccard": jaccard,
            })

    return pd.DataFrame(rows)


def topk_overlap(target_scores: pd.DataFrame, representation: str):
    metrics = [
        "weighted_raw_sum1",
        "herb_degree_normalized_sum1",
        "ubiquity_corrected_sum1",
    ]
    metrics = [m for m in metrics if m in target_scores.columns]
    pairs = [("PDOL", "PHOL"), ("PDOL", "WHIL"), ("PHOL", "WHIL")]
    rows = []

    for metric in metrics:
        for k in TOP_K_VALUES:
            sets = {}
            for s in SYNDROME_CODES:
                x = target_scores[
                    target_scores["syndrome_code"].eq(s)
                ].nlargest(k, metric)
                sets[s] = set(x["target_id"])

            for a, b in pairs:
                inter = sets[a] & sets[b]
                union = sets[a] | sets[b]
                rows.append({
                    "representation": representation,
                    "metric": metric,
                    "top_k": k,
                    "syndrome_A": a,
                    "syndrome_B": b,
                    "overlap_n": len(inter),
                    "jaccard": len(inter) / len(union) if union else np.nan,
                })

    return pd.DataFrame(rows)


def syndrome_target_summary(target_scores: pd.DataFrame, representation: str):
    rows = []
    for s in SYNDROME_CODES:
        x = target_scores[target_scores["syndrome_code"].eq(s)]
        rows.append({
            "representation": representation,
            "syndrome_code": s,
            "binary_union_targets": int((x["binary_union"] > 0).sum()),
            "weighted_raw_nonzero_targets": int((x["weighted_raw"] > 0).sum()),
            "degree_norm_nonzero_targets": int((x["herb_degree_normalized"] > 0).sum()),
            "ubiquity_nonzero_targets": int((x["ubiquity_corrected"] > 0).sum()),
            "sum_weighted_raw": float(x["weighted_raw"].sum()),
            "sum_degree_norm": float(x["herb_degree_normalized"].sum()),
            "sum_ubiquity_corrected": float(x["ubiquity_corrected"].sum()),
        })
    return pd.DataFrame(rows)


def top_targets(target_scores, representation, n=100):
    x = target_scores.copy()
    x["rank_ubiquity"] = (
        x.groupby("syndrome_code")["ubiquity_corrected_sum1"]
        .rank(method="min", ascending=False)
    )
    x = x[
        (x["rank_ubiquity"] <= n)
        & (x["ubiquity_corrected_sum1"] > 0)
    ].copy()
    x["representation"] = representation
    return x.sort_values(
        ["syndrome_code", "rank_ubiquity", "target_id"]
    )


# =============================================================================
# 3. MAIN
# =============================================================================

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("[1/8] Loading frozen FULL and LOPO syndrome-herb weights...")
    full_w = load_weights(FULL_WEIGHTS_FILE, "FULL")
    lopo_w = load_weights(LOPO_WEIGHTS_FILE, "LOPO")

    print("[2/8] Loading and canonicalizing herb-target associations...")
    ht, detected_herb_col, detected_target_col = load_herb_target()

    print(f"    Detected herb column:   {detected_herb_col}")
    print(f"    Detected target column: {detected_target_col}")
    print(f"    Unique target-DB herbs: {ht['canonical_name'].nunique()}")
    print(f"    Unique targets:         {ht['target_id'].nunique()}")
    print(f"    Unique herb-target pairs: {len(ht)}")

    print("[3/8] Computing herb promiscuity and target ubiquity...")
    herb_stats, target_stats = build_target_database_stats(ht)

    herb_stats.to_csv(
        OUTPUT_DIR / "target_database_herb_degree_v1_1.csv",
        index=False, encoding="utf-8-sig"
    )
    target_stats.to_csv(
        OUTPUT_DIR / "target_database_target_ubiquity_v1_1.csv",
        index=False, encoding="utf-8-sig"
    )
    ht.to_csv(
        OUTPUT_DIR / "canonical_herb_target_pairs_v1_1.csv",
        index=False, encoding="utf-8-sig"
    )

    # Freeze the exact prescription-herb -> target-database-herb mapping used.
    mapping_dictionary = pd.concat([
        build_weight_target_mapping(full_w, ht, "FULL"),
        build_weight_target_mapping(lopo_w, ht, "LOPO"),
    ], ignore_index=True)[
        ["representation", "syndrome_code", "canonical_name",
         "target_mapping_name", "mapping_method",
         "primary_weight_W", "article_support_m", "mapped_to_target_db"]
    ].drop_duplicates()
    mapping_dictionary.to_csv(
        OUTPUT_DIR / "FROZEN_herb_to_target_mapping_dictionary_v1_1.csv",
        index=False, encoding="utf-8-sig"
    )

    all_mapping = []
    all_unmapped = []
    all_scores = []
    all_contrib = []
    all_similarity = []
    all_topk = []
    all_summary = []
    all_top_targets = []

    for representation, weights in [("FULL", full_w), ("LOPO", lopo_w)]:
        print(f"[4/8] Mapping QA for {representation}...")
        qa = mapping_qa(weights, ht, representation)
        unmapped = unmapped_table(weights, ht, representation)
        all_mapping.append(qa)
        all_unmapped.append(unmapped)

        print(f"[5/8] Projecting {representation} syndrome weights to targets...")
        scores, contributions = project_one(
            weights, ht, herb_stats, target_stats, representation
        )
        all_scores.append(scores)
        all_contrib.append(contributions)

        print(f"[6/8] Quantifying {representation} target-space collapse...")
        all_similarity.append(pairwise_similarity(scores, representation))
        all_topk.append(topk_overlap(scores, representation))
        all_summary.append(syndrome_target_summary(scores, representation))
        all_top_targets.append(top_targets(scores, representation, n=100))

    mapping_qa_df = pd.concat(all_mapping, ignore_index=True)
    unmapped_df = pd.concat(all_unmapped, ignore_index=True)
    scores_df = pd.concat(all_scores, ignore_index=True)
    contributions_df = pd.concat(all_contrib, ignore_index=True)
    similarity_df = pd.concat(all_similarity, ignore_index=True)
    topk_df = pd.concat(all_topk, ignore_index=True)
    summary_df = pd.concat(all_summary, ignore_index=True)
    top_targets_df = pd.concat(all_top_targets, ignore_index=True)

    print("[7/8] Writing projection outputs...")

    mapping_qa_df.to_csv(
        OUTPUT_DIR / "herb_target_mapping_QA_v1_1.csv",
        index=False, encoding="utf-8-sig"
    )
    unmapped_df.to_csv(
        OUTPUT_DIR / "unmapped_frozen_herbs_v1_1.csv",
        index=False, encoding="utf-8-sig"
    )
    scores_df.to_csv(
        OUTPUT_DIR / "syndrome_target_projection_all_methods_v1_1.csv",
        index=False, encoding="utf-8-sig"
    )
    contributions_df.to_csv(
        OUTPUT_DIR / "herb_target_weight_contributions_v1_1.csv",
        index=False, encoding="utf-8-sig"
    )
    similarity_df.to_csv(
        OUTPUT_DIR / "target_space_pairwise_similarity_v1_1.csv",
        index=False, encoding="utf-8-sig"
    )
    topk_df.to_csv(
        OUTPUT_DIR / "target_space_topk_overlap_v1_1.csv",
        index=False, encoding="utf-8-sig"
    )
    summary_df.to_csv(
        OUTPUT_DIR / "syndrome_target_projection_summary_v1_1.csv",
        index=False, encoding="utf-8-sig"
    )
    top_targets_df.to_csv(
        OUTPUT_DIR / "top100_ubiquity_corrected_targets_v1_1.csv",
        index=False, encoding="utf-8-sig"
    )

    # Separate downstream-ready LOPO primary candidate.
    lopo_primary = scores_df[
        scores_df["representation"].eq("LOPO")
        & scores_df["ubiquity_corrected_sum1"].gt(0)
    ][
        [
            "syndrome_code", "target_id",
            "ubiquity_corrected", "ubiquity_corrected_sum1",
            "contributing_herbs", "target_herb_degree", "target_idf"
        ]
    ].copy()

    lopo_primary.to_csv(
        OUTPUT_DIR / "LOPO_FROZEN_primary_target_weights_for_RWR_v1_1.csv",
        index=False, encoding="utf-8-sig"
    )

    # Hard freeze checks for the downstream seed vector.
    frozen_sums = lopo_primary.groupby("syndrome_code")["ubiquity_corrected_sum1"].sum()
    for syndrome in SYNDROME_CODES:
        if syndrome not in frozen_sums.index:
            raise RuntimeError(f"Missing frozen LOPO seed vector for {syndrome}")
        if abs(float(frozen_sums.loc[syndrome]) - 1.0) > 1e-10:
            raise RuntimeError(
                f"Frozen LOPO seed vector for {syndrome} does not sum to 1: "
                f"{frozen_sums.loc[syndrome]}"
            )

    freeze_statement = f"""
FROZEN TARGET PROJECTION {FREEZE_VERSION}
Status: {FREEZE_STATUS}

Primary downstream representation:
  LOPO ubiquity_corrected_sum1

Frozen formula:
  W_hs = article-balanced LOPO syndrome-herb evidence weight
  degree-normalized contribution = W_hs / |T_h|
  IDF_t = log((N_H + 1)/(n_t + 1)) + 1
  target score = sum_h [W_hs / |T_h| * IDF_t]
  final syndrome target vector is normalized to sum 1

IDF reference universe:
  {IDF_REFERENCE_UNIVERSE}
  N_H is computed from unique herbs with >=1 valid target in the completed
  study target table; processed-form aliases do not create additional IDF herbs.

Frozen target-layer processed-form mappings:
  制陈皮 -> 陈皮
  炙百部 -> 百部
  焦山楂 -> 山楂

Prescription-layer herb identities remain unchanged.
No CAP transcriptomic outcome was used to select or tune this projection.
Any future change requires a new projection version.
""".strip() + "\n"

    (OUTPUT_DIR / "TARGET_PROJECTION_FREEZE_STATEMENT_v1_1.txt").write_text(
        freeze_statement, encoding="utf-8"
    )

    if WRITE_EXCEL_SUMMARY:
        with pd.ExcelWriter(
            OUTPUT_DIR / "weighted_herb_target_projection_v1_1_summary.xlsx",
            engine="openpyxl"
        ) as writer:
            mapping_qa_df.to_excel(writer, sheet_name="mapping_QA", index=False)
            summary_df.to_excel(writer, sheet_name="projection_summary", index=False)
            similarity_df.to_excel(writer, sheet_name="pairwise_similarity", index=False)
            topk_df.to_excel(writer, sheet_name="topk_overlap", index=False)
            unmapped_df.to_excel(writer, sheet_name="unmapped_herbs", index=False)
            top_targets_df.to_excel(writer, sheet_name="top100_targets", index=False)
            herb_stats.sort_values(
                "herb_target_degree", ascending=False
            ).to_excel(writer, sheet_name="herb_degree", index=False)
            target_stats.sort_values(
                "target_herb_degree", ascending=False
            ).to_excel(writer, sheet_name="target_ubiquity", index=False)

    print("[8/8] Writing frozen metadata and safety checks...")

    # Basic mass-conservation QA for degree normalization:
    # before target IDF, each mapped herb should distribute exactly W across targets.
    qa_mass = (
        contributions_df.groupby(
            ["representation", "syndrome_code", "canonical_name"],
            as_index=False
        )
        .agg(
            primary_weight_W=("primary_weight_W", "first"),
            distributed_degree_norm_mass=("degree_norm_contribution", "sum"),
        )
    )
    qa_mass["absolute_error"] = (
        qa_mass["primary_weight_W"] - qa_mass["distributed_degree_norm_mass"]
    ).abs()
    qa_mass.to_csv(
        OUTPUT_DIR / "degree_normalization_mass_conservation_QA_v1_1.csv",
        index=False, encoding="utf-8-sig"
    )

    max_mass_error = float(qa_mass["absolute_error"].max()) if len(qa_mass) else np.nan
    if np.isfinite(max_mass_error) and max_mass_error > 1e-10:
        raise RuntimeError(
            f"Degree-normalization mass conservation failed. "
            f"Max absolute error={max_mass_error}"
        )

    metadata = {
        "analysis_name": "Weighted Herb-Target Projection v1.1 FROZEN",
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "status": FREEZE_STATUS,
        "inputs": {
            "FULL_weights": str(FULL_WEIGHTS_FILE),
            "FULL_weights_sha256": sha256_file(FULL_WEIGHTS_FILE),
            "LOPO_weights": str(LOPO_WEIGHTS_FILE),
            "LOPO_weights_sha256": sha256_file(LOPO_WEIGHTS_FILE),
            "herb_target_file": str(HERB_TARGET_FILE),
            "herb_target_sha256": sha256_file(HERB_TARGET_FILE),
            "herb_target_sheet": HERB_TARGET_SHEET,
            "detected_herb_column": str(detected_herb_col),
            "detected_target_column": str(detected_target_col),
        },
        "projection_methods": {
            "binary_union": "1 if any mapped syndrome herb links to target",
            "weighted_raw": "sum_h W_hs * A_ht",
            "herb_degree_normalized": "sum_h W_hs * A_ht / |T_h|",
            "ubiquity_corrected": (
                "sum_h W_hs * (A_ht / |T_h|) * "
                "[log((N_H+1)/(n_t+1))+1]"
            ),
        },
        "frozen_primary_for_next_RWR": "LOPO ubiquity_corrected_sum1",
        "freeze_version": FREEZE_VERSION,
        "idf_reference_universe": IDF_REFERENCE_UNIVERSE,
        "processed_form_parent_map": PROCESSED_FORM_PARENT_MAP,
        "important_note": (
            "This projection is frozen upstream before CAP alignment. "
            "Do not alter mappings, IDF universe, projection formula, or primary "
            "method in response to downstream CAP results; create a new version "
            "for any future methodological change."
        ),
        "CAP_transcriptomics_used": False,
        "GAT_used": False,
        "target_database_unique_herbs": int(ht["canonical_name"].nunique()),
        "target_database_unique_targets": int(ht["target_id"].nunique()),
        "target_database_unique_pairs": int(len(ht)),
        "target_database_universe_definition": IDF_REFERENCE_UNIVERSE,
        "processed_form_parent_aliases_do_not_expand_idf_universe": True,
        "max_degree_normalization_mass_error": max_mass_error,
        "uppercase_target_ids": UPPERCASE_TARGET_IDS,
        "min_targets_per_herb": MIN_TARGETS_PER_HERB,
        "max_targets_per_herb": MAX_TARGETS_PER_HERB,
        "top_k_values": TOP_K_VALUES,
    }

    with open(
        OUTPUT_DIR / "weighted_herb_target_projection_run_metadata_v1_1.json",
        "w", encoding="utf-8"
    ) as f:
        json.dump(metadata, f, ensure_ascii=False, indent=2)

    print("\nMapping QA")
    print("----------")
    print(mapping_qa_df.to_string(index=False))

    print("\nProjection summary")
    print("------------------")
    print(summary_df.to_string(index=False))

    print("\nPairwise syndrome similarity")
    print("----------------------------")
    print(
        similarity_df[
            similarity_df["metric"].isin(
                ["binary_union", "herb_degree_normalized_sum1",
                 "ubiquity_corrected_sum1"]
            )
        ].to_string(index=False)
    )

    print(f"\nDegree-normalization max mass error: {max_mass_error:.3e}")
    print(f"\nOutputs written to:\n{OUTPUT_DIR}")
    print(
        "\nFROZEN: LOPO ubiquity_corrected_sum1 is now the prespecified target seed "
        "representation for downstream RWR. Do not retune it using CAP outcomes."
    )


if __name__ == "__main__":
    main()
