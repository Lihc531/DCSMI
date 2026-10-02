# -*- coding: utf-8 -*-
"""
audit_and_extend_herb_target_mapping_v1_0.py

Audit and deterministically extend herb -> target-database mapping BEFORE RWR.

Core principle
--------------
Do not invent molecular targets. This script only extends mapping when an
observed prescription herb can be connected to an EXISTING target-database herb
by a frozen, auditable name rule:

1) exact canonical-name match
2) frozen verified alias/correction
3) processed-form -> parent-herb mapping, using herb_dictionary when available
4) optional manually verified mapping table supplied by the researcher

Everything else remains UNMAPPED and is exported for review.

The script DOES NOT:
- scrape the web
- infer targets from similar herbs/formulas
- use CAP transcriptomics
- use downstream alignment results
- fuzzy-match Chinese herb names automatically
- collapse processed forms in the prescription database itself

Outputs include an extended herb-target long table that can be used to rerun
weighted_herb_target_projection_v1_0.py.
"""

from pathlib import Path
import hashlib
import json
import re
from datetime import datetime, timezone
import numpy as np
import pandas as pd


# =============================================================================
# 1. USER SETTINGS — EDIT PATHS HERE
# =============================================================================

CANONICAL_XLSX = Path(
    r"E:\00project\Tanre Yongfei\重构\canonical_prescription_tables_v1.1.xlsx"
)

FROZEN_WEIGHTS_DIR = Path(
    r"E:\00project\Tanre Yongfei\重构\04_frozen_syndrome_herb_weights_v1_0"
)

FULL_WEIGHTS_FILE = (
    FROZEN_WEIGHTS_DIR / "FULL_syndrome_herb_weights_observed_v1_0.csv"
)

LOPO_WEIGHTS_FILE = (
    FROZEN_WEIGHTS_DIR / "LOPO_syndrome_herb_weights_observed_v1_0.csv"
)

# Existing TCMID target workbook
HERB_TARGET_FILE = Path(
    r"E:\00project\Tanre Yongfei\重构\HERB_TARGET_LONG_TABLE.xlsx"
)
HERB_TARGET_SHEET = "targets_long"
HERB_TARGET_HERB_COLUMN = "Herb_ChineseName"
HERB_TARGET_GENE_COLUMN = "Gene_Symbol"

OUTPUT_DIR = Path(
    r"E:\00project\Tanre Yongfei\重构\05A_herb_target_mapping_audit_v1_0"
)

HERB_DICTIONARY_SHEET = "herb_dictionary"

# Optional researcher-verified manual mapping file.
# Leave as None for the first run.
#
# Accepted CSV/XLSX columns:
#   canonical_name,target_mapping_name
# Optional:
#   mapping_note
#
# Example:
#   蜜麻黄,麻黄,processed form manually verified
#
# IMPORTANT: only enter mappings you have independently verified.
MANUAL_MAPPING_FILE = None
MANUAL_MAPPING_SHEET = None

SYNDROME_CODES = ["PDOL", "PHOL", "WHIL"]


# =============================================================================
# 2. FROZEN VERIFIED NAME RULES
# =============================================================================

# These are deterministic corrections/aliases already adjudicated in the
# canonical reconstruction. They do not create biological equivalence.
VERIFIED_ALIAS_MAP = {
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

# A conservative fallback for processed forms when herb_dictionary lacks a
# usable parent_herb. This list is intentionally explicit rather than regex-
# stripping "炙/蜜/生/熟/法/焦" automatically.
#
# Add entries only when the parent identity is unambiguous.
VERIFIED_PROCESSED_PARENT_FALLBACK = {
    "蜜麻黄": "麻黄",
    "炙甘草": "甘草",
    "蜜甘草": "甘草",
    "法半夏": "半夏",
    "姜半夏": "半夏",
    "清半夏": "半夏",
    "蜜紫菀": "紫菀",
    "炙紫菀": "紫菀",
    "蜜款冬花": "款冬花",
    "炙款冬花": "款冬花",
    "炙百部": "百部",
    "焦山楂": "山楂",
    "生地黄": "地黄",
    "熟地黄": "地黄",
    "生地": "地黄",
    "熟地": "地黄",
    "炙黄芪": "黄芪",
    "生黄芪": "黄芪",
    "生大黄": "大黄",
    "酒大黄": "大黄",
    "生杜仲": "杜仲",
    "盐杜仲": "杜仲",
    "生牡蛎": "牡蛎",
    "煅牡蛎": "牡蛎",
    "生薏苡仁": "薏苡仁",
    "炒薏苡仁": "薏苡仁",
    "生麻黄": "麻黄",
    "炙麻黄": "麻黄",
    "蜜桑白皮": "桑白皮",
    "炙桑白皮": "桑白皮",
    "炙桑皮": "桑白皮",
    "姜竹茹": "竹茹",
}


# =============================================================================
# 3. HELPERS
# =============================================================================

def norm_text(x):
    if pd.isna(x):
        return ""
    return re.sub(r"\s+", "", str(x).strip())


def sha256_file(path: Path, chunk_size=1024 * 1024):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk_size)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def require_file(path: Path, label):
    if not path.exists():
        raise FileNotFoundError(f"{label} not found:\n{path}")


def require_columns(df, cols, label):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"{label} missing required columns: {missing}\n"
            f"Available columns: {list(df.columns)}"
        )


def read_optional_manual_mapping():
    if MANUAL_MAPPING_FILE is None:
        return pd.DataFrame(
            columns=["canonical_name", "target_mapping_name", "mapping_note"]
        )

    path = Path(MANUAL_MAPPING_FILE)
    require_file(path, "MANUAL_MAPPING_FILE")

    if path.suffix.lower() in {".xlsx", ".xlsm", ".xls"}:
        df = pd.read_excel(
            path, sheet_name=0 if MANUAL_MAPPING_SHEET is None else MANUAL_MAPPING_SHEET
        )
    else:
        df = pd.read_csv(path)

    require_columns(df, ["canonical_name", "target_mapping_name"], "manual mapping")

    if "mapping_note" not in df.columns:
        df["mapping_note"] = ""

    df = df[["canonical_name", "target_mapping_name", "mapping_note"]].copy()
    df["canonical_name"] = df["canonical_name"].map(norm_text)
    df["target_mapping_name"] = df["target_mapping_name"].map(norm_text)
    df["mapping_note"] = df["mapping_note"].fillna("").astype(str)

    df = df[
        df["canonical_name"].ne("") & df["target_mapping_name"].ne("")
    ].drop_duplicates("canonical_name")

    return df


def load_weights(path: Path, representation: str):
    require_file(path, f"{representation} weights")
    x = pd.read_csv(path)
    require_columns(
        x,
        ["canonical_name", "syndrome_code", "primary_weight_W", "article_support_m"],
        f"{representation} weights",
    )
    x = x.copy()
    x["canonical_name"] = x["canonical_name"].map(norm_text)
    x["syndrome_code"] = x["syndrome_code"].astype(str).str.strip()
    x["primary_weight_W"] = pd.to_numeric(x["primary_weight_W"], errors="coerce")
    x["article_support_m"] = pd.to_numeric(x["article_support_m"], errors="coerce")

    x = x[
        x["canonical_name"].ne("")
        & x["syndrome_code"].isin(SYNDROME_CODES)
        & x["article_support_m"].gt(0)
        & x["primary_weight_W"].notna()
    ].copy()

    x["representation"] = representation
    return x


def load_dictionary():
    require_file(CANONICAL_XLSX, "canonical workbook")
    d = pd.read_excel(CANONICAL_XLSX, sheet_name=HERB_DICTIONARY_SHEET)
    require_columns(d, ["canonical_name"], "herb_dictionary")

    d = d.copy()
    d["canonical_name"] = d["canonical_name"].map(norm_text)

    for col in ["parent_herb", "processed_form", "mapping_method",
                "mapping_note", "dictionary_note"]:
        if col not in d.columns:
            d[col] = ""
        d[col] = d[col].fillna("").map(norm_text)

    # one row per canonical herb expected; retain first deterministically if not
    d = d.sort_values("canonical_name").drop_duplicates("canonical_name", keep="first")
    return d


def load_target_db():
    require_file(HERB_TARGET_FILE, "herb-target workbook")

    if HERB_TARGET_FILE.suffix.lower() in {".xlsx", ".xlsm", ".xls"}:
        x = pd.read_excel(HERB_TARGET_FILE, sheet_name=HERB_TARGET_SHEET)
    elif HERB_TARGET_FILE.suffix.lower() == ".csv":
        x = pd.read_csv(HERB_TARGET_FILE)
    else:
        x = pd.read_csv(HERB_TARGET_FILE, sep="\t")

    require_columns(
        x,
        [HERB_TARGET_HERB_COLUMN, HERB_TARGET_GENE_COLUMN],
        "herb-target table",
    )

    # Preserve all original columns for provenance.
    x = x.copy()
    x["_db_herb"] = x[HERB_TARGET_HERB_COLUMN].map(norm_text)
    x["_gene"] = x[HERB_TARGET_GENE_COLUMN].map(norm_text).str.upper()

    x = x[x["_db_herb"].ne("") & x["_gene"].ne("")].copy()

    # A target-database herb name can itself contain an already-adjudicated alias.
    # Canonical DB lookup includes both raw DB name and verified canonical alias.
    x["_db_herb_canonical"] = x["_db_herb"].map(
        lambda s: VERIFIED_ALIAS_MAP.get(s, s)
    )

    return x


def build_parent_lookup(dictionary: pd.DataFrame):
    lookup = {}
    source = {}

    for _, row in dictionary.iterrows():
        herb = row["canonical_name"]
        parent = row.get("parent_herb", "")
        processed = row.get("processed_form", "")

        if herb and parent and parent != herb:
            lookup[herb] = parent
            source[herb] = "herb_dictionary_parent"

    for herb, parent in VERIFIED_PROCESSED_PARENT_FALLBACK.items():
        if herb not in lookup:
            lookup[herb] = parent
            source[herb] = "verified_processed_parent_fallback"

    return lookup, source


def build_db_name_index(target_db: pd.DataFrame):
    """
    Return mapping from normalized/canonicalized DB herb name to the actual raw
    DB herb name(s). Multiple raw names may collapse to one canonical key.
    """
    records = target_db[["_db_herb", "_db_herb_canonical"]].drop_duplicates()
    index = {}

    for _, r in records.iterrows():
        raw = r["_db_herb"]
        can = r["_db_herb_canonical"]
        index.setdefault(raw, set()).add(raw)
        index.setdefault(can, set()).add(raw)

    return index


def resolve_one_herb(herb, db_index, parent_lookup, parent_source, manual_lookup):
    """
    Strict deterministic priority:
      exact DB > verified alias > dictionary/frozen parent > manual verified > unresolved

    No fuzzy matching.
    """
    # 1. Exact name
    if herb in db_index:
        candidates = sorted(db_index[herb])
        return {
            "mapping_status": "mapped",
            "mapping_method": "exact_name",
            "target_mapping_name": herb,
            "db_herb_names": "|".join(candidates),
            "mapping_note": "Exact canonical herb name found in target database.",
        }

    # 2. Verified alias/correction
    alias = VERIFIED_ALIAS_MAP.get(herb, "")
    if alias and alias in db_index:
        candidates = sorted(db_index[alias])
        return {
            "mapping_status": "mapped",
            "mapping_method": "verified_alias",
            "target_mapping_name": alias,
            "db_herb_names": "|".join(candidates),
            "mapping_note": f"Frozen verified alias/correction: {herb}->{alias}.",
        }

    # 3. Processed form -> parent herb
    parent = parent_lookup.get(herb, "")
    if parent:
        parent2 = VERIFIED_ALIAS_MAP.get(parent, parent)
        if parent2 in db_index:
            candidates = sorted(db_index[parent2])
            src = parent_source.get(herb, "parent_mapping")
            return {
                "mapping_status": "mapped",
                "mapping_method": "processed_form_to_parent",
                "target_mapping_name": parent2,
                "db_herb_names": "|".join(candidates),
                "mapping_note": (
                    f"Target-layer parent mapping only: {herb}->{parent2}; "
                    f"source={src}. Prescription identity remains {herb}."
                ),
            }

    # 4. Explicit researcher-verified mapping
    if herb in manual_lookup:
        target_name, note = manual_lookup[herb]
        target_name = VERIFIED_ALIAS_MAP.get(target_name, target_name)
        if target_name in db_index:
            candidates = sorted(db_index[target_name])
            return {
                "mapping_status": "mapped",
                "mapping_method": "manual_verified",
                "target_mapping_name": target_name,
                "db_herb_names": "|".join(candidates),
                "mapping_note": note or (
                    f"Researcher-verified target-layer mapping: "
                    f"{herb}->{target_name}."
                ),
            }
        else:
            return {
                "mapping_status": "manual_target_not_in_database",
                "mapping_method": "manual_verified_but_db_absent",
                "target_mapping_name": target_name,
                "db_herb_names": "",
                "mapping_note": (
                    f"Manual mapping supplied, but target mapping name "
                    f"'{target_name}' is absent from target database. {note}"
                ).strip(),
            }

    # unresolved
    note = "No exact, frozen alias, verified parent, or manual verified mapping."
    if parent:
        note += f" Parent candidate '{parent}' exists but is absent from target database."

    return {
        "mapping_status": "unmapped",
        "mapping_method": "unresolved",
        "target_mapping_name": "",
        "db_herb_names": "",
        "mapping_note": note,
    }


def build_audit(full_w, lopo_w, dictionary, target_db, manual):
    herbs = sorted(
        set(full_w["canonical_name"])
        | set(lopo_w["canonical_name"])
    )

    db_index = build_db_name_index(target_db)
    parent_lookup, parent_source = build_parent_lookup(dictionary)

    manual_lookup = {
        r["canonical_name"]: (r["target_mapping_name"], r["mapping_note"])
        for _, r in manual.iterrows()
    }

    dict_idx = dictionary.set_index("canonical_name", drop=False)

    rows = []
    for herb in herbs:
        resolved = resolve_one_herb(
            herb, db_index, parent_lookup, parent_source, manual_lookup
        )

        if herb in dict_idx.index:
            drow = dict_idx.loc[herb]
            if isinstance(drow, pd.DataFrame):
                drow = drow.iloc[0]
            parent = norm_text(drow.get("parent_herb", ""))
            processed = norm_text(drow.get("processed_form", ""))
        else:
            parent = parent_lookup.get(herb, "")
            processed = ""

        rows.append({
            "canonical_name": herb,
            "parent_herb": parent,
            "processed_form": processed,
            **resolved,
        })

    audit = pd.DataFrame(rows)

    # Add FULL/LOPO evidence summaries by syndrome.
    for rep, weights in [("FULL", full_w), ("LOPO", lopo_w)]:
        for s in SYNDROME_CODES:
            sub = weights[weights["syndrome_code"].eq(s)].set_index("canonical_name")
            audit[f"{rep}_{s}_observed"] = audit["canonical_name"].map(
                lambda h: int(h in sub.index)
            )
            audit[f"{rep}_{s}_weight"] = audit["canonical_name"].map(
                lambda h: float(sub.loc[h, "primary_weight_W"])
                if h in sub.index else 0.0
            )
            audit[f"{rep}_{s}_article_support"] = audit["canonical_name"].map(
                lambda h: int(sub.loc[h, "article_support_m"])
                if h in sub.index else 0
            )

    return audit


def choose_db_rows_for_mapping(target_db, mapping_row):
    """
    Retrieve target DB rows matching the resolved target mapping name.
    Matches against canonicalized DB herb names, preserving original DB rows.
    """
    target_name = mapping_row["target_mapping_name"]
    if not target_name:
        return target_db.iloc[0:0].copy()

    return target_db[
        target_db["_db_herb_canonical"].eq(target_name)
        | target_db["_db_herb"].eq(target_name)
    ].copy()


def build_extended_long(audit, target_db):
    """
    Build one row per prescription-canonical-herb -> gene relation.

    The herb name in the output is the frozen prescription canonical_name.
    target_mapping_name records which database herb supplied targets.
    Thus processed forms remain distinct at the prescription representation layer.
    """
    parts = []

    mapped = audit[audit["mapping_status"].eq("mapped")].copy()

    for _, m in mapped.iterrows():
        db_rows = choose_db_rows_for_mapping(target_db, m)
        if db_rows.empty:
            continue

        z = db_rows.copy()
        z["canonical_name"] = m["canonical_name"]
        z["parent_herb"] = m["parent_herb"]
        z["processed_form"] = m["processed_form"]
        z["target_mapping_name"] = m["target_mapping_name"]
        z["mapping_method"] = m["mapping_method"]
        z["mapping_note"] = m["mapping_note"]
        z["Gene_Symbol"] = z["_gene"]

        # Preserve original target database herb explicitly.
        z["target_database_herb_raw"] = z["_db_herb"]
        parts.append(z)

    if not parts:
        return pd.DataFrame()

    out = pd.concat(parts, ignore_index=True)

    # Binary herb-gene association: database duplicate records do not multiply.
    out = out.drop_duplicates(["canonical_name", "Gene_Symbol"]).copy()

    preferred = [
        "canonical_name", "parent_herb", "processed_form",
        "target_mapping_name", "mapping_method", "mapping_note",
        "target_database_herb_raw", "Gene_Symbol"
    ]
    remaining = [
        c for c in out.columns
        if c not in preferred and not c.startswith("_")
    ]

    return out[preferred + remaining]


def coverage_table(weights, audit, representation):
    a = audit[
        ["canonical_name", "mapping_status", "mapping_method", "target_mapping_name"]
    ].copy()

    x = weights.merge(a, on="canonical_name", how="left")
    x["is_mapped"] = x["mapping_status"].eq("mapped")

    rows = []
    for s in SYNDROME_CODES:
        z = x[x["syndrome_code"].eq(s)].copy()
        n = z["canonical_name"].nunique()
        mapped_n = z.loc[z["is_mapped"], "canonical_name"].nunique()
        total_w = z["primary_weight_W"].sum()
        mapped_w = z.loc[z["is_mapped"], "primary_weight_W"].sum()

        method_counts = (
            z[z["is_mapped"]]
            .drop_duplicates("canonical_name")["mapping_method"]
            .value_counts()
            .to_dict()
        )

        rows.append({
            "representation": representation,
            "syndrome_code": s,
            "observed_herbs": int(n),
            "mapped_herbs_after_audit": int(mapped_n),
            "unmapped_herbs_after_audit": int(n - mapped_n),
            "herb_mapping_fraction_after_audit": mapped_n / n if n else np.nan,
            "total_frozen_weight": float(total_w),
            "mapped_frozen_weight_after_audit": float(mapped_w),
            "weight_mapping_fraction_after_audit": (
                mapped_w / total_w if total_w else np.nan
            ),
            "n_exact": int(method_counts.get("exact_name", 0)),
            "n_verified_alias": int(method_counts.get("verified_alias", 0)),
            "n_processed_to_parent": int(
                method_counts.get("processed_form_to_parent", 0)
            ),
            "n_manual_verified": int(method_counts.get("manual_verified", 0)),
        })

    return pd.DataFrame(rows)


def unresolved_priority(audit):
    x = audit[~audit["mapping_status"].eq("mapped")].copy()

    weight_cols = [
        c for c in x.columns
        if c.startswith("LOPO_") and c.endswith("_weight")
    ]
    support_cols = [
        c for c in x.columns
        if c.startswith("LOPO_") and c.endswith("_article_support")
    ]

    x["LOPO_max_weight"] = x[weight_cols].max(axis=1) if weight_cols else 0.0
    x["LOPO_total_weight"] = x[weight_cols].sum(axis=1) if weight_cols else 0.0
    x["LOPO_max_article_support"] = (
        x[support_cols].max(axis=1) if support_cols else 0
    )

    return x.sort_values(
        ["LOPO_max_weight", "LOPO_max_article_support", "canonical_name"],
        ascending=[False, False, True],
    )


def manual_review_template(unresolved):
    cols = [
        "canonical_name", "parent_herb", "processed_form",
        "mapping_status", "mapping_note",
        "LOPO_max_weight", "LOPO_total_weight", "LOPO_max_article_support"
    ]
    cols = [c for c in cols if c in unresolved.columns]
    x = unresolved[cols].copy()
    x["target_mapping_name"] = ""
    x["researcher_decision"] = ""
    x["mapping_note_manual"] = ""
    return x


# =============================================================================
# 4. MAIN
# =============================================================================

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("[1/7] Loading frozen FULL/LOPO weights...")
    full_w = load_weights(FULL_WEIGHTS_FILE, "FULL")
    lopo_w = load_weights(LOPO_WEIGHTS_FILE, "LOPO")

    print("[2/7] Loading canonical herb dictionary...")
    dictionary = load_dictionary()

    print("[3/7] Loading existing TCMID herb-target table...")
    target_db = load_target_db()
    print(f"    Target DB herb names: {target_db['_db_herb'].nunique()}")
    print(f"    Target genes:         {target_db['_gene'].nunique()}")
    print(f"    Herb-gene raw rows:   {len(target_db)}")

    print("[4/7] Applying deterministic mapping rules...")
    manual = read_optional_manual_mapping()
    audit = build_audit(full_w, lopo_w, dictionary, target_db, manual)

    print("[5/7] Building extended herb-target table...")
    extended = build_extended_long(audit, target_db)

    if extended.empty:
        raise RuntimeError("No herb-target mappings were resolved.")

    print("[6/7] Computing coverage and unresolved priorities...")
    coverage = pd.concat(
        [
            coverage_table(full_w, audit, "FULL"),
            coverage_table(lopo_w, audit, "LOPO"),
        ],
        ignore_index=True,
    )

    unresolved = unresolved_priority(audit)
    template = manual_review_template(unresolved)

    # Mapping method QA
    method_qa = (
        audit.groupby(["mapping_status", "mapping_method"], dropna=False)
        .size()
        .reset_index(name="n_herbs")
        .sort_values(["mapping_status", "n_herbs"], ascending=[True, False])
    )

    print("[7/7] Writing outputs...")

    audit.to_csv(
        OUTPUT_DIR / "herb_target_mapping_audit_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )
    coverage.to_csv(
        OUTPUT_DIR / "mapping_coverage_after_extension_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )
    unresolved.to_csv(
        OUTPUT_DIR / "unresolved_herbs_prioritized_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )
    template.to_csv(
        OUTPUT_DIR / "manual_mapping_review_template_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )
    method_qa.to_csv(
        OUTPUT_DIR / "mapping_method_QA_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )

    # This is the main file to feed back into weighted projection.
    extended.to_csv(
        OUTPUT_DIR / "extended_herb_target_long_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )

    # Also write a compact two-column version that the projection script can
    # consume with minimal configuration.
    compact = (
        extended[["canonical_name", "Gene_Symbol"]]
        .drop_duplicates()
        .sort_values(["canonical_name", "Gene_Symbol"])
    )
    compact.to_csv(
        OUTPUT_DIR / "extended_herb_target_pairs_compact_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )

    # Resolved mapping dictionary for provenance.
    resolved_dict = audit[
        [
            "canonical_name", "parent_herb", "processed_form",
            "mapping_status", "mapping_method",
            "target_mapping_name", "db_herb_names", "mapping_note"
        ]
    ].copy()
    resolved_dict.to_csv(
        OUTPUT_DIR / "resolved_target_mapping_dictionary_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )

    # Human-readable workbook.
    with pd.ExcelWriter(
        OUTPUT_DIR / "herb_target_mapping_audit_v1_0_summary.xlsx",
        engine="openpyxl"
    ) as writer:
        coverage.to_excel(writer, sheet_name="coverage", index=False)
        method_qa.to_excel(writer, sheet_name="mapping_methods", index=False)
        audit.to_excel(writer, sheet_name="all_herbs", index=False)
        unresolved.to_excel(writer, sheet_name="unresolved_priority", index=False)
        template.to_excel(writer, sheet_name="manual_review_template", index=False)
        resolved_dict.to_excel(writer, sheet_name="mapping_dictionary", index=False)

    # Compare coverage with exact-only mapping to quantify extension benefit.
    exact_audit = audit.copy()
    exact_audit["mapping_status"] = np.where(
        exact_audit["mapping_method"].eq("exact_name"), "mapped", "unmapped"
    )
    exact_cov = pd.concat(
        [
            coverage_table(full_w, exact_audit, "FULL_exact_only"),
            coverage_table(lopo_w, exact_audit, "LOPO_exact_only"),
        ],
        ignore_index=True,
    )
    exact_cov.to_csv(
        OUTPUT_DIR / "exact_only_coverage_baseline_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )

    metadata = {
        "analysis_name": "Herb-Target Mapping Audit and Deterministic Extension v1.0",
        "status": "UPSTREAM_MAPPING_AUDIT",
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "principle": (
            "Extend only by exact name, frozen verified alias/correction, "
            "verified processed-form-to-parent mapping, or explicit researcher-"
            "verified mapping. Never infer targets for unresolved herbs."
        ),
        "inputs": {
            "canonical_workbook": str(CANONICAL_XLSX),
            "canonical_workbook_sha256": sha256_file(CANONICAL_XLSX),
            "FULL_weights": str(FULL_WEIGHTS_FILE),
            "FULL_weights_sha256": sha256_file(FULL_WEIGHTS_FILE),
            "LOPO_weights": str(LOPO_WEIGHTS_FILE),
            "LOPO_weights_sha256": sha256_file(LOPO_WEIGHTS_FILE),
            "herb_target_file": str(HERB_TARGET_FILE),
            "herb_target_file_sha256": sha256_file(HERB_TARGET_FILE),
            "herb_target_sheet": HERB_TARGET_SHEET,
            "manual_mapping_file": (
                None if MANUAL_MAPPING_FILE is None else str(MANUAL_MAPPING_FILE)
            ),
        },
        "mapping_priority": [
            "exact_name",
            "verified_alias",
            "processed_form_to_parent",
            "manual_verified",
            "unresolved",
        ],
        "fuzzy_matching_used": False,
        "web_inference_used": False,
        "CAP_data_used": False,
        "processed_forms_preserved_in_prescription_identity": True,
        "processed_parent_targets_used_only_at_target_layer": True,
        "counts": {
            "audited_unique_frozen_herbs": int(len(audit)),
            "mapped_unique_frozen_herbs": int(
                audit["mapping_status"].eq("mapped").sum()
            ),
            "unresolved_unique_frozen_herbs": int(
                (~audit["mapping_status"].eq("mapped")).sum()
            ),
            "extended_unique_herbs": int(extended["canonical_name"].nunique()),
            "extended_unique_genes": int(extended["Gene_Symbol"].nunique()),
            "extended_unique_herb_gene_pairs": int(len(compact)),
        },
    }

    with open(
        OUTPUT_DIR / "herb_target_mapping_audit_metadata_v1_0.json",
        "w", encoding="utf-8"
    ) as f:
        json.dump(metadata, f, ensure_ascii=False, indent=2)

    print("\nMapping method QA")
    print("-----------------")
    print(method_qa.to_string(index=False))

    print("\nCoverage after deterministic extension")
    print("--------------------------------------")
    print(coverage.to_string(index=False))

    print("\nHighest-priority unresolved herbs (top 20)")
    print("------------------------------------------")
    show_cols = [
        "canonical_name", "parent_herb", "processed_form",
        "LOPO_max_weight", "LOPO_max_article_support",
        "mapping_status", "mapping_note"
    ]
    show_cols = [c for c in show_cols if c in unresolved.columns]
    print(unresolved[show_cols].head(20).to_string(index=False))

    print(f"\nExtended unique herbs: {extended['canonical_name'].nunique()}")
    print(f"Extended unique genes: {extended['Gene_Symbol'].nunique()}")
    print(f"Extended herb-gene pairs: {len(compact)}")

    print(
        "\nNEXT: inspect unresolved_herbs_prioritized_v1_0.csv before RWR. "
        "Do not manually map a herb merely to improve downstream coverage."
    )
    print(f"\nOutputs written to:\n{OUTPUT_DIR}")


if __name__ == "__main__":
    main()
