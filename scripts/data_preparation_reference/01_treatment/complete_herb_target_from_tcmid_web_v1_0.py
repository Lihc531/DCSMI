# -*- coding: utf-8 -*-
"""
complete_herb_target_from_tcmid_web_v1_0.py

Purpose
-------
Complete HERB_TARGET_LONG_TABLE.xlsx by retrieving missing herb-target records
directly from the public TCM-ID component pages.

Important methodological rule
-----------------------------
This script does NOT predict or infer targets.
It only appends target rows explicitly displayed by TCM-ID under:
"Targeted Human Proteins by the Ingredient of Component".

The original workbook is never overwritten.

Recommended workflow
--------------------
1. Put this script in PyCharm.
2. Edit the paths in USER SETTINGS.
3. Run once. A local cache is created, so interrupted runs can resume.
4. Upload the output folder for audit before downstream projection/RWR.

Dependencies
------------
pip install requests beautifulsoup4 pandas openpyxl lxml
"""

from __future__ import annotations

from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import re
import time
import hashlib
from io import StringIO

import pandas as pd
import requests
from bs4 import BeautifulSoup


# =============================================================================
# 1. USER SETTINGS
# =============================================================================

INPUT_XLSX = Path(
    r"E:\00project\Tanre Yongfei\重构\HERB_TARGET_LONG_TABLE.xlsx"
)

# Output from audit_and_extend_herb_target_mapping_v1_0.py.
# The script uses canonical_name values in this file as the herbs to retrieve.
UNRESOLVED_CSV = Path(
    r"E:\00project\Tanre Yongfei\重构\05A_herb_target_mapping_audit_v1_0"
    r"\unresolved_herbs_prioritized_v1_0.csv"
)

OUTPUT_DIR = Path(
    r"E:\00project\Tanre Yongfei\重构\05B_TCMID_web_completion_v1_0"
)

INPUT_TARGET_SHEET = "targets_long"
INPUT_MAPPING_SHEET = "name_to_tcmh_mapping"

# TCM-ID component IDs observed on the site are approximately TCMH1...TCMH2700+.
# A slightly wider range is used to avoid missing later IDs.
TCMH_ID_MIN = 1
TCMH_ID_MAX = 3000

# Network behavior.
MAX_WORKERS = 6
REQUEST_TIMEOUT = 25
RETRIES = 3
RETRY_SLEEP_SECONDS = 1.5
POLITE_DELAY_SECONDS = 0.08

BASE_URL = "https://www.bidd.group/TCMID/herb.php?herb={tcmh_id}"

# If False, only unresolved herbs from UNRESOLVED_CSV are fetched.
# Keep False for the primary run.
REFRESH_ALREADY_MAPPED_HERBS = False


# =============================================================================
# 2. FROZEN NAME NORMALIZATION
# =============================================================================

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

PROCESSED_PARENT_MAP = {
    "蜜麻黄": "麻黄",
    "炙麻黄": "麻黄",
    "法半夏": "半夏",
    "姜半夏": "半夏",
    "清半夏": "半夏",
    "制半夏": "半夏",
    "炙甘草": "甘草",
    "蜜甘草": "甘草",
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
    "蜜桑白皮": "桑白皮",
    "炙桑白皮": "桑白皮",
    "炙桑皮": "桑白皮",
    "姜竹茹": "竹茹",
}


# =============================================================================
# 3. HELPERS
# =============================================================================

def norm_text(x) -> str:
    if pd.isna(x):
        return ""
    return re.sub(r"\s+", "", str(x).strip())


def canonical_lookup_name(name: str) -> str:
    name = norm_text(name)
    name = VERIFIED_ALIAS_MAP.get(name, name)
    name = PROCESSED_PARENT_MAP.get(name, name)
    name = VERIFIED_ALIAS_MAP.get(name, name)
    return name


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def get_session() -> requests.Session:
    s = requests.Session()
    s.headers.update({
        "User-Agent": (
            "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
            "AppleWebKit/537.36 Chrome/128 Safari/537.36 "
            "AcademicResearch/1.0"
        ),
        "Accept-Language": "zh-CN,zh;q=0.9,en;q=0.8",
    })
    return s


def fetch_html(url: str, session: requests.Session) -> str | None:
    last_err = None
    for attempt in range(1, RETRIES + 1):
        try:
            r = session.get(url, timeout=REQUEST_TIMEOUT)
            if r.status_code == 200 and r.text:
                # TCM-ID pages may omit/incorrectly declare charset.
                r.encoding = r.apparent_encoding or "utf-8"
                return r.text
            last_err = RuntimeError(f"HTTP {r.status_code}")
        except Exception as e:
            last_err = e

        if attempt < RETRIES:
            time.sleep(RETRY_SLEEP_SECONDS * attempt)

    return None


def parse_component_identity(html: str, tcmh_id: str) -> dict | None:
    """
    Parse component identity. We require an explicit Chinese-name field or
    component heading; pages that do not resolve to a component are ignored.
    """
    soup = BeautifulSoup(html, "lxml")
    text = soup.get_text("\n", strip=True)

    # Component heading often: "Component: Bupleurum chinense 柴胡"
    heading = None
    for tag in soup.find_all(["h1", "h2", "h3", "h4"]):
        t = tag.get_text(" ", strip=True)
        if "Component:" in t:
            heading = t
            break

    # Extract labeled values using text sequence.
    def after_label(label):
        node = soup.find(string=lambda s: s and label in s)
        if node:
            parent = node.parent
            nxt = parent.find_next()
            # Move until a non-empty text not identical to label.
            for _ in range(8):
                if nxt is None:
                    break
                val = nxt.get_text(" ", strip=True)
                if val and label not in val and not val.startswith("#"):
                    return val
                nxt = nxt.find_next()
        return ""

    latin = after_label("Latin Name")
    pinyin = after_label("Chinese Pinyin Name")
    chinese = after_label("中文名")

    # More reliable fallback from heading: last CJK token(s).
    if not chinese and heading:
        m = re.search(r"([\u4e00-\u9fff（）()·\-]+)\s*$", heading)
        if m:
            chinese = m.group(1).strip()

    chinese = norm_text(chinese)
    if not chinese:
        return None

    return {
        "TCMH_ID": tcmh_id,
        "Chinese_Name": chinese,
        "Latin_Name": str(latin).strip(),
        "Pinyin_Name": str(pinyin).strip(),
        "Source_URL": BASE_URL.format(tcmh_id=tcmh_id),
    }


def parse_human_targets(html: str, identity: dict) -> pd.DataFrame:
    """
    Parse only the TCM-ID table titled:
    'Targeted Human Proteins by the Ingredient of Component'
    """
    try:
        tables = pd.read_html(StringIO(html))
    except Exception:
        return pd.DataFrame()

    chosen = None
    for t in tables:
        cols = [str(c).strip() for c in t.columns]
        norm_cols = {re.sub(r"\s+", " ", c).lower(): c for c in cols}
        if (
            any("target id" in c for c in norm_cols)
            and any("gene symbol" in c for c in norm_cols)
            and any("uniprot" in c for c in norm_cols)
        ):
            chosen = t.copy()
            break

    if chosen is None or chosen.empty:
        return pd.DataFrame()

    chosen.columns = [re.sub(r"\s+", " ", str(c)).strip() for c in chosen.columns]

    rename = {}
    for c in chosen.columns:
        lc = c.lower()
        if "target id" in lc:
            rename[c] = "Target_ID"
        elif "gene symbol" in lc:
            rename[c] = "Gene_Symbol"
        elif "target name" in lc:
            rename[c] = "Target_Name"
        elif "target class" in lc:
            rename[c] = "Target_Class"
        elif "uniprot" in lc:
            rename[c] = "Uniprot_ID"

    chosen = chosen.rename(columns=rename)
    required = ["Target_ID", "Gene_Symbol"]
    if not all(c in chosen.columns for c in required):
        return pd.DataFrame()

    for c in ["Target_Name", "Target_Class", "Uniprot_ID"]:
        if c not in chosen.columns:
            chosen[c] = ""

    chosen = chosen[
        ["Target_ID", "Gene_Symbol", "Target_Name", "Target_Class", "Uniprot_ID"]
    ].copy()

    for c in chosen.columns:
        chosen[c] = chosen[c].fillna("").astype(str).str.strip()

    chosen["Gene_Symbol"] = chosen["Gene_Symbol"].str.upper()
    chosen = chosen[
        chosen["Gene_Symbol"].ne("")
        & ~chosen["Gene_Symbol"].str.upper().isin(["NA", "NAN", "NONE"])
    ].copy()

    chosen.insert(0, "TCMH_ID", identity["TCMH_ID"])
    chosen["Source_URL"] = identity["Source_URL"]
    chosen["Herb_ChineseName"] = identity["Chinese_Name"]
    chosen["Retrieval_Method"] = "TCMID_public_component_page"
    chosen["Retrieval_Evidence"] = (
        "Targeted Human Proteins by the Ingredient of Component"
    )

    return chosen.drop_duplicates(["TCMH_ID", "Gene_Symbol"])


def load_unresolved_names() -> list[str]:
    if not UNRESOLVED_CSV.exists():
        raise FileNotFoundError(
            f"UNRESOLVED_CSV not found:\n{UNRESOLVED_CSV}"
        )
    x = pd.read_csv(UNRESOLVED_CSV)
    if "canonical_name" not in x.columns:
        raise ValueError("UNRESOLVED_CSV must contain canonical_name")
    names = [
        norm_text(v) for v in x["canonical_name"].tolist()
        if norm_text(v)
    ]
    return sorted(set(names))


def load_existing():
    if not INPUT_XLSX.exists():
        raise FileNotFoundError(f"INPUT_XLSX not found:\n{INPUT_XLSX}")

    targets = pd.read_excel(INPUT_XLSX, sheet_name=INPUT_TARGET_SHEET)
    mapping = pd.read_excel(INPUT_XLSX, sheet_name=INPUT_MAPPING_SHEET)

    required = [
        "TCMH_ID", "Target_ID", "Gene_Symbol", "Target_Name",
        "Target_Class", "Uniprot_ID", "Source_URL", "Herb_ChineseName"
    ]
    missing = [c for c in required if c not in targets.columns]
    if missing:
        raise ValueError(f"targets_long missing columns: {missing}")

    return targets, mapping


# =============================================================================
# 4. DISCOVERY CACHE
# =============================================================================

def discover_component(tcmh_num: int) -> dict:
    tcmh_id = f"TCMH{tcmh_num}"
    url = BASE_URL.format(tcmh_id=tcmh_id)
    session = get_session()
    html = fetch_html(url, session)
    time.sleep(POLITE_DELAY_SECONDS)

    if not html:
        return {"TCMH_ID": tcmh_id, "status": "fetch_failed"}

    ident = parse_component_identity(html, tcmh_id)
    if ident is None:
        return {"TCMH_ID": tcmh_id, "status": "not_component"}

    ident["status"] = "ok"
    return ident


def build_or_resume_component_index() -> pd.DataFrame:
    cache_csv = OUTPUT_DIR / "TCMID_component_index_cache_v1_0.csv"

    done = {}
    if cache_csv.exists():
        old = pd.read_csv(cache_csv)
        if "TCMH_ID" in old.columns:
            done = {
                str(r["TCMH_ID"]): r.to_dict()
                for _, r in old.iterrows()
            }
        print(f"Resuming component-index cache: {len(done)} IDs already checked.")

    todo = [
        n for n in range(TCMH_ID_MIN, TCMH_ID_MAX + 1)
        if f"TCMH{n}" not in done
    ]

    if todo:
        print(
            f"Scanning {len(todo)} TCM-ID component IDs "
            f"({TCMH_ID_MIN}..{TCMH_ID_MAX})..."
        )

        buffer = []
        completed = 0

        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
            futures = {ex.submit(discover_component, n): n for n in todo}
            for fut in as_completed(futures):
                try:
                    row = fut.result()
                except Exception as e:
                    n = futures[fut]
                    row = {
                        "TCMH_ID": f"TCMH{n}",
                        "status": "exception",
                        "error": repr(e),
                    }

                done[row["TCMH_ID"]] = row
                buffer.append(row)
                completed += 1

                if len(buffer) >= 100:
                    pd.DataFrame(list(done.values())).to_csv(
                        cache_csv, index=False, encoding="utf-8-sig"
                    )
                    buffer.clear()
                    print(f"  checked {completed}/{len(todo)} new IDs")

        pd.DataFrame(list(done.values())).to_csv(
            cache_csv, index=False, encoding="utf-8-sig"
        )

    index = pd.DataFrame(list(done.values()))
    if "status" in index.columns:
        index = index[index["status"].eq("ok")].copy()

    for c in ["Chinese_Name", "Latin_Name", "Pinyin_Name"]:
        if c not in index.columns:
            index[c] = ""
        index[c] = index[c].fillna("").astype(str).str.strip()

    index["lookup_name"] = index["Chinese_Name"].map(canonical_lookup_name)
    return index


# =============================================================================
# 5. MATCH UNRESOLVED HERBS TO TCM-ID COMPONENTS
# =============================================================================

def match_unresolved(unresolved: list[str], index: pd.DataFrame) -> pd.DataFrame:
    rows = []

    for herb in unresolved:
        lookup = canonical_lookup_name(herb)

        # Strict Chinese-name equality after frozen normalization.
        hits = index[index["lookup_name"].eq(lookup)].copy()

        if hits.empty:
            rows.append({
                "canonical_name": herb,
                "lookup_name": lookup,
                "match_status": "not_found_in_scanned_TCMID",
                "Matched_TCMH_IDs": "",
                "Chosen_TCMH_ID": "",
                "TCMID_Chinese_Name": "",
                "TCMID_Latin_Name": "",
                "Source_URL": "",
            })
            continue

        # If multiple IDs exist, do NOT silently choose based on target count.
        # Prefer exact raw Chinese-name equality. Otherwise mark ambiguous.
        exact = hits[hits["Chinese_Name"].map(norm_text).eq(norm_text(lookup))]
        candidates = exact if not exact.empty else hits

        ids = sorted(candidates["TCMH_ID"].astype(str).unique())

        if len(ids) == 1:
            chosen = candidates.iloc[0]
            status = "unique_exact_or_normalized_match"
            chosen_id = ids[0]
        else:
            # Multiple components can be variants/extracts. Exact herb name may
            # still uniquely identify one after removing obvious extract forms.
            plain = candidates[
                ~candidates["Chinese_Name"].str.contains(
                    "浸膏|提取|流浸膏|膏粉|油|炭|粉", regex=True, na=False
                )
            ]
            plain_ids = sorted(plain["TCMH_ID"].astype(str).unique())

            if len(plain_ids) == 1:
                chosen = plain.iloc[0]
                status = "unique_plain_component_among_multiple"
                chosen_id = plain_ids[0]
            else:
                chosen = candidates.iloc[0]
                status = "ambiguous_manual_review_required"
                chosen_id = ""

        rows.append({
            "canonical_name": herb,
            "lookup_name": lookup,
            "match_status": status,
            "Matched_TCMH_IDs": "|".join(ids),
            "Chosen_TCMH_ID": chosen_id,
            "TCMID_Chinese_Name": chosen.get("Chinese_Name", ""),
            "TCMID_Latin_Name": chosen.get("Latin_Name", ""),
            "Source_URL": (
                BASE_URL.format(tcmh_id=chosen_id) if chosen_id else ""
            ),
        })

    return pd.DataFrame(rows)


# =============================================================================
# 6. FETCH TARGETS FOR UNIQUE MATCHES
# =============================================================================

def fetch_targets_for_match(row: dict) -> tuple[dict, pd.DataFrame]:
    herb = row["canonical_name"]
    tcmh_id = row["Chosen_TCMH_ID"]

    if not tcmh_id:
        return row, pd.DataFrame()

    session = get_session()
    url = BASE_URL.format(tcmh_id=tcmh_id)
    html = fetch_html(url, session)
    time.sleep(POLITE_DELAY_SECONDS)

    if not html:
        row = dict(row)
        row["target_fetch_status"] = "fetch_failed"
        return row, pd.DataFrame()

    ident = parse_component_identity(html, tcmh_id)
    if ident is None:
        row = dict(row)
        row["target_fetch_status"] = "identity_parse_failed"
        return row, pd.DataFrame()

    targets = parse_human_targets(html, ident)

    row = dict(row)
    row["target_fetch_status"] = (
        "ok" if not targets.empty else "no_human_targets_displayed"
    )
    row["n_targets_retrieved"] = int(len(targets))

    if not targets.empty:
        # Keep the canonical prescription herb as the analysis-facing herb name,
        # while retaining the exact TCM-ID component name separately.
        targets["TCMID_Herb_ChineseName"] = targets["Herb_ChineseName"]
        targets["Herb_ChineseName"] = herb
        targets["Canonical_Name"] = herb
        targets["TCMID_Matched_Name"] = ident["Chinese_Name"]

    return row, targets


# =============================================================================
# 7. MAIN
# =============================================================================

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("[1/7] Loading current herb-target workbook...")
    existing_targets, existing_mapping = load_existing()

    print("[2/7] Loading unresolved herbs...")
    unresolved = load_unresolved_names()
    print(f"    unresolved herbs requested: {len(unresolved)}")

    if REFRESH_ALREADY_MAPPED_HERBS:
        existing_names = sorted(
            set(existing_targets["Herb_ChineseName"].dropna().map(norm_text))
        )
        unresolved = sorted(set(unresolved) | set(existing_names))

    print("[3/7] Building/resuming TCM-ID component index...")
    index = build_or_resume_component_index()
    index.to_csv(
        OUTPUT_DIR / "TCMID_component_index_resolved_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )
    print(f"    valid TCM-ID components indexed: {len(index)}")

    print("[4/7] Strictly matching unresolved herb names...")
    matches = match_unresolved(unresolved, index)
    matches.to_csv(
        OUTPUT_DIR / "unresolved_to_TCMID_mapping_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )

    retrievable = matches[matches["Chosen_TCMH_ID"].ne("")].copy()
    ambiguous = matches[
        matches["match_status"].eq("ambiguous_manual_review_required")
    ].copy()
    not_found = matches[
        matches["match_status"].eq("not_found_in_scanned_TCMID")
    ].copy()

    print(f"    unique retrievable: {len(retrievable)}")
    print(f"    ambiguous:          {len(ambiguous)}")
    print(f"    not found:          {len(not_found)}")

    print("[5/7] Fetching explicitly displayed human targets...")
    fetch_audit = []
    new_parts = []

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futures = {
            ex.submit(fetch_targets_for_match, r.to_dict()): i
            for i, (_, r) in enumerate(retrievable.iterrows())
        }

        done = 0
        for fut in as_completed(futures):
            audit_row, targets = fut.result()
            fetch_audit.append(audit_row)
            if not targets.empty:
                new_parts.append(targets)
            done += 1
            print(f"    fetched {done}/{len(retrievable)}")

    fetch_audit = pd.DataFrame(fetch_audit)
    fetch_audit.to_csv(
        OUTPUT_DIR / "TCMID_target_fetch_audit_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )

    if new_parts:
        new_targets = pd.concat(new_parts, ignore_index=True)
    else:
        new_targets = pd.DataFrame()

    print("[6/7] Merging with original workbook without overwriting provenance...")

    original = existing_targets.copy()
    original["Data_Source"] = "original_uploaded_TCMID_table"
    original["Retrieval_Method"] = original.get(
        "Retrieval_Method", "preexisting"
    )

    if not new_targets.empty:
        # Align to core schema plus provenance.
        for c in original.columns:
            if c not in new_targets.columns:
                new_targets[c] = ""

        for c in new_targets.columns:
            if c not in original.columns:
                original[c] = ""

        new_targets["Data_Source"] = "TCM-ID_public_web_completion"

        combined = pd.concat(
            [original[new_targets.columns], new_targets],
            ignore_index=True
        )

        combined["Gene_Symbol"] = (
            combined["Gene_Symbol"].fillna("").astype(str).str.upper().str.strip()
        )
        combined["Herb_ChineseName"] = (
            combined["Herb_ChineseName"].fillna("").map(norm_text)
        )

        # De-duplicate herb-gene relations, preferring the original row.
        combined["_source_priority"] = combined["Data_Source"].map({
            "original_uploaded_TCMID_table": 0,
            "TCM-ID_public_web_completion": 1,
        }).fillna(9)

        combined = (
            combined.sort_values(
                ["Herb_ChineseName", "Gene_Symbol", "_source_priority"]
            )
            .drop_duplicates(
                ["Herb_ChineseName", "Gene_Symbol"], keep="first"
            )
            .drop(columns="_source_priority")
        )
    else:
        combined = original

    # Updated mapping sheet: preserve old mappings and append new unique matches.
    map_add = matches[matches["Chosen_TCMH_ID"].ne("")][
        ["canonical_name", "Matched_TCMH_IDs", "Chosen_TCMH_ID"]
    ].copy()
    map_add = map_add.rename(columns={"canonical_name": "Herb_ChineseName"})

    mapping_combined = pd.concat(
        [existing_mapping, map_add],
        ignore_index=True
    ).drop_duplicates("Herb_ChineseName", keep="first")

    # QA summary
    summary = pd.DataFrame([
        ["Original target rows", len(existing_targets)],
        ["Original unique herbs", existing_targets["Herb_ChineseName"].nunique()],
        ["Original unique genes", existing_targets["Gene_Symbol"].nunique()],
        ["Unresolved herbs requested", len(unresolved)],
        ["Unique TCM-ID component matches", len(retrievable)],
        ["Ambiguous matches withheld", len(ambiguous)],
        ["Not found in scanned TCM-ID", len(not_found)],
        ["New web target rows before dedup", len(new_targets)],
        ["Completed unique herbs", combined["Herb_ChineseName"].nunique()],
        ["Completed unique genes", combined["Gene_Symbol"].nunique()],
        ["Completed herb-gene rows", len(combined)],
    ], columns=["Metric", "Value"])

    print("[7/7] Writing completed workbook and audit files...")

    output_xlsx = OUTPUT_DIR / "HERB_TARGET_LONG_TABLE_TCMID_web_completed_v1_0.xlsx"

    with pd.ExcelWriter(output_xlsx, engine="openpyxl") as writer:
        combined.to_excel(writer, sheet_name="targets_long", index=False)
        mapping_combined.to_excel(
            writer, sheet_name="name_to_tcmh_mapping", index=False
        )
        matches.to_excel(writer, sheet_name="web_mapping_audit", index=False)
        fetch_audit.to_excel(writer, sheet_name="target_fetch_audit", index=False)
        summary.to_excel(writer, sheet_name="completion_summary", index=False)

    combined.to_csv(
        OUTPUT_DIR / "HERB_TARGET_LONG_TABLE_TCMID_web_completed_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )

    if not new_targets.empty:
        new_targets.to_csv(
            OUTPUT_DIR / "newly_retrieved_TCMID_targets_only_v1_0.csv",
            index=False, encoding="utf-8-sig"
        )

    ambiguous.to_csv(
        OUTPUT_DIR / "ambiguous_TCMID_matches_manual_review_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )
    not_found.to_csv(
        OUTPUT_DIR / "not_found_in_TCMID_v1_0.csv",
        index=False, encoding="utf-8-sig"
    )

    metadata = {
        "analysis": "TCM-ID public web completion v1.0",
        "input_xlsx": str(INPUT_XLSX),
        "input_sha256": sha256_file(INPUT_XLSX),
        "unresolved_csv": str(UNRESOLVED_CSV),
        "source_database": "TCM-ID",
        "source_base_url": "https://www.bidd.group/TCMID/",
        "target_evidence": (
            "Targeted Human Proteins by the Ingredient of Component"
        ),
        "target_prediction_performed": False,
        "fuzzy_name_matching_performed": False,
        "ambiguous_components_automatically_used": False,
        "TCMH_scan_range": [TCMH_ID_MIN, TCMH_ID_MAX],
        "rules": [
            "Preserve original rows.",
            "Append only explicit human target rows displayed by TCM-ID.",
            "Use strict Chinese-name matching after frozen aliases.",
            "Do not use ambiguous component matches automatically.",
            "Do not infer targets for herbs absent from TCM-ID.",
        ],
    }
    with open(
        OUTPUT_DIR / "TCMID_web_completion_metadata_v1_0.json",
        "w", encoding="utf-8"
    ) as f:
        json.dump(metadata, f, ensure_ascii=False, indent=2)

    print("\nCompletion summary")
    print("------------------")
    print(summary.to_string(index=False))

    print("\nIMPORTANT:")
    print(
        "Do not run RWR yet. First upload this whole output folder so mapping "
        "coverage and database-induced target-space collapse can be re-audited."
    )
    print(f"\nCompleted workbook:\n{output_xlsx}")


if __name__ == "__main__":
    main()
