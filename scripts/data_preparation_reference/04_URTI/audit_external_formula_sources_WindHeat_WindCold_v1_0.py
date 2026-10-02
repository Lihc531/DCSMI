#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
audit_external_formula_sources_WindHeat_WindCold_v1_0.py

DCSMI v1.0 Phase 1A
External formula source audit for Wind-Heat / Wind-Cold syndrome priors.

Purpose:
Audit formula-level evidence BEFORE any herb-target projection.

This script DOES NOT:
- read patient transcriptomics
- read URTI RNA-seq
- construct molecular priors
- modify formulas based on DCSMI results

Input:
    DCSMI_external_WindHeat_WindCold_canonical_v1_0.xlsx

Output:
    external_formula_source_audit_WindHeat_WindCold_v1_0.xlsx
    external_formula_source_audit_summary_v1_0.json

"""

from pathlib import Path
import json
from datetime import datetime

import pandas as pd


# ============================================================
# USER CONFIGURATION
# ============================================================

INPUT_XLSX = r"E:\00project\Tanre Yongfei\重构\DCSMI_external_WindHeat_WindCold_canonical_v1_0.xlsx"

OUTPUT_DIR = r"E:\00project\Tanre Yongfei\重构\DCSMI_external_WindHeat_WindCold_source_audit_v1_0"


# ============================================================
# AUDIT RULES
# ============================================================

REQUIRED_FORMULA_FIELDS = [
    "formula_id",
    "formula_name_cn",
    "formula_name_en",
    "source_type",
    "source_reference"
]

REQUIRED_HERB_FIELDS = [
    "formula_id",
    "herb_cn"
]


# Formula evidence categories
ALLOWED_SOURCE_TYPES = {
    "classical_formula",
    "pharmacopoeia",
    "clinical_guideline",
    "peer_reviewed_literature"
}


# ============================================================
# FUNCTIONS
# ============================================================

def load_input(path):
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(
            f"Input file not found:\n{path}"
        )

    xls = pd.ExcelFile(path)

    formula = pd.read_excel(
        xls,
        sheet_name="formula_canonical"
    )

    herb = pd.read_excel(
        xls,
        sheet_name="formula_herb_long"
    )

    return formula, herb


def check_columns(df, required, name):
    missing = [
        x for x in required
        if x not in df.columns
    ]

    if missing:
        raise ValueError(
            f"{name} missing columns: {missing}\n"
            f"Available: {list(df.columns)}"
        )


def audit_formula_table(df):

    records=[]

    for _, row in df.iterrows():

        issues=[]

        for col in REQUIRED_FORMULA_FIELDS:
            if pd.isna(row.get(col)) or str(row.get(col)).strip()=="":
                issues.append(
                    f"missing_{col}"
                )

        source_type=str(
            row.get("source_type","")
        ).strip()

        if source_type not in ALLOWED_SOURCE_TYPES:
            issues.append(
                "source_type_not_verified"
            )

        records.append({

            "formula_id":
                row["formula_id"],

            "formula_name":
                row["formula_name_cn"],

            "source_type":
                source_type,

            "formula_audit_status":
                "PASS" if len(issues)==0 else "REVIEW",

            "issues":
                ";".join(issues)

        })

    return pd.DataFrame(records)


def audit_formula_herb(df):

    records=[]

    for fid, group in df.groupby("formula_id"):

        issues=[]

        herb_count = (
            group["herb_cn"]
            .dropna()
            .astype(str)
            .str.strip()
            .replace("", pd.NA)
            .dropna()
            .nunique()
        )

        if herb_count == 0:
            issues.append(
                "no_explicit_herb_composition"
            )

        if herb_count < 3:
            issues.append(
                "low_herb_count_review"
            )

        records.append({

            "formula_id":fid,

            "herb_count":
                herb_count,

            "herb_audit_status":
                "PASS" if len(issues)==0 else "REVIEW",

            "issues":
                ";".join(issues)

        })

    return pd.DataFrame(records)


def build_manual_review_table(formula_df):

    rows=[]

    for _, r in formula_df.iterrows():

        rows.append({

            "formula_id":
                r["formula_id"],

            "formula_name_cn":
                r["formula_name_cn"],

            "source_verified":
                "NO",

            "source_reference_final":
                "",

            "composition_verified":
                "NO",

            "dose_required":
                "YES",

            "modification_components_removed":
                "YES",

            "final_freeze_decision":
                "PENDING"

        })

    return pd.DataFrame(rows)


def main():

    print("="*70)
    print("DCSMI external formula source audit v1.0")
    print("Wind-Heat / Wind-Cold")
    print("="*70)


    out=Path(OUTPUT_DIR)
    out.mkdir(
        parents=True,
        exist_ok=True
    )


    formula, herb = load_input(
        INPUT_XLSX
    )


    check_columns(
        formula,
        REQUIRED_FORMULA_FIELDS,
        "formula_canonical"
    )

    check_columns(
        herb,
        REQUIRED_HERB_FIELDS,
        "formula_herb_long"
    )


    formula_audit = audit_formula_table(
        formula
    )

    herb_audit = audit_formula_herb(
        herb
    )

    review_table = build_manual_review_table(
        formula
    )


    output_xlsx = (
        out /
        "external_formula_source_audit_WindHeat_WindCold_v1_0.xlsx"
    )


    with pd.ExcelWriter(
        output_xlsx,
        engine="openpyxl"
    ) as writer:

        formula_audit.to_excel(
            writer,
            sheet_name="formula_audit",
            index=False
        )

        herb_audit.to_excel(
            writer,
            sheet_name="herb_composition_audit",
            index=False
        )

        review_table.to_excel(
            writer,
            sheet_name="manual_freeze_review",
            index=False
        )


    summary={

        "analysis":
            "external_formula_source_audit_WindHeat_WindCold_v1_0",

        "timestamp":
            datetime.now().isoformat(),

        "input":
            str(INPUT_XLSX),

        "formula_number":
            int(len(formula)),

        "formula_pass":
            int(
                (formula_audit.formula_audit_status=="PASS")
                .sum()
            ),

        "formula_review":
            int(
                (formula_audit.formula_audit_status=="REVIEW")
                .sum()
            ),

        "freeze_rule":
            "No target projection until manual source audit PASS"

    }


    with open(
        out /
        "external_formula_source_audit_summary_v1_0.json",
        "w",
        encoding="utf-8"
    ) as f:

        json.dump(
            summary,
            f,
            ensure_ascii=False,
            indent=2
        )


    print("\nAudit completed.")
    print(output_xlsx)
    print(json.dumps(
        summary,
        ensure_ascii=False,
        indent=2
    ))


if __name__=="__main__":
    main()
