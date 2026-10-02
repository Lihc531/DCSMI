#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
freeze_external_formula_composition_WindHeat_WindCold_v1_0.py

DCSMI v1.0 Phase 1A
Freeze audited formula compositions before herb-target projection.

Input:
    external_formula_source_audit_WindHeat_WindCold_v1_0.xlsx

Output:
    DCSMI_external_WindHeat_WindCold_canonical_v1_1.xlsx
    freeze_external_formula_composition_summary_v1_0.json

Important:
- This script does NOT use transcriptomic data.
- This script does NOT modify formulas according to DCSMI results.
- Manual verification fields are required before freeze.
"""

from pathlib import Path
import json
from datetime import datetime
import pandas as pd


# =========================
# USER CONFIGURATION
# =========================

INPUT_XLSX = r"E:\00project\Tanre Yongfei\重构\DCSMI_external_WindHeat_WindCold_source_audit_v1_0\external_formula_source_audit_WindHeat_WindCold_v1_0.xlsx"

OUTPUT_DIR = r"E:\00project\Tanre Yongfei\重构\DCSMI_external_WindHeat_WindCold_canonical_v1_1"


# =========================
# FREEZE RULES
# =========================

REQUIRED_COLUMNS = [
    "formula_id",
    "formula_name_cn",
    "source_verified",
    "composition_verified",
    "final_freeze_decision"
]


def require_columns(df, cols, name):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"{name} missing columns: {missing}\n"
            f"Available: {list(df.columns)}"
        )


def freeze_formula(review_df):

    freeze = review_df.copy()

    failed = []

    for _, row in freeze.iterrows():

        formula = row["formula_id"]

        checks = {
            "source_verified":
                str(row["source_verified"]).upper() == "YES",

            "composition_verified":
                str(row["composition_verified"]).upper() == "YES",

            "final_freeze_decision":
                str(row["final_freeze_decision"]).upper() == "PASS"
        }

        if not all(checks.values()):
            failed.append(formula)

    if failed:
        raise RuntimeError(
            "Freeze blocked. These formulas are not fully verified:\n"
            + ", ".join(failed)
        )

    freeze["freeze_status"] = "FROZEN"

    return freeze


def main():

    print("="*70)
    print("DCSMI external formula composition freeze v1.0")
    print("Wind-Heat / Wind-Cold")
    print("="*70)

    out = Path(OUTPUT_DIR)
    out.mkdir(parents=True, exist_ok=True)

    xls = pd.ExcelFile(INPUT_XLSX)

    review = pd.read_excel(
        xls,
        sheet_name="manual_freeze_review"
    )

    require_columns(
        review,
        REQUIRED_COLUMNS,
        "manual_freeze_review"
    )

    frozen_formula = freeze_formula(review)


    # create canonical workbook

    output_xlsx = (
        out /
        "DCSMI_external_WindHeat_WindCold_canonical_v1_1.xlsx"
    )

    with pd.ExcelWriter(
        output_xlsx,
        engine="openpyxl"
    ) as writer:

        frozen_formula.to_excel(
            writer,
            sheet_name="formula_master_frozen",
            index=False
        )

        review.to_excel(
            writer,
            sheet_name="audit_trace",
            index=False
        )


        pd.DataFrame(
            [{
                "version":
                    "DCSMI_external_WindHeat_WindCold_canonical_v1_1",

                "status":
                    "FROZEN",

                "freeze_time":
                    datetime.now().isoformat(),

                "next_step":
                    "external syndrome prior builder",

                "patient_omics_used":
                    "NO"
            }]
        ).to_excel(
            writer,
            sheet_name="metadata",
            index=False
        )


    summary = {

        "analysis":
            "freeze_external_formula_composition_WindHeat_WindCold_v1_0",

        "status":
            "FROZEN",

        "formula_number":
            int(len(frozen_formula)),

        "patient_omics_used":
            False,

        "output":
            str(output_xlsx)

    }


    with open(
        out /
        "freeze_external_formula_composition_summary_v1_0.json",
        "w",
        encoding="utf-8"
    ) as f:

        json.dump(
            summary,
            f,
            ensure_ascii=False,
            indent=2
        )


    print("Freeze completed:")
    print(output_xlsx)


if __name__ == "__main__":
    main()
