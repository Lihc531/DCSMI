#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
freeze_external_formula_herb_canonicalization_v1_0.py

DCSMI Phase 1A

Purpose:
Freeze the formula-herb canonicalization layer after:
1. formula composition verification
2. herb alias audit
3. canonical herb normalization

Input:
    DCSMI_external_WindHeat_WindCold_canonical_v1_2.xlsx

Output:
    DCSMI_external_WindHeat_WindCold_canonical_v1_2_FROZEN.xlsx

Frozen principles:
- No patient transcriptomic data.
- No target information.
- No DCSMI downstream result used.
- 荆芥穗 remains independent from 荆芥.

"""

from pathlib import Path
from datetime import datetime
import json
import pandas as pd


# =========================
# USER CONFIGURATION
# =========================

INPUT_XLSX = r"E:\00project\Tanre Yongfei\重构\DCSMI_external_WindHeat_WindCold_canonical_v1_2.xlsx"

OUTPUT_DIR = r"E:\00project\Tanre Yongfei\重构\DCSMI_external_WindHeat_WindCold_canonical_v1_2_FROZEN"


# =========================
# FREEZE VALIDATION
# =========================

REQUIRED_COLUMNS = [
    "formula_id",
    "formula_name_cn",
    "herb_raw",
    "herb_canonical"
]


PROTECTED_RULES = {
    "荆芥穗": "荆芥穗",
    "荆芥": "荆芥"
}


def validate(df):

    missing = [
        c for c in REQUIRED_COLUMNS
        if c not in df.columns
    ]

    if missing:
        raise ValueError(
            f"Missing required columns: {missing}"
        )

    errors = []

    for raw, canonical in PROTECTED_RULES.items():

        subset = df[df["herb_raw"] == raw]

        if len(subset) > 0:

            wrong = subset[
                subset["herb_canonical"] != canonical
            ]

            if len(wrong) > 0:
                errors.append(
                    f"{raw} incorrectly mapped"
                )

    if errors:
        raise RuntimeError(
            "Freeze blocked:\n" +
            "\n".join(errors)
        )


def main():

    print("=" * 70)
    print("Freeze external formula herb canonicalization v1.0")
    print("=" * 70)

    out = Path(OUTPUT_DIR)
    out.mkdir(
        parents=True,
        exist_ok=True
    )

    df = pd.read_excel(
        INPUT_XLSX,
        sheet_name="formula_herb_canonical"
    )

    validate(df)

    output = (
        out /
        "DCSMI_external_WindHeat_WindCold_canonical_v1_2_FROZEN.xlsx"
    )

    with pd.ExcelWriter(
        output,
        engine="openpyxl"
    ) as writer:

        df.to_excel(
            writer,
            sheet_name="frozen_formula_herb_mapping",
            index=False
        )

        pd.DataFrame([
            {
                "version":
                    "external_formula_herb_canonicalization_v1_0",

                "status":
                    "FROZEN",

                "freeze_time":
                    datetime.now().isoformat(),

                "patient_omics_used":
                    "NO",

                "target_information_used":
                    "NO",

                "critical_entity_rule":
                    "荆芥穗 != 荆芥"

            }
        ]).to_excel(
            writer,
            sheet_name="metadata",
            index=False
        )


    summary = {
        "status": "FROZEN",
        "rows": int(len(df)),
        "unique_formula": int(df.formula_id.nunique()),
        "unique_canonical_herbs": int(df.herb_canonical.nunique()),
        "patient_omics_used": False
    }

    with open(
        out / "freeze_summary_v1_0.json",
        "w",
        encoding="utf-8"
    ) as f:
        json.dump(
            summary,
            f,
            ensure_ascii=False,
            indent=2
        )

    print(output)


if __name__ == "__main__":
    main()
