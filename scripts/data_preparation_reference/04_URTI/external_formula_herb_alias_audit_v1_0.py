#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
external_formula_herb_alias_audit_v1_0.py

DCSMI Phase 1A
Audit raw herb aliases before canonical herb mapping.

Purpose:
Detect naming variants that may represent the same canonical herb.

Frozen principles:
1. Alias normalization is performed BEFORE herb-target projection.
2. Do not merge biologically/pharmacologically distinct herb entities.
3. 荆芥穗 must remain independent from 荆芥.

Current frozen alias candidates:
银花 -> 金银花
芥穗 -> 荆芥穗
生甘草 -> 甘草
苦桔梗 -> 桔梗

No patient transcriptomics.
No target data.
No DCSMI output used.
"""

from pathlib import Path
import pandas as pd
import json
from datetime import datetime


# =========================
# USER CONFIGURATION
# =========================

INPUT_XLSX = r"E:\00project\Tanre Yongfei\重构\DCSMI_external_WindHeat_WindCold_formula_composition_verified_candidate_v1_1.xlsx"

OUTPUT_DIR = r"E:\00project\Tanre Yongfei\重构\DCSMI_external_formula_alias_audit_v1_0"


# =========================
# FROZEN ALIAS TABLE
# =========================

ALIAS_RULES = {

    "银花": {
        "canonical": "金银花",
        "rule": "common_short_name"
    },

    "芥穗": {
        "canonical": "荆芥穗",
        "rule": "abbreviation_expansion"
    },

    "生甘草": {
        "canonical": "甘草",
        "rule": "processing_modifier_removed"
    },

    "苦桔梗": {
        "canonical": "桔梗",
        "rule": "descriptor_removed"
    },

    # Explicit protection rule
    "荆芥": {
        "canonical": "荆芥",
        "rule": "independent_entity"
    },

    "荆芥穗": {
        "canonical": "荆芥穗",
        "rule": "independent_entity"
    }
}


def audit_alias(name):

    if name in ALIAS_RULES:
        item = ALIAS_RULES[name]

        return (
            item["canonical"],
            item["rule"],
            "REVIEW" if name in ["荆芥", "荆芥穗"] else "PASS"
        )

    return (
        name,
        "no_alias_rule",
        "PASS"
    )


def main():

    print("=" * 70)
    print("External formula herb alias audit v1.0")
    print("=" * 70)

    out = Path(OUTPUT_DIR)
    out.mkdir(parents=True, exist_ok=True)

    df = pd.read_excel(
        INPUT_XLSX,
        sheet_name="formula_composition_verified"
    )

    required = [
        "formula_id",
        "formula_name_cn",
        "herb_cn"
    ]

    for c in required:
        if c not in df.columns:
            raise ValueError(
                f"Missing required column: {c}"
            )


    records = []

    for _, row in df.iterrows():

        canonical, rule, status = audit_alias(
            row["herb_cn"]
        )

        records.append({

            "formula_id":
                row["formula_id"],

            "formula_name_cn":
                row["formula_name_cn"],

            "herb_raw":
                row["herb_cn"],

            "alias_canonical_candidate":
                canonical,

            "alias_rule":
                rule,

            "audit_status":
                status

        })


    result = pd.DataFrame(records)

    output = (
        out /
        "external_formula_herb_alias_audit_v1_0.xlsx"
    )


    with pd.ExcelWriter(
        output,
        engine="openpyxl"
    ) as writer:

        result.to_excel(
            writer,
            sheet_name="alias_audit",
            index=False
        )

        pd.DataFrame([
            {
                "version":
                    "external_formula_herb_alias_audit_v1_0",

                "patient_omics_used":
                    "NO",

                "target_data_used":
                    "NO",

                "critical_rule":
                    "荆芥穗 != 荆芥"

            }
        ]).to_excel(
            writer,
            sheet_name="metadata",
            index=False
        )


    summary = {

        "version":
            "external_formula_herb_alias_audit_v1_0",

        "timestamp":
            datetime.now().isoformat(),

        "rows":
            len(result),

        "review_rows":
            int(
                (result.audit_status=="REVIEW")
                .sum()
            ),

        "alias_rules":
            len(ALIAS_RULES)

    }


    with open(
        out/"alias_audit_summary_v1_0.json",
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
