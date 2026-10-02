# -*- coding: utf-8 -*-
"""
external_formula_herb_canonical_mapping_v1_0.py

DCSMI Phase 1A external formula canonical herb mapping.

Frozen rules:
生甘草 -> 甘草
芍药 -> 白芍
苦桔梗 -> 桔梗
荆芥穗 -> 荆芥穗 (retain independently)

No patient omics or target information are used.
"""

from pathlib import Path
import pandas as pd
import json
from datetime import datetime


INPUT_XLSX = r"E:\00project\Tanre Yongfei\重构\DCSMI_external_WindHeat_WindCold_formula_composition_verified_candidate_v1_1.xlsx"

OUTPUT_DIR = r"E:\00project\Tanre Yongfei\重构\DCSMI_external_WindHeat_WindCold_canonical_v1_2"


RULES = {
    "生甘草": ("甘草", "processing_modifier_removed"),
    "芍药": ("白芍", "classical_formula_normalization"),
    "苦桔梗": ("桔梗", "processing_descriptor_removed"),
    "荆芥穗": ("荆芥穗", "independent_retention"),
}


def normalize(name):
    if name in RULES:
        return RULES[name]
    return name, "exact_name_retained"


def main():

    out = Path(OUTPUT_DIR)
    out.mkdir(parents=True, exist_ok=True)

    df = pd.read_excel(
        INPUT_XLSX,
        sheet_name="formula_composition_verified"
    )

    records = []

    for _, row in df.iterrows():

        canonical, rule = normalize(row["herb_cn"])

        records.append({
            "formula_id": row["formula_id"],
            "formula_name_cn": row["formula_name_cn"],
            "herb_raw": row["herb_cn"],
            "herb_canonical": canonical,
            "mapping_rule": rule,
            "mapping_status": "FROZEN_RULE_APPLIED"
        })

    result = pd.DataFrame(records)

    output = out / "DCSMI_external_WindHeat_WindCold_canonical_v1_2.xlsx"

    with pd.ExcelWriter(output, engine="openpyxl") as writer:
        result.to_excel(
            writer,
            sheet_name="formula_herb_canonical",
            index=False
        )

        pd.DataFrame([{
            "version": "external_formula_herb_canonicalization_v1_0",
            "patient_omics_used": "NO",
            "target_information_used": "NO",
            "荆芥穗_policy": "separate_from_荆芥"
        }]).to_excel(
            writer,
            sheet_name="metadata",
            index=False
        )

    with open(out/"mapping_summary_v1_0.json","w",encoding="utf-8") as f:
        json.dump({
            "version": "external_formula_herb_canonical_mapping_v1_0",
            "timestamp": datetime.now().isoformat(),
            "rows": len(result)
        }, f, ensure_ascii=False, indent=2)

    print(output)


if __name__ == "__main__":
    main()
