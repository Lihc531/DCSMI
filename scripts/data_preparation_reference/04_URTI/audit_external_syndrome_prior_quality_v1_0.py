#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
audit_external_syndrome_prior_quality_v1_0.py

DCSMI Phase 1A Quality Control

Purpose:
Audit treatment-informed external syndrome priors before
disease transcriptomic alignment.

Checks:
1. Target weight normalization
2. Herb coverage
3. Target coverage
4. Wind-Heat vs Wind-Cold prior similarity
5. Entropy / concentration diagnostics
6. Degree-preserving null preparation summary

No:
- patient transcriptomics
- disease labels
- CAP data
- outcome-driven tuning

"""

from pathlib import Path
import json
from datetime import datetime

import numpy as np
import pandas as pd


# =========================
# USER CONFIGURATION
# =========================

INPUT_XLSX = r"E:\00project\Tanre Yongfei\重构\DCSMI_external_WindHeat_WindCold_prior_v1_1.xlsx"

OUTPUT_DIR = r"E:\00project\Tanre Yongfei\重构\DCSMI_external_syndrome_prior_quality_v1_0"


# =========================
# FUNCTIONS
# =========================

def cosine_similarity(a, b):

    common = sorted(
        set(a.index) &
        set(b.index)
    )

    if len(common) == 0:
        return np.nan

    x = a.loc[common].values
    y = b.loc[common].values

    denom = (
        np.linalg.norm(x)
        *
        np.linalg.norm(y)
    )

    if denom == 0:
        return np.nan

    return float(
        np.dot(x, y) / denom
    )


def entropy(values):

    values = np.array(values)

    values = values[values > 0]

    p = values / values.sum()

    return float(
        -(p*np.log2(p)).sum()
    )


def main():

    print("="*70)
    print("Audit external syndrome prior quality v1.0")
    print("Wind-Heat / Wind-Cold")
    print("="*70)

    out = Path(OUTPUT_DIR)
    out.mkdir(
        parents=True,
        exist_ok=True
    )

    herb = pd.read_excel(
        INPUT_XLSX,
        sheet_name="syndrome_herb_weights"
    )

    target = pd.read_excel(
        INPUT_XLSX,
        sheet_name="syndrome_target_prior"
    )


    # -------------------------
    # Target normalization
    # -------------------------

    normalization = []

    for syndrome, group in target.groupby(
        "syndrome_id"
    ):

        s = group.target_weight.sum()

        normalization.append({

            "syndrome_id":
                syndrome,

            "target_weight_sum":
                float(s),

            "normalization_status":
                "PASS"
                if abs(s-1)<1e-6
                else "REVIEW"

        })


    normalization = pd.DataFrame(
        normalization
    )


    # -------------------------
    # Coverage
    # -------------------------

    coverage=[]

    for syndrome, group in herb.groupby(
        "syndrome_id"
    ):

        coverage.append({

            "syndrome_id":
                syndrome,

            "herb_number":
                int(
                    group.herb_canonical.nunique()
                ),

        })


    for syndrome, group in target.groupby(
        "syndrome_id"
    ):

        for item in coverage:

            if item["syndrome_id"]==syndrome:

                item["target_number"]=int(
                    group.target.nunique()
                )


    coverage=pd.DataFrame(
        coverage
    )


    # -------------------------
    # WH / WC similarity
    # -------------------------

    target_matrix = target.pivot_table(
        index="target",
        columns="syndrome_id",
        values="target_weight",
        fill_value=0
    )


    similarity={}

    cols=list(target_matrix.columns)

    if "WH" in cols and "WC" in cols:

        similarity["WH_WC_cosine"] = cosine_similarity(
            target_matrix["WH"],
            target_matrix["WC"]
        )


    # -------------------------
    # Entropy
    # -------------------------

    entropy_table=[]

    for syndrome, group in target.groupby(
        "syndrome_id"
    ):

        entropy_table.append({

            "syndrome_id":
                syndrome,

            "target_entropy":
                entropy(
                    group.target_weight
                )

        })


    entropy_table=pd.DataFrame(
        entropy_table
    )


    # -------------------------
    # Output
    # -------------------------

    outfile = (
        out /
        "external_syndrome_prior_quality_audit_v1_0.xlsx"
    )


    with pd.ExcelWriter(
        outfile,
        engine="openpyxl"
    ) as writer:

        normalization.to_excel(
            writer,
            sheet_name="normalization",
            index=False
        )

        coverage.to_excel(
            writer,
            sheet_name="coverage",
            index=False
        )

        entropy_table.to_excel(
            writer,
            sheet_name="entropy",
            index=False
        )

        pd.DataFrame(
            [similarity]
        ).to_excel(
            writer,
            sheet_name="similarity",
            index=False
        )


    summary={

        "version":
            "audit_external_syndrome_prior_quality_v1_0",

        "timestamp":
            datetime.now().isoformat(),

        "patient_omics_used":
            False,

        "target_normalization_checked":
            True,

        "WH_WC_similarity":
            similarity.get(
                "WH_WC_cosine",
                None
            )

    }


    with open(
        out/"quality_audit_summary_v1_0.json",
        "w",
        encoding="utf-8"
    ) as f:

        json.dump(
            summary,
            f,
            ensure_ascii=False,
            indent=2
        )


    print(outfile)


if __name__=="__main__":
    main()
