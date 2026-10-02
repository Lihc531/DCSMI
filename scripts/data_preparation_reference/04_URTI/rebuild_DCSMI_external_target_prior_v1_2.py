#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
rebuild_DCSMI_external_target_prior_v1_2.py

DCSMI Phase 1C -> Phase 2 preparation

Purpose:
    Rebuild frozen syndrome target prior from:

    syndrome herb weights
            +
    herb -> gene target dictionary

Result:
    syndrome_target_prior

    syndrome_id | target | target_weight


Concept:

Herb prior:
    syndrome -> herb weight

Herb target:
    herb -> gene symbol

Projection:

    syndrome
        |
        v
      herb
        |
        v
   gene target


Information boundary:
    - Uses frozen herb knowledge only
    - Does not use RNA-seq
    - Does not modify targets using disease results


Inputs:

1. Frozen syndrome herb prior Excel

Expected sheet:
    syndrome_herb_prior

Columns:
    syndrome_id
    herb_canonical
    herb_weight


2. Herb target long table CSV

Expected columns:
    Herb_ChineseName
    Gene_Symbol

Optional:
    Target_Name
    Uniprot_ID


Outputs:

DCSMI_external_WindHeat_WindCold_prior_v1_2.xlsx

Sheets:

syndrome_target_prior
syndrome_herb_prior
metadata

"""

from pathlib import Path
from datetime import datetime
import json

import pandas as pd


# =========================
# USER CONFIGURATION
# =========================

HERB_PRIOR_FILE = (
    r"E:\00project\Tanre Yongfei\重构"
    r"\DCSMI_external_WindHeat_WindCold_prior_v1_1.xlsx"
)

HERB_TARGET_FILE = (
    r"E:\00project\Tanre Yongfei\重构"
    r"\HERB_TARGET_LONG_TABLE.csv"
)

OUTPUT_DIR = (
    r"E:\00project\Tanre Yongfei\重构"
    r"\DCSMI_external_prior_v1_2"
)


HERB_PRIOR_SHEET = "syndrome_herb_weights"


# =========================
# FUNCTIONS
# =========================

def normalize_text(x):

    if pd.isna(x):
        return ""

    return str(x).strip()


def main():

    out = Path(OUTPUT_DIR)
    out.mkdir(
        parents=True,
        exist_ok=True
    )


    herb_prior = pd.read_excel(
        HERB_PRIOR_FILE,
        sheet_name=HERB_PRIOR_SHEET
    )


    herb_target = pd.read_csv(
        HERB_TARGET_FILE
    )


    required_prior = [
        "syndrome_id",
        "herb_canonical",
        "herb_weight"
    ]

    required_target = [
        "Herb_ChineseName",
        "Gene_Symbol"
    ]


    for c in required_prior:
        if c not in herb_prior.columns:
            raise ValueError(
                f"Missing herb prior column: {c}"
            )


    for c in required_target:
        if c not in herb_target.columns:
            raise ValueError(
                f"Missing herb target column: {c}"
            )


    herb_prior["herb_canonical"] = (
        herb_prior["herb_canonical"]
        .apply(normalize_text)
    )

    herb_target["Herb_ChineseName"] = (
        herb_target["Herb_ChineseName"]
        .apply(normalize_text)
    )

    herb_target["Gene_Symbol"] = (
        herb_target["Gene_Symbol"]
        .apply(normalize_text)
    )


    # herb -> gene projection

    projected = herb_prior.merge(
        herb_target,
        left_on="herb_canonical",
        right_on="Herb_ChineseName",
        how="inner"
    )


    if len(projected) == 0:
        raise ValueError(
            "No herb-target overlap. Check herb naming."
        )


    # combine herb weights with target mapping

    target_prior = (
        projected
        .assign(
            target_weight_product =
            projected["herb_weight"]
            .astype(float)
        )
        .groupby(
            ["syndrome_id", "Gene_Symbol"],
            as_index=False
        )["target_weight_product"]
        .sum()
        .rename(
            columns={
                "Gene_Symbol":"target",
                "target_weight_product":
                    "target_weight"
            }
        )
    )


    # normalize within syndrome

    target_prior["target_weight"] = (
        target_prior
        .groupby("syndrome_id")
        ["target_weight"]
        .transform(
            lambda x:
            x / x.sum()
        )
    )


    outfile = (
        out /
        "DCSMI_external_WindHeat_WindCold_prior_v1_2.xlsx"
    )


    with pd.ExcelWriter(
        outfile,
        engine="openpyxl"
    ) as writer:

        target_prior.to_excel(
            writer,
            sheet_name="syndrome_target_prior",
            index=False
        )

        herb_prior.to_excel(
            writer,
            sheet_name="syndrome_herb_prior",
            index=False
        )

        pd.DataFrame([{

            "version":
                "DCSMI_external_prior_v1_2",

            "RNAseq_used":
                "NO",

            "target_reweighted_by_disease":
                "NO",

            "timestamp":
                datetime.now().isoformat()

        }]).to_excel(
            writer,
            sheet_name="metadata",
            index=False
        )


    pd.DataFrame(
        {
            "target":
                sorted(
                    target_prior["target"]
                    .unique()
                )
        }
    ).to_csv(
        out/"syndrome_target_list_v1_2.csv",
        index=False
    )


    with open(
        out/"prior_rebuild_summary.json",
        "w",
        encoding="utf-8"
    ) as f:

        json.dump(
            {
                "version":
                    "v1_2",

                "target_number":
                    int(
                        target_prior["target"]
                        .nunique()
                    ),

                "record_number":
                    int(
                        len(target_prior)
                    )

            },
            f,
            ensure_ascii=False,
            indent=2
        )


    print(outfile)


if __name__ == "__main__":
    main()
