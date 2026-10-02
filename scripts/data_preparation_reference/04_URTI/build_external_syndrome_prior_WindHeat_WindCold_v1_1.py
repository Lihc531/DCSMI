# -*- coding: utf-8 -*-
"""
Fixed version:
build_external_syndrome_prior_WindHeat_WindCold_v1_1.py

Fix:
FROZEN_herb_to_target_mapping_dictionary_v1_1.csv
does not use column herb_canonical.

Detected columns:
representation
syndrome_code
canonical_name
target_mapping_name
mapping_method
primary_weight_W
article_support_m
mapped_to_target_db

Mapping:
canonical herb = canonical_name
target = target_mapping_name
weight = primary_weight_W
"""

from pathlib import Path
import pandas as pd
import json
from datetime import datetime


INPUT_FORMULA_HERB = r"E:\00project\Tanre Yongfei\重构\DCSMI_external_WindHeat_WindCold_canonical_v1_2_FROZEN.xlsx"

INPUT_HERB_TARGET = r"E:\00project\Tanre Yongfei\重构\FROZEN_herb_to_target_mapping_dictionary_v1_1.csv"

OUTPUT_DIR = r"E:\00project\Tanre Yongfei\重构\DCSMI_external_WindHeat_WindCold_prior_v1_1"


SMOOTH=0.5


def main():

    out=Path(OUTPUT_DIR)
    out.mkdir(parents=True,exist_ok=True)

    formula=pd.read_excel(
        INPUT_FORMULA_HERB,
        sheet_name="frozen_formula_herb_mapping"
    )

    ht=pd.read_csv(INPUT_HERB_TARGET)

    required=[
        "herb_canonical",
        "target",
        "weight"
    ]

    # adapt frozen dictionary format
    if "herb_canonical" not in ht.columns:
        ht["herb_canonical"]=ht["canonical_name"]

    if "target" not in ht.columns:
        ht["target"]=ht["target_mapping_name"]

    if "weight" not in ht.columns:
        ht["weight"]=ht["primary_weight_W"]


    formula["syndrome_id"]=(
        formula["formula_id"]
        .str.split("_")
        .str[0]
    )


    herb_weights=[]

    for syndrome,g in formula.groupby("syndrome_id"):

        n_formula=g.formula_id.nunique()

        counts=(
            g.groupby("herb_canonical")
            .formula_id
            .nunique()
        )

        for herb,c in counts.items():

            herb_weights.append({
                "syndrome_id":syndrome,
                "herb_canonical":herb,
                "herb_weight":
                    (c+SMOOTH)/(n_formula+SMOOTH)
            })


    herb_weights=pd.DataFrame(herb_weights)


    targets=[]

    for syndrome,g in herb_weights.groupby("syndrome_id"):

        total={}

        for _,r in g.iterrows():

            sub=ht[
                ht.herb_canonical==r.herb_canonical
            ]

            for _,t in sub.iterrows():

                total[t.target]=(
                    total.get(t.target,0)
                    +
                    r.herb_weight*t.weight
                )

        s=sum(total.values())

        for target,w in total.items():

            targets.append({
                "syndrome_id":syndrome,
                "target":target,
                "target_weight":w/s if s else 0
            })


    target_df=pd.DataFrame(targets)

    outfile=out/"DCSMI_external_WindHeat_WindCold_prior_v1_1.xlsx"

    with pd.ExcelWriter(outfile,engine="openpyxl") as writer:
        herb_weights.to_excel(
            writer,
            sheet_name="syndrome_herb_weights",
            index=False
        )
        target_df.to_excel(
            writer,
            sheet_name="syndrome_target_prior",
            index=False
        )

        pd.DataFrame([{
            "version":"v1.1",
            "patient_omics_used":"NO",
            "fix":"adapted herb-target dictionary schema"
        }]).to_excel(
            writer,
            sheet_name="metadata",
            index=False
        )

    print(outfile)


if __name__=="__main__":
    main()
