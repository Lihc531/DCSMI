#!/usr/bin/env python3
"""Prepare the exact composition CSV consumed by frozen script 01.
This helper does NOT infer herbs or verify sources. It only:
1) reads the existing 5 reference formulas from the workbook,
2) reads the user-filled missing-3 sheet,
3) requires non-empty evidence fields and verification_status=VERIFIED for used rows,
4) uses herb_standardized_cn when provided, otherwise herb_cn_source,
5) writes syndrome_id,formula,herb_cn,verification_status.
"""
import argparse, pandas as pd
from pathlib import Path

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--xlsx",required=True)
    ap.add_argument("--out",required=True)
    a=ap.parse_args()
    xls=pd.ExcelFile(a.xlsx)
    need={"TO_FILL_missing_3","REFERENCE_existing_5"}
    if not need.issubset(xls.sheet_names):
        raise ValueError(f"Workbook missing sheets: {need-set(xls.sheet_names)}")
    old=pd.read_excel(a.xlsx,sheet_name="REFERENCE_existing_5",dtype=str).fillna("")
    new=pd.read_excel(a.xlsx,sheet_name="TO_FILL_missing_3",dtype=str).fillna("")
    # Only nonblank herb rows are considered.
    new=new[new["herb_cn_source"].str.strip()!=""].copy()
    if new.empty: raise SystemExit("STOP: no herb rows entered in TO_FILL_missing_3")
    required=["source_title","source_type","source_identifier_or_URL","source_location","verification_status"]
    bad=[]
    for i,r in new.iterrows():
        miss=[c for c in required if not str(r.get(c,"")).strip()]
        if str(r.get("verification_status","")).strip()!="VERIFIED": miss.append("verification_status!=VERIFIED")
        if miss: bad.append((i+2,r.get("formula",""),r.get("herb_cn_source",""),miss))
    if bad:
        print("Rows not ready:")
        for b in bad: print(b)
        raise SystemExit("STOP: every used missing-formula herb row must have source evidence and VERIFIED status.")
    # Existing reference formulas are carried forward as verified source-layer inputs.
    old2=pd.DataFrame({
        "syndrome_id":old["syndrome_id"].str.strip(),
        "formula":old["formula"].str.strip(),
        "herb_cn":old["herb_cn"].str.strip(),
        "verification_status":"VERIFIED"
    })
    herb=new["herb_standardized_cn"].str.strip()
    herb=herb.where(herb!="",new["herb_cn_source"].str.strip())
    new2=pd.DataFrame({
        "syndrome_id":new["syndrome_id"].str.strip(),
        "formula":new["formula"].str.strip(),
        "herb_cn":herb,
        "verification_status":"VERIFIED"
    })
    out=pd.concat([old2,new2],ignore_index=True)
    out=out[(out.herb_cn!="")].drop_duplicates()
    Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    out.to_csv(a.out,index=False,encoding="utf-8-sig")
    print(out.groupby(["syndrome_id","formula"]).size().to_string())
    print(f"Wrote {len(out)} rows -> {a.out}")
if __name__=="__main__": main()
