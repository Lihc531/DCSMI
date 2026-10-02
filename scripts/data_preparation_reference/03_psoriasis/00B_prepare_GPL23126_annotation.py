#!/usr/bin/env python3
"""Convert GEO GPL23126 platform annotation ZIP/TXT to frozen two-column annotation.
Exact single-symbol rule:
- parse gene_assignment entries separated by ' /// '
- each entry is expected as accession // SYMBOL // description // cytoband // gene_id
- retain a probe only when all non-placeholder assignments resolve to exactly ONE unique gene symbol
- no synonym/fuzzy mapping.
"""
import argparse, zipfile, io, re, pandas as pd
from pathlib import Path

PLACE={"","---","NA","N/A","nan"}

def symbols_from_assignment(x):
    syms=[]
    for entry in str(x).split(" /// "):
        parts=[p.strip() for p in entry.split(" // ")]
        if len(parts)>=2 and parts[1] not in PLACE:
            s=parts[1].strip()
            # reject obvious multi-symbol fields
            if any(z in s for z in [";","/","|"]): return []
            syms.append(s)
    return sorted(set(syms))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--gpl",required=True,help="GPL23126-131.zip or extracted GPL23126-131.txt")
    ap.add_argument("--out",required=True)
    a=ap.parse_args()
    p=Path(a.gpl)
    if p.suffix.lower()==".zip":
        z=zipfile.ZipFile(p)
        names=[n for n in z.namelist() if n.lower().endswith(".txt")]
        if len(names)!=1: raise ValueError(f"Expected one TXT in ZIP, found {names}")
        fh=io.TextIOWrapper(z.open(names[0]),encoding="utf-8",errors="replace")
    else:
        fh=open(p,encoding="utf-8",errors="replace")
    # skip # metadata lines; pandas accepts the first non-comment row as header
    df=pd.read_csv(fh,sep="\t",comment="#",dtype=str,low_memory=False)
    if not {"ID","gene_assignment"}.issubset(df.columns):
        raise ValueError(f"Required GPL columns not found. First columns: {list(df.columns)[:20]}")
    rows=[]
    multi=zero=0
    for pid,ga in zip(df["ID"],df["gene_assignment"].fillna("")):
        sy=symbols_from_assignment(ga)
        if len(sy)==1: rows.append((pid,sy[0]))
        elif len(sy)>1: multi+=1
        else: zero+=1
    out=pd.DataFrame(rows,columns=["ID_REF","Gene.Symbol"]).drop_duplicates()
    if out.ID_REF.duplicated().any(): raise ValueError("Duplicate retained ID_REF after conversion")
    Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    out.to_csv(a.out,index=False,sep="\t")
    print({"platform_rows":len(df),"retained_exact_single_symbol":len(out),
           "excluded_zero_or_unparseable":zero,"excluded_multi_symbol":multi,
           "unique_symbols":out["Gene.Symbol"].nunique()})
    print("Wrote",a.out)
if __name__=="__main__": main()
