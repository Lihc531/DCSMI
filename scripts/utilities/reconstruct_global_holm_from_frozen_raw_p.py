#!/usr/bin/env python3
"""Reconstruct manuscript-facing global Holm adjustments from frozen raw endpoint P values.

This is a transparent post-processing utility reconstructed from the frozen multiplicity
rule. It is NOT claimed to be the original historical endpoint-execution script.
Family A: all 5 P1 tests (3 psoriasis + 2 URTI).
Family B: psoriasis P3 tests only; URTI P3 is NOT_IDENTIFIABLE / NOT_TESTED by frozen amendment.
"""
from pathlib import Path
import argparse
import pandas as pd
import numpy as np

def holm(p):
    p=np.asarray(p,float); m=len(p); order=np.argsort(p); adj=np.empty(m,float); running=0.0
    for rank,idx in enumerate(order):
        val=(m-rank)*p[idx]; running=max(running,val); adj[idx]=min(1.0,running)
    return adj

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--psoriasis',required=True)
    ap.add_argument('--urti',required=True)
    ap.add_argument('--out',required=True)
    a=ap.parse_args()
    ps=pd.read_csv(a.psoriasis); ur=pd.read_csv(a.urti)
    fam=[]
    for _,r in ps.iterrows(): fam.append({'context':'psoriasis','syndrome':r['syndrome'],'endpoint':'P1','raw_p':r['P1_empirical_p']})
    for _,r in ur.iterrows(): fam.append({'context':'URTI','syndrome':r['syndrome'],'endpoint':'P1','raw_p':r['P1_empirical_p']})
    p1=pd.DataFrame(fam); p1['Holm_p_global_family_A']=holm(p1.raw_p)
    p3=ps[['syndrome','P3_empirical_p']].copy(); p3.insert(0,'context','psoriasis'); p3['Holm_p_family_B']=holm(p3.P3_empirical_p)
    out=Path(a.out); out.mkdir(parents=True,exist_ok=True)
    p1.to_csv(out/'P1_global_Holm_family_A.csv',index=False)
    p3.to_csv(out/'P3_Holm_family_B_psoriasis_only.csv',index=False)
    print(p1.to_string(index=False)); print(); print(p3.to_string(index=False))
if __name__=='__main__': main()
