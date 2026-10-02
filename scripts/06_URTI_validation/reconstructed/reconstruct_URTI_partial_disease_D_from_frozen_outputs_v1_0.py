#!/usr/bin/env python3
"""
Forensic reconstruction of the identifiable portion of the URTI disease
Hallmark vector D from frozen outputs.

STATUS: RECONSTRUCTED_FORENSIC_UTILITY, not historical analysis code.

Given the reconstructed Aplus_s, frozen Mtilde_s, and frozen compatibility R_s:
    denom_s = R_s * ||Aplus_s||     (because ||D||_2 = 1)
    D_k = Mtilde_{s,k} / Aplus_{s,k} * denom_s, when Aplus_{s,k} > 0.

Hallmarks for which both WH and WC have Aplus=0 are not identifiable from these
outputs. They are emitted as NA. This utility quantifies exactly how much of D
can be recovered and why the historical P1-null generator cannot be fully
reconstructed from the currently supplied outputs alone.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd

# Import the replay implementation without duplicating the structural-null logic.
from reconstruct_DCSMI_URTI_endpoints_REPLAY_v1_0 import reconstruct_aplus

SYND = ["WH", "WC"]


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--prior',required=True)
    p.add_argument('--strata',required=True)
    p.add_argument('--hallmark',required=True)
    p.add_argument('--mtilde',required=True)
    p.add_argument('--endpoints',required=True)
    p.add_argument('--outdir',required=True)
    a=p.parse_args()
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    M=pd.read_csv(a.mtilde,index_col=0).loc[SYND]
    E=pd.read_csv(a.endpoints).set_index('syndrome').loc[SYND]
    Aplus,_=reconstruct_aplus(a.prior,a.strata,a.hallmark,list(M.columns))
    D=np.full(M.shape[1],np.nan,float)
    source=np.array(['UNIDENTIFIABLE_BOTH_APLUS_ZERO']*M.shape[1],object)
    conflict=[]
    for i,s in enumerate(SYND):
        den=float(E.loc[s,'compatibility_R'])*float(np.linalg.norm(Aplus[i]))
        mask=Aplus[i]>0
        vals=(M.loc[s].to_numpy(float)[mask]/Aplus[i,mask])*den
        for idx,val in zip(np.where(mask)[0],vals):
            if np.isfinite(D[idx]):
                conflict.append(abs(D[idx]-val))
                D[idx]=(D[idx]+val)/2.0
                source[idx]='IDENTIFIED_FROM_BOTH_SYNDROMES'
            else:
                D[idx]=val
                source[idx]=f'IDENTIFIED_FROM_{s}'
    known=np.isfinite(D)
    df=pd.DataFrame({'hallmark':M.columns,'D_reconstructed':D,'status':source})
    df.to_csv(out/'URTI_disease_Hallmark_D_PARTIAL_RECONSTRUCTED.csv',index=False)
    meta={
      'status':'FORENSIC_RECONSTRUCTION_NOT_HISTORICAL_CODE',
      'n_hallmarks_total':int(len(D)),
      'n_identifiable':int(known.sum()),
      'n_unidentifiable':int((~known).sum()),
      'known_D_l2_norm_squared':float(np.nansum(D**2)),
      'unidentified_D_l2_norm_squared_total':float(1.0-np.nansum(D**2)),
      'max_cross_syndrome_D_disagreement':float(max(conflict) if conflict else 0.0),
      'implication':'P1 observed/Mtilde can be audited, but the historical P1 null cannot be regenerated exactly without the missing full GSE63990-derived D vector because null Aplus may load Hallmarks that are zero in both observed syndrome Aplus profiles.'
    }
    (out/'URTI_partial_D_RECONSTRUCTION_NOTES.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
    print(json.dumps(meta,indent=2))

if __name__=='__main__': main()
