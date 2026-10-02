#!/usr/bin/env python3
"""Freeze independent psoriasis disease background from 1-2 standardized cohort gene-stat CSVs.
Each input must contain unique gene, logFC, pvalue from psoriasis-vs-control with NO DEG prefilter.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import norm,spearmanr

def sha(p):
 h=hashlib.sha256(); h.update(Path(p).read_bytes()); return h.hexdigest()
def load(p):
 x=pd.read_csv(p); req={'gene','logFC','pvalue'}
 if not req.issubset(x.columns): raise ValueError(f'{p}: missing {req-set(x.columns)}')
 if x.gene.duplicated().any(): raise ValueError(f'{p}: duplicated gene symbols')
 x=x.dropna(subset=['gene','logFC','pvalue']).copy(); x['pvalue']=x.pvalue.clip(1e-300,1); x['Z']=np.sign(x.logFC)*norm.ppf(1-x.pvalue/2); return x.set_index('gene')
def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--cohort',action='append',required=True); ap.add_argument('--outdir',required=True); a=ap.parse_args()
 if len(a.cohort)>2: raise ValueError('Frozen implementation accepts max 2 predeclared cohorts')
 out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True); xs=[load(p) for p in a.cohort]
 qa=[]
 if len(xs)==2:
  common=xs[0].index.intersection(xs[1].index); z1=xs[0].loc[common,'Z']; z2=xs[1].loc[common,'Z']; l1=xs[0].loc[common,'logFC']; l2=xs[1].loc[common,'logFC']
  rz=spearmanr(z1,z2).statistic; rl=spearmanr(l1,l2).statistic; dc=np.mean(np.sign(l1)==np.sign(l2)); qa=[len(common),rz,rl,dc]
  status='Tier_A' if (rz>0 and dc>0.50) else 'SECONDARY_STRESS_TEST'
  Z=(z1+z2)/np.sqrt(2); bg=pd.DataFrame({'gene':common,'Z_disease':Z.values})
 else:
  status='Tier_B'; bg=pd.DataFrame({'gene':xs[0].index,'Z_disease':xs[0].Z.values})
 bg.to_csv(out/'psoriasis_independent_disease_background_FROZEN.csv',index=False)
 pd.DataFrame([qa],columns=['shared_genes','signedZ_spearman','logFC_spearman','direction_concordance']).to_csv(out/'disease_background_QA.csv',index=False) if qa else None
 man={'status':status,'cohorts':[{ 'path':p,'sha256':sha(p)} for p in a.cohort],'rule':'signed Z; equal-weight Stouffer if two cohorts; no FDR/FC prefilter'}
 (out/'disease_background_FREEZE_MANIFEST.json').write_text(json.dumps(man,indent=2),encoding='utf-8'); print(man)
if __name__=='__main__': main()
