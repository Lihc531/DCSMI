#!/usr/bin/env python3
"""Prepare GSE192867 full-expression gene matrix and frozen within-disease O profiles.
Requires a COMPLETE GPL23126 annotation TSV/CSV with ID_REF and Gene.Symbol columns.
Implementation decision frozen before P1/P2/P3: exact symbol mapping; probes mapping to >1 symbols are excluded;
multiple probes per symbol are collapsed by median expression BEFORE differential testing; primary contrast is s vs all other psoriasis syndromes;
Welch two-sample t test, no covariates; logFC=mean(s)-mean(other); signed Z from two-sided p; no DEG prefilter.
"""
import argparse,gzip,csv,json,hashlib
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import ttest_ind,norm,rankdata

def sha(p):
 h=hashlib.sha256();
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()
def read_series(p):
 op=gzip.open if str(p).endswith('.gz') else open
 with op(p,'rt',errors='replace') as f:
  for line in f:
   if line.rstrip()=='!series_matrix_table_begin': break
  hdr=[x.strip('"') for x in next(f).rstrip().split('\t')]
  rows=[]
  for line in f:
   if line.rstrip()=='!series_matrix_table_end': break
   z=line.rstrip().split('\t'); rows.append([z[0].strip('"')]+[float(v) for v in z[1:]])
 return pd.DataFrame(rows,columns=hdr).set_index('ID_REF')
def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--series',required=True); ap.add_argument('--sample-manifest',required=True); ap.add_argument('--annotation',required=True); ap.add_argument('--hallmark-gmt',required=True); ap.add_argument('--outdir',required=True); a=ap.parse_args()
 out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
 sep='\t' if str(a.annotation).lower().endswith(('.tsv','.txt')) else ','
 an=pd.read_csv(a.annotation,sep=sep,dtype=str).fillna('')
 if not {'ID_REF','Gene.Symbol'}.issubset(an.columns): raise ValueError('Annotation needs ID_REF and Gene.Symbol')
 # exact single-symbol only
 an=an[(an['Gene.Symbol'].str.strip()!='') & ~an['Gene.Symbol'].str.contains(r'[;/|]')].copy(); an['Gene.Symbol']=an['Gene.Symbol'].str.strip()
 expr=read_series(a.series); sm=pd.read_csv(a.sample_manifest)
 if list(expr.columns)!=sm.sample_id.tolist(): raise ValueError('Sample order mismatch')
 common=expr.index.intersection(an.ID_REF); x=expr.loc[common].copy(); mp=an.drop_duplicates('ID_REF').set_index('ID_REF').loc[common,'Gene.Symbol']
 x['Gene.Symbol']=mp.values; gene=x.groupby('Gene.Symbol').median(numeric_only=True)
 groups=sm.set_index('sample_id').loc[gene.columns,'group']
 ps=groups[groups!='Healthy']; gene=gene[ps.index]; groups=ps
 # Hallmarks
 hall={}
 with open(a.hallmark_gmt) as f:
  for line in f:
   z=line.rstrip().split('\t'); hall[z[0]]=set(z[2:])
 allstats=[]; O=[]
 for s in ['Blood_Heat','Blood_Stasis','Blood_Dryness']:
  A=gene.loc[:,groups==s].to_numpy(); B=gene.loc[:,groups!=s].to_numpy()
  tt,p=ttest_ind(A,B,axis=1,equal_var=False,nan_policy='omit'); lfc=np.nanmean(A,axis=1)-np.nanmean(B,axis=1)
  p=np.clip(p,1e-300,1.0); Z=np.sign(lfc)*norm.ppf(1-p/2)
  st=pd.DataFrame({'gene':gene.index,'logFC':lfc,'pvalue':p,'Zobs':Z,'syndrome':s}); allstats.append(st)
  valid=np.isfinite(Z); genes=np.asarray(gene.index)[valid]; rz=rankdata(np.abs(Z[valid]),method='average')/(valid.sum()+1)
  rmap=dict(zip(genes,rz)); vals=[]
  for h,gs in hall.items():
   gg=[g for g in gs if g in rmap]; vals.append([s,h,len(gg),np.mean([rmap[g] for g in gg]) if gg else np.nan])
  od=pd.DataFrame(vals,columns=['syndrome','hallmark','n_measured','O_raw'])
  eligible=od.n_measured>=15
  if eligible.sum()<40: raise SystemExit(f'STOP: {s} has only {eligible.sum()} eligible Hallmarks')
  den=np.sqrt(np.nansum(od.loc[eligible,'O_raw']**2)); od['eligible']=eligible; od['O_l2']=np.where(eligible,od.O_raw/den,np.nan); O.append(od)
 pd.concat(allstats).to_csv(out/'GSE192867_within_disease_gene_stats_FROZEN.csv',index=False)
 pd.concat(O).to_csv(out/'GSE192867_observed_Hallmark_O_FROZEN.csv',index=False)
 pd.DataFrame({'gene':gene.index}).to_csv(out/'GSE192867_gene_universe_FROZEN.csv',index=False)
 manifest={'implementation':'median probe-to-gene collapse; Welch t test; syndrome vs all other psoriasis syndromes; no healthy controls; no DEG filter','input_sha256':{p:sha(p) for p in [a.series,a.sample_manifest,a.annotation,a.hallmark_gmt]},'n_genes':len(gene)}
 (out/'GSE192867_O_FREEZE_MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
 print(manifest)
if __name__=='__main__': main()
