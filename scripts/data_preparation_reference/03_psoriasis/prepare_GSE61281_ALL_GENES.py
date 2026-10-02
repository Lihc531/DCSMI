#!/usr/bin/env python3
import argparse,gzip,csv,pandas as pd,numpy as np
from scipy.stats import ttest_ind
def main():
 p=argparse.ArgumentParser()
 p.add_argument("--series",required=True); p.add_argument("--annotation",required=True)
 p.add_argument("--manifest",required=True); p.add_argument("--out",required=True)
 a=p.parse_args()
 ann=pd.read_csv(a.annotation,sep=None,engine="python",dtype=str)
 # Accept official GEO names from GPL4133/GPL6480-derived annotation.
 probe=next((c for c in ["NAME","ID","ID_REF"] if c in ann.columns),None)
 sym=next((c for c in ["GENE_SYMBOL","Gene Symbol","Gene.Symbol"] if c in ann.columns),None)
 if not probe or not sym: raise ValueError("Need probe-name and gene-symbol columns")
 ann=ann[[probe,sym]].dropna(); ann.columns=["ID_REF","gene"]
 ann["ID_REF"]=ann.ID_REF.str.strip(); ann["gene"]=ann.gene.str.strip()
 ann=ann[(ann.gene!="") & (~ann.gene.str.contains(r"[;/|]",regex=True))]
 ann=ann.drop_duplicates()
 # only unambiguous one-symbol probe mappings
 ann=ann[ann.groupby("ID_REF").gene.transform("nunique")==1].drop_duplicates("ID_REF")
 with gzip.open(a.series,"rt",errors="replace") as h:
  for line in h:
   if line.startswith("!series_matrix_table_begin"): break
  df=pd.read_csv(h,sep="\t",dtype={0:str})
 df=df.rename(columns={df.columns[0]:"ID_REF"})
 df["ID_REF"]=df.ID_REF.str.strip('"')
 df=df.merge(ann,on="ID_REF",how="inner")
 man=pd.read_csv(a.manifest)
 use=man[man.primary_use.isin(["Psoriasis","Control"])]
 cols=use.GSM.tolist()
 # multiple probes per gene collapsed by median BEFORE testing
 gx=df.groupby("gene",as_index=False)[cols].median()
 ps=use.loc[use.primary_use.eq("Psoriasis"),"GSM"].tolist()
 co=use.loc[use.primary_use.eq("Control"),"GSM"].tolist()
 A=gx[ps].to_numpy(float); B=gx[co].to_numpy(float)
 fc=np.nanmean(A,1)-np.nanmean(B,1)
 pv=ttest_ind(A,B,axis=1,equal_var=False,nan_policy="omit").pvalue
 out=pd.DataFrame({"gene":gx.gene,"logFC":fc,"pvalue":pv})
 out=out[np.isfinite(out.logFC)&np.isfinite(out.pvalue)&(out.pvalue>0)&(out.pvalue<=1)]
 out.to_csv(a.out,index=False)
 print({"mapped_genes":len(gx),"output_genes":len(out),"PsC":len(ps),"controls":len(co)})
if __name__=="__main__": main()
