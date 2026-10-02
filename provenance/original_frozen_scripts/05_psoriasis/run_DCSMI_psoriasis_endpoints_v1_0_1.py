#!/usr/bin/env python3
import argparse,gzip,json,hashlib
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import ttest_ind,norm,rankdata,spearmanr
from scipy.special import stdtr

SEED=20260922
N_NULL=10000
N_BOOT=2000

def sha(p):
 h=hashlib.sha256()
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
   z=line.rstrip().split('\t')
   rows.append([z[0].strip('"')]+[float(v) for v in z[1:]])
 return pd.DataFrame(rows,columns=hdr).set_index('ID_REF')

def hallmarks(path):
 h={}
 df=pd.read_csv(path,sep='\t')
 for _,r in df.iterrows(): h[r.hallmark]=set(str(r.genes).split(';'))
 return h


def welch_fast(A,B):
 n1=A.shape[1]; n2=B.shape[1]
 m1=np.nanmean(A,axis=1); m2=np.nanmean(B,axis=1)
 v1=np.nanvar(A,axis=1,ddof=1); v2=np.nanvar(B,axis=1,ddof=1)
 a=v1/n1; b=v2/n2; se=np.sqrt(a+b)
 t=np.divide(m1-m2,se,out=np.zeros_like(m1),where=se>0)
 df=np.divide((a+b)**2, a*a/(n1-1)+b*b/(n2-1), out=np.full_like(m1,np.inf), where=(a*a/(n1-1)+b*b/(n2-1))>0)
 p=2*stdtr(df,-np.abs(t))
 p=np.where(se==0,1.0,p)
 return m1-m2,p

def signed_z(fc,p):
 p=np.clip(np.asarray(p,float),1e-300,1.0)
 return np.sign(np.asarray(fc,float))*norm.ppf(1-p/2)

def percentile_module(z,genes,hall,eligible_names=None):
 ok=np.isfinite(z); zz=np.asarray(z)[ok]; gg=np.asarray(genes)[ok]
 rr=rankdata(np.abs(zz),method='average')/(len(zz)+1)
 mp=dict(zip(gg,rr))
 names=list(hall) if eligible_names is None else list(eligible_names)
 raw=[]; n=[]
 for k in names:
  x=[mp[g] for g in hall[k] if g in mp]
  raw.append(np.mean(x) if x else np.nan); n.append(len(x))
 raw=np.asarray(raw,float); n=np.asarray(n,int)
 elig=n>=15
 if elig.sum()<40: raise ValueError("fewer than 40 eligible Hallmarks")
 names=np.asarray(names)[elig]; raw=raw[elig]
 v=raw/np.sqrt(np.sum(raw**2))
 return names,v,n[elig]

def holm(p):
 p=np.asarray(p,float); m=len(p); order=np.argsort(p); adj=np.empty(m); running=0
 for rank,i in enumerate(order):
  val=(m-rank)*p[i]; running=max(running,val); adj[i]=min(1,running)
 return adj

def main():
 ap=argparse.ArgumentParser()
 for x in ['prior','strata','hallmark','disease','series','sample_manifest','annotation','outdir']:
  ap.add_argument('--'+x.replace('_','-'),required=True)
 a=ap.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
 rng=np.random.default_rng(SEED)

 # fixed universes
 H=hallmarks(a.hallmark)
 strata=pd.read_csv(a.strata)
 if len(strata)!=374 or strata.gene.nunique()!=374 or strata.matched_stratum.nunique()!=19:
  raise ValueError("strata universe invariant failed")
 genes374=strata.gene.astype(str).tolist()
 groups=[np.where(strata.matched_stratum.values==s)[0] for s in sorted(strata.matched_stratum.unique())]

 pr=pd.read_csv(a.prior)
 synd=['BH','BS','BD']; synd_long={'BH':'Blood_Heat','BS':'Blood_Stasis','BD':'Blood_Dryness'}
 T=np.zeros((3,374))
 for i,s in enumerate(synd):
  z=pr[pr.syndrome_id==s].set_index('target_id').target_weight
  T[i]=np.array([z.get(g,0.0) for g in genes374])
  if not np.isclose(T[i].sum(),1,atol=1e-8): raise ValueError("prior weight conservation failed")

 # disease D, Tier-B GSE147339
 dis=pd.read_csv(a.disease)
 Zd=signed_z(dis.logFC,dis.pvalue)
 dnames,D,_=percentile_module(Zd,dis.gene.astype(str).values,H)
 if len(dnames)<40: raise ValueError("disease Hallmark eligibility failed")
 # common primary Hallmarks are disease-eligible; O was already 50/50 eligible
 hidx={h:i for i,h in enumerate(dnames)}
 K=len(dnames)

 # target x Hallmark incidence
 Inc=np.zeros((374,K))
 for j,h in enumerate(dnames):
  gs=H[h]; Inc[:,j]=[g in gs for g in genes374]

 # SAME donor permutation across all syndromes in each replicate
 Qobs=T@Inc
 Qnull=np.empty((N_NULL,3,K))
 perms=np.empty((N_NULL,374),dtype=np.int16)
 for b in range(N_NULL):
  donor=np.arange(374)
  for ix in groups: donor[ix]=rng.permutation(ix)
  perms[b]=donor
  Qnull[b]=T[:,donor]@Inc

 mu=Qnull.mean(axis=0); sd=Qnull.std(axis=0,ddof=1)
 A=np.divide(Qobs-mu,sd,out=np.zeros_like(Qobs),where=sd>0)
 Aplus=np.maximum(A,0)
 nonzero=(Aplus>0).sum(axis=1)
 if np.any(nonzero<10): raise SystemExit("STOP: Aplus has fewer than 10 nonzero eligible Hallmarks")

 M=D[None,:]*Aplus
 Mtilde=M/M.sum(axis=1,keepdims=True)
 Anorm=Aplus/Aplus.sum(axis=1,keepdims=True)

 # observed O recomputed/frozen file equivalence source: patient expression
 an=pd.read_csv(a.annotation,sep='\t',dtype=str).fillna('')
 an=an[(an['Gene.Symbol'].str.strip()!='') & ~an['Gene.Symbol'].str.contains(r'[;/|]')].copy()
 an['Gene.Symbol']=an['Gene.Symbol'].str.strip()
 expr=read_series(a.series); sm=pd.read_csv(a.sample_manifest)
 if list(expr.columns)!=sm.sample_id.tolist(): raise ValueError("sample order mismatch")
 common=expr.index.intersection(an.ID_REF)
 x=expr.loc[common].copy()
 mp=an.drop_duplicates('ID_REF').set_index('ID_REF').loc[common,'Gene.Symbol']
 x['Gene.Symbol']=mp.values
 gene=x.groupby('Gene.Symbol').median(numeric_only=True)
 ggrp=sm.set_index('sample_id').loc[gene.columns,'group']
 ps=ggrp[ggrp!='Healthy']; gene=gene[ps.index]; ggrp=ps
 O=np.zeros((3,K))
 for i,s in enumerate(synd):
  sl=synd_long[s]; aa=gene.loc[:,ggrp==sl].to_numpy(); bb=gene.loc[:,ggrp!=sl].to_numpy()
  fc,pv=welch_fast(aa,bb)
  _,ov,_=percentile_module(signed_z(fc,pv),gene.index.astype(str).values,H,dnames)
  O[i]=ov

 # P1 observed and structural null
 p1obs=np.array([spearmanr(Mtilde[i],O[i]).statistic for i in range(3)])
 p1null=np.empty((N_NULL,3))
 Cnull=np.empty((N_NULL,3))
 # P3 observed
 corr_obs=np.array([[spearmanr(M[i],O[j]).statistic for j in range(3)] for i in range(3)])
 Cobs=np.array([corr_obs[i,i]-np.mean(np.delete(corr_obs[i],i)) for i in range(3)])
 for b in range(N_NULL):
  Ab=np.divide(Qnull[b]-mu,sd,out=np.zeros((3,K)),where=sd>0)
  Ap=np.maximum(Ab,0); Mb=D[None,:]*Ap
  for i in range(3):
   if Mb[i].sum()<=0:
    p1null[b,i]=-1.0; Cnull[b,i]=-2.0; continue
   mt=Mb[i]/Mb[i].sum()
   p1null[b,i]=spearmanr(mt,O[i]).statistic
   cc=np.array([spearmanr(Mb[i],O[j]).statistic for j in range(3)])
   Cnull[b,i]=cc[i]-np.mean(np.delete(cc,i))
 p1p=(1+(p1null>=p1obs).sum(axis=0))/(N_NULL+1)
 p3p=(1+(Cnull>=Cobs).sum(axis=0))/(N_NULL+1)
 p1holm=holm(p1p); p3holm=holm(p3p)

 # P2 patient bootstrap, stratified within the three psoriasis syndrome groups
 syndonly=np.array([spearmanr(Anorm[i],O[i]).statistic for i in range(3)])
 delta=p1obs-syndonly
 # Vectorized exact bootstrap implementation: resampling is represented by integer
 # multiplicities on the original patient columns. This is algebraically identical
 # to explicit repeated-column sampling and preserves group sizes.
 G=gene.to_numpy(float); gnames=gene.index.astype(str).values; ng=G.shape[0]
 Hmat=np.zeros((ng,K),float)
 for j,h in enumerate(dnames):
  gs=H[h]; Hmat[:,j]=[g in gs for g in gnames]
 Hcount=Hmat.sum(axis=0)
 if np.any(Hcount<15): raise ValueError('bootstrap Hallmark eligibility invariant failed')
 colpos={c:i for i,c in enumerate(gene.columns)}
 pos={s:np.array([colpos[c] for c in group_cols],int) for s,group_cols in
      {s:gene.columns[ggrp.values==synd_long[s]].tolist() for s in synd}.items()}
 # Freeze all bootstrap draws before computing statistics.
 draws={}
 for s in synd:
  n=len(pos[s]); C=np.zeros((N_BOOT,n),dtype=np.int16)
  for b in range(N_BOOT):
   ix=rng.choice(n,size=n,replace=True); C[b]=np.bincount(ix,minlength=n)
  draws[s]=C
 boot=np.empty((N_BOOT,3))
 # fixed ranks for fast Spearman against each bootstrap O
 rM=np.array([rankdata(Mtilde[i],method='average') for i in range(3)])
 rA=np.array([rankdata(Anorm[i],method='average') for i in range(3)])
 def rowcorr_fixed(rfixed,Y):
  ry=rankdata(Y,axis=1,method='average')
  xc=rfixed-rfixed.mean(); yc=ry-ry.mean(axis=1,keepdims=True)
  return (yc@xc)/(np.sqrt((yc*yc).sum(axis=1))*np.sqrt((xc*xc).sum()))
 batch=50
 for i,s in enumerate(synd):
  others=[t for t in synd if t!=s]
  n1=len(pos[s]); n2=sum(len(pos[t]) for t in others)
  for lo in range(0,N_BOOT,batch):
   hi=min(N_BOOT,lo+batch); C1=draws[s][lo:hi]
   S1=C1@G[:,pos[s]].T; SS1=C1@(G[:,pos[s]]**2).T
   C2=np.concatenate([draws[t][lo:hi] for t in others],axis=1)
   P2=np.concatenate([pos[t] for t in others])
   S2=C2@G[:,P2].T; SS2=C2@(G[:,P2]**2).T
   m1=S1/n1; m2=S2/n2
   v1=np.maximum((SS1-S1*S1/n1)/(n1-1),0); v2=np.maximum((SS2-S2*S2/n2)/(n2-1),0)
   aa=v1/n1; bb=v2/n2; se=np.sqrt(aa+bb)
   tt=np.divide(m1-m2,se,out=np.zeros_like(m1),where=se>0)
   den=aa*aa/(n1-1)+bb*bb/(n2-1)
   df=np.divide((aa+bb)**2,den,out=np.full_like(m1,np.inf),where=den>0)
   pv=2*stdtr(df,-np.abs(tt)); pv=np.where(se==0,1.0,pv)
   zz=signed_z(m1-m2,pv)
   rr=rankdata(np.abs(zz),axis=1,method='average')/(ng+1)
   raw=(rr@Hmat)/Hcount
   ob=raw/np.sqrt((raw*raw).sum(axis=1,keepdims=True))
   boot[lo:hi,i]=rowcorr_fixed(rM[i],ob)-rowcorr_fixed(rA[i],ob)
 ci=np.quantile(boot,[.025,.975],axis=0)

 rows=[]
 for i,s in enumerate(synd):
  rows.append(dict(syndrome=s,syndrome_label=synd_long[s],n_hallmarks=K,
    Aplus_nonzero=int(nonzero[i]),compatibility_R=float(np.dot(D,Aplus[i])/(np.linalg.norm(D)*np.linalg.norm(Aplus[i]))),
    P1_spearman=float(p1obs[i]),P1_empirical_p=float(p1p[i]),P1_Holm_p=float(p1holm[i]),
    syndrome_only_spearman=float(syndonly[i]),P2_delta=float(delta[i]),
    P2_boot_CI_low=float(ci[0,i]),P2_boot_CI_high=float(ci[1,i]),
    P2_support=bool(ci[0,i]>0),
    P3_selectivity=float(Cobs[i]),P3_empirical_p=float(p3p[i]),P3_Holm_p=float(p3holm[i])))
 pd.DataFrame(rows).to_csv(out/'DCSMI_psoriasis_primary_endpoints.csv',index=False)
 pd.DataFrame(Mtilde,index=synd,columns=dnames).to_csv(out/'DCSMI_psoriasis_Mtilde.csv')
 pd.DataFrame(O,index=synd,columns=dnames).to_csv(out/'DCSMI_psoriasis_O_recomputed.csv')
 pd.DataFrame({'syndrome':np.repeat(synd,N_BOOT),'delta_boot':boot.T.reshape(-1)}).to_csv(out/'P2_bootstrap_distribution.csv',index=False)

 manifest={'seed':SEED,'N_null':N_NULL,'N_boot':N_BOOT,'n_targets':374,'n_strata':19,'n_hallmarks':K,
  'inputs_sha256':{k:sha(getattr(a,k)) for k in ['prior','strata','hallmark','disease','series','sample_manifest','annotation']}}
 (out/'ENDPOINT_RUN_MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
 print(pd.DataFrame(rows).to_string(index=False))

if __name__=='__main__': main()
