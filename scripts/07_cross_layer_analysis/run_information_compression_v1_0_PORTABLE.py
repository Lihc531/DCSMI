from pathlib import Path
from collections import defaultdict
import numpy as np, pandas as pd, gzip, hashlib, json, time

import os
PROJECT_ROOT=Path(os.environ.get('DCSMI_PROJECT_ROOT', Path(__file__).resolve().parents[2]))
OUT=Path(os.environ.get('DCSMI_CROSSLAYER_OUTPUT', PROJECT_ROOT/'outputs'/'07_cross_layer_analysis')); OUT.mkdir(parents=True,exist_ok=True)
TP=PROJECT_ROOT/'results'/'full_frozen_outputs'/'01_target_projection'
P12=PROJECT_ROOT/'results'/'full_frozen_outputs'/'02_topology_null'
FUNC=PROJECT_ROOT/'data'/'frozen_inputs'/'functional'
GAF=Path(os.environ.get('DCSMI_GAF', PROJECT_ROOT/'external_data'/'HUMAN-uniprot.gaf.gz'))
OBO=Path(os.environ.get('DCSMI_GO_OBO', PROJECT_ROOT/'external_data'/'go-basic.obo'))
KLINK=Path(os.environ.get('DCSMI_KEGG_LINK', PROJECT_ROOT/'external_data'/'KEGG_hsa_gene_to_pathway_2026-09-27.txt'))
KNAMES=Path(os.environ.get('DCSMI_KEGG_NAMES', PROJECT_ROOT/'external_data'/'KEGG_hsa_pathway_names_2026-09-27.txt'))
MAP=PROJECT_ROOT/'data'/'frozen_inputs'/'functional'/'SYMBOL_ENTREZ_annotation.xlsx'
SEED=20260922; N=10000; CHAINS=4; PER=N//CHAINS; CODES=['PDOL','PHOL','WHIL']

def sha(p):
 h=hashlib.sha256();
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def norm(P):
 P=np.asarray(P,float); s=P.sum(1,keepdims=True)
 if np.any(s<=0): raise ValueError('zero mass')
 return P/s
def gjsd(P):
 P=norm(P); M=P.mean(0); z=0.
 for p in P:
  m=p>0; z+=np.sum(p[m]*np.log2(p[m]/M[m]))
 return z/len(P)
def trade(nb,rng):
 n=len(nb); a=int(rng.integers(n)); b=int(rng.integers(n-1)); b += (b>=a)
 A=nb[a]; B=nb[b]; oa=A-B; ob=B-A
 if not oa or not ob:return False
 pool=np.array(list(oa|ob),dtype=int); rng.shuffle(pool); na=len(oa); common=A&B
 nb[a]=common|set(pool[:na].tolist()); nb[b]=common|set(pool[na:].tolist()); return True
def trades(nb,rng,k):
 e=0
 for _ in range(k): e+=trade(nb,rng)
 return e
def profiles(nb,coef,idf,nt):
 S=np.zeros((3,nt)); active=np.any(coef!=0,axis=0)
 for h,ts in enumerate(nb):
  if not active[h] or not ts:continue
  jj=np.fromiter(ts,dtype=int); S[:,jj]+=coef[:,h,None]
 S*=idf[None,:]; return norm(S)

# frozen graph and universes
pairs=pd.read_csv(TP/'canonical_herb_target_pairs_v1_1.csv')
hcol='target_mapping_name' if 'target_mapping_name' in pairs else pairs.columns[0]
tcol='target_id' if 'target_id' in pairs else pairs.columns[1]
herbs=sorted(pairs[hcol].astype(str).unique()); targets=sorted(pairs[tcol].astype(str).unique())
hi={h:i for i,h in enumerate(herbs)}; ti={t:i for i,t in enumerate(targets)}
nb=[set() for _ in herbs]
for h,t in pairs[[hcol,tcol]].astype(str).itertuples(index=False):nb[hi[h]].add(ti[t])
rowdeg=np.array([len(x) for x in nb]); coldeg=np.zeros(len(targets),int)
for x in nb:
 for j in x: coldeg[j]+=1

# frozen herb weights + projection coefficients
hc=pd.read_csv(P12/'phase12A_frozen_LOPO_herb_coefficients_v1_0.csv')
H=np.zeros((3,len(herbs))); coef=np.zeros_like(H)
for si,s in enumerate(CODES):
 d=hc[hc.syndrome_code==s].groupby('target_mapping_name',as_index=False).agg(primary_weight_W=('primary_weight_W','sum'),degree_normalized_coefficient=('degree_normalized_coefficient','sum'))
 for _,r in d.iterrows():
  h=str(r.target_mapping_name)
  if h in hi: H[si,hi[h]]=r.primary_weight_W; coef[si,hi[h]]=r.degree_normalized_coefficient
H=norm(H)
D_H=gjsd(H)

ub=pd.read_csv(TP/'target_database_target_ubiquity_v1_1.csv').set_index('target_id').reindex(targets)
idf=ub.target_idf.to_numpy(float)
Tobs=profiles(nb,coef,idf,len(targets)); D_T=gjsd(Tobs)
# exact crosscheck frozen target
fr=pd.read_csv(TP/'LOPO_FROZEN_primary_target_weights_for_RWR_v1_1.csv')
Fm=fr.pivot(index='target_id',columns='syndrome_code',values='ubiquity_corrected_sum1').reindex(targets).fillna(0)
F=np.vstack([Fm[s].to_numpy() for s in CODES]); maxerr=float(np.max(np.abs(Tobs-F)))

# GO parse
names={}; ns={}; parents=defaultdict(set); obsolete=set(); cur=None
with open(OBO,errors='replace') as f:
 for line in f:
  line=line.rstrip()
  if line=='[Term]':cur=None
  elif line.startswith('id: GO:'):cur=line[4:].strip()
  elif cur and line.startswith('name: '):names[cur]=line[6:].strip()
  elif cur and line.startswith('namespace: '):ns[cur]=line[11:].strip()
  elif cur and line.startswith('is_a: '):parents[cur].add(line.split()[1])
  elif cur and line.startswith('relationship: part_of '):parents[cur].add(line.split()[2])
  elif cur and line=='is_obsolete: true':obsolete.add(cur)
memo={}
def anc(x):
 if x in memo:return memo[x]
 a={x}
 for p in parents.get(x,()):
  if p!=x:a|=anc(p)
 memo[x]=a;return a
sg=defaultdict(set); tset=set(targets)
with gzip.open(GAF,'rt',errors='replace') as f:
 for line in f:
  if line.startswith('!'):continue
  z=line.rstrip().split('\t')
  if len(z)<9:continue
  sym,qual,go,aspect=z[2],z[3],z[4],z[8]
  if sym in tset and aspect=='P' and 'NOT' not in qual.split('|') and go not in obsolete:sg[sym].add(go)
gosets=defaultdict(set)
for g,tt in sg.items():
 for t in tt:
  for a in anc(t):
   if ns.get(a)=='biological_process' and a not in obsolete:gosets[a].add(g)
go_aud=pd.read_csv(FUNC/'GO_BP_TERM_ELIGIBILITY_AUDIT.csv'); goterms=sorted(go_aud.loc[go_aud.primary_eligible.astype(str).str.lower().isin(['true','1']),'term_id'].astype(str))
Mgo=np.zeros((len(targets),len(goterms)),np.float32)
for j,t in enumerate(goterms):
 for g in gosets.get(t,set()):
  if g in ti:Mgo[ti[g],j]=1

# KEGG
mp=pd.read_excel(MAP,sheet_name='SYMBOL_ENTREZ',dtype=str).dropna(); mp.ENTREZID=mp.ENTREZID.str.replace(r'\.0$','',regex=True)
c=mp.groupby('SYMBOL').ENTREZID.nunique(); us=set(c[c==1].index); s2e=mp[mp.SYMBOL.isin(us)].drop_duplicates('SYMBOL').set_index('SYMBOL').ENTREZID.to_dict(); e2s={s2e[g]:g for g in targets if g in s2e}
ksets=defaultdict(set)
with open(KLINK) as f:
 for line in f:
  a,b=line.rstrip().split('\t'); e=a.split(':')[1]; p=b.replace('path:','')
  if e in e2s:ksets[p].add(e2s[e])
kaud=pd.read_csv(FUNC/'KEGG_TERM_ELIGIBILITY_AUDIT.csv'); kterms=sorted(kaud.loc[kaud.primary_eligible.astype(str).str.lower().isin(['true','1']),'term_id'].astype(str))
Mk=np.zeros((len(targets),len(kterms)),np.float32)
for j,t in enumerate(kterms):
 for g in ksets.get(t,set()):
  if g in ti:Mk[ti[g],j]=1

def func_repr(T,M):return norm(T@M)
GOobs=func_repr(Tobs,Mgo); Kobs=func_repr(Tobs,Mk); D_GO=gjsd(GOobs); D_K=gjsd(Kobs)
obs={'layer':'HERB','n_features':len(herbs),'gjsd':D_H,'retention_vs_previous':1.0,'retention_vs_herb':1.0}, {'layer':'TARGET','n_features':len(targets),'gjsd':D_T,'retention_vs_previous':D_T/D_H,'retention_vs_herb':D_T/D_H}, {'layer':'GO_BP','n_features':len(goterms),'gjsd':D_GO,'retention_vs_previous':D_GO/D_T,'retention_vs_herb':D_GO/D_H}, {'layer':'KEGG','n_features':len(kterms),'gjsd':D_K,'retention_vs_previous':D_K/D_T,'retention_vs_herb':D_K/D_H}
pd.DataFrame(obs).to_csv(OUT/'OBSERVED_CROSS_LAYER_DISCRIMINATION.csv',index=False)

# exact curveball null
burn=100*len(herbs); thin=5*len(herbs); ss=np.random.SeedSequence(SEED).spawn(CHAINS); rows=[]
for ch,cs in enumerate(ss,1):
 rng=np.random.default_rng(cs); x=[set(q) for q in nb]; trades(x,rng,burn)
 for it in range(1,PER+1):
  trades(x,rng,thin); T=profiles(x,coef,idf,len(targets)); dt=gjsd(T); dg=gjsd(func_repr(T,Mgo)); dk=gjsd(func_repr(T,Mk))
  rows.append((ch,it,(ch-1)*PER+it,dt,dg,dk,dt/D_H,dg/dt,dk/dt,dg/D_H,dk/D_H))
 print('chain',ch,'done',flush=True)
null=pd.DataFrame(rows,columns=['chain_id','iteration','global_index','target_gjsd','go_gjsd','kegg_gjsd','target_retention_vs_herb','go_retention_vs_target','kegg_retention_vs_target','go_retention_vs_herb','kegg_retention_vs_herb'])
null.to_csv(OUT/'EXACT_CURVEBALL_CROSS_LAYER_NULL_10000.csv.gz',index=False,compression='gzip')

def p_up(o,x):return (1+np.sum(np.asarray(x)>=o))/(len(x)+1)
def p_lo(o,x):return (1+np.sum(np.asarray(x)<=o))/(len(x)+1)
summary=[]
for layer,o,col in [('TARGET',D_T,'target_gjsd'),('GO_BP',D_GO,'go_gjsd'),('KEGG',D_K,'kegg_gjsd')]:
 x=null[col].to_numpy(); summary.append({'layer':layer,'observed_gjsd':o,'null_mean':x.mean(),'null_sd':x.std(ddof=1),'z':(o-x.mean())/x.std(ddof=1),'empirical_upper_p':p_up(o,x),'empirical_lower_p':p_lo(o,x)})
pd.DataFrame(summary).to_csv(OUT/'STRUCTURAL_NULL_CROSS_LAYER_SUMMARY.csv',index=False)
mon_go=D_H>D_T>D_GO; mon_k=D_H>D_T>D_K
# structurally calibrated functional specificity = function obs > null (upper p<.05); expected false from prior but prespecified
pgo=p_up(D_GO,null.go_gjsd); pk=p_up(D_K,null.kegg_gjsd)
if mon_go and mon_k and pgo>=.05 and pk>=.05: decision='GO_INFORMATION_COMPRESSION_FRAMEWORK'
elif (mon_go or mon_k): decision='PIVOT_INCONSISTENT_OR_PARTIAL_CONTRACTION'
else: decision='STOP_NO_COHERENT_CROSS_LAYER_CONTRACTION'
ver={'decision':decision,'observed':{'herb_gjsd':D_H,'target_gjsd':D_T,'go_gjsd':D_GO,'kegg_gjsd':D_K,'target_retention_vs_herb':D_T/D_H,'go_retention_vs_target':D_GO/D_T,'kegg_retention_vs_target':D_K/D_T},'monotonic':{'GO':bool(mon_go),'KEGG':bool(mon_k)},'structural_null_upper_p':{'target':p_up(D_T,null.target_gjsd),'GO':pgo,'KEGG':pk},'target_reconstruction_max_abs_error':maxerr,'interpretation_limit':'projection-associated contraction of syndrome discrimination; not mutual information and not proof of absent patient biology'}
with open(OUT/'MANUSCRIPT_DECISION.json','w') as f:json.dump(ver,f,indent=2)
qa=pd.DataFrame([['target_reconstruction_max_abs_error',maxerr,'PASS' if maxerr<1e-12 else 'FAIL'],['herb_count',len(herbs),'PASS'],['target_count',len(targets),'PASS' if len(targets)==374 else 'FAIL'],['go_terms',len(goterms),'PASS' if len(goterms)==1560 else 'FAIL'],['kegg_terms',len(kterms),'PASS' if len(kterms)==207 else 'FAIL'],['null_n',len(null),'PASS' if len(null)==10000 else 'FAIL']],columns=['metric','value','status']); qa.to_csv(OUT/'EXECUTION_QA.csv',index=False)
print(json.dumps(ver,indent=2))
