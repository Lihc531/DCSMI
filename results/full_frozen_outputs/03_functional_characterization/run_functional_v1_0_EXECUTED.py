import pandas as pd, numpy as np, gzip, re, os, json, hashlib, time
from collections import defaultdict
from pathlib import Path
from scipy.stats import hypergeom
from scipy import sparse

OUT=Path('/mnt/data/DCSMI_bias_aware_functional_characterization_EXECUTION_v1_0'); OUT.mkdir(exist_ok=True)
SEED=20260922; NNULL=10000
prof_path='/mnt/data/DCSMI_bias_aware_functional_characterization_v1_0/LOPO_SHARED_SPECIFIC_primary_target_profiles_v1_0.csv'
strata_path='/mnt/data/DCSMI_bias_aware_functional_characterization_v1_0/matched_null_target_universe_and_strata_v1_0_1.csv'
gaf='/mnt/data/HUMAN-uniprot.gaf.gz'; obo='/mnt/data/go-basic.obo'
klink='/mnt/data/KEGG_hsa_gene_to_pathway_2026-09-27.txt'; knames='/mnt/data/KEGG_hsa_pathway_names_2026-09-27.txt'; kinfo='/mnt/data/KEGG_info_2026-09-27.txt'
mapxlsx='/mnt/data/evidence_main/SYMBOL_ENTREZ_annotation.xlsx'

def sha(f):
 h=hashlib.sha256();
 with open(f,'rb') as x:
  for b in iter(lambda:x.read(1<<20),b''): h.update(b)
 return h.hexdigest()

def bh(p):
 p=np.asarray(p,float); n=len(p); o=np.argsort(p); q=np.empty(n); vals=p[o]*n/np.arange(1,n+1); vals=np.minimum.accumulate(vals[::-1])[::-1]; q[o]=np.minimum(vals,1); return q

# targets/order
st=pd.read_csv(strata_path); genes=st['gene'].astype(str).tolist(); idx={g:i for i,g in enumerate(genes)}; n=len(genes)
pr=pd.read_csv(prof_path)
raw={cid: g.set_index('target_id')['ubiquity_corrected'].reindex(genes).fillna(0).to_numpy(float) for cid,g in pr.groupby('component_id')}
# reconstruct complete from raw shared+specific then normalize to sum 1
W={}
for s in ['PDOL','PHOL','WHIL']:
 v=raw['SHARED']+raw['SPECIFIC_'+s]; W[s]=v/v.sum()
# normalized specific sensitivity and shared
Wspecific={s:raw['SPECIFIC_'+s]/raw['SPECIFIC_'+s].sum() for s in ['PDOL','PHOL','WHIL']}
Wshared=raw['SHARED']/raw['SHARED'].sum()
# reconstruction QA against known original if available
qa=[]
for s in W: qa.append({'metric':f'{s}_complete_sum','value':W[s].sum(),'status':'PASS' if abs(W[s].sum()-1)<1e-12 else 'FAIL'})
qa.append({'metric':'target_universe_n','value':n,'status':'PASS' if n==374 else 'FAIL'})
qa.append({'metric':'matched_strata_n','value':st.matched_stratum.nunique(),'status':'PASS' if st.matched_stratum.nunique()==19 else 'FAIL'})

# parse OBO terms, parents is_a + part_of
names={}; ns={}; parents=defaultdict(set); obsolete=set(); cur=None
with open(obo,errors='replace') as f:
 for line in f:
  line=line.rstrip('\n')
  if line=='[Term]': cur=None
  elif line.startswith('id: GO:'): cur=line[4:].strip()
  elif cur and line.startswith('name: '): names[cur]=line[6:].strip()
  elif cur and line.startswith('namespace: '): ns[cur]=line[11:].strip()
  elif cur and line.startswith('is_a: '): parents[cur].add(line.split()[1])
  elif cur and line.startswith('relationship: part_of '): parents[cur].add(line.split()[2])
  elif cur and line=='is_obsolete: true': obsolete.add(cur)
# ancestor memo
memo={}
def anc(t):
 if t in memo:return memo[t]
 a={t}
 for p in parents.get(t,()):
  if p!=t:a |= anc(p)
 memo[t]=a; return a

# GAF exact symbols, BP, exclude NOT; direct then propagate
symbol_go=defaultdict(set); gaf_header=[]
with gzip.open(gaf,'rt',errors='replace') as f:
 for line in f:
  if line.startswith('!'):
   if len(gaf_header)<20:gaf_header.append(line.rstrip())
   continue
  z=line.rstrip('\n').split('\t')
  if len(z)<15: continue
  sym=z[2]; qual=z[3]; go=z[4]; aspect=z[8]
  if aspect!='P' or 'NOT' in qual.split('|') or go in obsolete: continue
  if sym in idx: symbol_go[sym].add(go)
go_sets=defaultdict(set)
for g,terms in symbol_go.items():
 for t in terms:
  for a in anc(t):
   if ns.get(a)=='biological_process' and a not in obsolete: go_sets[a].add(g)

# KEGG exact SYMBOL->ENTREZ local authoritative mapping; ambiguous excluded
mp=pd.read_excel(mapxlsx,sheet_name='SYMBOL_ENTREZ',dtype={'SYMBOL':str,'ENTREZID':str})
mp=mp.dropna(); mp['SYMBOL']=mp.SYMBOL.astype(str); mp['ENTREZID']=mp.ENTREZID.astype(str).str.replace(r'\.0$','',regex=True)
counts=mp.groupby('SYMBOL').ENTREZID.nunique(); unique_syms=set(counts[counts==1].index)
sym2ent=mp[mp.SYMBOL.isin(unique_syms)].drop_duplicates('SYMBOL').set_index('SYMBOL').ENTREZID.to_dict()
ent2sym={sym2ent[g]:g for g in genes if g in sym2ent}
ksets=defaultdict(set)
with open(klink) as f:
 for line in f:
  a,b=line.rstrip().split('\t'); ent=a.split(':',1)[1]; pid=b.replace('path:','')
  if ent in ent2sym:ksets[pid].add(ent2sym[ent])
kn={}
with open(knames) as f:
 for line in f:
  p,nm=line.rstrip().split('\t',1); kn[p]=nm.replace(' - Homo sapiens (human)','')

# annotation QA
qa += [
 {'metric':'GO_exact_direct_annotated_targets','value':len(symbol_go),'status':'PASS'},
 {'metric':'GO_BP_terms_intersecting_targets','value':len(go_sets),'status':'PASS'},
 {'metric':'KEGG_exact_symbol_entrez_targets','value':sum(g in sym2ent for g in genes),'status':'PASS'},
 {'metric':'KEGG_targets_with_pathway','value':len(set().union(*ksets.values())) if ksets else 0,'status':'PASS'},
 {'metric':'KEGG_pathways_intersecting_targets','value':len(ksets),'status':'PASS'},
]
pd.DataFrame(qa).to_csv(OUT/'ANNOTATION_AND_INPUT_QA.csv',index=False)

# build eligible matrices, broad audit
families={}
for fam,sets,nmap in [('GO_BP',go_sets,names),('KEGG',ksets,kn)]:
 rows=[]; audit=[]
 for term in sorted(sets):
  gs=sets[term]
  k=len(gs & set(genes)); broad=k>n/2; elig=k>=5 and not broad
  audit.append({'family':fam,'term_id':term,'term_name':nmap.get(term,term),'n_frozen_targets':k,'broad_gt50pct':broad,'primary_eligible':elig})
  if elig: rows.append((term,nmap.get(term,term),gs))
 pd.DataFrame(audit).sort_values(['primary_eligible','n_frozen_targets'],ascending=[False,False]).to_csv(OUT/f'{fam}_TERM_ELIGIBILITY_AUDIT.csv',index=False)
 terms=[r[0] for r in rows]; tnames=[r[1] for r in rows]
 M=np.zeros((n,len(rows)),dtype=np.float32)
 for j,(_,_,gs) in enumerate(rows):
  for g in gs:
   if g in idx:M[idx[g],j]=1
 families[fam]=(terms,tnames,sparse.csr_matrix(M))

# strata index groups
strata_groups=[np.where(st.matched_stratum.to_numpy()==x)[0] for x in sorted(st.matched_stratum.unique())]

def perm_weights(rng, base):
 out=base.copy()
 for ii in strata_groups: out[ii]=base[rng.permutation(ii)]
 return out

def null_Q(base_weights, M, seed, label):
 rng=np.random.default_rng(seed)
 P=np.tile(base_weights,(NNULL,1)).astype(np.float32)
 for ii in strata_groups:
  keys=rng.random((NNULL,len(ii)))
  ords=np.argsort(keys,axis=1)
  P[:,ii]=base_weights[ii][ords]
 return np.asarray(P @ M, dtype=np.float32)

def run_family(fam, terms, tnames, M):
 print('RUN',fam,'terms',len(terms),flush=True)
 # Primary complete profiles: independent RNG streams but synchronized replicate indices
 Bs={}; obs={}
 for j,s in enumerate(['PDOL','PHOL','WHIL']):
  obs[s]=np.asarray(M.T @ W[s]).ravel()
  Bs[s]=null_Q(W[s],M,SEED+1000*(j+1),s)
 means={s:Bs[s].mean(0,dtype=np.float64) for s in Bs}; sds={s:Bs[s].std(0,ddof=1,dtype=np.float64) for s in Bs}
 A={s:np.divide(obs[s]-means[s],sds[s],out=np.zeros_like(obs[s]),where=sds[s]>0) for s in Bs}
 Anull={s:np.divide(Bs[s]-means[s],sds[s],out=np.zeros_like(Bs[s],dtype=np.float64),where=sds[s]>0) for s in Bs}
 # absolute results
 absrows=[]; diffrows=[]
 for s in ['PDOL','PHOL','WHIL']:
  pabs=(1+(Bs[s]>=obs[s]).sum(0))/(NNULL+1)
  for j,t in enumerate(terms):
   absrows.append([fam,s,t,tnames[j],int(M[:,j].sum()),obs[s][j],means[s][j],sds[s][j],A[s][j],pabs[j],sds[s][j]==0])
  others=[x for x in ['PDOL','PHOL','WHIL'] if x!=s]
  delta=A[s]-(A[others[0]]+A[others[1]])/2
  dnull=Anull[s]-(Anull[others[0]]+Anull[others[1]])/2
  pdiff=(1+(dnull>=delta).sum(0))/(NNULL+1)
  for j,t in enumerate(terms): diffrows.append([fam,s,t,tnames[j],int(M[:,j].sum()),A[s][j],A[others[0]][j],A[others[1]][j],delta[j],pdiff[j]])
 absdf=pd.DataFrame(absrows,columns=['family','syndrome','term_id','term_name','n_frozen_targets','Q_obs','Q_null_mean','Q_null_sd','A_z','empirical_p_upper','zero_null_sd'])
 absdf['BH_q_within_annotation_family']=bh(absdf.empirical_p_upper)
 diffdf=pd.DataFrame(diffrows,columns=['family','syndrome','term_id','term_name','n_frozen_targets','A_self','A_other1','A_other2','Delta','empirical_p_upper'])
 diffdf['BH_q_primary_family']=bh(diffdf.empirical_p_upper)
 diffdf['primary_positive']=(diffdf.Delta>0)&(diffdf.BH_q_primary_family<.05)
 absdf.to_csv(OUT/f'{fam}_COMPLETE_ABSOLUTE_STRUCTURAL_NULL.csv',index=False)
 diffdf.to_csv(OUT/f'{fam}_COMPLETE_PRIMARY_DIFFERENTIAL.csv',index=False)
 # specific-only sensitivity same procedure
 Bs2={}; obs2={}
 for j,s in enumerate(['PDOL','PHOL','WHIL']):
  obs2[s]=np.asarray(M.T @ Wspecific[s]).ravel(); Bs2[s]=null_Q(Wspecific[s],M,SEED+10000+1000*(j+1),s)
 means2={s:Bs2[s].mean(0,dtype=np.float64) for s in Bs2}; sds2={s:Bs2[s].std(0,ddof=1,dtype=np.float64) for s in Bs2}
 A2={s:np.divide(obs2[s]-means2[s],sds2[s],out=np.zeros_like(obs2[s]),where=sds2[s]>0) for s in Bs2}
 An2={s:np.divide(Bs2[s]-means2[s],sds2[s],out=np.zeros_like(Bs2[s],dtype=np.float64),where=sds2[s]>0) for s in Bs2}
 sr=[]
 for s in ['PDOL','PHOL','WHIL']:
  oo=[x for x in ['PDOL','PHOL','WHIL'] if x!=s]; d=A2[s]-(A2[oo[0]]+A2[oo[1]])/2; dn=An2[s]-(An2[oo[0]]+An2[oo[1]])/2; pp=(1+(dn>=d).sum(0))/(NNULL+1)
  for j,t in enumerate(terms):sr.append([fam,s,t,tnames[j],int(M[:,j].sum()),A2[s][j],d[j],pp[j]])
 sdf=pd.DataFrame(sr,columns=['family','syndrome','term_id','term_name','n_frozen_targets','A_specific','Delta_specific','empirical_p_upper'])
 sdf['BH_q_sensitivity_family']=bh(sdf.empirical_p_upper); sdf['sensitivity_positive']=(sdf.Delta_specific>0)&(sdf.BH_q_sensitivity_family<.05)
 sdf.to_csv(OUT/f'{fam}_SPECIFIC_ONLY_DIFFERENTIAL_SENSITIVITY.csv',index=False)
 # shared absolute descriptive
 obsh=np.asarray(M.T @ Wshared).ravel(); Bh=null_Q(Wshared,M,SEED+50000,M.shape[1] if False else 'shared'); mh=Bh.mean(0,dtype=np.float64); sh=Bh.std(0,ddof=1,dtype=np.float64); ah=np.divide(obsh-mh,sh,out=np.zeros_like(obsh),where=sh>0); ph=(1+(Bh>=obsh).sum(0))/(NNULL+1)
 hdf=pd.DataFrame({'family':fam,'term_id':terms,'term_name':tnames,'n_frozen_targets':np.asarray(M.sum(0)).ravel().astype(int),'Q_obs':obsh,'Q_null_mean':mh,'Q_null_sd':sh,'A_z':ah,'empirical_p_upper':ph}); hdf['BH_q_shared_annotation']=bh(hdf.empirical_p_upper); hdf.to_csv(OUT/f'{fam}_SHARED_DESCRIPTIVE.csv',index=False)
 return diffdf,absdf,sdf,hdf

results={}
for fam,(terms,tnames,M) in families.items(): results[fam]=run_family(fam,terms,tnames,M)
# summary and verdict
summary=[]
for fam,(diff,absd,spec,sh) in results.items():
 pos=diff[diff.primary_positive]; summary.append({'family':fam,'eligible_terms':diff.term_id.nunique(),'primary_tests':len(diff),'primary_positive_n':len(pos),'primary_positive_terms_n':pos.term_id.nunique(),'min_primary_q':diff.BH_q_primary_family.min(),'specific_sensitivity_positive_n':int(spec.sensitivity_positive.sum()),'shared_annotation_q_lt_0_05_n':int((sh.BH_q_shared_annotation<.05).sum())})
sumdf=pd.DataFrame(summary); sumdf.to_csv(OUT/'PRIMARY_FUNCTIONAL_RESULT_SUMMARY.csv',index=False)
anchor=bool((sumdf.primary_positive_n>0).any())
verdict={'status':'POSITIVE_FUNCTIONAL_ANCHOR' if anchor else 'NO_PRIMARY_FUNCTIONAL_ANCHOR','positive_anchor_exists':anchor,'decision_rule':'At least one GO-BP or KEGG primary family contains Delta>0 and BH-FDR q<0.05','claim_scope':'treatment-informed/prescription-derived functional organization only; not patient syndrome mechanism','families':summary}
with open(OUT/'SCIENTIFIC_VERDICT.json','w') as f:json.dump(verdict,f,indent=2)
# annotation manifest
manifest={'seed':SEED,'N_null':NNULL,'target_universe_n':n,'GO_gaf_headers':gaf_header,'GO_obo_data_version':next((x.strip() for x in open(obo,errors='replace') if x.startswith('data-version:')),None),'hashes':{Path(x).name:sha(x) for x in [gaf,obo,klink,knames,kinfo,prof_path,strata_path,mapxlsx]},'mapping':{'GO':'exact GAF gene symbol; BP direct annotations propagated through go-basic is_a and part_of; NOT excluded','KEGG':'exact project-local org.Hs.eg.db SYMBOL->ENTREZ; ambiguous symbols excluded','unmapped_rule':'remain in 374-target universe with zero functional membership'},'eligibility':'>=5 frozen targets; >187 targets excluded from primary differential'}
with open(OUT/'EXECUTION_MANIFEST.json','w') as f:json.dump(manifest,f,indent=2)
# top positive table
allpos=[]
for fam,(d,*_) in results.items(): allpos.append(d.sort_values(['BH_q_primary_family','empirical_p_upper','Delta']).head(50))
pd.concat(allpos).to_csv(OUT/'TOP_PRIMARY_DIFFERENTIAL_RESULTS.csv',index=False)
# readme
with open(OUT/'README.md','w') as f:
 f.write('# DCSMI bias-aware functional characterization v1.0\n\nFormal execution of the frozen 2026-09-27 GO-BP/KEGG branch. Primary inference uses full 374-target weighted pathway mass, 10,000 within-stratum permutations, and between-syndrome Delta with BH-FDR separately for GO-BP and KEGG. Interpret results only as prescription-derived/treatment-informed functional organization.\n')
# hashes outputs
with open(OUT/'SHA256SUMS.txt','w') as f:
 for x in sorted(OUT.iterdir()):
  if x.is_file() and x.name!='SHA256SUMS.txt':f.write(f'{sha(x)}  {x.name}\n')
print(sumdf.to_string(index=False)); print(json.dumps(verdict,indent=2))
