#!/usr/bin/env python3
"""Build exact DCSMI v1.0 external psoriasis priors from a VERIFIED composition CSV.
No fuzzy mapping; no specificity weighting; no patient omics.
"""
import argparse, hashlib, json
from pathlib import Path
import numpy as np, pandas as pd

FORMULAS={
 'BH':['Liangxue Jiedu Tang','Liangxue Huoxue Fang','Xijiao Dihuang Tang'],
 'BS':['Huoxue Sanyu Xiaoyin Tang','Taohong Siwu Tang'],
 'BD':['Yangxue Jiedu Tang','Danggui Yinzi','Yangxue Runfu Yin']}

def sha(p):
 h=hashlib.sha256(); h.update(Path(p).read_bytes()); return h.hexdigest()
def norm_formula(x):
 x=str(x).lower()
 for z in [' (犀角地黄汤)',' (桃红四物汤加减)',' (养血润肤饮)',' modification']: x=x.replace(z,'')
 return ' '.join(x.split())
def main():
 ap=argparse.ArgumentParser()
 ap.add_argument('--composition-csv',required=True,help='Columns: syndrome_id,formula,herb_cn,verification_status')
 ap.add_argument('--herb-target-pairs',required=True)
 ap.add_argument('--target-idf',required=True)
 ap.add_argument('--outdir',required=True)
 a=ap.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
 comp=pd.read_csv(a.composition_csv,dtype=str).fillna('')
 req={'syndrome_id','formula','herb_cn','verification_status'}
 if not req.issubset(comp.columns): raise ValueError(f'Missing columns {req-set(comp.columns)}')
 comp=comp[(comp.herb_cn.str.strip()!='') & (comp.verification_status=='VERIFIED')].copy()
 comp['formula_norm']=comp.formula.map(norm_formula)
 expected={s:{norm_formula(f) for f in fs} for s,fs in FORMULAS.items()}
 for s,efs in expected.items():
  got=set(comp.loc[comp.syndrome_id==s,'formula_norm'])
  if got!=efs: raise ValueError(f'{s}: verified formula set mismatch. expected={sorted(efs)}, got={sorted(got)}')
 comp=comp.drop_duplicates(['syndrome_id','formula_norm','herb_cn'])
 pairs=pd.read_csv(a.herb_target_pairs,dtype=str)
 pairs=pairs[['canonical_name','target_id']].drop_duplicates()
 idf=pd.read_csv(a.target_idf)
 idf=idf[['target_id','target_herb_degree','N_database_herbs','target_idf']].copy()
 if len(idf)!=374 or idf.target_id.nunique()!=374: raise ValueError('Target IDF file is not frozen 374-target universe')
 dbherbs=set(pairs.canonical_name)
 qa=[]; pri=[]; herbrows=[]
 for s,fs in FORMULAS.items():
  M=len(fs); x=comp[comp.syndrome_id==s].copy()
  m=x.groupby('herb_cn').formula_norm.nunique()
  E=((m+0.5)/(M+1)).rename('E_hs').reset_index()
  E['mapped_exact']=E.herb_cn.isin(dbherbs)
  retained=E.loc[E.mapped_exact,'E_hs'].sum()/E.E_hs.sum()
  mapped=E[E.mapped_exact].copy()
  herbrows.append(E.assign(syndrome_id=s,M_s=M))
  z=mapped.merge(pairs,left_on='herb_cn',right_on='canonical_name',how='left')
  deg=z.groupby('herb_cn').target_id.nunique().rename('herb_target_degree')
  z=z.merge(deg,on='herb_cn').merge(idf,on='target_id',how='left')
  z['contribution']=z.E_hs/z.herb_target_degree*z.target_idf
  tw=z.groupby('target_id').contribution.sum().reindex(idf.target_id,fill_value=0.0)
  tw=tw/tw.sum()
  positive=int((tw>0).sum())
  eligible=(positive>=5 and retained>=0.70)
  qa.append([s,M,len(E),int(E.mapped_exact.sum()),float(E.E_hs.sum()),float(E.loc[E.mapped_exact,'E_hs'].sum()),float(retained),positive,eligible])
  pri.append(pd.DataFrame({'syndrome_id':s,'target_id':tw.index,'target_weight':tw.values}))
 qa=pd.DataFrame(qa,columns=['syndrome_id','M_s','unique_herbs','mapped_herbs','preprojection_weight_mass','mapped_weight_mass','retained_weight_fraction','positive_targets','eligible_primary'])
 if not qa.eligible_primary.all():
  qa.to_csv(out/'psoriasis_prior_QA_FAILED.csv',index=False); raise SystemExit('STOP: >=1 syndrome prior failed frozen eligibility QA')
 pd.concat(herbrows).to_csv(out/'psoriasis_formula_evidence_weights_FROZEN.csv',index=False)
 pd.concat(pri).to_csv(out/'psoriasis_target_prior_374_FROZEN.csv',index=False)
 qa.to_csv(out/'psoriasis_prior_QA_FROZEN.csv',index=False)
 manifest={'rule':'E_hs=(m_hs+0.5)/(M_s+1); contribution=E_hs/|T_h|*IDF_g; sum1 target vector','no_fuzzy_mapping':True,'specificity_weighting':False,'patient_omics_used':False,'inputs':{p:sha(p) for p in [a.composition_csv,a.herb_target_pairs,a.target_idf]}}
 (out/'psoriasis_prior_FREEZE_MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
 print(qa.to_string(index=False))
if __name__=='__main__': main()
