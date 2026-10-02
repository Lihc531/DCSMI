#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DCSMI Master Asset Audit and Final Freeze v1.0
Audits extracted project assets, computes SHA256, assigns provenance/freeze roles,
and exports JSON/CSV audit tables. Workbook formatting is produced by the companion
artifact build in the audited release. No biological scores are recomputed.
"""
from pathlib import Path
import hashlib, csv, json, re, argparse

def sha256_file(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1024*1024), b''): h.update(b)
    return h.hexdigest()

def classify(group, rel):
    low=rel.lower(); role='shared_support'; status='RETAIN_SUPPORT'; reason=''
    if group=='main_results':
        if re.match(r'^(0[2-9]|1[0-2][ab]?)_', low):
            role='WH_WC_formal_chain'; status='FREEZE_CANDIDATE'; reason='formal 02–12B chain; consult final freeze/audit statements'
        elif 'psoriasis' in low or 'gse192867' in low:
            role='Psoriasis_chain'; status='RETAIN_CURRENT_OR_PROVENANCE'
        elif 'dcsmi_' in low: role='WH_WC_legacy_or_support'
    elif group=='new_scripts': role='reproducibility_script_or_bundle'
    elif group.startswith('legacy_'): role='legacy_code_or_raw_asset'
    elif group=='manuscript_old': role='old_manuscript_asset'; status='HISTORICAL_MANUSCRIPT_NOT_CURRENT_FINAL'
    if any(x in low for x in ['psoriasis_target_pathway_mapping_v1_0','psoriasis_target_pathway_mapping_v1_1','psoriasis_target_pathway_mapping_v1_3_offline','pathway_annotation_database']):
        status='SUPERSEDED_DO_NOT_USE_FOR_FINAL_PATHWAY_CLAIMS'; reason='historical target→pathway scaffold route superseded by GO-based v2 chain'
    if 'external_prior_pathway_bridge_v1_' in low:
        status='SUPERSEDED_EXPLORATORY'; reason='superseded by syndrome-specific GO bridge v2.x'
    if '07_cap_transcriptomic_metaz_v1_0/' in low and 'v1_0_1' not in low: status='SUPERSEDED_BY_v1_0_1'
    if '08_syndrome_disease_alignment_matched_null_v1_0/' in low and 'v1_0_1' not in low: status='SUPERSEDED_BY_v1_0_1'
    return role,status,reason

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--root', required=True, help='directory containing main,new,r1,py,cbc extracted folders')
    ap.add_argument('--outdir', default='DCSMI_Master_Asset_Audit_v1_0')
    args=ap.parse_args(); base=Path(args.root); out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)
    groups={'main_results':'main','new_scripts':'new','legacy_R':'r1','legacy_python':'py','manuscript_old':'cbc'}
    rows=[]
    for group,sub in groups.items():
        root=base/sub
        if not root.exists(): continue
        for p in root.rglob('*'):
            if not p.is_file(): continue
            rel=str(p.relative_to(root)); role,status,reason=classify(group,rel)
            rows.append({'group':group,'relative_path':rel,'filename':p.name,'extension':p.suffix.lower(),'size_bytes':p.stat().st_size,'sha256':sha256_file(p),'role':role,'freeze_status':status,'status_reason':reason})
    fields=list(rows[0]) if rows else []
    with open(out/'asset_inventory_v1_0.csv','w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)
    summary={'version':'DCSMI_Master_Asset_Audit_and_Final_Freeze_v1.0','asset_rows':len(rows),'groups':{g:sum(r['group']==g for r in rows) for g in groups},'note':'Classification is provenance/freeze policy; it does not recompute biological results.'}
    (out/'master_audit_summary_v1_0.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps(summary,indent=2,ensure_ascii=False))
if __name__=='__main__': main()
