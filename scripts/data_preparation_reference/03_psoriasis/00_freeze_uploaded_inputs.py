#!/usr/bin/env python3
import argparse, gzip, hashlib, json, csv, re
from pathlib import Path
from datetime import datetime, timezone


def sha256(p, chunk=1024*1024):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(chunk),b''): h.update(b)
    return h.hexdigest()

def parse_gmt(p):
    sets=[]
    with open(p,encoding='utf-8') as f:
        for line in f:
            x=line.rstrip('\n').split('\t')
            if len(x)<3: raise ValueError('Malformed GMT line')
            sets.append((x[0],x[1],x[2:]))
    if len(sets)!=50: raise ValueError(f'Expected 50 Hallmarks, found {len(sets)}')
    if len({x[0] for x in sets})!=50: raise ValueError('Duplicate Hallmark names')
    return sets

def parse_series(p):
    meta={}; header=None; nprobe=0
    opener=gzip.open if str(p).endswith('.gz') else open
    with opener(p,'rt',encoding='utf-8',errors='replace') as f:
        in_table=False
        for line in f:
            line=line.rstrip('\n')
            if line.startswith('!Sample_title\t'):
                meta['title']=[x.strip('"') for x in line.split('\t')[1:]]
            elif line.startswith('!Sample_geo_accession\t'):
                meta['gsm']=[x.strip('"') for x in line.split('\t')[1:]]
            elif line=='!series_matrix_table_begin': in_table=True
            elif in_table and header is None:
                header=[x.strip('"') for x in line.split('\t')]
            elif in_table and line=='!series_matrix_table_end': break
            elif in_table and header is not None: nprobe+=1
    if not header or header[0] != 'ID_REF': raise ValueError('Series matrix table not found')
    if len(header)-1 != 72: raise ValueError(f'Expected 72 samples, found {len(header)-1}')
    if nprobe != 135750: raise ValueError(f'Expected 135750 probes, found {nprobe}')
    if meta.get('gsm') != header[1:]: raise ValueError('Metadata GSM order differs from matrix header')
    rows=[]
    for gsm,title in zip(meta['gsm'],meta['title']):
        tl=title.lower()
        if 'blood dryness' in tl: group='Blood_Dryness'
        elif 'blood heat' in tl: group='Blood_Heat'
        elif 'blood stasis' in tl: group='Blood_Stasis'
        elif 'healthy' in tl or 'normal' in tl: group='Healthy'
        else: group='UNRESOLVED'
        rows.append((gsm,title,group))
    counts={g:sum(r[2]==g for r in rows) for g in sorted({r[2] for r in rows})}
    expected={'Blood_Dryness':23,'Blood_Heat':16,'Blood_Stasis':23,'Healthy':10}
    if counts != expected: raise ValueError(f'Unexpected sample groups: {counts}')
    return rows,nprobe

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--series',required=True); ap.add_argument('--gmt',required=True); ap.add_argument('--outdir',required=True)
    a=ap.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    series=Path(a.series); gmt=Path(a.gmt)
    sets=parse_gmt(gmt); samples,nprobe=parse_series(series)
    with open(out/'GSE192867_sample_manifest_FROZEN.csv','w',newline='',encoding='utf-8') as f:
        w=csv.writer(f); w.writerow(['sample_id','title','group','use_primary_within_disease'])
        for gsm,title,g in samples: w.writerow([gsm,title,g,g!='Healthy'])
    with open(out/'Hallmark_2026_1_Hs_FROZEN.tsv','w',newline='',encoding='utf-8') as f:
        w=csv.writer(f,delimiter='\t'); w.writerow(['hallmark','source','n_genes','genes'])
        for name,src,genes in sets: w.writerow([name,src,len(genes),';'.join(genes)])
    manifest={
      'freeze_version':'DCSMI_validation_input_freeze_v1_0',
      'timestamp_utc':datetime.now(timezone.utc).isoformat(),
      'series_matrix':{'path':str(series),'sha256':sha256(series),'n_samples':72,'n_probes':nprobe},
      'hallmark_gmt':{'path':str(gmt),'sha256':sha256(gmt),'n_gene_sets':50},
      'sample_counts':{g:sum(r[2]==g for r in samples) for g in ['Blood_Heat','Blood_Stasis','Blood_Dryness','Healthy']},
      'guardrails':['No syndrome omics used to construct T_s or D_d','Primary O uses syndrome-vs-all-other-syndromes within psoriasis','Healthy controls excluded from primary O','No Hallmark removed based on results']
    }
    (out/'INPUT_FREEZE_MANIFEST_v1_0.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps(manifest,indent=2,ensure_ascii=False))
if __name__=='__main__': main()
