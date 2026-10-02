#!/usr/bin/env python3
from pathlib import Path
import csv, hashlib, sys
ROOT=Path(__file__).resolve().parents[2]
def sha(p):
 h=hashlib.sha256();
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()
manifest=ROOT/'provenance'/'SHA256SUMS.csv'
if not manifest.exists(): raise SystemExit('Missing provenance/SHA256SUMS.csv')
rows=list(csv.DictReader(open(manifest,encoding='utf-8-sig')))
bad=[]
for r in rows:
 p=ROOT/r['path']
 if not p.exists(): bad.append((r['path'],'MISSING')); continue
 if sha(p)!=r['sha256']: bad.append((r['path'],'HASH_MISMATCH'))
print(f'Checked {len(rows)} repository files.')
if bad:
 for x in bad: print(*x)
 raise SystemExit(1)
print('PASS: all packaged-file SHA256 checks match.')
print('Note: external third-party inputs and the lost historical URTI endpoint generator/raw inputs are outside this packaged-file check; reconstructed code is included and hashed.')
