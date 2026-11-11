"""Check the published file allowlist and hash manifest without opening clinical data."""
from pathlib import Path
import ast,csv,hashlib,json,re,sys

ROOT=Path(__file__).resolve().parents[1]

def main():
 manifest=json.loads((ROOT/'provenance/public_files.json').read_text(encoding='utf-8'))
 for name,digest in manifest['sha256'].items():
  p=ROOT/name
  assert p.is_file(),name
  assert hashlib.sha256(p.read_bytes()).hexdigest()==digest,name
  if p.suffix=='.py':ast.parse(p.read_text(encoding='utf-8'))
 for name in manifest['sha256']:
  text=(ROOT/name).read_text(encoding='utf-8-sig')
  assert not re.search(r'(?i)(gh[pousr]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{20,}|-----BEGIN (?:RSA |OPENSSH )?PRIVATE KEY-----)',text),name
  assert not re.search(r'[A-Za-z]:[\\/](?:Users|restricted_data)[\\/]',text),name
 allowed_csv={'data/figure_source/cohort_flow.csv','data/figure_source/time_windows.csv'}
 assert {n for n in manifest['sha256']if n.endswith('.csv')}==allowed_csv
 for name in allowed_csv:
  with open(ROOT/name,encoding='utf-8-sig',newline='')as f:
   reader=csv.DictReader(f);assert not {'subject_id','hadm_id','stay_id','specimen_id','labevent_id','charttime','storetime','comments'}&set(reader.fieldnames)
   rows=list(reader);assert len(rows)<=16
 print(json.dumps({'status':'passed','hashed_files':len(manifest['sha256']),'public_tables':'two aggregate CSVs only; no patient rows','scope':'file integrity and allowlist checks, not clinical validation'}))

if __name__=='__main__':main()
