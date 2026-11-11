"""Run the preserved v1.0 implementations against independently authorized data."""
from pathlib import Path
import argparse,gc,os,subprocess,sys

ROOT=Path(__file__).resolve().parents[1]

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--data-root',type=Path,required=True,help='Authorized MIMIC-IV 3.1 root; not downloaded by this script')
 parser.add_argument('--stage',choices=['main','sensitivity','tables','duration','all'],default='all')
 args=parser.parse_args();data=args.data_root.expanduser().resolve()
 if sys.flags.optimize:parser.error('Run without -O: the frozen source uses scientific integrity assertions.')
 os.environ['MIMIC_IV_ROOT']=str(data)
 sys.path.insert(0,str(ROOT/'work'))
 import hyperkalemia_formal_v1_core as m
 paths={
  'main':ROOT/'outputs/hyperkalemia_v1_main_analysis/formal_v1_main_aggregate.json',
  'sensitivity':ROOT/'outputs/hyperkalemia_v1_sensitivities/sensitivity_v1_aggregate.json',
  'tables':ROOT/'outputs/hyperkalemia_v1_manuscript_tables/tables_v1_aggregate.json',
  'duration':ROOT/'outputs/hyperkalemia_v1_observation_duration/observation_duration_v1_aggregate.json',
 }
 stages=['main','sensitivity','tables','duration']if args.stage=='all'else[args.stage]
 for stage in stages:
  if paths[stage].parent.exists()and any(paths[stage].parent.iterdir()):
   parser.error(f'Output directory is not empty: {paths[stage].parent}. Use a separate working copy; existing outputs are preserved.')
 for name in m.FILES+['SHA256SUMS.txt']:
  if not(data/name).is_file():parser.error(f'Missing required source: {name}')
 if args.stage not in ['main','all']and not paths['main'].is_file():parser.error('Run main first using this same code copy.')
 if args.stage=='duration'and not paths['tables'].is_file():parser.error('Run tables before duration.')
 if any(x!='duration'for x in stages):
  ctx=m.load_data()
  for stage in stages:
   if stage=='main':m.run(ctx)
   elif stage=='sensitivity':
    import hyperkalemia_sensitivity_v1_core as s
    s.run(ctx)
   elif stage=='tables':
    import hyperkalemia_tables_v1_core as a
    a.run(ctx)
  del ctx;gc.collect()
 if 'duration'in stages:
  subprocess.run([sys.executable,str(ROOT/'work/describe_hyperkalemia_observation_duration_v1.py'),'--data-root',str(data)],check=True,cwd=ROOT)

if __name__=='__main__':main()
