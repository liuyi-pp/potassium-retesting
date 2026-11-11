"""Actual observation duration only; immutable v1.0 selector and aggregate checks.
Usage: bundled-python this_script.py [--data-root AUTHORIZED_LOCAL_ROOT]
Patient-level records stay in RAM. Original result files are never overwritten.
"""
from pathlib import Path
import argparse,csv,hashlib,json,sys,time
import numpy as np
import pandas as pd
import hyperkalemia_formal_v1_core as m
import hyperkalemia_tables_v1_core as a

parser=argparse.ArgumentParser();parser.add_argument('--data-root',type=Path)
args=parser.parse_args()
if args.data_root is not None:m.ROOT=args.data_root
ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/'outputs/hyperkalemia_v1_observation_duration'
OUT.mkdir(parents=True,exist_ok=True)
contract=ROOT/'work/hyperkalemia_observation_duration_contract_v1.json'
main_path=ROOT/'outputs/hyperkalemia_v1_main_analysis/formal_v1_main_aggregate.json'
tables_path=ROOT/'outputs/hyperkalemia_v1_manuscript_tables/tables_v1_aggregate.json'
prior=json.loads(main_path.read_text(encoding='utf-8'))
oldtabs=json.loads(tables_path.read_text(encoding='utf-8'))
old_hashes={str(p.relative_to(ROOT)):m.sha(p)for p in [main_path,tables_path,m.PROTOCOL,Path(m.__file__),ROOT/'work/potassium_validation_core.py',Path(a.__file__)]}
assert old_hashes[str(Path(m.__file__).relative_to(ROOT))]==prior['analysis_code_sha256']
assert m.sha(m.PROTOCOL)==prior['protocol_sha256']
assert m.sha(ROOT/'work/potassium_validation_core.py')==prior['comment_rules_code_sha256']
start=time.monotonic();m.say('observation_duration_addendum_started')
ctx=m.load_data()
k=ctx['k'];c=ctx['c'];ce=ctx['ce'];flow=[dict(x)for x in ctx['flow']];prev=flow[-1]['remaining_admissions']
x=k[k.source.eq('chemistry')&k.numeric_ok&k.v.ge(5.5)&k.original_window].copy();prev=m.add_flow(flow,'05_numeric_chemistry_high_in_first_ICU_record_window',x,prev)
x=x[x.specimen_id.notna()&~x.spec_conflict].copy();prev=m.add_flow(flow,'06_identifiable_nonconflicting_index_specimen',x,prev)
x=x[x.valid_available&x.storetime.lt(x.stop)].copy();prev=m.add_flow(flow,'07_valid_index_availability_within_record_window',x,prev)
x=x[~x.refined_flag].copy();prev=m.add_flow(flow,'08_negation_aware_index_comment_rule',x,prev)
m.say('observation_reconstructing_ICU_intervals',candidate_admissions=prev)
ep,bounds,t=m.fast_intervals(c,ctx['tr'],set(x.hadm_id))
x=m.attach_location(x,ep,bounds,t);prev=m.add_flow(flow,'09_sample_and_result_in_same_verified_continuous_ICU',x,prev)
idx=x.sort_values(['hadm_id','storetime','charttime','labevent_id']).drop_duplicates('hadm_id').copy()
idx=m.earlier_flags(idx,k,ce)
idx=idx[~idx.earlier_lab_high].copy();prev=m.add_flow(flow,'10_no_earlier_numeric_lab_high_regardless_comment',idx,prev)
idx=idx[~idx.earlier_chart_high].copy();prev=m.add_flow(flow,'11_no_earlier_charted_high',idx,prev)
idx=idx[~(idx.uncertain_earlier_lab_high|idx.uncertain_earlier_chart_high)].copy();prev=m.add_flow(flow,'12_no_unresolvable_prior_high_timing',idx,prev)
assert len(idx)==5577
idx=idx.sort_values(['subject_id','storetime','hadm_id','labevent_id']).drop_duplicates('subject_id').copy();prev=m.add_flow(flow,'13_first_qualifying_event_per_patient',idx,prev)
assert len(idx)==5230 and idx.subject_id.is_unique and idx.hadm_id.is_unique and flow==prior['flow']
idx=idx.merge(c.reset_index()[['hadm_id','gender','age','dischtime','deathtime']],on='hadm_id',validate='one_to_one')
m.say('observation_cohort_reproduced',patients=len(idx),flow_exact=True)
groups={h:g for h,g in k[k.hadm_id.isin(set(idx.hadm_id))].sort_values(['charttime','source_rank','labevent_id']).groupby('hadm_id',sort=False)}
rows=[];direct=[]
for i,r in enumerate(idx.itertuples(),1):
    g=groups[r.hadm_id]
    event=m.analyze_event(r,g)
    event.update(a.states_at_horizons(r,g,event))
    clinical_end=min(r.episode_stop,r.dischtime,r.deathtime if pd.notna(r.deathtime)else r.episode_stop)
    observation_end=min(clinical_end,r.storetime+pd.Timedelta(hours=24))
    duration=(observation_end-r.storetime).total_seconds()/3600
    assert bool(np.isclose(duration,min(event['follow_h'],24),rtol=0,atol=1e-10))
    direct.append(duration);event['observation_duration_h']=duration;rows.append(event)
    if i%1500==0:m.say('observation_verification_progress',patients=i)
d=pd.DataFrame(rows)
summaries=[m.summarize(d,'overall')]+[m.summarize(d[d.stratum.eq(g)],g)for g in a.LABELS[1:]]
assert summaries==prior['summaries'],'Unchanged v1.0 main results failed reproduction'
baseline=[a.baseline_summary(idx,g)for g in a.LABELS];assert baseline==oldtabs['baseline']
windows=[a.window_summary(d,g,h)for g in a.LABELS for h in a.HORIZONS];assert windows==oldtabs['time_windows']
values=d.observation_duration_h.to_numpy();assert np.isfinite(values).all()and(values>0).all()and(values<=24).all()
result=[]
for group in a.LABELS:
    part=d if group=='overall'else d[d.stratum.eq(group)]
    v=part.observation_duration_h.to_numpy(dtype=float);N=len(v)
    q=np.quantile(v,[.25,.5,.75],method='linear');full=int((v==24).sum());early=int((v<24).sum())
    assert full+early==N
    old=next(x for x in prior['summaries']if x['group']==group)
    beyond4=int(part.observation_duration_h.gt(4).sum())
    assert beyond4==old['supplementary_no_repeat_with_full4h']['denominator']
    result.append({'group':group,'N':N,'observation_duration_h':{'n':N,'median':float(q[1]),'q25':float(q[0]),'q75':float(q[2]),'total_patient_hours':float(v.sum())},'duration_lt24_n':early,'duration_lt24_pct':100*early/N,'duration_24h_n':full,'duration_24h_pct':100*full/N,'duration_gt4h_n_for_existing_QC':beyond4,'missing_duration_n':0})
assert [x['N']for x in result]==[5230,2895,1079,1256]
for key in ['N','duration_lt24_n','duration_24h_n','duration_gt4h_n_for_existing_QC']:assert sum(r[key]for r in result[1:])==result[0][key]
assert np.isclose(sum(r['observation_duration_h']['total_patient_hours']for r in result[1:]),result[0]['observation_duration_h']['total_patient_hours'],atol=1e-8)
m.say('observation_all_frozen_results_reproduced',main_groups=4,baseline_groups=4,time_window_rows=16)
cohort_fingerprint=hashlib.sha256(idx.sort_values('subject_id')[['subject_id','hadm_id','labevent_id','storetime']].to_csv(index=False).encode('utf-8')).hexdigest()
integrity=[];before={x['file']:x for x in prior['source_integrity']}
for name in m.FILES:
    f=m.ROOT/name;st=f.stat();digest=m.sha(f)
    integrity.append({'file':name,'sha256':digest,'matches_formal_main':digest==before[name]['sha256'],'size_mtime_unchanged':ctx['stamps'][name]==(st.st_size,st.st_mtime_ns)})
assert all(x['matches_formal_main']and x['size_mtime_unchanged']for x in integrity)
assert all(m.sha(ROOT/p)==h for p,h in old_hashes.items())
payload={'analysis':'Actual observation-duration reporting addendum only','protocol_version':'1.0','date':'2026-11-10','contract_sha256':m.sha(contract),'supplement_script_sha256':m.sha(Path(__file__)),'provenance':'Additional descriptive reporting requested after manuscript review; no claim of prospective prespecification','definition':'min(verified ICU interval end, hospital discharge, death if available, t0+24h) minus index storetime; not stopped at first repeat','duration_cap_h':24,'groups':result,'verification':{'frozen_flow_exact':True,'main_summaries_exact':True,'baseline_exact':True,'all_16_windows_exact':True,'direct_timestamp_vs_existing_follow_h_exact':True,'existing_gt4h_denominator_exact':True,'no_additional_exclusions':True,'no_imputation':True,'no_new_clinical_variables_or_models':True,'cohort_reconstruction_fingerprint_sha256':cohort_fingerprint},'original_artifact_sha256':old_hashes,'source_integrity':integrity,'software':{'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__},'elapsed_s':round(time.monotonic()-start,2)}
(OUT/'observation_duration_v1_aggregate.json').write_text(json.dumps(payload,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
flat=[]
for r in result:
    flat.append({'group':r['group'],'N':r['N'],'duration_median_h':r['observation_duration_h']['median'],'duration_q25_h':r['observation_duration_h']['q25'],'duration_q75_h':r['observation_duration_h']['q75'],'total_patient_hours':r['observation_duration_h']['total_patient_hours'],'duration_lt24_n':r['duration_lt24_n'],'duration_lt24_pct':r['duration_lt24_pct'],'duration_24h_n':r['duration_24h_n'],'duration_24h_pct':r['duration_24h_pct'],'missing_duration_n':r['missing_duration_n']})
with open(OUT/'observation_duration_v1.csv','w',encoding='utf-8-sig',newline='')as f:
    w=csv.DictWriter(f,fieldnames=list(flat[0]));w.writeheader();w.writerows(flat)
m.say('observation_duration_complete',patients=5230,groups=4,elapsed_s=round(time.monotonic()-start,1))
