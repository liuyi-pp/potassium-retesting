"""Five v1.0 sensitivities. Patient/admission rows remain in RAM only."""
from pathlib import Path
import json,csv,time,sys
import numpy as np
import pandas as pd
import hyperkalemia_formal_v1_core as m

HERE=Path(__file__).resolve().parent
OUT=HERE.parent/'outputs'/'hyperkalemia_v1_sensitivities'
MAIN=HERE.parent/'outputs'/'hyperkalemia_v1_main_analysis'/'formal_v1_main_aggregate.json'
EXECUTION=HERE/'hyperkalemia_sensitivity_v1_execution.json'
STRATA=['5.5_to_lt6','6_to_lt6.5','ge6.5']
NAMES={'MAIN':'主分析','S1':'起点前有较低血钾记录','S2':'允许仅溶血警示的起点结果','S3':'仅计生化复查，固定主队列','S4':'排除预采集后报告标本','S5':'全部符合条件住院，患者整群bootstrap'}
GROUPS={'overall':'总体','5.5_to_lt6':'5.5–<6.0','6_to_lt6.5':'6.0–<6.5','ge6.5':'≥6.5'}

def earliest_patient(x):
 return x.sort_values(['subject_id','storetime','hadm_id','labevent_id']).drop_duplicates('subject_id').copy()

def select_admissions(ctx,ep,bounds,t,policy):
 k=ctx['k'];flow=[]
 x=k[k.source.eq('chemistry')&k.numeric_ok&k.v.ge(5.5)&k.original_window].copy()
 prev=m.add_flow(flow,'chemistry_high_in_first_record',x)
 x=x[x.specimen_id.notna()&~x.spec_conflict&x.valid_available&x.storetime.lt(x.stop)].copy()
 prev=m.add_flow(flow,'identifiable_index_and_valid_availability',x,prev)
 x=x[~x[policy]].copy();prev=m.add_flow(flow,'index_comment_policy',x,prev)
 x=m.attach_location(x,ep,bounds,t);prev=m.add_flow(flow,'verified_continuous_ICU',x,prev)
 x=x.sort_values(['hadm_id','storetime','charttime','labevent_id']).drop_duplicates('hadm_id').copy()
 x=m.earlier_flags(x,k,ctx['ce'])
 x=x[~x.earlier_lab_high].copy();prev=m.add_flow(flow,'no_earlier_numeric_lab_high',x,prev)
 x=x[~x.earlier_chart_high].copy();prev=m.add_flow(flow,'no_earlier_chart_high',x,prev)
 x=x[~(x.uncertain_earlier_lab_high|x.uncertain_earlier_chart_high)].copy()
 m.add_flow(flow,'no_uncertain_prior_high',x,prev)
 x=x.merge(ctx['c'].reset_index()[['hadm_id','age','gender','dischtime','deathtime']],on='hadm_id',validate='one_to_one')
 assert x.hadm_id.is_unique
 return x,flow

def prior_lower_mask(index,k):
 g=k[k.numeric_ok&k.v.lt(5.5)&k.valid_available].groupby('hadm_id',sort=False).storetime.min()
 return index.hadm_id.map(g).lt(index.storetime)if len(g)else pd.Series(False,index=index.index)

def pending_mask(r,g):
 return (g.numeric_ok&g.specimen_id.notna()&g.specimen_id.ne(r.specimen_id)&g.charttime.ge(r.charttime)&g.charttime.le(r.storetime)&g.valid_available&g.storetime.gt(r.storetime))

def cluster_bootstrap(d,reps=2000,seed=20261107):
 """Percentile CIs for admission-weighted ratios with full-patient resampling."""
 ids=np.sort(d.subject_id.unique());lookup=pd.Index(ids);codes=lookup.get_indexer(d.subject_id)
 groups=['overall']+STRATA;mat=np.zeros((len(ids),4*len(groups)),dtype=np.int64)
 for j,label in enumerate(groups):
  a=np.ones(len(d),dtype=bool)if label=='overall'else d.stratum.eq(label).to_numpy()
  columns=[a,a&d.state4.eq('repeat_recorded').to_numpy(),a&d.first_value_status.isin(['high','below_high_threshold']).to_numpy(),a&d.first_value_status.eq('high').to_numpy()]
  for q,v in enumerate(columns):np.add.at(mat[:,4*j+q],codes,v.astype(int))
 rng=np.random.default_rng(seed);values=np.full((reps,len(groups),2),np.nan)
 for b in range(reps):
  totals=mat[rng.integers(0,len(ids),size=len(ids))].sum(axis=0).reshape(-1,4)
  for j in range(len(groups)):
   if totals[j,0]:values[b,j,0]=totals[j,1]/totals[j,0]
   if totals[j,2]:values[b,j,1]=totals[j,3]/totals[j,2]
 out={}
 for j,label in enumerate(groups):
  out[label]={}
  for q,name in enumerate(['retest4','first_repeat_still_high']):
   a=values[:,j,q];a=a[np.isfinite(a)]
   out[label][name]={'ci95_lower':float(np.quantile(a,.025,method='linear'))if len(a)else None,'ci95_upper':float(np.quantile(a,.975,method='linear'))if len(a)else None,'valid_replicates':len(a),'undefined_replicates':int(reps-len(a))}
 return out

def summarize_variant(d,analysis,boot=None):
 out=[]
 for label in ['overall']+STRATA:
  part=d if label=='overall'else d[d.stratum.eq(label)]
  z=m.summarize(part,label);z['n_events']=z.pop('n_patients');z['n_unique_patients']=int(part.subject_id.nunique())
  z['analysis']=analysis;z['analysis_unit']='admission'if analysis=='S5'else 'patient'
  for name in ['retest4','first_repeat_still_high']:
   z[name]['ci_method']='patient-cluster percentile bootstrap (2000 replicates)'if boot is not None else 'Wilson'
   if boot is not None:z[name].update(boot[label][name])
  out.append(z)
 for name in ['retest4','first_repeat_still_high']:
  for key in ['numerator','denominator']:assert sum(a[name][key]for a in out[1:])==out[0][name][key]
 assert sum(a['n_events']for a in out[1:])==out[0]['n_events']
 for state in out[0]['states4']:assert sum(a['states4'][state]for a in out[1:])==out[0]['states4'][state]
 return out

def export_tables(result):
 OUT.mkdir(parents=True,exist_ok=True)
 rows=[];readable=[]
 for variant,groups in result['results'].items():
  for z in groups:
   a=z['retest4'];b=z['first_repeat_still_high'];base=next(x for x in result['results']['MAIN']if x['group']==z['group'])
   row={'analysis':variant,'group':z['group'],'analysis_unit':z['analysis_unit'],'N_events':z['n_events'],'N_unique_patients':z['n_unique_patients'],'repeat4h_n':a['numerator'],'repeat4h_N':a['denominator'],'repeat4h_pct':None if a['proportion']is None else 100*a['proportion'],'repeat4h_CI95_lower_pct':None if a['ci95_lower']is None else 100*a['ci95_lower'],'repeat4h_CI95_upper_pct':None if a['ci95_upper']is None else 100*a['ci95_upper'],'first_high_n':b['numerator'],'first_interpretable_N':b['denominator'],'first_high_pct':None if b['proportion']is None else 100*b['proportion'],'first_high_CI95_lower_pct':None if b['ci95_lower']is None else 100*b['ci95_lower'],'first_high_CI95_upper_pct':None if b['ci95_upper']is None else 100*b['ci95_upper'],'first_value_indeterminate':z['first_repeat_value_status']['indeterminate'],'no_numeric_repeat':z['first_repeat_value_status']['no_numeric_repeat'],'repeat4h_indeterminate':z['states4']['indeterminate'],'repeat4h_delta_pp':None if a['proportion']is None else 100*(a['proportion']-base['retest4']['proportion']),'first_high_delta_pp':None if b['proportion']is None else 100*(b['proportion']-base['first_repeat_still_high']['proportion']),'ci_method':a['ci_method']}
   rows.append(row)
   def fmt(v):
    if v['proportion']is None:return '无可用分母'
    return f"{v['numerator']}/{v['denominator']}；{100*v['proportion']:.2f}%（{100*v['ci95_lower']:.2f}–{100*v['ci95_upper']:.2f}）"
   readable.append({'分析':NAMES[variant],'血钾分层（mmol/L）':GROUPS[z['group']],'分析单位':z['analysis_unit'],'事件数':z['n_events'],'独立患者数':z['n_unique_patients'],'4小时确认复测：n/N；%（95%CI）':fmt(a),'首次复查仍高钾：n/可判读N；%（95%CI）':fmt(b),'4小时复测较主分析变化（百分点）':row['repeat4h_delta_pp'],'首次复查高钾较主分析变化（百分点）':row['first_high_delta_pp'],'4小时状态不确定':z['states4']['indeterminate'],'首次复查结果不确定':z['first_repeat_value_status']['indeterminate'],'观察期内无数值复查':z['first_repeat_value_status']['no_numeric_repeat'],'区间方法':a['ci_method']})
 for file,data in [('sensitivity_results.csv',rows),('五项敏感性分析_总体及分层.csv',readable)]:
  with open(OUT/file,'w',encoding='utf-8-sig',newline='')as f:
   w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
 temp=OUT/'sensitivity_v1_aggregate.tmp';temp.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8');temp.replace(OUT/'sensitivity_v1_aggregate.json')

def run(ctx):
 start=time.monotonic();old=json.loads(MAIN.read_text(encoding='utf-8'))
 assert m.sha(m.PROTOCOL)==ctx['protocol_sha256']==old['protocol_sha256']
 assert m.sha(Path(m.__file__))==old['analysis_code_sha256'],'Main-analysis engine changed'
 m.say('sensitivity_cohort_reconstruction')
 k=ctx['k'];c=ctx['c']
 broad=k[k.source.eq('chemistry')&k.numeric_ok&k.v.ge(5.5)&k.original_window&k.specimen_id.notna()&~k.spec_conflict&k.valid_available&k.storetime.lt(k.stop)&~k.explicit_invalid]
 ep,bounds,t=m.fast_intervals(c,ctx['tr'],set(broad.hadm_id))
 all_idx,mainflow=select_admissions(ctx,ep,bounds,t,'refined_flag')
 main_idx=earliest_patient(all_idx)
 s2_all,s2flow=select_admissions(ctx,ep,bounds,t,'explicit_invalid');s2_idx=earliest_patient(s2_all)
 assert len(all_idx)==old['eligible_admissions_before_patient_selection']==5577
 assert len(main_idx)==old['final_patients']==5230
 assert set(all_idx.hadm_id)<=set(s2_all.hadm_id)
 assert set(main_idx.subject_id)<=set(s2_idx.subject_id)
 needed=set(all_idx.hadm_id)|set(s2_idx.hadm_id)
 groups={h:g for h,g in k[k.hadm_id.isin(needed)].sort_values(['charttime','source_rank','labevent_id']).groupby('hadm_id',sort=False)}
 cache={}
 def outcomes(index,chemistry_only=False):
  rows=[]
  for r in index.itertuples():
   key=(int(r.labevent_id),chemistry_only)
   if key not in cache:
    g=groups[r.hadm_id]
    if chemistry_only:g=g[g.source.eq('chemistry')]
    cache[key]=m.analyze_event(r,g)
   rows.append(cache[key])
  return pd.DataFrame(rows)
 all_d=outcomes(all_idx);main_d=all_d[all_d.index_id.isin(set(main_idx.labevent_id))].copy()
 # Reproduce all previously published main summaries before running sensitivities.
 check=[m.summarize(main_d,'overall')]+[m.summarize(main_d[main_d.stratum.eq(g)],g)for g in STRATA]
 assert check==old['summaries'],'Main-analysis results not exactly reproduced'
 m.say('main_reproduced',patients=len(main_d),admissions=len(all_d))
 keep1=prior_lower_mask(main_idx,k);s1_ids=set(main_idx.loc[keep1,'labevent_id'])
 d1=main_d[main_d.index_id.isin(s1_ids)].copy()
 d2=outcomes(s2_idx);d3=outcomes(main_idx,True)
 pending_ids=set();pending_at_left=set();low_same_only=0;low_flag_only=0
 for r in main_idx.itertuples():
  g=groups[r.hadm_id];pending=pending_mask(r,g)
  if pending.any():pending_ids.add(int(r.labevent_id))
  if (pending&g.charttime.eq(r.charttime)).any():pending_at_left.add(int(r.labevent_id))
  if int(r.labevent_id)in s1_ids:
   lower=g[g.numeric_ok&g.v.lt(5.5)&g.valid_available&g.storetime.lt(r.storetime)]
   low_same_only+=int(not lower.specimen_id.ne(r.specimen_id).any())
   low_flag_only+=int(lower.refined_flag.all())
 d4=main_d[~main_d.index_id.isin(pending_ids)].copy()
 assert len(d1)==int(keep1.sum())and len(d4)+len(pending_ids)==len(main_d)
 assert set(d3.index_id)==set(main_d.index_id)
 pair=main_d[['index_id','state4','event_h','first_value_status']].merge(d3[['index_id','state4','event_h','first_value_status']],on='index_id',suffixes=('_main','_chem'),validate='one_to_one')
 assert not((pair.state4_chem=='repeat_recorded')&(pair.state4_main!='repeat_recorded')).any()
 observed=pair.event_h_main.notna()&pair.event_h_chem.notna()
 assert (pair.loc[observed,'event_h_chem']+1e-12>=pair.loc[observed,'event_h_main']).all()
 joined=main_idx[['subject_id','hadm_id','labevent_id','v']].merge(s2_idx[['subject_id','hadm_id','labevent_id','v']],on='subject_id',suffixes=('_main','_S2'),validate='one_to_one')
 paired_qc={}
 for group in ['overall']+STRATA:
  ids=set(main_d.index_id if group=='overall'else main_d.loc[main_d.stratum.eq(group),'index_id'])
  p=pair[pair.index_id.isin(ids)]
  both=p.first_value_status_main.isin(['high','below_high_threshold'])&p.first_value_status_chem.isin(['high','below_high_threshold'])
  paired_qc[group]={'fixed_patients':len(p),'main_repeat4_lost_with_chemistry_only':int((p.state4_main.eq('repeat_recorded')&~p.state4_chem.eq('repeat_recorded')).sum()),'repeat4_gained_with_chemistry_only':int((~p.state4_main.eq('repeat_recorded')&p.state4_chem.eq('repeat_recorded')).sum()),'both_first_values_interpretable':int(both.sum()),'first_high_classification_disagrees_among_both_interpretable':int((both&p.first_value_status_main.ne(p.first_value_status_chem)).sum())}
 m.say('sensitivity_outcomes_complete',S1=len(d1),S2=len(d2),S3=len(d3),S4=len(d4),S5_admissions=len(all_d))
 boot=cluster_bootstrap(all_d,2000,20261107)
 results={label:summarize_variant(data,label,boot if label=='S5'else None)for label,data in [('MAIN',main_d),('S1',d1),('S2',d2),('S3',d3),('S4',d4),('S5',all_d)]}
 integrity=[];prior={x['file']:x for x in old['source_integrity']}
 for name in m.FILES:
  path=m.ROOT/name;stat=path.stat();digest=m.sha(path)
  integrity.append({'file':name,'sha256':digest,'matches_previous_main':digest==prior[name]['sha256'],'size_mtime_unchanged':ctx['stamps'][name]==(stat.st_size,stat.st_mtime_ns)})
 assert all(z['matches_previous_main']and z['size_mtime_unchanged']for z in integrity)
 qc={'S1_removed':len(main_d)-len(d1),'S1_only_prior_lower_on_index_specimen':low_same_only,'S1_only_flagged_prior_lower':low_flag_only,'S2_primary_patients_preserved':len(joined),'S2_new_patients':len(d2)-len(joined),'S2_existing_patients_with_reselected_event':int(joined.labevent_id_main.ne(joined.labevent_id_S2).sum()),'S2_index_with_hemolysis_warning':int(s2_idx.refined_flag.sum()),'S2_explicitly_invalid_indexes':int(s2_idx.explicit_invalid.sum()),'S2_eligible_admissions':len(s2_all),'S3_paired_checks':paired_qc,'S4_removed_precollected_later_reported':len(pending_ids),'S4_patients_with_lower_bound_timestamp_tie':len(pending_at_left),'S5_unique_patients':int(all_d.subject_id.nunique()),'S5_repeated_admissions':len(all_d)-int(all_d.subject_id.nunique()),'S5_bootstrap_replicates':2000,'S5_bootstrap_seed':20261107}
 result={'protocol_version':'1.0','protocol_sha256':m.sha(m.PROTOCOL),'execution_specification':json.loads(EXECUTION.read_text(encoding='utf-8')),'execution_specification_sha256':m.sha(EXECUTION),'main_results_exactly_reproduced':True,'main_engine_sha256':m.sha(Path(m.__file__)),'sensitivity_engine_sha256':m.sha(Path(__file__)),'comment_rules_sha256':m.sha(HERE/'potassium_validation_core.py'),'results':results,'qc':qc,'cohort_flows':{'primary_admissions':mainflow,'S2_comment_tolerance_admissions':s2flow},'source_integrity':integrity,'rows_scanned':{'labevents':ctx['nlab'],'chartevents':ctx['nchart']},'software':{'python':sys.version,'pandas':pd.__version__,'numpy':np.__version__},'elapsed_analysis_seconds':round(time.monotonic()-start,2),'provenance':'v1.0 five planned sensitivities following feasibility-informed protocol and formal main analysis; no new clinical variables or effect models; no claim of prospective preregistration'}
 result['severity_pattern_checks']={label:{'repeat4h_proportions':[x['retest4']['proportion']for x in groups[1:]],'repeat4h_increases_across_three_strata':bool(all(a['retest4']['proportion']<b['retest4']['proportion']for a,b in zip(groups[1:],groups[2:]))),'first_high_proportions':[x['first_repeat_still_high']['proportion']for x in groups[1:]],'minimum_stratum_events':min(x['n_events']for x in groups[1:])}for label,groups in results.items()}
 export_tables(result)
 m.say('sensitivity_complete',output=str(OUT),elapsed_s=round(time.monotonic()-start,1))
