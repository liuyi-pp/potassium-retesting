"""Complete v1.0 baseline and fixed time windows. Aggregate-only exports."""
from pathlib import Path
import json,csv,sys,time
import numpy as np
import pandas as pd
import hyperkalemia_formal_v1_core as m
import hyperkalemia_sensitivity_v1_core as s

HERE=Path(__file__).resolve().parent
OUT=HERE.parent/'outputs'/'hyperkalemia_v1_manuscript_tables'
EXECUTION=HERE/'hyperkalemia_tables_v1_execution.json'
HORIZONS=[2,4,6,24]
LABELS=['overall']+s.STRATA
UNIT_ORDER=['Medical Intensive Care Unit (MICU)','Medical/Surgical Intensive Care Unit (MICU/SICU)','Coronary Care Unit (CCU)','Cardiac Vascular Intensive Care Unit (CVICU)','Surgical Intensive Care Unit (SICU)','Trauma SICU (TSICU)','Neuro Surgical Intensive Care Unit (Neuro SICU)','Intensive Care Unit (ICU)']
UNIT_SHORT=['Medical ICU','Medical/surgical ICU','Coronary care unit','Cardiac vascular ICU','Surgical ICU','Trauma surgical ICU','Neurosurgical ICU','General ICU']
STATES=['repeat_recorded','death_before_repeat','interval_end_before_repeat','still_in_interval_without_repeat','indeterminate']

def states_at_horizons(r,g,event):
 """Extend the exact formal 4h rule to prespecified 2/6/24h horizons."""
 t0=r.storetime;end_time=min(r.episode_stop,r.dischtime,r.deathtime if pd.notna(r.deathtime)else r.episode_stop)
 zero=t0.value;end=end_time.value;ct=m.ns(g.charttime);st=m.ns(g.storetime)
 spec=g.specimen_id.to_numpy(dtype=float);has=np.isfinite(spec);distinct=(~has)|(spec!=r.specimen_id)
 numeric=g.numeric_ok.to_numpy(dtype=bool);explicit=g.explicit_invalid.to_numpy(dtype=bool);conflict=g.spec_conflict.to_numpy(dtype=bool)
 time_conf=(st!=m.NAT)&(ct!=m.NAT)&(st<ct);unresolved_time=(ct==m.NAT)|time_conf
 rows={}
 for h in HORIZONS:
  before=distinct&(ct>zero)&(ct<=zero+h*m.HOUR)&(ct<end)
  uncertain=before&(~numeric|~has|explicit|time_conf|conflict)
  evidence=distinct&unresolved_time&(st>zero)&(st<=zero+h*m.HOUR)&(st<end)
  unknown=bool(np.any(uncertain|evidence));hh=event['event_h']
  if hh is not None and hh<=h:state='repeat_recorded'
  elif unknown:state='indeterminate'
  elif event['follow_h']<=h:state='death_before_repeat'if event['end_reason']=='death'else 'interval_end_before_repeat'
  else:state='still_in_interval_without_repeat'
  rows[f'state_{h}h']=state;rows[f'uncertain_evidence_{h}h']=unknown
 assert rows['state_4h']==event['state4']and rows['uncertain_evidence_4h']==event['unknown_before4']
 return rows

def quantiles(v):
 x=np.asarray(v,dtype=float);x=x[np.isfinite(x)]
 return {'n':len(x),'median':float(np.median(x))if len(x)else None,'q25':float(np.quantile(x,.25))if len(x)else None,'q75':float(np.quantile(x,.75))if len(x)else None}

def baseline_summary(idx,label):
 part=idx if label=='overall'else idx[idx.v.ge(5.5)&(np.where(idx.v.lt(6),'5.5_to_lt6',np.where(idx.v.lt(6.5),'6_to_lt6.5','ge6.5'))==label)]
 n=len(part);age=part.age.clip(upper=90)
 sex=part.gender.fillna('Unknown').where(part.gender.isin(['F','M']),'Unknown')
 assert part.actual_unit.isin(UNIT_ORDER).all()
 result={'group':label,'N':n,'age_capped90':quantiles(age),'age_ge90_n':int(part.age.ge(90).sum()),'sex':{code:int(sex.eq(code).sum())for code in ['F','M','Unknown']},'ICU':{name:int(part.actual_unit.eq(name).sum())for name in UNIT_ORDER},'index_potassium':quantiles(part.v),'missing':{'age':int(part.age.isna().sum()),'sex':int(sex.eq('Unknown').sum()),'ICU':int(part.actual_unit.isna().sum()),'index_potassium':int(part.v.isna().sum())}}
 assert sum(result['sex'].values())==sum(result['ICU'].values())==n
 assert result['age_capped90']['n']+result['missing']['age']==n
 return result

def window_summary(d,label,h):
 part=d if label=='overall'else d[d.stratum.eq(label)];n=len(part);states=part[f'state_{h}h'].value_counts()
 counts={state:int(states.get(state,0))for state in STATES};k=counts['repeat_recorded'];u=counts['indeterminate']
 assert sum(counts.values())==n
 return {'group':label,'horizon_h':h,'N':n,'confirmed_repeat':m.wilson(k,n),'states':counts,'repeat_identification_bounds':[k/n,(k+u)/n]if n else [None,None]}

def write_csv(path,rows):
 with open(path,'w',encoding='utf-8-sig',newline='')as f:
  writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

def baseline_rows(b):
 keys=[('patients','Patients',lambda z:str(z['N'])),('age','Age, years',lambda z:fmt_quantile(z['age_capped90'],1,True)),('age90','Age ≥90 years, n (%)',lambda z:fmt_count(z['age_ge90_n'],z['N'])),('sex_header','Recorded sex, n (%)',lambda z:''),('female','  Female',lambda z:fmt_count(z['sex']['F'],z['N'])),('male','  Male',lambda z:fmt_count(z['sex']['M'],z['N'])),('sex_unknown','  Unknown',lambda z:fmt_count(z['sex']['Unknown'],z['N'])),('ICU_header','Actual ICU type at index, n (%)',lambda z:'')]
 for name,short in zip(UNIT_ORDER,UNIT_SHORT):keys.append(('ICU_'+short,'  '+short,lambda z,name=name:fmt_count(z['ICU'][name],z['N'])))
 keys.append(('index_K','Index potassium, mmol/L',lambda z:fmt_quantile(z['index_potassium'],2)))
 return [{'variable_key':key,'Characteristic':text,**{label:func(z)for label,z in zip(LABELS,b)}}for key,text,func in keys]

def fmt_count(n,N):return f'{n:,} ({100*n/N:.1f})'if N else '—'

def fmt_quantile(q,decimals=2,age=False):
 def num(v):return '≥90'if age and v>=90 else f'{v:.{decimals}f}'
 return f"{num(q['median'])} ({num(q['q25'])}–{num(q['q75'])})"if q['n']else '—'

def export(result):
 OUT.mkdir(parents=True,exist_ok=True)
 write_csv(OUT/'Table_1_baseline.csv',baseline_rows(result['baseline']))
 rows=[]
 for a in result['time_windows']:
  q=a['confirmed_repeat'];rows.append({'group':a['group'],'horizon_h':a['horizon_h'],'N':a['N'],'repeat_n':q['numerator'],'repeat_pct':100*q['proportion'],'CI95_lower_pct':100*q['ci95_lower'],'CI95_upper_pct':100*q['ci95_upper'],**a['states'],'identification_lower_pct':100*a['repeat_identification_bounds'][0],'identification_upper_pct':100*a['repeat_identification_bounds'][1]})
 write_csv(OUT/'Table_2_and_S1_time_windows.csv',rows)
 readable=[]
 for a in rows:
  readable.append({'血钾分层（mmol/L）':s.GROUPS[a['group']],'时间窗（小时）':a['horizon_h'],'总人数':a['N'],'确认复测人数':a['repeat_n'],'确认复测比例（%）':a['repeat_pct'],'Wilson95%CI下限（%）':a['CI95_lower_pct'],'Wilson95%CI上限（%）':a['CI95_upper_pct'],'复测前死亡':a['death_before_repeat'],'复测前ICU观察区间结束':a['interval_end_before_repeat'],'仍在区间无复测记录':a['still_in_interval_without_repeat'],'无法判定':a['indeterminate'],'未知状态识别范围下限（%）':a['identification_lower_pct'],'未知状态识别范围上限（%）':a['identification_upper_pct']})
 write_csv(OUT/'复测时间窗_总体及分层.csv',readable)
 rows3=[]
 for a in result['first_repeat']:
  q=a['first_repeat_still_high'];t=a['observed_repeat_interval_h'];v=a['first_repeat_value_status']
  rows3.append({'group':a['group'],'N':a['n_patients'],'repeat_observed_n':a['numeric_repeat_within_followup'],'interval_median_h':t['median'],'interval_q25_h':t['q25'],'interval_q75_h':t['q75'],'first_value_interpretable_n':q['denominator'],'first_high_n':q['numerator'],'first_high_pct':100*q['proportion'],'first_high_CI95_lower_pct':100*q['ci95_lower'],'first_high_CI95_upper_pct':100*q['ci95_upper'],'first_value_indeterminate_n':v['indeterminate'],'no_numeric_repeat_n':v['no_numeric_repeat']})
 write_csv(OUT/'Table_3_first_repeat.csv',rows3)
 (OUT/'tables_v1_aggregate.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')

def run(ctx):
 start=time.monotonic();prior=json.loads(s.MAIN.read_text(encoding='utf-8'))
 assert m.sha(m.PROTOCOL)==ctx['protocol_sha256']==prior['protocol_sha256']
 assert m.sha(Path(m.__file__))==prior['analysis_code_sha256']
 k=ctx['k'];broad=k[k.source.eq('chemistry')&k.numeric_ok&k.v.ge(5.5)&k.original_window&k.specimen_id.notna()&~k.spec_conflict&k.valid_available&k.storetime.lt(k.stop)&~k.refined_flag]
 m.say('table_cohort_reconstruction')
 ep,bounds,t=m.fast_intervals(ctx['c'],ctx['tr'],set(broad.hadm_id))
 admissions,flow=s.select_admissions(ctx,ep,bounds,t,'refined_flag');idx=s.earliest_patient(admissions)
 assert len(admissions)==5577 and len(idx)==5230 and idx.subject_id.is_unique and idx.hadm_id.is_unique
 groups={h:g for h,g in k[k.hadm_id.isin(set(idx.hadm_id))].sort_values(['charttime','source_rank','labevent_id']).groupby('hadm_id',sort=False)}
 rows=[]
 for i,r in enumerate(idx.itertuples(),1):
  event=m.analyze_event(r,groups[r.hadm_id]);event.update(states_at_horizons(r,groups[r.hadm_id],event));rows.append(event)
  if i%1500==0:m.say('table_outcomes',patients=i)
 d=pd.DataFrame(rows)
 first=[m.summarize(d,'overall')]+[m.summarize(d[d.stratum.eq(label)],label)for label in s.STRATA]
 assert first==prior['summaries'],'Previous main summaries did not reproduce exactly'
 baseline=[baseline_summary(idx,label)for label in LABELS]
 windows=[window_summary(d,label,h)for label in LABELS for h in HORIZONS]
 for label in LABELS:
  w=[x for x in windows if x['group']==label]
  assert all(a['confirmed_repeat']['numerator']<=b['confirmed_repeat']['numerator']for a,b in zip(w,w[1:]))
  original=next(x for x in first if x['group']==label)
  assert next(x for x in w if x['horizon_h']==4)['confirmed_repeat']==original['retest4']
  assert next(x for x in w if x['horizon_h']==24)['confirmed_repeat']['numerator']==original['numeric_repeat_within_followup']
 for h in HORIZONS:
  all_w=next(x for x in windows if x['group']=='overall'and x['horizon_h']==h)
  subs=[x for x in windows if x['group']!='overall'and x['horizon_h']==h]
  assert sum(x['N']for x in subs)==all_w['N']
  assert all(sum(x['states'][state]for x in subs)==all_w['states'][state]for state in STATES)
 assert sum(x['N']for x in baseline[1:])==baseline[0]['N']
 assert all(sum(x['sex'][code]for x in baseline[1:])==baseline[0]['sex'][code]for code in ['M','F','Unknown'])
 assert all(sum(x['ICU'][unit]for x in baseline[1:])==baseline[0]['ICU'][unit]for unit in UNIT_ORDER)
 old_age=ctx['c'].anchor_age.ge(90)&ctx['c'].age.lt(90)
 assert not old_age.any(),'Older-age reporting needs an explicit interpretation review'
 before={x['file']:x for x in prior['source_integrity']};integrity=[]
 for name in m.FILES:
  file=m.ROOT/name;stat=file.stat();digest=m.sha(file)
  integrity.append({'file':name,'sha256':digest,'matches_formal_main':digest==before[name]['sha256'],'size_mtime_unchanged':ctx['stamps'][name]==(stat.st_size,stat.st_mtime_ns)})
 assert all(x['matches_formal_main']and x['size_mtime_unchanged']for x in integrity)
 result={'protocol_version':'1.0','protocol_sha256':m.sha(m.PROTOCOL),'execution_specification':json.loads(EXECUTION.read_text(encoding='utf-8')),'execution_specification_sha256':m.sha(EXECUTION),'analysis_code_sha256':m.sha(Path(__file__)),'formal_engine_sha256':m.sha(Path(m.__file__)),'cohort_selector_sha256':m.sha(Path(s.__file__)),'main_results_exactly_reproduced':True,'final_patients':len(idx),'eligible_admissions':len(admissions),'baseline':baseline,'time_windows':windows,'first_repeat':first,'flow':flow,'source_integrity':integrity,'quality':{'linked_lab_K_rows':len(k),'linked_charted_high_rows':len(ctx['ce']),'anchor_older_group_derived_age_lt90':int(old_age.sum()),'no_extra_clinical_variables':True,'no_models_or_P_values':True},'software':{'python':sys.version,'pandas':pd.__version__,'numpy':np.__version__},'provenance':'v1.0 was defined after feasibility analysis; this is not prospective preregistration','elapsed_analysis_seconds':round(time.monotonic()-start,2)}
 export(result)
 m.say('tables_analysis_complete',patients=len(idx),output=str(OUT),elapsed_s=round(time.monotonic()-start,1))
