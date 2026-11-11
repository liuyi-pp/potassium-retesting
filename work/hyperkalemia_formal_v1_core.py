"""Formal v1.0 primary analysis. Restricted rows remain in memory; exports aggregate only."""
from pathlib import Path
import os
from concurrent.futures import ThreadPoolExecutor
from collections import Counter
import json,hashlib,time,sys,math,csv
import numpy as np
import pandas as pd
from potassium_validation_core import classify_comments,unique_intervals

HERE=Path(__file__).resolve().parent
ROOT=Path(os.environ.get('MIMIC_IV_ROOT', 'data/mimic-iv-3.1'))
OUT=HERE.parent/'outputs'/'hyperkalemia_v1_main_analysis'
PROTOCOL=HERE/'hyperkalemia_protocol_sap_v1_0.json'
FILES=['hosp/patients.csv.gz','hosp/admissions.csv.gz','icu/icustays.csv.gz','hosp/transfers.csv.gz','hosp/d_labitems.csv.gz','icu/d_items.csv.gz','hosp/labevents.csv.gz','icu/chartevents.csv.gz']
ICU={'Medical Intensive Care Unit (MICU)','Medical/Surgical Intensive Care Unit (MICU/SICU)','Cardiac Vascular Intensive Care Unit (CVICU)','Coronary Care Unit (CCU)','Surgical Intensive Care Unit (SICU)','Trauma SICU (TSICU)','Neuro Surgical Intensive Care Unit (Neuro SICU)','Intensive Care Unit (ICU)'}
CHEM={50971,52610};BG={50822,52452};CHART={220640,226535,227442,227464}
HOUR=3600000000000
NAT=np.iinfo(np.int64).min
def say(stage,**v):print(json.dumps({'stage':stage,**v}),flush=True)
def dt(x):return pd.to_datetime(x,errors='coerce')
def ns(x):return pd.to_datetime(x).to_numpy(dtype='datetime64[ns]').view('int64')
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb')as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def wilson(k,n):
 if n==0:return {'numerator':int(k),'denominator':0,'proportion':None,'ci95_lower':None,'ci95_upper':None}
 assert 0<=k<=n
 z=1.959963984540054;p=k/n;den=1+z*z/n
 mid=(p+z*z/(2*n))/den;half=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
 return {'numerator':int(k),'denominator':int(n),'proportion':p,'ci95_lower':max(0.,mid-half),'ci95_upper':min(1.,mid+half)}
def add_flow(flow,stage,x,previous=None):
 n=int(x.hadm_id.nunique()) if 'hadm_id'in x.columns else int(x.index.nunique())
 people=int(x.subject_id.nunique())
 row={'stage':stage,'remaining_admissions':n,'remaining_patients':people,'excluded_admissions_this_step':None if previous is None else int(previous-n)}
 assert row['excluded_admissions_this_step']is None or row['excluded_admissions_this_step']>=0
 flow.append(row);return n

def load_data():
 start=time.monotonic();qc=Counter();stamps={n:(ROOT.joinpath(n).stat().st_size,ROOT.joinpath(n).stat().st_mtime_ns)for n in FILES}
 metadata={};flow=[]
 pa=pd.read_csv(ROOT/FILES[0],usecols=['subject_id','gender','anchor_age','anchor_year'])
 ad=pd.read_csv(ROOT/FILES[1],usecols=['subject_id','hadm_id','admittime','dischtime','deathtime'])
 for s in ['admittime','dischtime','deathtime']:ad[s]=dt(ad[s])
 all_icu=pd.read_csv(ROOT/FILES[2]);all_icu['intime']=dt(all_icu.intime);all_icu['outtime']=dt(all_icu.outtime)
 metadata.update(raw_ICU_records=len(all_icu),raw_ICU_admissions=int(all_icu.hadm_id.nunique()),raw_ICU_patients=int(all_icu.subject_id.nunique()))
 first=all_icu.sort_values(['hadm_id','intime','stay_id']).drop_duplicates('hadm_id')
 prev=add_flow(flow,'01_first_ICU_record_per_admission',first)
 c=first.merge(ad,on=['subject_id','hadm_id'],how='left',validate='one_to_one').merge(pa,on='subject_id',how='left',validate='many_to_one')
 c=c[c.admittime.notna()&c.anchor_age.notna()&c.anchor_year.notna()].copy();prev=add_flow(flow,'02_linked_patient_admission_metadata',c,prev)
 c['age']=c.anchor_age+c.admittime.dt.year-c.anchor_year
 c=c[c.age.ge(18)].copy();prev=add_flow(flow,'03_adult_age_ge18',c,prev)
 c=c[c.intime.notna()&c.outtime.gt(c.intime)&c.dischtime.gt(c.intime)&(c.deathtime.isna()|c.deathtime.gt(c.intime))].copy()
 c['stop']=c[['outtime','dischtime','deathtime']].min(axis=1);c=c[c.stop.gt(c.intime)].copy();prev=add_flow(flow,'04_valid_positive_time_and_alive_at_ICU_entry',c,prev)
 c=c.set_index('hadm_id');assert c.index.is_unique
 H=set(c.index);S=set(c.subject_id)
 tr=pd.read_csv(ROOT/FILES[3]);tr=tr[tr.hadm_id.isin(H)].copy();tr['intime']=dt(tr.intime);tr['outtime']=dt(tr.outtime)
 ld=pd.read_csv(ROOT/FILES[4]);di=pd.read_csv(ROOT/FILES[5])
 assert set(ld.loc[ld.label.eq('Potassium')&ld.fluid.eq('Blood')&ld.category.eq('Chemistry'),'itemid'])==CHEM
 assert set(ld.loc[ld.label.eq('Potassium, Whole Blood')&ld.fluid.eq('Blood')&ld.category.eq('Blood Gas'),'itemid'])==BG
 assert set(di.loc[di.linksto.eq('chartevents')&di.label.isin(['Potassium (serum)','Potassium (whole blood)','ZPotassium (serum)','ZPotassium (whole blood)']),'itemid'])==CHART
 stay_map=all_icu[['stay_id','subject_id','hadm_id']].rename(columns={'subject_id':'expected_subject','hadm_id':'expected_hadm'})
 target_stays=set(stay_map.loc[stay_map.expected_hadm.isin(H),'stay_id'])
 say('population',first_adult_admissions=len(c),patients=int(c.subject_id.nunique()),protocol_sha256=sha(PROTOCOL))
 def scan_lab():
  a=[];n=0;cols=['labevent_id','subject_id','hadm_id','specimen_id','itemid','charttime','storetime','valuenum','valueuom','comments']
  dtype={'labevent_id':'int64','subject_id':'int32','hadm_id':'float64','specimen_id':'float64','itemid':'int32','valuenum':'float64'}
  for z in pd.read_csv(ROOT/FILES[6],usecols=cols,dtype=dtype,chunksize=2000000,low_memory=False):
   n+=len(z);z=z[z.itemid.isin(CHEM|BG)&(z.hadm_id.isin(H)|(z.hadm_id.isna()&z.subject_id.isin(S)))].copy()
   if len(z):a.append(z)
   if n%40000000==0:say('lab_scan',rows=n,elapsed_s=round(time.monotonic()-start))
  say('lab_scan_complete',rows=n);return pd.concat(a,ignore_index=True),n
 def scan_chart():
  a=[];n=0;cols=['subject_id','hadm_id','stay_id','charttime','storetime','itemid','valuenum','valueuom']
  dtype={'subject_id':'int32','hadm_id':'float64','stay_id':'int64','itemid':'int32','valuenum':'float64'}
  for z in pd.read_csv(ROOT/FILES[7],usecols=cols,dtype=dtype,chunksize=2000000,low_memory=False):
   n+=len(z);z=z[z.itemid.isin(CHART)&z.valuenum.ge(5.5)&z.valuenum.le(30)&(z.hadm_id.isin(H)|z.stay_id.isin(target_stays))].copy()
   if len(z):a.append(z)
   if n%60000000==0:say('chart_scan',rows=n,elapsed_s=round(time.monotonic()-start))
  say('chart_scan_complete',rows=n);return pd.concat(a,ignore_index=True),n
 with ThreadPoolExecutor(max_workers=2)as ex:
  lf=ex.submit(scan_lab);cf=ex.submit(scan_chart);k,nlab=lf.result();ce,nchart=cf.result()
 k['charttime']=dt(k.charttime);k['storetime']=dt(k.storetime)
 amap={s:g.sort_values('admittime')for s,g in ad[ad.subject_id.isin(S)].groupby('subject_id',sort=False)}
 for subject,g in k[k.hadm_id.isna()].groupby('subject_id',sort=False):
  m=amap[subject];ids,amb=unique_intervals(m.admittime.values,m.dischtime.values,m.hadm_id.values,g.charttime.values)
  k.loc[g.index,'hadm_id']=ids;qc['lab_missing_hadm_repaired_rows']+=int(np.isfinite(ids).sum());qc['lab_ambiguous_admission_rows']+=int(amb.sum())
 k=k[k.hadm_id.isin(H)].copy();k['hadm_id']=k.hadm_id.astype('int64')
 m=c.reset_index()[['hadm_id','subject_id','intime','stop']]
 kn=len(k);k=k.merge(m,on=['hadm_id','subject_id'],validate='many_to_one');qc['lab_identity_mismatch_rows']=kn-len(k)
 k['source']=np.where(k.itemid.isin(CHEM),'chemistry','bloodgas');k['source_rank']=np.where(k.source.eq('chemistry'),0,1)
 u=k.valueuom.fillna('').str.lower().str.replace(' ','',regex=False)
 k['unit_ok']=u.isin(['mmol/l','meq/l']);k['numeric_ok']=k.valuenum.gt(0)&k.valuenum.le(30)&k.unit_ok
 k['v']=k.valuenum;k['valid_available']=k.charttime.notna()&k.storetime.notna()&k.storetime.ge(k.charttime)
 k['original_window']=k.charttime.ge(k.intime)&k.charttime.lt(k.stop)
 k=k.join(classify_comments(k.comments));k=k.drop(columns=['comments'])
 keys=['subject_id','hadm_id','specimen_id','source']
 present=k[k.specimen_id.notna()]
 st=present.groupby(keys,sort=False).agg(n=('labevent_id','size'),v_n=('v','nunique'),time_n=('charttime','nunique'),any_refined=('refined_flag','max'),any_explicit=('explicit_invalid','max')).reset_index()
 st['spec_conflict']=st.n.gt(1)&(st.v_n.gt(1)|st.time_n.gt(1));qc['duplicate_specimen_source_groups']=int(st.n.gt(1).sum());qc['conflicting_specimen_source_groups']=int(st.spec_conflict.sum())
 k=k.merge(st[keys+['spec_conflict','any_refined','any_explicit']],on=keys,how='left',validate='many_to_one')
 k['spec_conflict']=k.spec_conflict.fillna(False).astype(bool)
 k['refined_flag']=k.refined_flag|k.any_refined.fillna(False).astype(bool);k['explicit_invalid']=k.explicit_invalid|k.any_explicit.fillna(False).astype(bool)
 # Keep missing-specimen rows separate. Earliest valid availability wins only for identical duplicates.
 has=k.specimen_id.notna();ded=k[has].sort_values(['valid_available','storetime','labevent_id'],ascending=[False,True,True]).drop_duplicates(keys)
 k=pd.concat([ded,k[~has]],ignore_index=True)
 ce=ce.merge(stay_map,on='stay_id',how='left',validate='many_to_one')
 mismatch=ce.expected_hadm.isna()|ce.subject_id.ne(ce.expected_subject)|(ce.hadm_id.notna()&ce.hadm_id.ne(ce.expected_hadm))
 qc['chart_identity_mismatch_rows']=int(mismatch.sum());ce=ce[~mismatch].copy();ce['hadm_id']=ce.expected_hadm.astype('int64');ce=ce[ce.hadm_id.isin(H)].copy()
 ce['charttime']=dt(ce.charttime);ce['storetime']=dt(ce.storetime)
 cu=ce.valueuom.fillna('').str.lower().str.replace(' ','',regex=False);ce['unit_ok']=cu.isin(['mmol/l','meq/l'])
 qc['charted_high_unknown_unit_rows']=int((~ce.unit_ok).sum());ce=ce[ce.unit_ok].copy()
 ce['valid_available']=ce.charttime.notna()&ce.storetime.notna()&ce.storetime.ge(ce.charttime)
 ce=ce.drop_duplicates(['subject_id','hadm_id','itemid','charttime','storetime','valuenum'])
 qc['retained_linked_lab_K_rows']=len(k);qc['retained_linked_chart_high_rows']=len(ce)
 say('data_loaded',lab_rows=nlab,chart_rows=nchart,retained_lab=len(k),retained_chart_high=len(ce),elapsed_s=round(time.monotonic()-start))
 return dict(k=k,ce=ce,c=c,tr=tr,flow=flow,metadata=metadata,qc=qc,start=start,stamps=stamps,nlab=nlab,nchart=nchart,load_code_sha256=sha(Path(__file__)),protocol_sha256=sha(PROTOCOL))

def fast_intervals(c,tr,hadms):
 z=tr[tr.hadm_id.isin(hadms)&tr.intime.notna()&tr.outtime.notna()&tr.outtime.gt(tr.intime)].copy()
 z=z.merge(c.reset_index()[['hadm_id','subject_id','intime','stop']].rename(columns={'intime':'first_start'}),on=['hadm_id','subject_id'],validate='many_to_one')
 z=z[z.outtime.gt(z.first_start)&z.intime.lt(z.stop)].copy()
 z['lo']=z[['intime','first_start']].max(axis=1);z['hi']=z[['outtime','stop']].min(axis=1);z=z[z.hi.gt(z.lo)].copy()
 bounds={h:(np.sort(ns(g.lo)),np.sort(ns(g.hi)))for h,g in z.groupby('hadm_id',sort=False)}
 t=z[z.careunit.isin(ICU)].sort_values(['hadm_id','lo','hi','transfer_id']).copy();t['hadm_id']=t.hadm_id.astype('int64')
 # A location overlap interrupts verified continuity, including after t0.
 # Sweep all transfers so that a later ambiguous location cannot extend follow-up.
 episodes=[]
 for h,g in z.groupby('hadm_id',sort=False):
  lo=ns(g.lo);hi=ns(g.hi);true=g.careunit.isin(ICU).to_numpy(dtype=int)
  points,inv=np.unique(np.r_[lo,hi],return_inverse=True);n=len(lo)
  delta=np.zeros(len(points),dtype=int);icudelta=delta.copy()
  np.add.at(delta,inv[:n],1);np.add.at(delta,inv[n:],-1)
  np.add.at(icudelta,inv[:n],true);np.add.at(icudelta,inv[n:],-true)
  good=(np.cumsum(delta)[:-1]==1)&(np.cumsum(icudelta)[:-1]==1)
  starts=np.flatnonzero(good&~np.r_[False,good[:-1]])
  ends=np.flatnonzero(good&~np.r_[good[1:],False])+1
  for a,b in zip(starts,ends):episodes.append((int(h),pd.Timestamp(points[a]),pd.Timestamp(points[b])))
 ep=pd.DataFrame(episodes,columns=['hadm_id','episode_start','episode_stop'])
 return ep,bounds,t

def attach_location(x,ep,bounds,t):
 x=x.copy();x['candidate_row']=np.arange(len(x));x['hadm_id']=x.hadm_id.astype('int64')
 if not len(ep):return x.iloc[0:0].assign(episode_start=pd.NaT,episode_stop=pd.NaT,actual_unit='')
 e=ep[['hadm_id','episode_start','episode_stop']].sort_values('episode_start')
 e=e.copy();e['episode_start']=e.episode_start.astype(x.storetime.dtype)
 x=pd.merge_asof(x.sort_values('storetime'),e,left_on='storetime',right_on='episode_start',by='hadm_id',direction='backward')
 valid=x.episode_start.notna()&x.storetime.lt(x.episode_stop)&x.charttime.ge(x.episode_start)&x.charttime.lt(x.episode_stop)
 x=x[valid].copy();unclear=pd.Series(False,index=x.index)
 for h,g in x.groupby('hadm_id',sort=False):
  a,b=bounds[h];ct=ns(g.charttime);st=ns(g.storetime)
  nc=np.searchsorted(a,ct,side='right')-np.searchsorted(b,ct,side='right')
  nr=np.searchsorted(a,st,side='right')-np.searchsorted(b,st,side='right')
  unclear.loc[g.index]=(nc!=1)|(nr!=1)
 x=x[~unclear].copy()
 units=t[['hadm_id','lo','hi','careunit']].rename(columns={'lo':'unit_start','hi':'unit_end','careunit':'actual_unit'}).sort_values('unit_start')
 x=pd.merge_asof(x.sort_values('storetime'),units,left_on='storetime',right_on='unit_start',by='hadm_id',direction='backward')
 x=x[x.storetime.lt(x.unit_end)&x.actual_unit.isin(ICU)]
 return x.sort_values(['hadm_id','storetime','charttime','labevent_id'])

def earlier_flags(index,k,ce):
 h=k[k.numeric_ok&k.v.ge(5.5)]
 known=h[h.valid_available].groupby('hadm_id',sort=False).storetime.min()
 unknown=h[~h.valid_available].copy();unknown['possible_time']=unknown[['charttime','storetime']].min(axis=1)
 unmin=unknown.groupby('hadm_id',sort=False).possible_time.min()
 cknown=ce[ce.valid_available].groupby('hadm_id',sort=False).storetime.min()
 cu=ce[~ce.valid_available].copy();cu['possible_time']=cu[['charttime','storetime']].min(axis=1);cumin=cu.groupby('hadm_id',sort=False).possible_time.min()
 out=index.copy()
 for name,series in [('earlier_lab_high',known),('earlier_chart_high',cknown),('uncertain_earlier_lab_high',unmin),('uncertain_earlier_chart_high',cumin)]:
  out[name]=out.hadm_id.map(series).lt(out.storetime)if len(series)else False
 return out

def four_hour_state(event_h,unknown_before4,follow_h,end_reason):
 if event_h is not None and event_h<=4:return 'repeat_recorded'
 if unknown_before4:return 'indeterminate'
 if follow_h<=4:return 'death_before_repeat'if end_reason=='death'else 'interval_end_before_repeat'
 return 'still_in_interval_without_repeat'

def analyze_event(r,g):
 t0=r.storetime;clinical_end=min(r.episode_stop,r.dischtime,r.deathtime if pd.notna(r.deathtime)else r.episode_stop)
 end_reason='interval_end'
 if pd.notna(r.deathtime)and r.deathtime==clinical_end:end_reason='death'
 elif r.dischtime==clinical_end:end_reason='hospital_discharge'
 cap=t0+pd.Timedelta(hours=24)
 ct=ns(g.charttime);st=ns(g.storetime);zero=t0.value;end=clinical_end.value;capns=cap.value
 spec=g.specimen_id.to_numpy(dtype=float);has_spec=np.isfinite(spec);distinct=(~has_spec)|(spec!=r.specimen_id)
 post=(ct>zero)&(ct<=capns)&(ct<end)&distinct
 time_conf=(st!=NAT)&(ct!=NAT)&(st<ct)
 numeric=g.numeric_ok.to_numpy(dtype=bool);explicit=g.explicit_invalid.to_numpy(dtype=bool);conflict=g.spec_conflict.to_numpy(dtype=bool)
 good=post&numeric&has_spec&~explicit&~time_conf&~conflict
 gi=np.flatnonzero(good)
 event_h=float((ct[gi[0]]-zero)/HOUR)if len(gi)else None
 firstcand=np.flatnonzero(post&numeric)
 value_status='no_numeric_repeat';first_value=None;repeat_source=None
 firsttime=None
 if len(firstcand):
  q=firstcand[0];firsttime=ct[q];tie=np.flatnonzero(post&numeric&(ct==firsttime))
  vv=g.v.to_numpy(dtype=float);diff=float(np.max(vv[tie])-np.min(vv[tie]))>.5
  unsupported=(not has_spec[q])or explicit[q]or time_conf[q]or conflict[q]or bool(g.refined_flag.iloc[q])or diff
  repeat_source=g.source.iloc[q]
  if unsupported:value_status='indeterminate_first_value'
  else:first_value=float(vv[q]);value_status='high'if first_value>=5.5 else 'below_high_threshold'
 # Plausible failed/ambiguous measurement before the horizon remains an unknown category.
 horizon=min(zero+4*HOUR,end)
 post4=distinct&(ct>zero)&(ct<=zero+4*HOUR)&(ct<end)
 uncertain_rows=post4&(~numeric|~has_spec|explicit|time_conf|conflict)
 missing_or_conflict_time=(ct==NAT)|time_conf
 evidence_by_store=distinct&missing_or_conflict_time&(st>zero)&(st<=zero+4*HOUR)&(st<end)
 unknown_before4=bool(np.any(uncertain_rows|evidence_by_store))
 state=four_hour_state(event_h,unknown_before4,(clinical_end-t0).total_seconds()/3600,end_reason)
 # A confirmed later event cannot retrospectively identify an earlier ambiguous first value.
 if value_status in ['high','below_high_threshold'] and firsttime is not None:
  earlier_uncertain=(distinct&(ct>zero)&(ct<firsttime)&(ct<end)&(explicit|time_conf|conflict|(numeric&~has_spec)))
  earlier_unknown_time=distinct&missing_or_conflict_time&(st>zero)&(st<firsttime)&(st<end)
  if np.any(earlier_uncertain|earlier_unknown_time):value_status='indeterminate_first_value';first_value=None
 return {'subject_id':int(r.subject_id),'hadm_id':int(r.hadm_id),'index_id':int(r.labevent_id),'index_K':float(r.v),'stratum':'5.5_to_lt6'if r.v<6 else '6_to_lt6.5'if r.v<6.5 else 'ge6.5','event_h':event_h,'follow_h':(clinical_end-t0).total_seconds()/3600,'state4':state,'end_reason':end_reason,'first_value_status':value_status,'first_value':first_value,'repeat_source':repeat_source,'unknown_before4':unknown_before4,'nonnumeric_or_failed_post_records':int(np.sum(post&(~numeric|explicit))),'same_t0_distinct_records':int(np.sum(distinct&(ct==zero))),'numeric_repeat_time_conflicts':int(np.sum(post&numeric&time_conf))}

def summarize(d,label):
 n=len(d);counts=d.state4.value_counts().to_dict();tested=int(counts.get('repeat_recorded',0));unk=int(counts.get('indeterminate',0))
 supported=d.first_value_status.isin(['high','below_high_threshold']);pos=int(d.first_value_status.eq('high').sum());den=int(supported.sum())
 assert sum(counts.values())==n and pos<=den<=int(d.event_h.notna().sum())<=n
 out={'group':label,'n_patients':n,'retest4':wilson(tested,n),'states4':{key:int(counts.get(key,0))for key in ['repeat_recorded','death_before_repeat','interval_end_before_repeat','still_in_interval_without_repeat','indeterminate']},'retest4_identification_bounds':[tested/n,(tested+unk)/n]if n else [None,None],'first_repeat_still_high':wilson(pos,den),'first_repeat_value_status':{'supported_high':pos,'supported_below_threshold':int(d.first_value_status.eq('below_high_threshold').sum()),'indeterminate':int(d.first_value_status.eq('indeterminate_first_value').sum()),'no_numeric_repeat':int(d.first_value_status.eq('no_numeric_repeat').sum())},'numeric_repeat_within_followup':int(d.event_h.notna().sum())}
 assert sum(out['first_repeat_value_status'].values())==n
 times=d.event_h.dropna().to_numpy(dtype=float)
 out['observed_repeat_interval_h']={'n':len(times),'median':float(np.median(times))if len(times)else None,'q25':float(np.quantile(times,.25))if len(times)else None,'q75':float(np.quantile(times,.75))if len(times)else None}
 full=d.follow_h.gt(4);none=full&(~d.event_h.le(4));out['supplementary_no_repeat_with_full4h']={'numerator':int(none.sum()),'denominator':int(full.sum()),'unknown_evidence_in_numerator':int((none&d.unknown_before4).sum())}
 return out

def run(ctx):
 assert sha(PROTOCOL)==ctx['protocol_sha256'],'Protocol changed during extraction'
 k=ctx['k'];ce=ctx['ce'];c=ctx['c'];flow=[dict(z)for z in ctx['flow']];qc=ctx['qc'].copy()
 prev=flow[-1]['remaining_admissions']
 x=k[k.source.eq('chemistry')&k.numeric_ok&k.v.ge(5.5)&k.original_window].copy();prev=add_flow(flow,'05_numeric_chemistry_high_in_first_ICU_record_window',x,prev)
 x=x[x.specimen_id.notna()&~x.spec_conflict].copy();prev=add_flow(flow,'06_identifiable_nonconflicting_index_specimen',x,prev)
 x=x[x.valid_available&x.storetime.lt(x.stop)].copy();prev=add_flow(flow,'07_valid_index_availability_within_record_window',x,prev)
 x=x[~x.refined_flag].copy();prev=add_flow(flow,'08_negation_aware_index_comment_rule',x,prev)
 say('constructing_ICU_intervals',candidate_admissions=prev,elapsed_s=round(time.monotonic()-ctx['start']))
 ep,bounds,t=fast_intervals(c,ctx['tr'],set(x.hadm_id))
 x=attach_location(x,ep,bounds,t);prev=add_flow(flow,'09_sample_and_result_in_same_verified_continuous_ICU',x,prev)
 idx=x.sort_values(['hadm_id','storetime','charttime','labevent_id']).drop_duplicates('hadm_id').copy()
 idx=earlier_flags(idx,k,ce)
 idx=idx[~idx.earlier_lab_high].copy();prev=add_flow(flow,'10_no_earlier_numeric_lab_high_regardless_comment',idx,prev)
 idx=idx[~idx.earlier_chart_high].copy();prev=add_flow(flow,'11_no_earlier_charted_high',idx,prev)
 idx=idx[~(idx.uncertain_earlier_lab_high|idx.uncertain_earlier_chart_high)].copy();prev=add_flow(flow,'12_no_unresolvable_prior_high_timing',idx,prev)
 admission_n=len(idx);patient_n=int(idx.subject_id.nunique())
 idx=idx.sort_values(['subject_id','storetime','hadm_id','labevent_id']).drop_duplicates('subject_id').copy();prev=add_flow(flow,'13_first_qualifying_event_per_patient',idx,prev)
 assert len(idx)==patient_n and idx.subject_id.is_unique and idx.hadm_id.is_unique
 idx=idx.merge(c.reset_index()[['hadm_id','gender','age','dischtime','deathtime']],on='hadm_id',validate='one_to_one')
 say('cohort_final',eligible_admissions=admission_n,unique_patients=len(idx),elapsed_s=round(time.monotonic()-ctx['start']))
 groups={h:g for h,g in k[k.hadm_id.isin(set(idx.hadm_id))].sort_values(['charttime','source_rank','labevent_id']).groupby('hadm_id',sort=False)}
 transfer_groups={h:g for h,g in ctx['tr'][ctx['tr'].hadm_id.isin(set(idx.hadm_id))].groupby('hadm_id',sort=False)}
 rows=[]
 for i,r in enumerate(idx.itertuples(),1):
  item=analyze_event(r,groups[r.hadm_id]);tg=transfer_groups[r.hadm_id]
  at_end=tg[tg.intime.le(r.episode_stop)&tg.outtime.gt(r.episode_stop)&tg.subject_id.eq(r.subject_id)]
  item['interval_end_detail']='overlapping_locations'if len(at_end)>1 else 'unverified_location'if len(at_end)==0 else 'first_ICU_record_boundary'if at_end.careunit.iloc[0]in ICU else 'documented_non_ICU_transfer'
  rows.append(item)
  if i%1500==0:say('outcomes_progress',patients=i,elapsed_s=round(time.monotonic()-ctx['start']))
 d=pd.DataFrame(rows);assert d.subject_id.is_unique and len(d)==len(idx)
 summaries=[summarize(d,'overall')]+[summarize(d[d.stratum.eq(z)],z)for z in ['5.5_to_lt6','6_to_lt6.5','ge6.5']]
 for key in ['n_patients','numeric_repeat_within_followup']:
  assert sum(a[key]for a in summaries[1:])==summaries[0][key]
 for key in ['retest4','first_repeat_still_high']:
  for count in ['numerator','denominator']:assert sum(a[key][count]for a in summaries[1:])==summaries[0][key][count]
 for key in summaries[0]['states4']:assert sum(a['states4'][key]for a in summaries[1:])==summaries[0]['states4'][key]
 qc.update({'primary_index_ge10':int(idx.v.ge(10).sum()),'nonnumeric_or_failed_post_records':int(d.nonnumeric_or_failed_post_records.sum()),'patients_with_nonnumeric_or_failed_post_record':int(d.nonnumeric_or_failed_post_records.gt(0).sum()),'same_t0_distinct_records':int(d.same_t0_distinct_records.sum()),'numeric_repeat_time_conflicts':int(d.numeric_repeat_time_conflicts.sum())})
 qc['interval_end_before_repeat_4h_detail']={str(a):int(b)for a,b in d.loc[d.state4.eq('interval_end_before_repeat'),'interval_end_detail'].value_counts().items()}
 manifest={}
 for line in(ROOT/'SHA256SUMS.txt').read_text().splitlines():
  try:h,n=line.split(maxsplit=1);manifest[n.lstrip('*').replace('\\','/').lstrip('./')]=h
  except ValueError:pass
 integrity=[]
 for name in FILES:
  dig=sha(ROOT/name);s=ROOT.joinpath(name).stat();integrity.append({'file':name,'sha256':dig,'release_match':manifest.get(name)==dig,'size_mtime_unchanged':ctx['stamps'][name]==(s.st_size,s.st_mtime_ns)})
 assert all(z['release_match']and z['size_mtime_unchanged']for z in integrity)
 result={'analysis':'Formal v1.0 main analysis only, no additional clinical variables or models','protocol_version':'1.0','protocol_sha256':sha(PROTOCOL),'load_code_sha256':ctx['load_code_sha256'],'analysis_code_sha256':sha(Path(__file__)),'comment_rules_code_sha256':sha(HERE/'potassium_validation_core.py'),'source_inventory':ctx['metadata'],'rows_scanned':{'labevents':ctx['nlab'],'chartevents':ctx['nchart']},'flow':flow,'eligible_admissions_before_patient_selection':admission_n,'final_patients':len(idx),'summaries':summaries,'qc':dict(qc),'source_integrity':integrity,'elapsed_seconds':round(time.monotonic()-ctx['start'],1),'primary_denominator':'all enrolled unique patients; early deaths and verified ICU interval ends retained as separate states','secondary_denominator':'patients with an interpretable first new-specimen numeric repeat during <=24h or earlier interval end; unavailable/flagged first values not replaced by later values','provenance':'Protocol developed after feasibility work; these estimates are descriptive, not prospective preregistered confirmation','execution_notes':['Earlier high from the same specimen but another assay is included in prior-history exclusion.','No prior high with suspect comments is bypassed to select a later clean index.','Earlier chart high is checked across all stays linked to the admission; valid timestamps need not lie inside the first ICU interval.','Missing/inverted earlier high timestamps are conservatively classified using the earliest available time stamp as potential pre-index evidence; not imputed.','Explicitly invalid/nonnumeric or unresolved repeat evidence cannot silently be labelled no sampling. Unknown states and identification bounds are reported.','Sampling at an exact death or ICU-exit time is not counted before that end.','No sensitivity analyses or clinical-effect models were added in this main-analysis run.']}
 result['software']={'python':sys.version,'pandas':pd.__version__,'numpy':np.__version__}
 result['execution_notes'].append('Any transfer-location overlap interrupts verified ICU continuity, including an overlap starting after the index.')
 OUT.mkdir(parents=True,exist_ok=True)
 dest=OUT/'formal_v1_main_aggregate.json';tmp=OUT/'formal_v1_main_aggregate.tmp'
 tmp.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8');tmp.replace(dest)
 with open(OUT/'cohort_flow.csv','w',newline='',encoding='utf-8-sig')as f:
  w=csv.DictWriter(f,fieldnames=['stage','excluded_admissions_this_step','remaining_admissions','remaining_patients']);w.writeheader();w.writerows(flow)
 flat=[]
 for z in summaries:
  a=z['retest4'];b=z['first_repeat_still_high'];flat.append({'group':z['group'],'N_patients':z['n_patients'],'repeat4h_n':a['numerator'],'repeat4h_denominator':a['denominator'],'repeat4h_pct':100*a['proportion']if a['proportion']is not None else None,'repeat4h_CI95_lower_pct':100*a['ci95_lower']if a['ci95_lower']is not None else None,'repeat4h_CI95_upper_pct':100*a['ci95_upper']if a['ci95_upper']is not None else None,'first_repeat_high_n':b['numerator'],'first_repeat_interpretable_n':b['denominator'],'first_repeat_high_pct':100*b['proportion']if b['proportion']is not None else None,'first_repeat_high_CI95_lower_pct':100*b['ci95_lower']if b['ci95_lower']is not None else None,'first_repeat_high_CI95_upper_pct':100*b['ci95_upper']if b['ci95_upper']is not None else None,**z['states4']})
 with open(OUT/'primary_results.csv','w',newline='',encoding='utf-8-sig')as f:
  w=csv.DictWriter(f,fieldnames=list(flat[0]));w.writeheader();w.writerows(flat)
 say('formal_complete',final_patients=len(idx),eligible_admissions=admission_n,output=str(OUT),elapsed_s=round(time.monotonic()-ctx['start']))
