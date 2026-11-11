"""All examples below are artificial software inputs, not MIMIC patient records."""
from pathlib import Path
from types import SimpleNamespace
import math,sys,unittest
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'work'))
import hyperkalemia_formal_v1_core as m
import hyperkalemia_tables_v1_core as a
from potassium_validation_core import classify_comments,unique_intervals

T=pd.Timestamp('2000-01-01')  # Deliberately invented test clock.

def index(end=12,death=False):
 return SimpleNamespace(storetime=T,episode_stop=T+pd.Timedelta(hours=end),dischtime=T+pd.Timedelta(hours=48),deathtime=T+pd.Timedelta(hours=end)if death else pd.NaT,specimen_id=1,subject_id=1,hadm_id=1,labevent_id=1,v=6.0)

def record(hours,value=5.8,specimen=2,source='chemistry',flag=False,invalid=False,conflict=False,store_hours=None):
 return dict(charttime=T+pd.Timedelta(hours=hours)if hours is not None else pd.NaT,storetime=T+pd.Timedelta(hours=store_hours if store_hours is not None else(hours+.1 if hours is not None else 1)),specimen_id=specimen,numeric_ok=math.isfinite(value),explicit_invalid=invalid,spec_conflict=conflict,refined_flag=flag,source=source,source_rank=0 if source=='chemistry'else 1,labevent_id=2,v=value)

def event(records,r=None):
 g=pd.DataFrame(records).sort_values(['charttime','source_rank','labevent_id']).reset_index(drop=True)
 r=r or index();return m.analyze_event(r,g),g,r

class ArtificialRules(unittest.TestCase):
 def test_negation_and_explicit_invalidity(self):
  z=classify_comments(pd.Series(['no evidence of hemolysis','hemolyzed','specimen rejected','no contamination','___',None]))
  self.assertEqual(z.refined_flag.tolist(),[False,True,True,False,False,False])
 def test_historical_hemolysis_remains_primary_warning(self):
  z=classify_comments(pd.Series(['prior specimen hemolyzed']))
  self.assertTrue(z.refined_flag.iloc[0]);self.assertFalse(z.explicit_invalid.iloc[0])
 def test_exact_horizon_counts_if_observation_continues(self):
  e,_,_=event([record(4)])
  self.assertEqual(e['state4'],'repeat_recorded')
 def test_exact_interval_end_does_not_count(self):
  e,_,_=event([record(4)],index(end=4))
  self.assertEqual(e['state4'],'interval_end_before_repeat');self.assertIsNone(e['event_h'])
 def test_death_tie_ends_observation(self):
  e,_,_=event([record(4)],index(end=4,death=True))
  self.assertEqual(e['state4'],'death_before_repeat')
 def test_same_t0_is_not_post_result_repeat(self):
  e,_,_=event([record(0)])
  self.assertEqual(e['same_t0_distinct_records'],1);self.assertIsNone(e['event_h'])
 def test_same_index_specimen_does_not_count(self):
  e,_,_=event([record(1,specimen=1)])
  self.assertIsNone(e['event_h'])
 def test_precollected_later_report_is_not_new_sampling(self):
  e,_,_=event([record(-1,store_hours=1)])
  self.assertIsNone(e['event_h'])
 def test_missing_specimen_remains_indeterminate(self):
  e,_,_=event([record(1,specimen=np.nan)])
  self.assertEqual(e['state4'],'indeterminate');self.assertEqual(e['first_value_status'],'indeterminate_first_value')
 def test_flagged_numeric_sampling_differs_from_value_interpretability(self):
  e,_,_=event([record(1,flag=True)])
  self.assertEqual(e['state4'],'repeat_recorded');self.assertEqual(e['first_value_status'],'indeterminate_first_value')
 def test_later_clean_value_does_not_replace_first_warning(self):
  e,_,_=event([record(1,flag=True),record(2,specimen=3,value=4.5)])
  self.assertEqual(e['first_value_status'],'indeterminate_first_value');self.assertIsNone(e['first_value'])
 def test_earlier_failed_record_prevents_later_first_value_replacement(self):
  e,_,_=event([record(1,invalid=True),record(2,specimen=3,value=4.5)])
  self.assertEqual(e['state4'],'repeat_recorded');self.assertEqual(e['first_value_status'],'indeterminate_first_value')
 def test_simultaneous_discordant_assays_leave_value_indeterminate(self):
  e,_,_=event([record(1,value=6.2),record(1,value=4.8,specimen=3,source='bloodgas')])
  self.assertEqual(e['state4'],'repeat_recorded');self.assertEqual(e['first_value_status'],'indeterminate_first_value')
 def test_bloodgas_can_establish_sampling(self):
  e,_,_=event([record(1,source='bloodgas')])
  self.assertEqual(e['state4'],'repeat_recorded');self.assertEqual(e['repeat_source'],'bloodgas')
 def test_missing_collection_time_can_be_unknown(self):
  e,_,_=event([record(None,store_hours=1)])
  self.assertEqual(e['state4'],'indeterminate')
 def test_cumulative_horizons_remain_distinct(self):
  e,g,r=event([record(5)])
  states=a.states_at_horizons(r,g,e)
  self.assertEqual(states['state_2h'],'still_in_interval_without_repeat')
  self.assertEqual(states['state_4h'],'still_in_interval_without_repeat')
  self.assertEqual(states['state_6h'],'repeat_recorded');self.assertEqual(states['state_24h'],'repeat_recorded')
 def test_wilson_boundaries(self):
  for k,n in [(0,10),(10,10),(3,10)]:
   q=m.wilson(k,n);self.assertLessEqual(q['ci95_lower'],k/n+1e-14);self.assertGreaterEqual(q['ci95_upper'],k/n-1e-14)
  self.assertIsNone(m.wilson(0,0)['proportion'])
 def test_ambiguous_admission_is_not_assigned(self):
  starts=np.array(['2000-01-01','2000-01-02'],dtype='datetime64[ns]')
  ends=np.array(['2000-01-04','2000-01-05'],dtype='datetime64[ns]')
  values,amb=unique_intervals(starts,ends,np.array([1,2]),np.array(['2000-01-03'],dtype='datetime64[ns]'))
  self.assertTrue(amb[0]);self.assertTrue(np.isnan(values[0]))

if __name__=='__main__':unittest.main()
