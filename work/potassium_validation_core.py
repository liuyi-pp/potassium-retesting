"""Exact formal dependencies extracted from the historical validation module.
No data loader or exploratory analysis is included.
"""
import numpy as np
import pandas as pd

def unique_intervals(starts,ends,ids,times):
 lo=np.searchsorted(starts,times,side='right')-1;safe=np.clip(lo,0,len(starts)-1);valid=~np.isnat(times)
 hit=(lo>=0)&(times<=ends[safe])&valid
 prefix=np.maximum.accumulate(np.where(np.isnat(ends),np.datetime64('1678-01-01'),ends))
 earlier=(lo>0)&(times<=prefix[np.clip(safe-1,0,len(starts)-1)])&valid
 assigned=np.full(len(times),np.nan);assigned[hit&~earlier]=ids[safe[hit&~earlier]];amb=np.zeros(len(times),dtype=bool)
 for j in np.flatnonzero(earlier):
  z=np.flatnonzero((starts<=times[j])&(ends>=times[j]))
  if len(z)==1:assigned[j]=ids[z[0]]
  elif len(z)>1:amb[j]=True
 return assigned,amb

def classify_comments(s):
 s=s.fillna('').astype(str);low=s.str.lower()
 hem=low.str.contains('hemol|haemol',regex=True)
 old=low.str.contains('hemol|haemol|contam|cancel|disregard|clott|error',regex=True)
 neghem=low.str.contains(r'\b(?:no|without)\s+(?:evidence\s+of\s+|significant\s+|visible\s+)?h[ae]?emol|\bnot\s+h[ae]?emoly[sz]ed|h[ae]?emol\w*\s*[:=-]?\s*(?:negative|none|absent)',regex=True)
 # Explicit invalid/cancel context is deliberately narrower than any mention of 'error'.
 explicit=low.str.contains(r'\b(?:cancel\w*|disregard\w*|contaminat\w*|clott\w*|erroneous)\b|unable\s+to\s+(?:report|perform)|not\s+report(?:ed|able)|invalid\s+(?:sample|specimen|result)|specimen\s+(?:was\s+)?reject',regex=True)
 neginvalid=low.str.contains(r'\b(?:no|without)\s+(?:evidence\s+of\s+)?(?:contaminat\w*|clott\w*|error)|\bnot\s+(?:contaminat\w*|clott\w*)',regex=True)
 history=hem&low.str.contains(r'\b(?:previous|prior|earlier|original)\s+(?:sample|specimen|result)|(?:sample|specimen)\s+(?:previously|earlier)',regex=True)
 unaffected=hem&low.str.contains(r'not\s+(?:significantly\s+)?affected|no\s+(?:significant\s+)?interference|does\s+not\s+(?:affect|interfere)',regex=True)
 positive=hem&~neghem
 refined=positive|(explicit&~neginvalid)
 narrow=explicit&~neginvalid&~history
 return pd.DataFrame({'old_flag':old,'hemolysis_mention':hem,'negated_hemolysis':neghem,'historical_hemolysis_context':history,'hemolysis_unaffected_context':unaffected,'positive_hemolysis':positive,'explicit_invalid':explicit&~neginvalid,'refined_flag':refined,'narrow_flag':narrow,'comment_present':s.ne(''),'comment_deidentified':s.eq('___')},index=s.index)
