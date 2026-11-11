"""Supplementary File 1. Exact v1.0 comment classifier and item/unit definitions.
This excerpt does not load data or contain patient comments. Requires pandas.
The clinical selection and outcome algorithms are described in Supplementary Methods.
"""
import pandas as pd

CHEMISTRY_ITEMS = [50971, 52610]
BLOOD_GAS_ITEMS = [50822, 52452]
CHARTED_ITEMS = [220640, 226535, 227442, 227464]
TRUE_ICU_UNITS = ['Cardiac Vascular Intensive Care Unit (CVICU)', 'Coronary Care Unit (CCU)', 'Intensive Care Unit (ICU)', 'Medical Intensive Care Unit (MICU)', 'Medical/Surgical Intensive Care Unit (MICU/SICU)', 'Neuro Surgical Intensive Care Unit (Neuro SICU)', 'Surgical Intensive Care Unit (SICU)', 'Trauma SICU (TSICU)']

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
