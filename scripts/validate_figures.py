from pathlib import Path
import csv,json,hashlib,math,xml.etree.ElementTree as ET
from PIL import Image
from pypdf import PdfReader
import argparse

parser=argparse.ArgumentParser();parser.add_argument('--figure-dir',type=Path,default=Path(__file__).resolve().parents[1]/'reproduced_figures')
parser.add_argument('--source-dir',type=Path,default=Path(__file__).resolve().parents[1]/'data/figure_source')
args=parser.parse_args();p=args.figure_dir;src=args.source_dir
qc=json.loads((p/'figure_qc.json').read_text(encoding='utf-8'))
assert qc['plotted_points']==12 and qc['groups_N']==[2895,1079,1256]
assert all(qc[k]for k in ['all_input_proportions_and_CIs_independently_checked','flow_counts_conserve','flow_box_text_fits','fixed_group_denominators','categorical_time_windows_without_interpolation'])
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
for name,value in qc['source_sha256'].items():assert sha(src/name)==value
for name,value in qc['exports_sha256'].items():assert sha(p/name)==value
with open(p/'Figure_2_source_data.csv',encoding='utf-8-sig',newline='')as f:plotted=list(csv.DictReader(f))
with open(src/'time_windows.csv',encoding='utf-8-sig',newline='')as f:original=list(csv.DictReader(f))
assert len(plotted)==12
for a in plotted:
 b=next(x for x in original if x['group']==a['group']and x['horizon_h']==a['horizon_h'])
 for key in a:
  if key=='group':assert a[key]==b[key]
  else:assert math.isclose(float(a[key]),float(b[key]),abs_tol=1e-10)
report={}
for stem,width,height in [('Figure_1_cohort_flow',183,235),('Figure_2_retesting_by_potassium',183,118)]:
 for ext in ['png','tiff']:
  with Image.open(p/(stem+'.'+ext))as im:
   dpi=im.info.get('dpi');assert dpi and abs(dpi[0]-600)<.1 and abs(dpi[1]-600)<.1
   assert abs(im.width-width/25.4*600)<2 and abs(im.height-height/25.4*600)<2
   if ext=='tiff':assert im.mode=='RGB'
   report[stem+'.'+ext]={'pixels':im.size,'dpi':[float(v)for v in dpi],'mode':im.mode}
 reader=PdfReader(p/(stem+'.pdf'));assert len(reader.pages)==1
 page=reader.pages[0];assert abs(float(page.mediabox.width)*25.4/72-width)<.01 and abs(float(page.mediabox.height)*25.4/72-height)<.01
 text=page.extract_text();assert len(text)>100
 fonts=page['/Resources']['/Font'].get_object();embedded=[]
 for key,value in fonts.items():
  font=value.get_object()
  if '/DescendantFonts'in font:
   descriptor=font['/DescendantFonts'][0].get_object().get('/FontDescriptor')
  else:descriptor=font.get('/FontDescriptor')
  if descriptor:embedded.append('/FontFile2'in descriptor.get_object())
 assert embedded and all(embedded)
 xml=ET.parse(p/(stem+'.svg'));nodes=xml.findall('.//{http://www.w3.org/2000/svg}text');assert nodes
 texts=[''.join(x.itertext())for x in nodes]
 if 'Figure_2'in stem:
  assert all(any(f'N = {N:,} patients'in value for value in texts)for N in [2895,1079,1256])
  for row in plotted:assert f"{float(row['repeat_pct']):.1f}"in texts
 report[stem]={'pdf_single_page':True,'embedded_TrueType_fonts':True,'pdf_text_extractable':True,'svg_text_nodes':len(nodes)}
(p/'export_validation.json').write_text(json.dumps({'data_checks':'passed','exports':report},indent=2),encoding='utf-8')
print('PASS: 12 plotted points exactly match frozen data; fixed N and CI checks; SVG editable text; PDF selectable embedded fonts; exact dimensions; 600dpi PNG/RGB TIFF; source and export hashes')
