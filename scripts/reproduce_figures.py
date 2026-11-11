"""Render manuscript figures using frozen aggregate source data only."""
from pathlib import Path
import argparse,csv,json,math,hashlib,sys,platform
from statistics import NormalDist
import importlib.metadata as metadata
import xml.etree.ElementTree as ET
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
from matplotlib import font_manager
import numpy as np
from PIL import Image

GROUPS=['5.5_to_lt6','6_to_lt6.5','ge6.5']
LABELS=['5.5–<6.0 mmol/L','6.0–<6.5 mmol/L','≥6.5 mmol/L']
COLORS=['#0072B2','#D55E00','#7B3294']
MARKERS=['o','s','^']
HOURS=[2,4,6,24]
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':7.5,'axes.titlesize':8,'axes.labelsize':8,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'axes.spines.right':False,'axes.spines.top':False,'axes.linewidth':.7,'legend.frameon':False,'svg.fonttype':'none','pdf.fonttype':42,'ps.fonttype':42,'figure.facecolor':'white','savefig.facecolor':'white','text.color':'#202124','axes.labelcolor':'#202124','axes.edgecolor':'#444444','svg.hashsalt':'hyperkalemia-v1-figures'})

def sha(file):return hashlib.sha256(file.read_bytes()).hexdigest()
def read_sources(src):
 with open(src/'cohort_flow.csv',encoding='utf-8-sig',newline='')as f:flow=list(csv.DictReader(f))
 with open(src/'time_windows.csv',encoding='utf-8-sig',newline='')as f:windows=list(csv.DictReader(f))
 totals=json.loads((src/'cohort_totals.json').read_text(encoding='utf-8'))
 assert totals['final_patients']==5230 and totals['eligible_admissions']==5577
 assert totals['raw_ICU_records']==94458 and totals['raw_ICU_admissions']==85242
 for row in flow:
  for key in ['remaining_admissions','remaining_patients']:row[key]=int(row[key])
  row['excluded_admissions_this_step']=int(row['excluded_admissions_this_step'])if row['excluded_admissions_this_step']else None
 for prev,row in zip(flow,flow[1:]):assert prev['remaining_admissions']-row['remaining_admissions']==row['excluded_admissions_this_step']
 z=NormalDist().inv_cdf(.975)
 for row in windows:
  for key in ['horizon_h','N','repeat_n','repeat_recorded','death_before_repeat','interval_end_before_repeat','still_in_interval_without_repeat','indeterminate']:row[key]=int(row[key])
  for key in ['repeat_pct','CI95_lower_pct','CI95_upper_pct','identification_lower_pct','identification_upper_pct']:row[key]=float(row[key])
  N=row['N'];k=row['repeat_n'];A=N*N+z*z*N;B=-2*k*N-z*z*N;C=k*k;disc=math.sqrt(max(0,B*B-4*A*C))
  assert math.isclose(k/N*100,row['repeat_pct'],abs_tol=1e-10)
  assert abs(100*(-B-disc)/(2*A)-row['CI95_lower_pct'])<1e-9
  assert abs(100*(-B+disc)/(2*A)-row['CI95_upper_pct'])<1e-9
  assert sum(row[key]for key in ['repeat_recorded','death_before_repeat','interval_end_before_repeat','still_in_interval_without_repeat','indeterminate'])==N
 assert len(windows)==16
 for h in HOURS:
  sub=[x for x in windows if x['group']in GROUPS and x['horizon_h']==h]
  whole=next(x for x in windows if x['group']=='overall'and x['horizon_h']==h)
  assert sum(x['N']for x in sub)==whole['N']==5230
  assert sum(x['repeat_n']for x in sub)==whole['repeat_n']
 return flow,windows,totals

def export(fig,out,stem):
 paths=[]
 for ext in ['svg','pdf','png','tiff']:
  dest=out/(stem+'.'+ext)
  kwargs={'dpi':600,'bbox_inches':None,'facecolor':'white'}
  if ext=='tiff':kwargs['pil_kwargs']={'compression':'tiff_lzw'}
  if ext=='pdf':kwargs['metadata']={'Title':stem,'Author':'','Creator':'Python matplotlib','CreationDate':None,'ModDate':None}
  if ext=='svg':kwargs['metadata']={'Date':None,'Title':stem,'Creator':'Python matplotlib'}
  fig.savefig(dest,**kwargs);paths.append(dest)
  if ext=='tiff':
   with Image.open(dest)as raster:rgb=raster.convert('RGB')
   rgb.save(dest,format='TIFF',dpi=(600,600),compression='tiff_lzw')
 preview=out/'previews';preview.mkdir(exist_ok=True)
 fig.savefig(preview/(stem+'.png'),dpi=180,bbox_inches=None,facecolor='white')
 return paths

def box(ax,x,y,w,h,text,*,fill='#F3F6F8',edge='#67737C',fs=7.5,lw=.75):
 patch=FancyBboxPatch((x,y-h/2),w,h,boxstyle='round,pad=0.003,rounding_size=0.003',facecolor=fill,edgecolor=edge,linewidth=lw,zorder=3)
 ax.add_patch(patch)
 artist=ax.text(x+w/2,y,text,ha='center',va='center',fontsize=fs,linespacing=1.27,zorder=4)
 return patch,artist

def arrow(ax,start,end,color='#67737C',style='-|>'):
 patch=FancyArrowPatch(start,end,arrowstyle=style,mutation_scale=7.5,linewidth=.75,color=color,zorder=1,shrinkA=0,shrinkB=0)
 ax.add_patch(patch)

def flow_figure(flow,totals):
 by={row['stage'][:2]:row for row in flow}
 nodes=[
  f"{totals['raw_ICU_records']:,} ICU records in {totals['raw_ICU_admissions']:,} admissions\nFirst ICU record retained per admission",
  f"Adult admissions with a valid first ICU record\n{by['04']['remaining_admissions']:,} admissions",
  f"Recorded chemistry potassium ≥5.5 mmol/L\n{by['05']['remaining_admissions']:,} admissions",
  f"Identifiable specimen and valid index availability\n{by['07']['remaining_admissions']:,} admissions",
  f"Index specimen comment criteria met\n{by['08']['remaining_admissions']:,} admissions",
  f"Sample and result within the same\nverified continuous ICU interval\n{by['09']['remaining_admissions']:,} admissions",
  f"No earlier numeric laboratory potassium ≥5.5\n{by['10']['remaining_admissions']:,} admissions",
  f"Eligible after prior high potassium checks\n{totals['eligible_admissions']:,} admissions in {totals['final_patients']:,} patients",
  f"Earliest eligible event retained per patient\n{totals['final_patients']:,} unique patients"]
 exclusions=[
  f"Timing or survival-at-entry criteria not met\nn = {by['04']['excluded_admissions_this_step']:,} admissions",
  f"No qualifying chemistry high result\nin the first ICU record window\nn = {by['05']['excluded_admissions_this_step']:,} admissions",
  f"Index availability time criteria not met\nn = {by['07']['excluded_admissions_this_step']:,} admissions",
  f"Index specimen comment rule not met\nn = {by['08']['excluded_admissions_this_step']:,} admissions",
  f"Continuous ICU interval criteria not met\nn = {by['09']['excluded_admissions_this_step']:,} admissions",
  f"Earlier numeric laboratory high result\nn = {by['10']['excluded_admissions_this_step']:,} admissions",
  f"Earlier charted high result\nn = {by['11']['excluded_admissions_this_step']:,} admissions",
  f"Additional eligible admissions\nfrom repeated patients\nn = {by['13']['excluded_admissions_this_step']:,} admissions"]
 selected=['01','04','05','07','08','09','10','12','13']
 for a,b,excluded in zip(selected,selected[1:],[126,75794,80,907,92,2556,110,347]):assert by[a]['remaining_admissions']-by[b]['remaining_admissions']==excluded
 fig=plt.figure(figsize=(183/25.4,235/25.4));ax=fig.add_axes([0,0,1,1]);ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis('off')
 ys=np.linspace(.95,.19,len(nodes));x=.02;w=.55;h=.061;cx=x+w/2;checks=[]
 for i,text in enumerate(nodes):
  final=i==len(nodes)-1
  checks.append(box(ax,x,ys[i],w,h,text,fill='#E3EDF5'if final else '#F3F6F8',edge='#31556D'if final else '#67737C',lw=1.1 if final else .75,fs=7.5))
  if i:
   upper=ys[i-1]-h/2-.003;lower=ys[i]+h/2+.003
   arrow(ax,(cx,upper),(cx,lower))
   mid=(upper+lower)/2;arrow(ax,(cx,mid),(.64,mid))
   checks.append(box(ax,.64,mid,.34,.060,exclusions[i-1],fill='white',edge='#9AA1A6',fs=7.0,lw=.6))
 # The final branching counts refer to unique patients and sum to 5,230.
 ax.plot([cx,cx],[ys[-1]-h/2-.003,.126],color='#67737C',lw=.75,zorder=1)
 ax.plot([.17,.83],[.126,.126],color='#67737C',lw=.75,zorder=1)
 for label,group,color,c in zip(LABELS,GROUPS,COLORS,[.17,.50,.83]):
  arrow(ax,(c,.126),(c,.107))
  checks.append(box(ax,c-.15,.071,.30,.068,f"{label}\nN = {totals['strata_N'][group]:,} patients",fill='#F8FAFC',edge=color,fs=7.5,lw=.9))
 ax.text(.5,.012,'Final cohort: one event per patient. Earlier eligibility steps count admissions.',ha='center',va='bottom',fontsize=6.8)
 fig.canvas.draw();renderer=fig.canvas.get_renderer()
 for patch,text in checks:
  pb=patch.get_window_extent(renderer);tb=text.get_window_extent(renderer)
  assert pb.contains(tb.x0,tb.y0)and pb.contains(tb.x1,tb.y1),'Flow text exceeds its box'
 return fig,{'main_nodes':nodes,'exclusion_boxes':exclusions,'unique_patient_strata':[totals['strata_N'][g]for g in GROUPS],'all_box_text_fits':True}

def repeat_figure(windows,totals):
 fig,axes=plt.subplots(1,3,figsize=(183/25.4,118/25.4),sharey=True)
 fig.subplots_adjust(left=.105,right=.985,bottom=.23,top=.79,wspace=.28)
 plotted=[]
 for j,(ax,g,label,color,marker)in enumerate(zip(axes,GROUPS,LABELS,COLORS,MARKERS)):
  rows=sorted([x for x in windows if x['group']==g],key=lambda x:x['horizon_h'])
  assert [x['horizon_h']for x in rows]==HOURS and all(x['N']==totals['strata_N'][g]for x in rows)
  y=np.array([x['repeat_pct']for x in rows]);lo=np.array([x['CI95_lower_pct']for x in rows]);hi=np.array([x['CI95_upper_pct']for x in rows]);xx=np.arange(4)
  assert np.all(np.diff(y)>=0)
  ax.set_axisbelow(True);ax.grid(axis='y',color='#DADDE0',linewidth=.45)
  ax.axvspan(.78,1.22,color='#EEEEEE',zorder=0)
  ax.errorbar(xx,y,yerr=np.vstack([y-lo,hi-y]),fmt=marker,color=color,ecolor=color,markersize=4.2,markeredgewidth=.8,elinewidth=.9,capsize=2.7,capthick=.8,zorder=3)
  for xpos,point,lower in zip(xx,y,lo):ax.text(xpos,lower-3.0,f'{point:.1f}',ha='center',va='top',fontsize=7.5,color='#333333')
  ax.set_xlim(-.42,3.42);ax.set_ylim(0,100);ax.set_xticks(xx,[str(h)for h in HOURS]);ax.set_yticks([0,25,50,75,100]);ax.tick_params(length=3,width=.6)
  ax.set_title(f'K {label}\nN = {totals["strata_N"][g]:,} patients',fontsize=8,pad=12,linespacing=1.45)
  ax.text(-.04,1.24,chr(97+j),transform=ax.transAxes,fontweight='bold',fontsize=9,ha='left',va='top')
  plotted.extend(rows)
 axes[0].set_ylabel('Patients with a documented repeat (%)',labelpad=7)
 fig.supxlabel('Prespecified time window after index result (h)',x=.55,y=.115,fontsize=8)
 fig.text(.55,.035,'Points: observed proportions; bars: Wilson 95% CI. Grey band: primary 4-h window.',ha='center',fontsize=6.8)
 return fig,plotted

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--source-dir',type=Path,default=Path(__file__).resolve().parent/'source_data');parser.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent)
 args=parser.parse_args();src=args.source_dir;out=args.output_dir;out.mkdir(parents=True,exist_ok=True)
 before={f.name:sha(f)for f in src.iterdir()if f.is_file()}
 flow,windows,totals=read_sources(src)
 f1,diagram=flow_figure(flow,totals);exports=export(f1,out,'Figure_1_cohort_flow');plt.close(f1)
 f2,plotted=repeat_figure(windows,totals);exports+=export(f2,out,'Figure_2_retesting_by_potassium');plt.close(f2)
 with open(out/'Figure_2_source_data.csv','w',encoding='utf-8-sig',newline='')as f:
  writer=csv.DictWriter(f,fieldnames=list(plotted[0]));writer.writeheader();writer.writerows(plotted)
 (out/'Figure_1_source_data.json').write_text(json.dumps(diagram,ensure_ascii=False,indent=2),encoding='utf-8')
 versions={name:metadata.version(name)for name in ['matplotlib','numpy','pillow','contourpy','cycler','fonttools','kiwisolver','packaging','pyparsing','python-dateutil','tzdata']}
 (out/'requirements-figures.txt').write_text('\n'.join(name+'=='+version for name,version in versions.items())+'\n',encoding='utf-8')
 editable={}
 for file in out.glob('Figure_*.svg'):
  tree=ET.parse(file);texts=tree.findall('.//{http://www.w3.org/2000/svg}text');assert texts,'Missing editable SVG text';editable[file.name]=len(texts)
 assert all(sha(src/name)==value for name,value in before.items())
 report={'analysis_version':'1.0','Python':platform.python_version(),'backend':'matplotlib only','versions':versions,'font':font_manager.FontProperties(family=plt.rcParams['font.sans-serif']).get_name(),'source_sha256':before,'plotted_points':len(plotted),'groups_N':[totals['strata_N'][g]for g in GROUPS],'CI_method':'Wilson patient-level 95% CI; no fitted model','all_input_proportions_and_CIs_independently_checked':True,'flow_counts_conserve':True,'flow_box_text_fits':diagram['all_box_text_fits'],'fixed_group_denominators':True,'categorical_time_windows_without_interpolation':True,'editable_svg_text_counts':editable,'exports_sha256':{file.name:sha(file)for file in exports}}
 (out/'figure_qc.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
 print(json.dumps({'status':'figures_created','points':len(plotted),'font':report['font'],'versions':versions}),flush=True)

if __name__=='__main__':main()
