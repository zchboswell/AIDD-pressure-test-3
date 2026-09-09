"""Assemble isolated replay from explicit inputs/code; never overwrite delivery."""
import argparse,subprocess,shutil,json,os,time,hashlib,resource
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--out',required=True);a=p.parse_args();repo=Path.cwd();dest=Path(a.out).resolve();assert not dest.exists(),'Use a new empty destination';dest.mkdir(parents=True)
t0=time.time();base=dest/'CEDAR-SELECT';base.mkdir();shutil.copytree(repo/'CEDAR-SELECT/inputs',base/'inputs')
for folder in ['analysis','chemistry_work','coverage_work','structure_work']:
 (base/folder).mkdir()
 for source in (repo/'CEDAR-SELECT'/folder).glob('*.py'):shutil.copy2(source,base/folder/source.name)
(base/'delivery').mkdir();(dest/'inputs/data').mkdir(parents=True);shutil.copy2(repo/'inputs/data/inventory.csv',dest/'inputs/data/inventory.csv')
env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',MPLCONFIGDIR=str(dest/'.mplconfig'))
python='/home/zbos/.local/share/cedar-runtime/20260908/cedar-python'
scripts=['analysis/predictive.py','analysis/censor_sensitivity.py','chemistry_work/design_check.py','analysis/design_applicability.py','coverage_work/reconcile.py','structure_work/analyze_structures.py','analysis/plan.py']
runs=[]
for s in scripts:
 start=time.time();r=subprocess.run([python,'CEDAR-SELECT/'+s],cwd=dest,env=env,text=True,capture_output=True,timeout=180);runs.append({'script':s,'exit_code':r.returncode,'wall_seconds':time.time()-start,'stderr':r.stderr});assert r.returncode==0,(s,r.stderr)
compared=[]
for sub in ['delivery/model','delivery','structure_work']:
 for f in (repo/'CEDAR-SELECT'/sub).glob('*.csv'):
  q=base/sub/f.name;assert q.exists();assert f.read_bytes()==q.read_bytes(),str(f);compared.append(sub+'/'+f.name)
for sub in ['delivery/designs.json','delivery/coverage_reconciliation.json','delivery/plan_validation.json']:
 assert json.loads((repo/'CEDAR-SELECT'/sub).read_text())==json.loads((base/sub).read_text());compared.append(sub)
src=json.load(open(repo/'CEDAR-SELECT/structure_work/results.json'));new=json.load(open(base/'structure_work/results.json'))
# Runtime counters differ; all scientific keys must match.
for k in src:
 if k not in ['resources','runtime','execution','elapsed_seconds','resource_usage']:
  if src[k]!=new[k]:
   assert any(z in k.lower() for z in ['time','cpu','rss','resource','wall','elapsed']),k
  else:compared.append('structure_work/results.json:'+k)
result={'status':'PASS','destination':str(dest),'runs':runs,'byte_identical_csv_and_equal_json':compared,'wall_seconds':time.time()-t0,'children_cpu_seconds':resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime+resource.getrusage(resource.RUSAGE_CHILDREN).ru_stime,'children_peak_rss_kib':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,'note':'Runtime counters excluded from equality; no original outputs overwritten.'}
(dest/'replay_check.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
