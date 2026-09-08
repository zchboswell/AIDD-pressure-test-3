"""Sequential replay into a NEW directory; supplied input bytes verified before/after."""
import argparse,hashlib,json,os,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
EXPECTED_MANIFEST='77934c089aa941e191e08c2459d12d6c522596ba4614f81c7fb7850aba6a7094'

def verify(inputs):
    b=(inputs/'MANIFEST.json').read_bytes();assert hashlib.sha256(b).hexdigest()==EXPECTED_MANIFEST,'Unexpected input manifest'
    files=json.loads(b)['files']
    for name,h in files.items():
        p=(inputs/name).resolve();assert p.is_relative_to(inputs.resolve());assert hashlib.sha256(p.read_bytes()).hexdigest()==h,name
    return {'manifest_sha256':EXPECTED_MANIFEST,'file_count':len(files),'all_hashes_match':True}

def main():
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--runtime',default='/home/zbos/.local/share/cedar-runtime/20260908/cedar-python');a=p.parse_args();inputs=a.inputs.resolve();out=a.out.resolve()
    assert not out.exists(),'Use a new output directory; originals must not be overwritten.'
    before=verify(inputs);out.mkdir(parents=True);(out/'support').mkdir()
    env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',MPLCONFIGDIR=str(out/'.cache/matplotlib'))
    commands=[['test_assay_math.py'],['curate.py','--inputs',str(inputs),'--out',str(out)],['model_audit.py','--inputs',str(inputs),'--curation',str(out/'curation.csv'),'--out',str(out)],['structure_audit.py','--inputs',str(inputs),'--out',str(out)],['decision.py','--inputs',str(inputs),'--out',str(out)],['render_pocket.py','--inputs',str(inputs),'--out',str(out)],['validate_delivery.py','--inputs',str(inputs),'--out',str(out)]]
    runs=[]
    for script,*rest in commands:
        start=time.monotonic();r=subprocess.run([a.runtime,str(ROOT/'analysis'/script),*rest],env=env,capture_output=True,text=True)
        (out/'support'/f'{script}.stdout.txt').write_text(r.stdout);(out/'support'/f'{script}.stderr.txt').write_text(r.stderr)
        runs.append({'script':script,'exit_code':r.returncode,'wall_seconds':time.monotonic()-start})
        if r.returncode:
            (out/'support/replay_status.json').write_text(json.dumps({'status':'failed','runs':runs},indent=2));raise RuntimeError(script+' failed; see retained logs')
    after=verify(inputs);report={'status':'pass','input_verification_before':before,'input_verification_after':after,'runs':runs,'cpu_affinity':sorted(os.sched_getaffinity(0)),'threads_per_process':1,'gpu_used':False,'notes':'Generated numerical/figure artifacts; hand-authored report and handoff are maintained separately. Figure SVG metadata may differ between runs.'}
    (out/'support/replay_status.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
