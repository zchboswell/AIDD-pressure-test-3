"""One-sided squared-hinge loss sensitivity; bounds remain inequalities."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import argparse,json,time,resource
from pathlib import Path
import numpy as np,pandas as pd
from scipy.optimize import minimize
from rdkit import Chem,DataStructs
from rdkit.Chem import rdFingerprintGenerator
from sklearn.model_selection import GroupKFold
p=argparse.ArgumentParser();p.add_argument('--base',default='CEDAR-SELECT/delivery/model');p.add_argument('--out',default='CEDAR-SELECT/delivery/model');a=p.parse_args();out=Path(a.out);out.mkdir(parents=True,exist_ok=True);t0=time.time()
w=pd.read_csv(Path(a.base)/'paired.csv');gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=True);fps=[gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in w.canonical];K=np.array([DataStructs.BulkTanimotoSimilarity(f,fps) for f in fps]);tr=np.where(w.partition=='train')[0];te=np.where(w.partition=='development')[0]
splits=[('development',-1,tr,te)]+[('train_group_cv',f,tr[x],tr[z]) for f,(x,z) in enumerate(GroupKFold(5).split(tr,groups=w.iloc[tr].validation_group))]
rows=[];fits=[]
for scheme,fold,tr,te in splits:
 A=K[np.ix_(tr,tr)];B=K[np.ix_(te,tr)]; y=w.CDK1.values[tr];exact=np.isfinite(y);upper=(w.CDK2-w.delta_lower_bound).values[tr];assert np.all(exact|np.isfinite(upper))
 for strength in [.25,1,4]:
  def objective(theta):
   alpha=theta[:-1];f=A@alpha+theta[-1];res=np.where(exact,f-np.nan_to_num(y),np.maximum(f-np.nan_to_num(upper),0));weights=np.where(exact,1,strength);v=res*weights
   return .5*np.sum(res*v)+.5*alpha@A@alpha, np.r_[A@v+A@alpha,np.sum(v)]
  fit=minimize(objective,np.r_[np.zeros(len(tr)),5.5],jac=True,method='L-BFGS-B',options={'maxiter':2000,'ftol':1e-12,'gtol':1e-7});assert fit.success,fit.message
  fits.append({'scheme':scheme,'fold':fold,'strength':strength,'iterations':fit.nit,'success':bool(fit.success)})
  estimate=B@fit.x[:-1]+fit.x[-1]; yy=w.CDK2.values[tr];p2=B@np.linalg.solve(A+np.eye(len(tr)),yy-yy.mean())+yy.mean()
  for i,p1,p2i in zip(te,estimate,p2):rows.append(dict(scheme=scheme,fold=fold,bound_weight=strength,compound_id=w.iloc[i].compound_id,pred_CDK1=p1,pred_delta=p2i-p1,actual_CDK1=w.iloc[i].CDK1,actual_delta=w.iloc[i].delta_p,CDK1_upper=w.iloc[i].CDK2-w.iloc[i].delta_lower_bound))
r=pd.DataFrame(rows);r.to_csv(out/'censor_predictions.csv',index=False);metrics=[]
for (scheme,weight),g in r.groupby(['scheme','bound_weight']):
 v=g.dropna(subset=['actual_CDK1']);metrics.append(dict(scheme=scheme,bound_weight=weight,exact_n=len(v),CDK1_mae=float(abs(v.pred_CDK1-v.actual_CDK1).mean()),delta_mae=float(abs(v.pred_delta-v.actual_delta).mean()),censored_n=int(g.CDK1_upper.notna().sum()),bound_violations=int((g.pred_CDK1>g.CDK1_upper).sum())))
pd.DataFrame(metrics).to_csv(out/'censor_metrics.csv',index=False)
(out/'censor_execution.json').write_text(json.dumps({'fits':fits,'seconds':time.time()-t0,'cpu_seconds':resource.getrusage(resource.RUSAGE_SELF).ru_utime+resource.getrusage(resource.RUSAGE_SELF).ru_stime,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'interpretation':'Sensitivity analysis, not validated censored likelihood or confidence interval. Bound weight varies relative penalty; no artificial exact threshold labels.'},indent=2)+'\n');print(pd.DataFrame(metrics).to_string(index=False))
