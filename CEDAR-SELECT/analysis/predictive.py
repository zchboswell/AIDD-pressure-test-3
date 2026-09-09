"""Deterministic within-study evaluation; no hidden labels or pooled assays."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: os.environ[k]='1'
import argparse,json,time,resource,hashlib
from pathlib import Path
import numpy as np,pandas as pd
from rdkit import Chem,DataStructs,rdBase
from rdkit.Chem import rdFingerprintGenerator,Descriptors,Crippen
from sklearn.model_selection import GroupKFold
from sklearn.metrics import mean_absolute_error,r2_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('--out',default='CEDAR-SELECT/delivery/model');a=p.parse_args();out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
t0=time.time(); src=json.load(open('CEDAR-SELECT/inputs/public_assay_context.json'))
rows=[]
for part in ['train','development','bounds_context_only']:
 for r in src[part]:
  r=dict(r);m=Chem.MolFromSmiles(r['smiles']);assert m is not None
  r['canonical']=Chem.MolToSmiles(m);r['connectivity']=Chem.MolToSmiles(m,isomericSmiles=False)
  r['curation_status']='context_bound' if part=='bounds_context_only' else 'exact_training' if part=='train' else 'exposed_development'
  rows.append(r)
df=pd.DataFrame(rows);assert df.row_id.is_unique
# Union source leakage relationships and connectivity, including isotopologues.
parent={s:s for s in df.canonical}
def root(s):
 while parent[s]!=s:s=parent[s]
 return s
for field in ['connectivity','leakage_group_id']:
 for _,g in df.groupby(field):
  members=list(g.canonical.unique()); base=root(members[0])
  for member in members[1:]:parent[root(member)]=base
df['validation_group']=df.canonical.map(lambda s:root(s))
assert not set(df[df.partition=='train'].validation_group)&set(df[df.partition=='development'].validation_group)
# Canonical stereo graphs retain source identities. Group by connectivity to keep stereoisomers together.
assert not set(df[df.partition=='train'].connectivity)&set(df[df.partition=='development'].connectivity)
df.to_csv(out/'curation.csv',index=False)
wide=df.pivot(index='canonical',columns='target',values='pIC50').reset_index()
meta=df.drop_duplicates('canonical').set_index('canonical')
for c in ['compound_id','connectivity','validation_group','partition']:wide[c]=wide.canonical.map(meta[c])
wide['delta_p']=wide.CDK2-wide.CDK1
b=df[df.curation_status=='context_bound'].set_index('canonical').bound_pIC50
wide['delta_lower_bound']=wide.CDK2-wide.canonical.map(b)
wide['ratio']=10**wide.delta_p;wide['ratio_lower_bound']=10**wide.delta_lower_bound
wide.to_csv(out/'paired.csv',index=False)
gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=True)
fps=[gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in wide.canonical]
K=np.array([DataStructs.BulkTanimotoSimilarity(f,fps) for f in fps])
train=np.where(wide.partition=='train')[0];dev=np.where(wide.partition=='development')[0]
y=wide[['CDK2','CDK1','delta_p']].values
# Hyperparameters fixed before fitting: Tanimoto kernel ridge alpha1; nearest3 similarity weights; median.
def pred(tr,te,col,method):
 tr=np.array([i for i in tr if np.isfinite(y[i,col])]); yy=y[tr,col]
 if method=='median':return np.repeat(np.median(yy),len(te))
 if method=='kernel_ridge':return K[np.ix_(te,tr)]@np.linalg.solve(K[np.ix_(tr,tr)]+np.eye(len(tr)),yy-yy.mean())+yy.mean()
 ans=[]
 for i in te:
  j=np.argsort(-K[i,tr],kind='stable')[:3];w=K[i,tr[j]];ans.append(np.average(yy[j],weights=w) if w.sum() else yy.mean())
 return np.array(ans)
records=[]
splits=[('development',-1,train,dev)]
for fold,(tr,te) in enumerate(GroupKFold(n_splits=5).split(train,groups=wide.iloc[train].validation_group)):
 splits.append(('train_group_cv',fold,train[tr],train[te]))
for scheme,fold,tr,te in splits:
 assert not set(wide.iloc[tr].validation_group)&set(wide.iloc[te].validation_group)
 for col,target in enumerate(['CDK2','CDK1','delta_p']):
  for method in ['median','nearest3','kernel_ridge']:
   estimates=pred(tr,te,col,method)
   validtr=[i for i in tr if np.isfinite(y[i,col])]
   for i,estimate in zip(te,estimates):
    records.append(dict(scheme=scheme,fold=fold,method=method,target=target,compound_id=wide.iloc[i].compound_id,canonical=wide.iloc[i].canonical,actual=y[i,col],predicted=estimate,max_similarity=K[i,validtr].max(),bound=wide.iloc[i].delta_lower_bound if col==2 else (wide.iloc[i].CDK2-wide.iloc[i].delta_lower_bound if col==1 else np.nan)))
r=pd.DataFrame(records);r.to_csv(out/'held_out_predictions.csv',index=False)
metrics=[]
for keys,g in r.groupby(['scheme','method','target']):
 v=g.dropna(subset=['actual']);err=np.abs(v.actual-v.predicted).values
 rng=np.random.default_rng(709); boot=np.mean(rng.choice(err,(2000,len(err)),replace=True),axis=1)
 metrics.append(dict(zip(['scheme','method','target'],keys))|dict(n=len(v),mae=float(err.mean()),mae_bootstrap_low=float(np.quantile(boot,.025)),mae_bootstrap_high=float(np.quantile(boot,.975)),r2=r2_score(v.actual,v.predicted),spearman=float(v[['actual','predicted']].corr(method='spearman').iloc[0,1]),median_similarity=float(v.max_similarity.median()),bound_n=int(g.bound.notna().sum()),bound_violations=int(((g.predicted<g.bound) if keys[2]=='delta_p' else (g.predicted>g.bound)).sum()) if keys[2]!='CDK2' else 0))
pd.DataFrame(metrics).to_csv(out/'metrics.csv',index=False)
# Source crosswalk is identity + target + assay interpretation, never count analogue exports as replicates.
a8=json.load(open('CEDAR-SELECT/inputs/coverage/analogue_evidence.json'));cross=[];pairs=[]
for c in a8['compounds']:
 s=Chem.MolToSmiles(Chem.MolFromSmiles(c['smiles']));z=df[df.canonical==s];assert len(z)==2
 obs=[o for o in a8['observations'] if o['compound_id']==c['compound_id']];v={o['target']:float(o['original_value'])*1000 for o in obs}
 for o in obs:
  rr=z[z.target==o['target']].iloc[0];assert abs(rr.pIC50-o['pIC50'])<1e-10
  assert abs(9-np.log10(float(o['original_value'])*1000)-o['pIC50'])<1e-10
  cross.append(dict(compound_id=c['compound_id'],HC=rr.compound_id,row_id=rr.row_id,observation_id=o['observation_id'],canonical=s,status='same_public_measurement_not_replicate',example=c['public_source_identifiers']['example']))
 pairs.append(dict(example=c['public_source_identifiers']['example'],compound_id=c['compound_id'],canonical=s,CDK2_nM=v['CDK2'],CDK1_nM=v['CDK1'],ratio=v['CDK1']/v['CDK2'],delta_p=np.log10(v['CDK1']/v['CDK2'])))
pd.DataFrame(cross).to_csv(out/'crosswalk.csv',index=False);pair=pd.DataFrame(pairs);pair.to_csv(out/'analogue_pairs.csv',index=False)
fig,axs=plt.subplots(1,2,figsize=(10,4));axs[0].scatter(pair.CDK2_nM,pair.ratio)
for _,z in pair.iterrows():axs[0].annotate(str(z['example']),(z.CDK2_nM,z.ratio))
axs[0].set(xscale='log',xlabel='CDK2 IC50 (nM); lower preferred',ylabel='CDK1 / CDK2 assay ratio',title='Historical analogues (example IDs)')
g=r[(r.scheme=='development')&(r.method=='kernel_ridge')&(r.target=='delta_p')].dropna(subset=['actual']);axs[1].scatter(g.actual,g.predicted);axs[1].plot([0,3],[0,3],c='gray');axs[1].set(xlabel='Observed delta pIC50',ylabel='Predicted delta pIC50',title='Exposed within-study development; n='+str(len(g)))
fig.tight_layout();fig.savefig(out/'diagnostic.png',dpi=180);plt.close(fig)
summary={'compounds':len(wide),'exact_rows':int(df.pIC50.notna().sum()),'bound_rows':len(b),'paired_exact':int(wide.delta_p.notna().sum()),'partition_counts':wide.partition.value_counts().to_dict(),'within_study_only':True,'hyperparameters':'fixed Morgan radius2 2048 chiral; kernel ridge alpha1; nearest3 similarity weighted; no tuning','bounds':'CDK1 pIC50 <=5 gives delta >= CDK2-5; excluded from exact regression, inequality violation diagnostics retained','bootstrap':'2000 resamples of per-compound absolute errors, fixed predictions; descriptive not external calibration','runtime_seconds':time.time()-t0,'cpu_seconds':resource.getrusage(resource.RUSAGE_SELF).ru_utime+resource.getrusage(resource.RUSAGE_SELF).ru_stime,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'rdkit':rdBase.rdkitVersion}
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary));print(pd.DataFrame(metrics).to_string(index=False))
