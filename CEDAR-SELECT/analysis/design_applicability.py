import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import json,csv,argparse
from pathlib import Path
from rdkit import Chem,DataStructs
from rdkit.Chem import rdFingerprintGenerator,Draw
import pandas as pd
p=argparse.ArgumentParser();p.add_argument('--out',default='CEDAR-SELECT/delivery');a=p.parse_args();out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
d=json.load(open('CEDAR-SELECT/delivery/designs.json'));w=pd.read_csv('CEDAR-SELECT/delivery/model/paired.csv');gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=True);fps=[gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in w.canonical]
old=[]
for r in csv.DictReader(open('inputs/data/inventory.csv')):
 m=Chem.MolFromSmiles(r['smiles'])
 if m:old.append(Chem.MolToSmiles(m))
rows=[];mols=[];labels=[]
for x in d['designs']:
 m=Chem.MolFromSmiles(x['canonical_isomeric_smiles']);s=Chem.MolToSmiles(m);assert s not in old and s not in set(w.canonical)
 sim=DataStructs.BulkTanimotoSimilarity(gen.GetFingerprint(m),fps);i=max(range(len(sim)),key=sim.__getitem__)
 rows.append(dict(design_id=x['design_id'],status='abstained',endpoint='prospective paired CDK2/CDK1 and cellular engagement',point_prediction='',lower='',upper='',max_public_similarity=sim[i],nearest_public_compound=w.iloc[i].compound_id,new_to_original_inventory=True,reason='Same-study selectivity ranking failed; no prospective cellular training labels. Similarity is descriptive, not validated applicability.'))
 mols.append(m);labels.append(x['design_id']+': '+('conditional nomination' if x['design_id'] in ['D1','D2'] else 'reserve; not ordered'))
for c in d['diagnostic_comparators']:mols.append(Chem.MolFromSmiles(c['canonical_isomeric_smiles']));labels.append('Known example '+str(c['example'])+' diagnostic')
pd.DataFrame(rows).to_csv(out/'predictions.csv',index=False)
Draw.MolsToGridImage(mols,molsPerRow=2,subImgSize=(550,330),legends=labels).save(str(out/'designs.png'))
print(pd.DataFrame(rows).to_string(index=False))
