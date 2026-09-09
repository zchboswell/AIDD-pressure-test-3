import json
from rdkit import Chem
x=json.load(open('CEDAR-SELECT/inputs/public_assay_context.json'))
d={}
for kind in ['train','development','bounds_context_only']:
 for r in x[kind]:
  k=Chem.MolToSmiles(Chem.MolFromSmiles(r['smiles']))
  d.setdefault(k,{})[r['target']]=r.get('pIC50',r.get('value',r.get('bound_pIC50')))
print(x['train'][0])
for s,v in d.items(): print(s,v)
