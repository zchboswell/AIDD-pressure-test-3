import json,datetime
from rdkit import Chem,rdBase
from rdkit.Chem import Descriptors,Crippen
from pathlib import Path
x=json.load(open('CEDAR-SELECT/inputs/public_assay_context.json'))
a=json.load(open('CEDAR-SELECT/inputs/coverage/analogue_evidence.json'))
known={Chem.MolToSmiles(Chem.MolFromSmiles(r['smiles'])) for k in ['train','development','bounds_context_only'] for r in x[k]}
known|={Chem.MolToSmiles(Chem.MolFromSmiles(r['smiles'])) for r in a['compounds']}
specs=[('D1','ethoxy/cycloheptyl','CCOc1n[nH]cc1Nc1ncc2c(n1)N([C@@H]1CCCC[C@@H](O)C1)C(=O)C21CC1','ring expansion'),('D2','2-fluoroethoxy/cycloheptyl','O=C1N([C@@H]2CCCC[C@@H](O)C2)c2nc(Nc3c[nH]nc3OCCF)ncc2C12CC2','ring expansion + fluoroethoxy'),('D3','1-fluoropropan-2-yloxy stereoisomer A','O=C1N([C@@H]2CCC[C@@H](O)C2)c2nc(Nc3c[nH]nc3O[C@H](C)CF)ncc2C12CC2','alpha-methyl fluoroethoxy A'),('D4','1-fluoropropan-2-yloxy stereoisomer B','O=C1N([C@@H]2CCC[C@@H](O)C2)c2nc(Nc3c[nH]nc3O[C@@H](C)CF)ncc2C12CC2','alpha-methyl fluoroethoxy B')]
out=[]
for id,name,s,edit in specs:
 m=Chem.MolFromSmiles(s); assert m
 cs=Chem.MolToSmiles(m,isomericSmiles=True); assert cs not in known
 out.append(dict(design_id=id,name=name,canonical_isomeric_smiles=cs,edit=edit,new_to_supplied_labeled_packet=True,worldwide_novelty_claim=False,measured_potency=None,predicted_potency=None,MW=Descriptors.MolWt(m),cLogP=Crippen.MolLogP(m),TPSA=Descriptors.TPSA(m),stereocenters=Chem.FindMolChiralCenters(m,includeUnassigned=True),priority='conditional first pair' if id in ['D1','D2'] else 'reserve paired diagnostic chemistry'))
comps=[]
for r in a['compounds']:
 if r['public_source_identifiers']['example'] in [23,39]:
  comps.append(dict(example=r['public_source_identifiers']['example'],compound_id=r['compound_id'],canonical_isomeric_smiles=Chem.MolToSmiles(Chem.MolFromSmiles(r['smiles'])),role='existing diagnostic comparator; not new design'))
payload=dict(schema='cedar-designs-v1',rdkit_version=rdBase.rdkitVersion,checked_scope='Canonical isomeric graph equality against public_assay_context train/development/bounds and coverage analogue compounds only; no global novelty claim',designs=out,diagnostic_comparators=comps)
Path('CEDAR-SELECT/delivery/designs.json').write_text(json.dumps(payload,indent=2)+'\n')
print(json.dumps(payload,indent=2))
