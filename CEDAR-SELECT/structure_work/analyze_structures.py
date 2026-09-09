"""Coordinate-only structural qualification; no docking/energy claims.
Run from repository root with pinned cedar-python, single-thread environment.
"""
from pathlib import Path
import json, csv, time, resource
import numpy as np
from Bio.PDB import MMCIFParser, Superimposer
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils import seq1
from Bio.Align import PairwiseAligner
start=time.time(); root=Path('CEDAR-SELECT'); src=root/'inputs/structural_sources'; out=root/'structure_work'
models={s:MMCIFParser(QUIET=True).get_structure(s,src/(s+'.cif'))[0] for s in ['5NEV','5LQF']}
def heavy(res):return [a for a in res if a.element not in ['H','D']]
def protein(chain):return [r for r in chain if 'CA' in r and r.id[0]==' ']
contacts=[]; instances=[]
for sid, model in models.items():
 d=MMCIF2Dict(str(src/(sid+'.cif')))
 for ch in model:
  for lig in ch:
   if lig.resname not in ['72L','4SP']: continue
   la=heavy(lig); local=[]
   for r in protein(ch):
    pa=heavy(r); distances=np.linalg.norm(np.array([a.coord for a in la])[:,None]-np.array([a.coord for a in pa])[None,:],axis=2)
    if distances.min() <= 4.0:
     ix=np.unravel_index(distances.argmin(), distances.shape)
     local+=pa
     contacts.append(dict(pdb=sid,chain=ch.id,ligand=lig.resname,residue=r.resname,resnum=r.id[1],min_distance_A=float(distances.min()),ligand_atom=la[ix[0]].name,protein_atom=pa[ix[1]].name))
   instances.append(dict(pdb=sid,chain=ch.id,ligand=lig.resname,n_heavy=len(la),mean_ligand_B=float(np.mean([a.bfactor for a in la])),mean_contact_residue_B=float(np.mean([a.bfactor for a in local])),resolution_A=float(d['_refine.ls_d_res_high'][0]),rfree=float(d['_refine.ls_R_factor_R_free'][0])))
with (out/'contacts.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(contacts[0]));w.writeheader();w.writerows(contacts)
aligner=PairwiseAligner(mode='global',match_score=2,mismatch_score=-1,open_gap_score=-10,extend_gap_score=-0.5)
comparisons=[]; clashes=[]
for c2 in ['A','C']:
 for c1,fit_region in [(ch,region) for ch in ['A','D'] for region in ['domain','pocket']]:
  p2=protein(models['5NEV'][c2]);p1=protein(models['5LQF'][c1])
  aln=aligner.align(''.join(seq1(r.resname) for r in p2),''.join(seq1(r.resname) for r in p1))[0]
  pairs=[(p2[i],p1[j]) for i,j in zip(*aln.indices) if i>=0 and j>=0]
  anchors=[(a,b) for a,b in pairs if a.resname==b.resname and a.id[1]>25 and a.id[1]<285 and (fit_region=='domain' or 30<=a.id[1]<=90 or 120<=a.id[1]<=145)]
  sup=Superimposer();sup.set_atoms([a['CA'] for a,b in anchors],[b['CA'] for a,b in anchors]);rot,tr=sup.rotran
  loop=[(a,b) for a,b in pairs if 9<=a.id[1]<=19]
  looprms=np.sqrt(np.mean([np.sum((a['CA'].coord-(b['CA'].coord@rot+tr))**2) for a,b in loop]))
  a,b=next((a,b) for a,b in pairs if a.id[1]==15)
  lig=next(r for r in models['5NEV'][c2] if r.resname=='72L'); la=heavy(lig)
  pa=[a for r in p1 for a in heavy(r)];coords=np.array([a.coord@rot+tr for a in pa]); ds=np.linalg.norm(np.array([a.coord for a in la])[:,None]-coords[None,:],axis=2)
  for i,j in zip(*np.where(ds<2.2)):
   clashes.append(dict(cdk2_chain=c2,cdk1_chain=c1,fit_region=fit_region,ligand_atom=la[i].name,cdk1_resnum=pa[j].parent.id[1],cdk1_residue=pa[j].parent.resname,cdk1_atom=pa[j].name,distance_A=float(ds[i,j])))
  if c2=='A' and c1=='A' and fit_region=='pocket':
   import matplotlib
   matplotlib.use('Agg')
   import matplotlib.pyplot as plt
   fig=plt.figure(figsize=(8,6));ax=fig.add_subplot(111,projection='3d')
   for protein_pairs,transform,color,label in [(loop,False,'#1565c0','CDK2 loop'),(loop,True,'#c62828','CDK1 loop')]:
    xyz=np.array([(r2['CA'].coord@rot+tr) if transform else r1['CA'].coord for r1,r2 in protein_pairs]);ax.plot(*xyz.T,color=color,lw=3,label=label)
   for tyr,transform,color in [(a,False,'#1565c0'),(b,True,'#c62828')]:
    names=['CA','CB','CG','CD1','CE1','CZ','OH','CZ','CE2','CD2','CG']
    xyz=np.array([(tyr[n].coord@rot+tr) if transform else tyr[n].coord for n in names]);ax.plot(*xyz.T,color=color,lw=2);ax.text(*xyz[6],' Tyr15 OH',fontsize=8,color=color)
   ccd=MMCIF2Dict(str(src/'72L.cif'));lookup={atom.name:atom.coord for atom in la}
   for n1,n2 in zip(ccd['_chem_comp_bond.atom_id_1'],ccd['_chem_comp_bond.atom_id_2']):
    if n1 in lookup and n2 in lookup:
     xyz=np.array([lookup[n1],lookup[n2]]);ax.plot(*xyz.T,color='#444444',lw=2)
   ax.scatter(*np.array([atom.coord for atom in la]).T,s=7,c='#444444',label='73 observed in CDK2')
   ax.set_title('Observed Tyr15 conformations differ; ligand transfer is a hypothesis',fontsize=11)
   ax.set_xlabel('X (Å)');ax.set_ylabel('Y (Å)');ax.set_zlabel('Z (Å)');ax.legend(fontsize=8,loc='upper left');ax.set_box_aspect((1,1,1));ax.view_init(elev=20,azim=125)
   fig.text(.04,.02,'5NEV A / 5LQF A; 65 conserved pocket Cα atoms fitted. Different ligands and cyclin complexes.\n73 is shown only at its CDK2 coordinates; no CDK1 pose or binding energy is established.',fontsize=8)
   fig.savefig(out/'structure_comparison.png',dpi=160,bbox_inches='tight');plt.close(fig)
  # Native control checks same overlap threshold within deposited kinase chain.
  pn=[a for r in p2 for a in heavy(r)];dn=np.linalg.norm(np.array([a.coord for a in la])[:,None]-np.array([a.coord for a in pn])[None,:],axis=2)
  comparisons.append(dict(cdk2_chain=c2,cdk1_chain=c1,fit_region=fit_region,conserved_anchor_CA_count=len(anchors),anchor_CA_rmsd_A=float(sup.rms),loop_CA_count=len(loop),loop_CA_rmsd_A=float(looprms),tyr15_OH_displacement_A=float(np.linalg.norm(a['OH'].coord-(b['OH'].coord@rot+tr))),transplanted_72L_min_heavy_distance_A=float(ds.min()),transplanted_pairs_lt2_2_A=int((ds<2.2).sum()),native_72L_min_heavy_distance_A=float(dn.min()),native_pairs_lt2_2_A=int((dn<2.2).sum())))
with (out/'transplant_clashes.csv').open('w') as f:
 if clashes:
  w=csv.DictWriter(f,fieldnames=list(clashes[0]));w.writeheader();w.writerows(clashes)
result=dict(instances=instances,comparisons=comparisons,resource=dict(wall_seconds=time.time()-start,maxrss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,cpu_seconds=resource.getrusage(resource.RUSAGE_SELF).ru_utime+resource.getrusage(resource.RUSAGE_SELF).ru_stime),method='Author chain IDs; observed standard residues; sequence aligned; fit conserved CA 26-284 (domain), or 30-90 plus 120-145 (pocket); no energy or relaxation. Heavy-atom distance only; contacts <=4 A; severe overlap <2.2 A. Crystallographic copies are not independent experiments.')
(out/'results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
