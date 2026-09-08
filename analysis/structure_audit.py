"""Audit supplied exports without fitting; render independently of commercial software.
Usage: cedar-python analysis/structure_audit.py --out delivery
"""
from pathlib import Path
import argparse, json, csv
import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import rdMolAlign
from rdkit.Chem.MolStandardize import rdMolStandardize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',default='delivery');p.add_argument('--inputs',type=Path,default=ROOT/'inputs');args=p.parse_args()
    out=Path(args.out); support=out/'structure_support';fig=out/'figures'
    support.mkdir(parents=True,exist_ok=True);fig.mkdir(parents=True,exist_ok=True)
    src=args.inputs/'structures'
    docs='RDKit '+rdBase.rdkitVersion+'\n\n'+str(rdMolAlign.CalcRMS.__doc__)+'\n\n'+str(rdMolAlign.GetBestRMS.__doc__)+'\n\n'+str(Chem.RenumberAtoms.__doc__)
    (support/'installed_interface.txt').write_text(docs)
    ref=Chem.SDMolSupplier(str(src/'1h1q_2a6_A.sdf'),removeHs=False)[0]
    xyz=ref.GetConformer().GetPositions()
    atoms=[]
    for line in (src/'1h1q.pdb').read_text().splitlines():
        if line[:6].strip() in ('ATOM','HETATM'):
            atoms.append(dict(record=line[:6].strip(),name=line[12:16].strip(),resname=line[17:20].strip(),chain=line[21],resid=line[22:26].strip(),element=line[76:78].strip(),xyz=[float(line[30:38]),float(line[38:46]),float(line[46:54])]))
    deposited=[a for a in atoms if a['resname']=='2A6' and a['chain']=='A']
    assert len(deposited)==len(xyz)==24
    assert np.max(np.abs(np.array([a['xyz'] for a in deposited])-xyz))<1e-6
    checks={}; mols={};maps={}; rows=[]
    for label,file in [('reference','vendor_reference_export.sdf'),('priority','vendor_pose_priority.sdf')]:
        mol=Chem.SDMolSupplier(str(src/file),removeHs=False)[0];mols[label]=mol
        assert mol.GetNumAtoms()==ref.GetNumAtoms()
        matches=mol.GetSubstructMatches(ref,uniquify=False,useChirality=True,maxMatches=100000)
        assert matches,'No full molecular graph isomorphism'
        coords=mol.GetConformer().GetPositions().copy()
        ranked=sorted((float(np.sqrt(np.mean(np.sum((coords[list(m)]-xyz)**2,axis=1)))),m) for m in matches)
        rms,mapping=ranked[0];maps[label]=mapping
        atommap=[(int(mapping[i]),i) for i in range(len(mapping))]
        calc=float(rdMolAlign.CalcRMS(mol,ref,map=[atommap]))
        assert abs(calc-rms)<1e-8
        assert np.array_equal(coords,mol.GetConformer().GetPositions()),'CalcRMS changed coordinates'
        reordered=Chem.RenumberAtoms(mol,list(mapping))
        with Chem.SDWriter(str(support/f'{label}_reordered_original_frame.sdf')) as writer:writer.write(reordered)
        # A fit is explicitly a separate demonstration, never the export integrity metric.
        fitted=Chem.Mol(mol)
        fitrms=float(rdMolAlign.AlignMol(fitted,ref,atomMap=atommap))
        with Chem.SDWriter(str(support/f'{label}_fitted_demonstration.sdf')) as writer:writer.write(fitted)
        delta=coords[list(mapping)]-xyz
        checks[label]=dict(fixed_frame_rmsd_A=rms,row_order_rmsd_A=float(np.sqrt(np.mean(np.sum((coords-xyz)**2,axis=1)))),graph_isomorphisms=len(matches),separately_fitted_rmsd_A=fitrms,mean_displacement_A=delta.mean(axis=0).tolist(),displacement_scatter_rms_A=float(np.sqrt(np.mean(np.sum((delta-delta.mean(axis=0))**2,axis=1)))) )
        for i,j in enumerate(mapping):rows.append(dict(export=label,source_index_zero_based=i,export_index_zero_based=j,source_atom_name=deposited[i]['name'],element=ref.GetAtomWithIdx(i).GetSymbol(),distance_A=float(np.linalg.norm(delta[i]))))
    with (support/'atom_maps.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    protein=[a for a in atoms if a['record']=='ATOM' and a['chain']=='A' and a['element']!='H']
    pp=np.array([a['xyz'] for a in protein]);contacts=[]
    for label,coords in [('deposited',xyz),('priority_original',mols['priority'].GetConformer().GetPositions()[list(maps['priority'])])]:
        d=np.linalg.norm(coords[:,None,:]-pp[None,:,:],axis=2)
        for i,j in zip(*np.where(d<=4.0)):
            contacts.append(dict(ligand=label,ligand_atom=deposited[i]['name'],ligand_element=ref.GetAtomWithIdx(int(i)).GetSymbol(),protein_atom=protein[j]['name'],protein_element=protein[j]['element'],residue=protein[j]['resname']+protein[j]['resid'],chain='A',distance_A=float(d[i,j])))
    with (support/'distance_contacts.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(contacts[0]));w.writeheader();w.writerows(contacts)
    result=dict(source_accession='1H1Q',ligand_identity_chain={'component':'2A6','name':'NU6094; 2-anilino-6-cyclohexylmethoxypurine','author_chain':'A','author_residue':1298,'receptor':'Thr160-phosphorylated CDK2 chain A; cyclin A2 chain B'},atom_mapping_method='All complete element/bond/stereochemistry-preserving RDKit graph isomorphisms; choose minimum heavy-atom RMSD in original coordinate frame. 24 heavy atoms, no fitting. Explicit map passed to CalcRMS; independently verified NumPy metric and coordinate immutability; export atoms reordered using RenumberAtoms.',reference_export_fixed_frame_rmsd_A=checks['reference']['fixed_frame_rmsd_A'],priority_export_fixed_frame_rmsd_A=checks['priority']['fixed_frame_rmsd_A'],priority_usable_for_contacts=False,export_details=checks,evidence_paths=['analysis/structure_audit.py','delivery/structure_support/atom_maps.csv','delivery/structure_support/distance_contacts.csv','delivery/structure_support/installed_interface.txt','delivery/figures/structure_audit.png'],interpretation='Original priority export does not preserve deposited binding-site frame. Fitted demonstration only recovers the reference molecule, not new evidence of a priority pose. Distances are geometric proximity, not proof of hydrogen bonds or selectivity.',deposited_source_resolution_A=2.5)
    inventory=list(csv.DictReader((args.inputs/'data/inventory.csv').open()))
    ced23=Chem.MolFromSmiles(next(r['smiles'] for r in inventory if r['compound_id']=='CED-23'))
    tautomer=rdMolStandardize.TautomerEnumerator()
    result['inventory_reference_identity']={
        'compound_id':'CED-23',
        'deposited_smiles':Chem.MolToSmiles(ref),
        'inventory_smiles':Chem.MolToSmiles(ced23),
        'exact_graph_match':Chem.MolToSmiles(ref)==Chem.MolToSmiles(ced23),
        'canonical_tautomer_match':Chem.MolToSmiles(tautomer.Canonicalize(ref))==Chem.MolToSmiles(tautomer.Canonicalize(ced23)),
        'interpretation':'Same tautomer family, different encoded purine proton/bond state. Canonicalization is an identity aid, not evidence of the bound protonation state.'}
    (out/'structure_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    # Clear two-dimensional orthographic views of actual 3D molecular coordinates.
    # Choose PCA plane from reference ligand; same projection/frame in both panels.
    center=xyz.mean(axis=0);_,_,vt=np.linalg.svd(xyz-center,full_matrices=False);basis=vt[:2].T
    project=lambda x:(np.asarray(x)-center)@basis
    colors={'C':'#287f9b','N':'#2854cc','O':'#d23d41'}
    plt.rcParams.update({'font.size':11})
    f,axes=plt.subplots(1,2,figsize=(15,7))
    def ligand(ax,coords,alpha=1,color=None):
        xy=project(coords)
        for bond in ref.GetBonds():
            a,b=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx();ax.plot(xy[[a,b],0],xy[[a,b],1],color=color or '#526774',lw=2.3,alpha=alpha,zorder=3)
        for i,atom in enumerate(ref.GetAtoms()):ax.scatter(*xy[i],s=65,c=color or colors.get(atom.GetSymbol(),'gray'),edgecolor='white',linewidth=.5,alpha=alpha,zorder=4)
        return xy
    ax=axes[0];ligand(ax,xyz)
    polar=[c for c in contacts if c['ligand']=='deposited' and c['ligand_element'] in ('N','O') and c['protein_element'] in ('N','O') and c['distance_A']<=3.5]
    for k,c in enumerate(polar):
        i=next(i for i,a in enumerate(deposited) if a['name']==c['ligand_atom']);j=next(j for j,a in enumerate(protein) if a['name']==c['protein_atom'] and a['resname']+a['resid']==c['residue'])
        a,b=project([xyz[i],pp[j]])
        ax.plot([a[0],b[0]],[a[1],b[1]],'--',color='#bd8425',lw=1.5)
        ax.scatter(*b,c='#bd8425',s=55,zorder=5)
        ax.annotate(f"{c['residue']} {c['protein_atom']}\n{c['distance_A']:.2f} Å to {c['ligand_atom']}",b,xytext=(12,16 if k%2 else -32),textcoords='offset points',fontsize=9,arrowprops=dict(arrowstyle='-',color='#bd8425'))
    ax.set_title('Deposited NU6094 · 1H1Q / 2A6 chain A\nPolar-atom proximity in the CDK2 binding site',loc='left',fontweight='bold');ax.margins(.30)
    ax=axes[1];ligand(ax,xyz,alpha=.7);ligand(ax,mols['priority'].GetConformer().GetPositions()[list(maps['priority'])],color='#c54f86')
    ax.set_title('Original export frame audit · no fitting\nBlue: deposited / reference; pink: priority',loc='left',fontweight='bold')
    ax.text(.02,.98,f"Mapped reference RMSD: {checks['reference']['fixed_frame_rmsd_A']:.3f} Å\nMapped priority RMSD: {checks['priority']['fixed_frame_rmsd_A']:.3f} Å\nPriority contacts are not interpretable",transform=ax.transAxes,va='top',bbox=dict(facecolor='white',alpha=.9,edgecolor='#ddd'))
    axes[1].set_ylim(-5.8,10.2)
    for ax in axes:ax.set_aspect('equal');ax.set_xlabel('Projection coordinate (Å)');ax.set_ylabel('Projection coordinate (Å)');ax.spines[['top','right']].set_visible(False);ax.grid(alpha=.12)
    f.text(.02,.02,'Same orthographic projection of supplied 3D coordinates. Dashed lines show distances, not assigned hydrogen bonds.\nThe source is NU6094, the supplied structural reference for CED-23; no pose transfer to the main series. No CDK1 selectivity follows from this structure.',fontsize=10)
    f.tight_layout(rect=[0,.08,1,1]);f.savefig(fig/'structure_audit.png',dpi=180,bbox_inches='tight');f.savefig(fig/'structure_audit.svg',bbox_inches='tight');plt.close(f)
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
