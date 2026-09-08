"""Headless PyMOL depiction of deposited receptor/ligand; no pose generation or fitting."""
import argparse,json
from pathlib import Path
import pymol
from pymol import cmd

def main(inputs,out):
    figures=out/'figures';figures.mkdir(parents=True,exist_ok=True)
    cmd.reinitialize();cmd.set('max_threads',1);cmd.set('ray_shadows',0);cmd.set('antialias',2)
    cmd.set('ray_opaque_background',1);cmd.bg_color('white')
    cmd.load(str(inputs/'structures/1h1q.pdb'),'deposited')
    cmd.remove('not chain A');cmd.remove('solvent');cmd.hide('everything')
    cmd.select('ligand','resn 2A6 and chain A and resi 1298')
    cmd.select('hinge','polymer and chain A and resi 81+82+83')
    cmd.show('cartoon','polymer');cmd.color('gray80','polymer');cmd.set('cartoon_transparency',.35)
    cmd.show('sticks','ligand or hinge');cmd.set('stick_radius',.18)
    cmd.color('cyan','ligand and elem C');cmd.color('orange','hinge and elem C')
    cmd.color('blue','(ligand or hinge) and elem N');cmd.color('red','(ligand or hinge) and elem O')
    cmd.orient('ligand');cmd.zoom('ligand',7)
    cmd.set('label_color','black');cmd.set('label_size',18);cmd.set('label_outline_color','white')
    cmd.label('chain A and resi 81 and name CA','"Glu81"');cmd.label('chain A and resi 83 and name CA','"Leu83"')
    cmd.png(str(figures/'deposited_pocket.png'),width=1400,height=1000,dpi=180,ray=1)
    info={'source':'inputs/structures/1h1q.pdb','selection':'deposited chain A, ligand 2A6 author residue1298; hinge81-83','ligand_atoms':cmd.count_atoms('ligand'),'hinge_atoms':cmd.count_atoms('hinge'),'pymol_version':cmd.get_version()[0],'interpretation':'Source crystal pose only. Cartoon CDK2, cyan ligand carbon, orange hinge carbon, blue nitrogen/red oxygen. Displayed proximity does not prove hydrogen bonds or CDK1 selectivity. No modeled candidate pose; no fitting.'}
    (out/'support/pocket_render.json').write_text(json.dumps(info,indent=2)+'\n');print(json.dumps(info))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,default=Path('inputs'));p.add_argument('--out',type=Path,default=Path('delivery'));a=p.parse_args();main(a.inputs,a.out)
