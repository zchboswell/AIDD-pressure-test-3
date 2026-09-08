"""Explicit decision policy, evidence accounting and menu-constrained plan; no TE model."""
import argparse,csv,json,math
from decimal import Decimal
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw
from curate import read_csv,write_csv

ADVANCE=['CED-02','CED-03','CED-08','CED-13','CED-17','CED-22']
DIAGNOSTIC=['CED-23']
REASONS={
'CED-01':'Reserve: TE65nM but CDK1 90nM and HLM14min; known narrow cross-assay profile adds less decision value than repeats of unknowns.',
'CED-02':'Advance: TE130nM/HLM22min; lactam analog adds chemistry contrast. Confirm borderline HLM and matched CDK2/CDK1 panel.',
'CED-03':'Advance provisionally: TE and CDK1 unmeasured; SOL35uM/HLM38min. Use VIAL-03 only; IDENTITY first due separate conflicting export; duplicate prospective ENGAGE.',
'CED-04':'Reserve outside this cycle: arrival day10 exceeds day8; TE850nM also weak despite HLM55min.',
'CED-05':'Reserve: purity72% requires IDENTITY before any new biology; TE35/CDK1 45nM and HLM8min do not justify spending on this stock.',
'CED-06':'Reserve: TE65nM but CDK1 90nM/HLM14min and SOL3uM; potency alone insufficient.',
'CED-07':'Reserve behind CED02/17/22: same TE130/CDK1 400/HLM22 profile but SOL8uM; no additional named decision question.',
'CED-08':'Advance provisionally: missing TE; CDK1 1800nM/HLM38min/SOL18uM support high information value, not potency prediction. Duplicate prospective ENGAGE.',
'CED-09':'Reserve: TE850nM weak despite HLM55min; missing CDK1 does not establish a useful window.',
'CED-10':'Reserve: TE35nM but CDK1 45nM/HLM8min; deprioritize potency-only shortcut.',
'CED-11':'Reserve: missing TE but CDK1 90nM/HLM14min reduce decision value versus CED08/03.',
'CED-12':'Reserve outside this cycle: arrival day10 exceeds day8; SOL3uM also unfavorable.',
'CED-13':'Conditional advance: TE320nM misses preference but CDK1 1800nM/HLM38min offers explicit potency/window tradeoff; SOL8uM warrants concentration QC.',
'CED-14':'Reserve: TE850nM weak and CDK1 missing despite HLM55min; not preferred over stronger current profiles.',
'CED-15':'Reserve: TE35nM but CDK1 45nM/HLM8min; not rescued by SOL35uM.',
'CED-16':'Reserve: missing TE but CDK1 90nM/HLM14min; lower information value than CED08/03.',
'CED-17':'Advance: TE130/HLM22/SOL120; cyclopropyl counterpart to tertbutyl CED22. Historical IC50 700nM vs54.9nM tests applicability of old ranking, not assay contradiction.',
'CED-18':'Reserve: historical IC50 22.5nM does not overcome current TE320nM and SOL3uM; prefer CED13 with same TE/CDK1/HLM and higher solubility.',
'CED-19':'Reserve fallback: HLM55min/CDK1 6000nM but TE850nM misses current preference; apparent cross-assay separation is not selectivity.',
'CED-20':'Reserve: TE35nM but CDK1 45nM/HLM8min; historical IC50143nM is contextual only.',
'CED-21':'Reserve: TE65nM but CDK1 90nM/HLM14min; matched window unproven.',
'CED-22':'Advance leading observed anchor: TE130/HLM22/SOL70; confirm borderline stability and matched window. Compare CED17 substitution pair.',
'CED-23':'Diagnostic: deposited-ligand tautomer family/different chemotype, SOL120/HLM38; TE and CDK1 unknown. Source structure is not a biological positive control.',
'CED-24':'Reserve structural analog: both endpoints unknown and SOL3uM versus CED23 SOL120; no supplied CED24 pose or current-window evidence.'}

def main(inputs,out):
    out.mkdir(parents=True,exist_ok=True);(out/'support').mkdir(exist_ok=True);(out/'figures').mkdir(exist_ok=True)
    inventory=read_csv(out/'support/inventory_identity.csv');reads={r['compound_id']:r for r in read_csv(inputs/'data/program_readouts.csv')}
    menu={r['experiment_id']:r for r in read_csv(inputs/'data/experiment_menu.csv')}
    curated=read_csv(out/'curation.csv');selection=[];predictions=[];valid={}
    for r in inventory:
        cid=r['compound_id'];status=r['identity_status'];reason=r['identity_reason']
        if status=='valid':
            status='advance' if cid in ADVANCE else 'diagnostic' if cid in DIAGNOSTIC else 'reserve'
            reason=REASONS[cid];valid[cid]=r
        selection.append(dict(request_id=r['request_id'],compound_id=cid,status=status,reason=reason))
        p=dict(request_id=r['request_id'],compound_id=cid,status='abstained',endpoint='cellular CDK2 engagement IC50',assay_id='SIM-TE-01',p_activity='',lower='',upper='',reason='Unmeasured in supplied program; no validated endpoint-specific model or uncertainty coverage.')
        if r['identity_status'] in ['duplicate','quarantined']:
            p['status']=r['identity_status'];p['reason']=r['identity_reason'];p['assay_id']=''
        elif reads[cid]['cell_target_engagement_ic50_nM']:
            p['status']='observed';p['p_activity']=9-math.log10(float(reads[cid]['cell_target_engagement_ic50_nM']));p['reason']='Supplied simulated program observation, nominal compound ID; not a public measurement or forecast. No replicate-derived interval available.'
            if cid=='CED-05':p['reason']+=' Current stock purity72%; identity/purity check required before interpreting new biological results.'
        predictions.append(p)
    write_csv(out/'selection.csv',selection,['request_id','compound_id','status','reason'])
    write_csv(out/'predictions.csv',predictions,['request_id','status','endpoint','assay_id','p_activity','lower','upper','reason','compound_id'])
    exp=[]
    for cid in ADVANCE+DIAGNOSTIC:
        for assay,purpose in [('CDK2-E1','Establish matched 1mM ATP CDK2 arm; historical assay is not a substitute.'),('CDK1-B','Establish matched 1mM ATP CDK1 arm; calculate biochemical window jointly with CDK2-E1.'),('ENGAGE','Measure current cellular CDK2 endpoint; review full curves and concentration feasibility.')]:
            n=2 if cid in ['CED-03','CED-08'] and assay=='ENGAGE' else 1
            if n==2:purpose+=' Two separately charged prospective experiments scheduled in parallel; assess consistency, not assumed calibrated intervals.'
            exp.append(dict(compound_id=cid,experiment_id=assay,replicates=n,purpose=purpose))
    exp.append(dict(compound_id='CED-03',experiment_id='IDENTITY',replicates=1,purpose='Check VIAL-03 identity/purity before biology because STOCK-28 has a conflicting graph; STOCK-28 remains quarantined.'))
    for cid in ['CED-02','CED-22']:
        exp.append(dict(compound_id=cid,experiment_id='MIC',replicates=1,purpose='Prospective confirmation of HLM half-life near20min preference; supplied value22min has no uncertainty.'))
    write_csv(out/'experiments.csv',exp,['compound_id','experiment_id','replicates','purpose'])
    budget=[];matrix=[]
    for cid in ADVANCE+DIAGNOSTIC:
        r=valid[cid];ee=[e for e in exp if e['compound_id']==cid]
        cost=sum(Decimal(menu[e['experiment_id']]['cost_usd'])*e['replicates'] for e in ee)
        mass=sum(Decimal(menu[e['experiment_id']]['mass_mg'])*e['replicates'] for e in ee)
        identity=sum(int(menu[e['experiment_id']]['turnaround_days'])*e['replicates'] for e in ee if e['experiment_id']=='IDENTITY')
        days=int(r['delivery_days'])+identity+max(int(menu[e['experiment_id']]['turnaround_days']) for e in ee if e['experiment_id']!='IDENTITY')
        b=dict(compound_id=cid,request_id=r['request_id'],vial_id=r['vial_id'],handling_usd=r['handling_usd'],experiments_usd=str(cost),total_usd=str(cost+Decimal(r['handling_usd'])),mass_used_mg=str(mass),stock_mass_mg=r['mass_mg'],remaining_mass_mg=str(Decimal(r['mass_mg'])-mass),arrival_day=r['delivery_days'],identity_days=identity,completion_day=days)
        assert mass<=Decimal(r['mass_mg']) and days<=8
        budget.append(b)
        history=[f"{x['relation_nm']}{float(x['value_nm']):g}" for x in curated if x['chemical_id']==r['source_id'] and x['assay_id']=='CHEMBL5736732' and x['status'] not in ['duplicate','quarantined']]
        matrix.append(dict(compound_id=cid,role='advance' if cid in ADVANCE else 'diagnostic',**{k:reads[cid][k] for k in ['cell_target_engagement_ic50_nM','cdk1_cyclinB_ic50_nM','kinetic_solubility_uM','microsomal_t_half_min']},historical_CHEMBL5736732_IC50_nM=';'.join(history),total_usd=b['total_usd'],completion_day=days))
    write_csv(out/'support/budget_material_schedule.csv',budget,list(budget[0]))
    write_csv(out/'support/decision_matrix.csv',matrix,list(matrix[0]))
    total=sum(Decimal(b['total_usd']) for b in budget)
    summary=dict(total_usd=str(total),unspent_usd=str(Decimal(18000)-total),unique_compounds=len(budget),advance_count=len(ADVANCE),diagnostic_count=len(DIAGNOSTIC),latest_completion_day=max(b['completion_day'] for b in budget),schedule_assumption='Arrival day plus identity first (when ordered) plus longest parallel assay. All repeats prebooked in parallel, not after seeing results.',new_experimental_runs=sum(e['replicates'] for e in exp),assay_combinations=len(exp),all_mass_checks_pass=True,all_deadline_checks_pass=True)
    assert total<=18000 and len(budget)<=8 and len(ADVANCE)<=6 and len(DIAGNOSTIC)<=2
    (out/'support/plan_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    fig,ax=plt.subplots(figsize=(13.4,5.2));ax.axis('off')
    labels=['Compound / role','Cell CDK2\nIC50 nM','CDK1-B\nIC50 nM','Solubility\nµM','HLM t½\nmin','Historical\nIC50 nM','Cycle\nUSD','Finish\nday']
    data=[]
    for m in matrix:
        data.append([m['compound_id']+(' / reference' if m['role']=='diagnostic' else ' / advance'),m['cell_target_engagement_ic50_nM'] or 'unmeasured',m['cdk1_cyclinB_ic50_nM'] or 'unmeasured',m['kinetic_solubility_uM'],m['microsomal_t_half_min'],m['historical_CHEMBL5736732_IC50_nM'] or 'not supplied',f"{int(m['total_usd']):,}",str(m['completion_day'])])
    t=ax.table(cellText=data,colLabels=labels,cellLoc='center',colWidths=[.20,.13,.13,.10,.09,.13,.10,.07],loc='center');t.auto_set_font_size(False);t.set_fontsize(10);t.scale(1,2.15)
    for (rr,cc),cell in t.get_celld().items():
        cell.set_edgecolor('#d9e1e8')
        if rr==0:
            cell.set_facecolor('#203c52');cell.set_text_props(color='white',weight='bold');cell.set_height(cell.get_height()*1.5)
        else:
            cell.set_facecolor('#eef5f8' if rr%2 else '#ffffff')
            if cc in [1,2] and data[rr-1][cc]=='unmeasured':cell.set_facecolor('#fff0cd')
    ax.set_title('CEDAR: measure the gaps that can change the decision',fontsize=17,weight='bold',loc='left',pad=18)
    fig.text(.06,.095,'Current cellular / CDK1 / ADME values are simulated program observations. Historical IC50: public CHEMBL5736732.\nCross-assay values do not establish selectivity. Proposed work: $17,500, seven vials, latest day 7; no new assay results yet.',fontsize=10,color='#34495e')
    fig.subplots_adjust(top=.83,bottom=.22,left=.035,right=.995)
    fig.savefig(out/'figures/decision_comparison.png',dpi=180);fig.savefig(out/'figures/decision_comparison.svg');plt.close(fig)
    mols=[Chem.MolFromSmiles(valid[c]['smiles']) for c in ADVANCE+DIAGNOSTIC]
    legends=[c+' | '+('structural reference family' if c=='CED-23' else 'conditional advance' if c=='CED-13' else 'prospective gap' if c in ['CED-03','CED-08'] else 'observed anchor') for c in ADVANCE+DIAGNOSTIC]
    image=Draw.MolsToGridImage(mols,molsPerRow=3,subImgSize=(480,320),legends=legends,useSVG=False)
    image.save(str(out/'figures/selected_chemical_structures.png'))
    print(json.dumps(summary,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,default=Path('inputs'));p.add_argument('--out',type=Path,default=Path('delivery'));a=p.parse_args();main(a.inputs,a.out)
