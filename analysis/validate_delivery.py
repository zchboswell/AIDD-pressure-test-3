"""Independent accounting/contract checks; does not certify scientific truth."""
import argparse,csv,json,math
from collections import Counter,defaultdict
from decimal import Decimal
from pathlib import Path

def rows(p):
    with p.open() as f:return list(csv.DictReader(f))
def exact_ids(actual,expected,key):
    assert Counter(r[key] for r in actual)==Counter(r[key] for r in expected)
    assert len({r[key] for r in actual})==len(actual)
def finite(r,keys):
    for k in keys:
        if r[k]!='':assert math.isfinite(float(r[k])),(r,k)

def main(inputs,out):
    inv=rows(inputs/'data/inventory.csv');acts=rows(inputs/'data/activity_archive.csv');cur=rows(out/'curation.csv');sel=rows(out/'selection.csv');pred=rows(out/'predictions.csv');exp=rows(out/'experiments.csv');menu={r['experiment_id']:r for r in rows(inputs/'data/experiment_menu.csv')}
    exact_ids(cur,acts,'row_id');exact_ids(sel,inv,'request_id');exact_ids(pred,inv,'request_id')
    for r in cur:
        assert r['status'] in ['eligible','context_only','duplicate','quarantined'];finite(r,['value_nm','p_activity'])
        if r['p_activity']!='':
            assert float(r['value_nm'])>0
            assert abs(float(r['p_activity'])-(9-math.log10(float(r['value_nm']))))<1e-10
            assert r['p_relation']=={'=':'=','<':'>','<=':'>=','>':'<','>=':'<='}[r['relation_nm']]
        if r['status']=='eligible':assert r['endpoint']=='IC50' and r['relation_nm']=='=' and r['assay_id']=='CHEMBL5736732'
    c={r['row_id']:r for r in cur}
    expected={'ACT-9003':('context_only','IC50',10000,5,'<'),'ACT-9004':('context_only','IC50',10,8,'>='),'ACT-9008':('context_only','Ki',300,6.522878745280337,'<='),'ACT-9009':('context_only','IC50',100,7,'=')}
    for rid,(status,endpoint,nm,p,relation) in expected.items():
        r=c[rid];assert (r['status'],r['endpoint'],r['p_relation'])==(status,endpoint,relation);assert abs(float(r['value_nm'])-nm)<1e-8 and abs(float(r['p_activity'])-p)<1e-10
    for rid in ['ACT-9005','ACT-9006','ACT-9007','ACT-9010']:assert c[rid]['status']=='quarantined' and c[rid]['p_activity']==''
    groups=defaultdict(list)
    for r in acts:
        if r['source_activity_id']:groups[r['source_activity_id']].append(c[r['row_id']])
    for g in groups.values():
        valid=[r for r in g if r['status']!='quarantined']
        if valid:assert sum(r['status']!='duplicate' for r in valid)==1
    for r in sel:assert r['status'] in ['advance','diagnostic','reserve','duplicate','quarantined']
    sm={r['request_id']:r for r in sel}
    for k in ['STOCK-25','STOCK-26']:assert sm[k]['status']=='duplicate'
    for k in ['STOCK-27','STOCK-28']:assert sm[k]['status']=='quarantined'
    for k in ['STOCK-04','STOCK-12']:assert sm[k]['status']=='reserve'
    current={r['compound_id']:r for r in rows(inputs/'data/program_readouts.csv')}
    for r in pred:
        assert r['status'] in ['predicted','observed','abstained','duplicate','quarantined'];finite(r,['p_activity','lower','upper'])
        if r['status'] in ['predicted','observed']:assert any(r[k]!='' for k in ['p_activity','lower','upper'])
        if r['lower']!='' and r['upper']!='':assert float(r['lower'])<=float(r['upper'])
        if r['status']=='observed':
            v=float(current[r['compound_id']]['cell_target_engagement_ic50_nM']);assert abs(float(r['p_activity'])-(9-math.log10(v)))<1e-10 and r['assay_id']=='SIM-TE-01'
        if r['status']=='abstained':assert current[r['compound_id']]['cell_target_engagement_ic50_nM']==''
    chosen=[r for r in sel if r['status'] in ['advance','diagnostic']];assert sum(r['status']=='advance' for r in chosen)<=6 and sum(r['status']=='diagnostic' for r in chosen)<=2
    im={r['request_id']:r for r in inv};vials=[im[r['request_id']]['vial_id'] for r in chosen];assert len(vials)==len(set(vials))<=8
    chosen_comp={r['compound_id']:im[r['request_id']] for r in chosen};assert set(e['compound_id'] for e in exp)==set(chosen_comp)
    assert len(set((e['compound_id'],e['experiment_id']) for e in exp))==len(exp)
    total=Decimal(0);checks=[]
    for cid,v in chosen_comp.items():
        ee=[r for r in exp if r['compound_id']==cid];mass=Decimal(0);cost=Decimal(v['handling_usd']);identity=0;parallel=0
        for e in ee:
            n=Decimal(e['replicates']);assert n.is_finite() and n>0 and n==n.to_integral_value();assert e['purpose'];m=menu[e['experiment_id']]
            cost+=n*Decimal(m['cost_usd']);mass+=n*Decimal(m['mass_mg'])
            if e['experiment_id']=='IDENTITY':identity+=int(n)*int(m['turnaround_days'])
            else:parallel=max(parallel,int(m['turnaround_days']))
        if Decimal(v['purity_pct'])<90:assert identity>0
        end=int(v['delivery_days'])+identity+parallel;assert end<=8 and mass<=Decimal(v['mass_mg']);total+=cost
        checks.append({'compound_id':cid,'cost_usd':str(cost),'mass_mg':str(mass),'completion_day':end})
    assert total==Decimal('17500') and total<=18000
    st=json.loads((out/'structure_checks.json').read_text());assert st['source_accession']=='1H1Q';assert type(st['priority_usable_for_contacts']) is bool and st['priority_usable_for_contacts'] is False
    assert abs(float(st['reference_export_fixed_frame_rmsd_A']))<.001;assert abs(float(st['priority_export_fixed_frame_rmsd_A'])-math.sqrt(77))<.001
    assert st['evidence_paths'] and all(isinstance(x,str) and x for x in st['evidence_paths'])
    result={'status':'pass','activity_rows':len(cur),'inventory_rows':len(sel),'prediction_status_counts':dict(Counter(r['status'] for r in pred)),'unique_compounds':len(chosen),'total_usd':str(total),'compound_checks':checks,'scope':'Accounting, declared endpoint, known converter sentinels, source-duplicate accounting, output numerics and fixed-frame reference values. Does not certify biological efficacy or future assay success.'}
    (out/'support/validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,default=Path('inputs'));p.add_argument('--out',type=Path,default=Path('delivery'));a=p.parse_args();main(a.inputs,a.out)
