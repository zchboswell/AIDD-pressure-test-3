"""Deterministic row-level curation with provenance and exact chemical identity checks."""
import argparse, csv, hashlib, json, math
from collections import Counter, defaultdict
from pathlib import Path
from rdkit import Chem
from assay_math import normalize


def read_csv(p):
    with p.open() as f:return list(csv.DictReader(f))

def write_csv(p,rows,fields):
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader()
        for row in rows:
            row=dict(row)
            for key in ('p_activity','lower','upper'):
                if key in row and row[key]!='':
                    value=str(row[key])
                    if 'e' not in value.lower():
                        if '.' not in value:value+='.'
                        value+='0'*max(0,3-len(value.split('.')[1]))
                    row[key]=value
            w.writerow(row)

def canonical(s):
    m=Chem.MolFromSmiles(s)
    return Chem.MolToSmiles(m,isomericSmiles=True) if m else None

def main(root,out):
    out.mkdir(parents=True,exist_ok=True);(out/'support').mkdir(exist_ok=True)
    rows=read_csv(root/'data/activity_archive.csv')
    sources={str(r['activity_id']):r for r in json.loads((root/'references/chembl_selected_records.json').read_text())}
    seen={}; result=[]; identities=[]; mismatch=[]
    for r in rows:
        n=normalize(r['value'],r['unit'],r['relation'],r['endpoint']);can=canonical(r['smiles'])
        group='CHEM-'+hashlib.sha256(can.encode()).hexdigest()[:16] if can else ''
        d=dict(row_id=r['row_id'],status='context_only',endpoint=r['endpoint'],value_nm='',relation_nm='',p_activity='',p_relation='',chemical_group=group,assay_group=r['assay_id']+'|'+r['endpoint'],reason='',chemical_id=r['chemical_id'],assay_id=r['assay_id'],source_activity_id=r['source_activity_id'],canonical_smiles=can or '')
        if n['status']=='ok':
            d.update({k:n[k] for k in ('value_nm','relation_nm','p_activity','p_relation')})
        issues=[];source=sources.get(r['source_activity_id'])
        if n['status']!='ok':issues.append(n['reason'])
        if not can: issues.append('invalid molecular graph')
        if r['data_validity']:issues.append('source validity flag: '+r['data_validity'])
        if source:
            sn=normalize(source['standard_value'],source['standard_units'],source['standard_relation'],source['standard_type'])
            checks={'chemical_id':r['chemical_id']==source['molecule_chembl_id'],'graph':can==canonical(source['canonical_smiles']),'assay':r['assay_id']==source['assay_chembl_id'],'document':r['source_document']==source['document_chembl_id'],'endpoint':r['endpoint']==source['standard_type'],'normalized_value':n['status']=='ok' and sn['status']=='ok' and math.isclose(n['value_nm'],sn['value_nm'],rel_tol=1e-9),'relation':n.get('relation_nm')==sn.get('relation_nm')}
            if not all(checks.values()):issues.append('source record mismatch');mismatch.append({'row_id':r['row_id'],'checks':checks})
        elif r['source_activity_id']:issues.append('source_activity_id absent from supplied reference snapshot')
        if issues:
            d['status']='quarantined';d['reason']='; '.join(issues)
            # Invalid chemistry cannot support an activity number assigned to a known chemical.
            if not can:
                for k in ['value_nm','relation_nm','p_activity','p_relation']:d[k]=''
        elif r['source_activity_id'] in seen:
            d['status']='duplicate';d['reason']='same source_activity_id as '+seen[r['source_activity_id']]+'; no independent observation'
        else:
            if r['source_activity_id']:seen[r['source_activity_id']]=r['row_id']
            if r['assay_id']=='CHEMBL5736732' and r['endpoint']=='IC50' and n['relation_nm']=='=' and source:
                d['status']='eligible';d['reason']='source-verified exact IC50 in declared historical assay; not current engagement'
            elif r['assay_id']=='SIM-QC-CHECK':
                d['reason']='fictional export/QC example; preserve endpoint and censoring only; not training or program biology'
            elif n['relation_nm']!='=':d['reason']='censored assay-specific threshold; not an exact regression label'
            else:d['reason']='different historical assay context; not pooled into primary-assay regression'
            if r['potential_duplicate']=='1':d['reason']+='; potential_duplicate flag retained, not proof of reingestion'
        result.append(d)
    write_csv(out/'curation.csv',result,['row_id','status','endpoint','value_nm','relation_nm','p_activity','p_relation','chemical_group','assay_group','reason','chemical_id','assay_id','source_activity_id','canonical_smiles'])
    inv=read_csv(root/'data/inventory.csv'); first_vial={};first_id={}; inventory=[]
    source_chem=defaultdict(set)
    for r in result:
        if r['status']!='quarantined':source_chem[r['chemical_id']].add(r['canonical_smiles'])
    for r in inv:
        can=canonical(r['smiles']);state='valid';why='valid distinct vial and chemical graph'
        if not can:state='quarantined';why='invalid SMILES'
        elif source_chem[r['source_id']] and can not in source_chem[r['source_id']]:state='quarantined';why='source identifier/graph conflicts with supplied activity reference'
        elif r['compound_id'] in first_id and can!=first_id[r['compound_id']]['canonical_smiles']:state='quarantined';why='compound identifier has conflicting graph; differs from established inventory representative '+first_id[r['compound_id']]['request_id']
        elif r['vial_id'] in first_vial:
            previous=first_vial[r['vial_id']]
            if can==previous['canonical_smiles'] and r['compound_id']==previous['compound_id']:state='duplicate';why='same physical vial as '+previous['request_id']+'; stock not additive'
            else:state='quarantined';why='vial identifier graph/compound conflict'
        x=dict(r,canonical_smiles=can or '',identity_status=state,identity_reason=why)
        if state=='valid':first_vial[r['vial_id']]=x;first_id[r['compound_id']]=x
        inventory.append(x)
    write_csv(out/'support/inventory_identity.csv',inventory,list(inventory[0]))
    summary={'declared_model_endpoint':'Exact pIC50 for CHEMBL5736732 only; assay-specific historical biochemical endpoint, not cellular engagement','chemical_group_method':'RDKit canonical isomeric SMILES SHA256 prefix; preserve fragments, charge, stereochemistry; no tautomer/salt normalization','row_counts':dict(Counter(r['status'] for r in result)),'assay_status_counts':{a:dict(Counter(r['status'] for r in result if r['assay_id']==a)) for a in sorted({r['assay_id'] for r in result})},'source_mismatches':mismatch,'source_verified_rows':sum(bool(sources.get(r['source_activity_id'])) for r in rows),'inventory_counts':dict(Counter(r['identity_status'] for r in inventory)),'eligible_chemical_groups':len({r['chemical_group'] for r in result if r['status']=='eligible'})}
    (out/'support/curation_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,default=Path('inputs'));p.add_argument('--out',type=Path,default=Path('delivery'));a=p.parse_args();main(a.inputs,a.out)
