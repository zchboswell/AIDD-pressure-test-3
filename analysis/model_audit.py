#!/usr/bin/env python3
"""Historical-only, fixed-specification model diagnostic. No prospective TE model.
Run only within an approved compute lease. All source inputs remain unchanged.
"""
import argparse, csv, json, pathlib, subprocess, sys
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupKFold, train_test_split
from sklearn.metrics import mean_absolute_error, r2_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def rows(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f))


def save_csv(path, data):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(data[0]))
        w.writeheader(); w.writerows(data)


def score(y, p):
    return {'n': len(y), 'mae_log10': float(mean_absolute_error(y, p)),
            'r2': float(r2_score(y, p)) if len(y) > 1 else None}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--inputs', type=pathlib.Path, default=pathlib.Path('inputs'))
    ap.add_argument('--curation', type=pathlib.Path, default=pathlib.Path('delivery/curation.csv'))
    ap.add_argument('--output', '--out', type=pathlib.Path, default=pathlib.Path('delivery'))
    a = ap.parse_args(); out = a.output; art = out/'model_artifacts'
    art.mkdir(parents=True, exist_ok=True); (out/'figures').mkdir(parents=True, exist_ok=True)
    # Reproduce the supplied implementation without repairs, record stdout AND stderr.
    original = subprocess.run([sys.executable, str(a.inputs/'legacy/vendor_model.py'), str(a.inputs)],
                              capture_output=True, text=True, check=False)
    (art/'vendor_original_stdout.txt').write_text(original.stdout)
    (art/'vendor_original_stderr.txt').write_text(original.stderr)
    archive = rows(a.inputs/'data/activity_archive.csv')
    acts = {r['row_id']:r for r in archive}
    feats = {r['row_id']:r for r in rows(a.inputs/'data/vendor_features.csv')}
    eligible = [r for r in rows(a.curation) if r['status']=='eligible']
    assert eligible, 'No eligible historical labels; no model fit performed.'
    assert all(acts[r['row_id']]['assay_id']=='CHEMBL5736732' and r['endpoint']=='IC50'
               and r['relation_nm']=='=' for r in eligible), 'Endpoint must be fixed before fitting'
    assert len({acts[r['row_id']]['source_activity_id'] for r in eligible}) == len(eligible)
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    mols = [Chem.MolFromSmiles(acts[r['row_id']]['smiles']) for r in eligible]
    assert all(m is not None for m in mols)
    chem = np.array([Chem.MolToSmiles(m, isomericSmiles=True) for m in mols])
    scaffold = np.array([MurckoScaffold.MurckoScaffoldSmiles(mol=m, includeChirality=False) or '<acyclic>' for m in mols])
    fps = [generator.GetFingerprint(m) for m in mols]
    X = np.array([generator.GetFingerprintAsNumPy(m) for m in mols])
    desc = np.array([[float(feats[r['row_id']][k]) for k in ['molecular_weight','logp','heavy_atoms']] for r in eligible])
    leak = np.array([float(feats[r['row_id']]['vendor_activity_index']) for r in eligible])
    y = np.array([float(r['p_activity']) for r in eligible])
    spec = {'n_estimators':80, 'random_state':7, 'n_jobs':1}
    summary = {'endpoint':'historical exact IC50; CHEMBL5736732 only', 'original_exit_code':original.returncode,
      'n_eligible_rows':len(y), 'n_canonical_chemicals':len(set(chem)), 'n_murcko_scaffolds':len(set(scaffold)),
      'model_spec':spec, 'fingerprint':'Morgan radius=2, 2048 binary bits, RDKit default nonchiral fingerprint',
      'leak_feature_mae_to_curated_label':float(np.mean(abs(leak-y))),
      'leak_feature_pearson_r':float(np.corrcoef(leak,y)[0,1]), 'evaluations':{}}
    # Same test rows across ablations; score is diagnostic, not model selection.
    tr, te = train_test_split(np.arange(len(y)), test_size=.25, random_state=7)
    summary['curated_random_split_chemical_overlap'] = len(set(chem[tr]) & set(chem[te]))
    for name, x in [('retrospective_feature_plus_descriptors',np.column_stack([desc,leak])),
                    ('descriptors_only',desc),('fingerprints_only',X)]:
        p = RandomForestRegressor(**spec).fit(x[tr], y[tr]).predict(x[te])
        summary['evaluations']['random_'+name] = score(y[te],p)
    summary['evaluations']['random_train_median'] = score(y[te],np.repeat(np.median(y[tr]),len(te)))
    pred_rows=[]
    for split_name, groups in [('chemical_group',chem), ('murcko_scaffold',scaffold)]:
        ng = len(set(groups))
        if ng<2:
            summary['evaluations'][split_name]={'status':'unavailable: fewer than two groups'}
            continue
        cv = GroupKFold(n_splits=min(5,ng)); preds=np.empty(len(y)); baseline=np.empty(len(y)); sims=np.empty(len(y)); folds=[]
        for fold,(tr,te) in enumerate(cv.split(X,y,groups),1):
            assert not (set(groups[tr]) & set(groups[te]))
            preds[te]=RandomForestRegressor(**spec).fit(X[tr],y[tr]).predict(X[te])
            baseline[te]=np.median(y[tr])
            for j in te:
                sims[j]=max(DataStructs.BulkTanimotoSimilarity(fps[j],[fps[k] for k in tr]))
                pred_rows.append({'row_id':eligible[j]['row_id'],'split':split_name,'fold':fold,
                                  'chemical_group':chem[j],'scaffold':scaffold[j], 'observed_pIC50':y[j],
                                  'predicted_pIC50':preds[j], 'training_median_pIC50':baseline[j],
                                  'nearest_train_Tanimoto':sims[j]})
            folds.append({'fold':fold,'n_train':len(tr),'n_test':len(te),'fingerprint':score(y[te],preds[te]),
                          'train_median':score(y[te],baseline[te]), 'nearest_train_Tanimoto_median':float(np.median(sims[te]))})
        summary['evaluations'][split_name]={'fingerprint':score(y,preds),'train_median':score(y,baseline),
            'folds':folds,'nearest_train_Tanimoto_min':float(min(sims)), 'nearest_train_Tanimoto_median':float(np.median(sims))}
    save_csv(art/'historical_out_of_fold.csv', pred_rows)
    inventory=[]
    for r in rows(a.inputs/'data/inventory.csv'):
        m=Chem.MolFromSmiles(r['smiles']); result={'request_id':r['request_id'], 'compound_id':r['compound_id'],
            'source_id':r['source_id'],'valid_smiles':m is not None, 'canonical_match_rows':'',
            'nearest_training_row':'','nearest_Tanimoto':'','endpoint_scope':'historical IC50 only; not cellular engagement'}
        if m is not None:
            c=Chem.MolToSmiles(m,isomericSmiles=True); sim=np.array(DataStructs.BulkTanimotoSimilarity(generator.GetFingerprint(m),fps))
            j=int(sim.argmax()); result.update(canonical_match_rows=';'.join(eligible[k]['row_id'] for k in np.where(chem==c)[0]),
                 nearest_training_row=eligible[j]['row_id'],nearest_Tanimoto=float(sim[j]))
        inventory.append(result)
    save_csv(art/'inventory_applicability.csv',inventory)
    # Audit original random rows separately using the original normalizer import.
    import importlib.util
    sp=importlib.util.spec_from_file_location('untouched_normalizer',a.inputs/'legacy/assay_math.py')
    legacy=importlib.util.module_from_spec(sp); sp.loader.exec_module(legacy)
    original_rows=[]
    for rid,f in feats.items():
        r=acts[rid]
        if f['vendor_activity_index'] and legacy.normalize(r['value'],r['unit'],r['relation'],r['endpoint'])['status']=='ok':
            m=Chem.MolFromSmiles(r['smiles'])
            original_rows.append((rid,r['source_activity_id'], Chem.MolToSmiles(m) if m else '<invalid>'))
    tr,te=train_test_split(np.arange(len(original_rows)),test_size=.25,random_state=7)
    summary['original_random_split_audit']={'rows':len(original_rows),'train_rows':len(tr),'test_rows':len(te),
       'shared_source_activity_ids':sorted(({original_rows[i][1] for i in tr}&{original_rows[i][1] for i in te})-{''}),
       'shared_canonical_chemicals':len({original_rows[i][2] for i in tr}&{original_rows[i][2] for i in te}),
       'assays':sorted({acts[r[0]]['assay_id'] for r in original_rows}),
       'endpoints':sorted({acts[r[0]]['endpoint'] for r in original_rows}),
       'nonexact_rows':sum(acts[r[0]]['relation']!='=' for r in original_rows)}
    (art/'model_summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    fig,ax=plt.subplots(1,2,figsize=(11,4.4),layout='constrained')
    ax[0].scatter(y,leak,s=17,alpha=.7); lo=min(y.min(),leak.min())-.1; hi=max(y.max(),leak.max())+.1
    ax[0].plot([lo,hi],[lo,hi],color='black',lw=1); ax[0].set(xlabel='Curated historical pIC50',ylabel='Retrospective vendor activity index',title='Post-assay information mirrors the label')
    labels=[]; rf=[]; med=[]
    for key,label in [('chemical_group','Chemical-group CV'),('murcko_scaffold','Scaffold CV')]:
        e=summary['evaluations'][key]
        if 'fingerprint' in e:
            labels.append(label); rf.append(e['fingerprint']['mae_log10']); med.append(e['train_median']['mae_log10'])
    xx=np.arange(len(labels)); ax[1].bar(xx-.18,rf,.36,label='Structure-only RF',color='#246d98'); ax[1].bar(xx+.18,med,.36,label='Training median',color='#92979e')
    ax[1].set(xticks=xx,xticklabels=labels,ylabel='Out-of-fold MAE (log10 nM)',title='Honest historical holdouts');ax[1].legend(frameon=False)
    fig.suptitle('Historical biochemical assay diagnostic — no cellular engagement model',fontsize=12)
    fig.savefig(out/'figures/model_diagnostic.png',dpi=190); fig.savefig(out/'figures/model_diagnostic.svg'); plt.close(fig)
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
