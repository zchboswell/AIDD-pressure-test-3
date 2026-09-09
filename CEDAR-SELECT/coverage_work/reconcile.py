#!/usr/bin/env python3
"""Reconcile supplied inert scenario records; never infer numerical biology from status."""
import json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / 'inputs' / 'coverage'
def read(name):
    return json.loads((INPUT / name).read_text())
def unique(rows, key, fields):
    out = {}
    for row in rows:
        ident = row[key]
        normalized = {f: row[f] for f in fields}
        if ident in out:
            assert out[ident] == normalized, f'Conflicting duplicate {ident}'
        out[ident] = normalized
    return out
w = read('work_orders.json')
a = read('analogue_evidence.json')
b = read('entity_bindings.json')
r = read('retrieval_archive.json')
financial = unique(w['financial_records'], 'obligation_id', ['state','amount'])
stock = unique(w['inventory'], 'vial_id', ['compound_slot','usable_mg','reserved_mg','available_business_day'])
paid = sum(x['amount'] for x in financial.values() if x['state']=='paid')
reserved = sum(x['amount'] for x in financial.values() if x['state']=='reserved')
selected = ['O2','O1']
options = {x['option_id']:x for x in w['new_options']}
new_cost = sum(options[x]['incremental_usd'] for x in selected)
left = w['budget']['earmarked_within_tranche']-paid-reserved-new_cost
assert (paid,reserved,new_cost,left)==(1800,2200,3400,2600)
schedule = [dict(work='W2/CDK1/context-B',start=0,end=3,lane=1,material_mg=2,vial_id='V2',cost_already_reserved=True),dict(work='O2',start=0,end=3,lane=2,material_mg=2,vial_id='V3'),dict(work='O1',start=3,end=5,lane=1,material_mg=1,vial_id='V1')]
for day in range(w['contract']['completion_deadline_business_day']):
    active=[x for x in schedule if x['start']<=day<x['end']]
    assert len(active)<=w['contract']['parallel_capacity']
    assert len({x['lane'] for x in active})==len(active)
for x in schedule:
    assert x['start']>=stock[x['vial_id']]['available_business_day']
    assert x['end']<=w['contract']['completion_deadline_business_day']
    if x['work'] in options:
        assert x['end']-x['start']==options[x['work']]['duration_business_days']
remaining_stock={}
for vial,s in stock.items():
    new_used=sum(options[o]['material_mg'] for o in selected if options[o]['vial_id']==vial)
    assert new_used<=s['usable_mg']-s['reserved_mg']
    planned=sum(x['material_mg'] for x in schedule if x['vial_id']==vial)
    remaining_stock[vial]=dict(total_at_cutoff_mg=s['usable_mg'],reserved_at_cutoff_mg=s['reserved_mg'],free_at_cutoff_mg=s['usable_mg']-s['reserved_mg'],after_planned_work_mg=s['usable_mg']-planned)
    assert remaining_stock[vial]['after_planned_work_mg']>=0
assays={}
for o in a['observations']:
    assert o['record_kind']=='public_measurement' and o['original_relation']=='=' and o['original_unit']=='uM'
    assays.setdefault(o['compound_id'],{})[o['target']]=float(o['original_value'])*1000
slotmap={x['compound_id']:x['slot'] for x in b['bindings']}
analogues=[]
for c in a['compounds']:
    val=assays[c['compound_id']]
    analogues.append(dict(compound_id=c['compound_id'],slot=slotmap.get(c['compound_id']),example=c['public_source_identifiers']['example'],CDK2_nM=val['CDK2'],CDK1_nM=val['CDK1'],ratio_CDK1_over_CDK2=val['CDK1']/val['CDK2'],scope='historical public paired assay; no uncertainty supplied'))
classifications={
'R1':('valid_index_hit','Scope-matched complete fictional index; no primary numerical measurement. Retrieve primary activity and context.'),
'R2':('unsupported_filter_wrong_target','Executed target key unsupported; P06493 is CDK1, not requested P24941. Repeat with target_accession and verify returned target.'),
'R3':('wrong_entity_namespace','Executed assay:771 differs from requested compound:771. Repeat compound namespace; numeric suffix is not identity.'),
'R4':('incomplete_pagination','Page 2 missing; one of two advertised records. Retrieve and reconcile all pages before coverage claim.'),
'R5':('valid_scoped_empty','Complete zero within executed compound:999 index query only. Not inactivity or universal absence.'),
'R6':('service_error_unknown','HTTP503 records null is failed retrieval; retry or check an authorized independent source. Not a negative result.'),
'R7':('form_index_empty','Mechanism operation valid, but annotations indexed on parent; salt-form emptiness is not absence of parent mechanism.'),
'R8':('parent_annotation_only','Linked parent has simulated mechanism annotation; no potency or salt-to-parent assay transfer.'),
'R9':('endpoint_mismatch','Returned percent inhibition at10uM fails requested IC50. Retrieve a concentration-response IC50; no conversion from this row.')}
assert set(classifications)=={x['id'] for x in r['responses']}
result=dict(schema='cedar-coverage-reconciliation-v1',evidence_separation='All response and logistics records are simulated. Only analogue_evidence observations supply historical numerical biology.',retrieval=[dict(id=k,status=v[0],interpretation_and_next_action=v[1],eligible_numerical_IC50=False) for k,v in classifications.items()],financial=dict(unique_obligations=financial,paid_usd=paid,reserved_usd=reserved,free_before_new_usd=6000,selected_options=selected,new_proposed_usd=new_cost,total_paid_reserved_proposed_usd=paid+reserved+new_cost,unallocated_usd=left,earmarked_usd=10000,parent_tranche_usd=120000,parent_unearmarked_usd=110000,added_to_parent=False),option_decisions=dict(O1='Propose: deliberate independent preparation/day reproducibility control; not an accidental duplicate.',O2='Propose: recover QC-failed S3 CDK1 endpoint; absence of eligible value does not imply inactivity.',O3='Do not propose: overlaps authorized pending S2 endpoint; wait through day3 and inspect QC before considering a separately authorized repeat.',O4='Defer: null all-in cost/duration, unverified catalogue quantity and lead time; route proposal is not experimental execution.'),endpoint_state=dict(W1='completed-valid attestation only; numerical CDK2 value absent',W2_CDK2='completed-valid attestation only; numerical CDK2 value absent',W2_CDK1='pending due day3; preserve reservation',W3_CDK1='failed QC; no eligible potency estimate',W4='intentional repeat proposed, not yet authorized'),schedule_assumption='Conservatively occupy one lane for pending W2 from day0 through day3; new orders start only after proposed authorization. Scenario timing is not provider verification.',schedule=schedule,material=remaining_stock,material_claims=dict(M1='Unexecuted two-step retrosynthetic proposal for S4; no procedure or yield.',M2='Listing only; stock, price and lead time unknown.',M3='Simulated received S3 4mg is V3 inventory already counted, cost in OB-17.',M4='Simulated historical S1 preparation with41% yield and identity attestations; not current stock, S4 feasibility or potency.'),public_analogue_summary=analogues,biological_decision='Prioritize S3/S2 comparative validation over acquiring S4 solely for potency. S3 versus S2 ratio difference171vs162 is small and unreplicated; S4 ratio146 despite improved potency. Test new matched paired contexts before cellular-window claims.')
out=ROOT/'delivery'/'coverage_reconciliation.json'
out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(checks='PASS',output=str(out.relative_to(ROOT)),total_usd=paid+reserved+new_cost,remaining_usd=left,last_completion_day=max(x['end'] for x in schedule),unique_vials=len(stock))))
