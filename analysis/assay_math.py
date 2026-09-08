"""Validated concentration conversion; derived repair of inputs/legacy/assay_math.py.

Concentrations and their logarithms preserve censoring as thresholds, not point labels.
No affinity/functional endpoint conversion, mass-unit conversion, or assay correction.
"""
import math

_FACTORS = {'nM': 1.0, 'uM': 1e3, 'µM': 1e3, 'μM': 1e3, 'mM': 1e6, 'M': 1e9}
_REVERSE = {'=': '=', '<': '>', '<=': '>=', '>': '<', '>=': '<='}

def normalize(value, unit, relation='=', endpoint='IC50'):
    """Return concentration in nM and p=-log10(molar), or explicit invalid status."""
    if not isinstance(endpoint, str) or endpoint not in ('IC50', 'Ki', 'Kd'):
        return {'status': 'invalid', 'reason': 'unsupported endpoint'}
    if not isinstance(unit, str) or unit not in _FACTORS:
        return {'status': 'invalid', 'reason': 'unsupported concentration unit'}
    if not isinstance(relation, str):
        return {'status': 'invalid', 'reason': 'unsupported relation'}
    relation = {'≤': '<=', '≥': '>='}.get(relation, relation)
    if relation not in _REVERSE:
        return {'status': 'invalid', 'reason': 'unsupported relation'}
    if isinstance(value, bool):
        return {'status': 'invalid', 'reason': 'boolean is not a concentration'}
    try:
        value_nm = float(value) * _FACTORS[unit]
    except (ValueError, TypeError, OverflowError):
        return {'status': 'invalid', 'reason': 'missing or nonnumeric concentration'}
    if not math.isfinite(value_nm) or value_nm <= 0:
        return {'status': 'invalid', 'reason': 'concentration must be positive and finite'}
    return {'status': 'ok', 'endpoint': endpoint, 'value_nm': value_nm,
            'relation_nm': relation, 'p_activity': 9.0 - math.log10(value_nm),
            'p_relation': _REVERSE[relation]}
