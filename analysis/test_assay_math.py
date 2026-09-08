"""Meaningful regression checks for the supplied converter's failure modes."""
import unittest
from assay_math import normalize

class AssayMathTests(unittest.TestCase):
    def test_equivalent_units(self):
        for value, unit in [(100,'nM'),(.1,'uM'),(.1,'µM'),(.1,'μM'),(.0001,'mM'),(1e-7,'M')]:
            r=normalize(value,unit); self.assertEqual(r['status'],'ok')
            self.assertAlmostEqual(r['value_nm'],100); self.assertAlmostEqual(r['p_activity'],7)
    def test_inequality_reversal(self):
        for src,dst in [('=','='),('<','>'),('<=','>='),('>','<'),('>=','<='),('≤','>='),('≥','<=')]:
            r=normalize(.01,'µM',src); self.assertEqual(r['p_relation'],dst);self.assertEqual(r['p_activity'],8)
    def test_endpoint_preservation(self):
        for endpoint in ['IC50','Ki','Kd']:
            self.assertEqual(normalize(.3,'uM','>=',endpoint)['endpoint'],endpoint)
        self.assertEqual(normalize(1,'nM',endpoint='EC50')['status'],'invalid')
    def test_invalid(self):
        for value in [None,'', 'bad',0,-3,float('nan'),float('inf'),'-inf',True,1e309]:
            self.assertEqual(normalize(value,'nM')['status'],'invalid',repr(value))
        for unit in ['ng/mL','NM','',None,[]]: self.assertEqual(normalize(20,unit)['status'],'invalid')
        for relation in ['~','',None,[]]: self.assertEqual(normalize(20,'nM',relation)['status'],'invalid')
    def test_no_clipping(self):
        self.assertAlmostEqual(normalize(1e-15,'nM')['p_activity'],24)
        self.assertEqual(normalize(1e308,'M')['status'],'invalid')

if __name__=='__main__': unittest.main()
