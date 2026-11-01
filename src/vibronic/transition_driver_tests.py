"""Run with the configured scientific Python: python transition_driver_tests.py.
Optional integration: set BEAVYR_TEST_PYSCF_PYTHON / BEAVYR_TEST_PSI4_PYTHON.
The test process must be limited to four CPU cores, e.g. taskset -c 0-3.
"""
import sys
sys.dont_write_bytecode = True
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
import numpy as np
import transition_driver as driver

class StateTracking(unittest.TestCase):
    def reference(self):
        c=np.eye(4)
        return dict(orbitals=[(c[:,:2],c[:,2:])],vectors=[[np.array([[1.,0.],[0.,0.]])],[np.array([[0.,1.],[0.,0.]])]],dipoles=np.array([[1.,2.,3.],[2.,0.,0.]]))

    def test_root_permutation_and_orbital_and_electronic_phases(self):
        ref=self.reference();c=np.eye(4)
        o=np.diag([-1.,1.]);v=np.diag([1.,-1.])
        vectors=[o@x[0]@v for x in ref['vectors']]
        cur=dict(orbitals=[(c[:,:2]@o,c[:,2:]@v)],vectors=[[vectors[1]],[-vectors[0]]],dipoles=np.array([[2.,0.,0.],[-1.,-2.,-3.]]))
        dipole,root,overlap=driver.match(ref,cur,np.eye(4),0)
        self.assertEqual(root,1);self.assertAlmostEqual(overlap,-1.)
        np.testing.assert_allclose(dipole,[1.,2.,3.])

    def test_ambiguous_roots_are_rejected(self):
        ref=self.reference();cur=self.reference()
        x,y=[v[0] for v in cur['vectors']]
        cur['vectors']=[[(x+y)/np.sqrt(2)],[(x-y)/np.sqrt(2)]]
        with self.assertRaisesRegex(ValueError,'Ambiguous'):driver.match(ref,cur,np.eye(4),0)

    def test_singular_occupied_overlap_is_rejected(self):
        ref=self.reference()
        with self.assertRaisesRegex(ValueError,'Singular'):driver.match(ref,ref,np.zeros((4,4)),0)

class BackendIntegration(unittest.TestCase):
    def calculate(self,module,open_shell=False,tda=True):
        binary=os.environ.get('BEAVYR_TEST_'+module.upper()+'_PYTHON')
        if not binary:self.skipTest('Set BEAVYR_TEST_'+module.upper()+'_PYTHON for integration tests')
        with tempfile.TemporaryDirectory(prefix='beavyr-properties-test-')as folder:
            directory=Path(folder)
            (directory/'geometry.bohr').write_text('O 0 0 0\nH 0 0 1.8\n' if open_shell else 'O 0 0 0\nH 1.43 1.1 0\nH -1.43 1.1 0\n')
            env=os.environ.copy()
            env.update({f'BEAVYR_{k}':v for k,v in dict(PROPERTY_BACKEND=module,ROOT='1',ROOTS='2' if open_shell else '3',CHARGE='0',MULTIPLICITY='2' if open_shell else '1',BASIS='STO-3G',FUNCTIONAL='pbe',TDA=str(tda).lower(),HF='true',HT='true',STEP='0.001',CORES='1',MEMORY='512').items()})
            env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
            result=subprocess.run([binary,str(Path(__file__).with_name('transition_driver.py').resolve())],cwd=directory,env=env,capture_output=True,text=True,timeout=300)
            self.assertEqual(result.returncode,0,result.stdout[-3000:]+result.stderr[-3000:])
            self.assertIn('BEAVYR_TRANSITION_SUCCESS',result.stdout)
            count=2 if open_shell else 3
            data=np.loadtxt(directory/'transition.properties',skiprows=2+count)
            self.assertEqual(data.shape,(3*count+1,3));self.assertTrue(np.isfinite(data).all())
            tracking=np.loadtxt(directory/'tracking.csv',delimiter=',',skiprows=1)
            self.assertTrue(np.all(np.abs(tracking[:,3])>0.7))
            return data

    def test_pyscf_full_response(self):self.calculate('pyscf',tda=False)
    def test_psi4_full_response(self):self.calculate('psi4',tda=False)
    def test_pyscf_open_shell(self):self.calculate('pyscf',open_shell=True)
    def test_psi4_open_shell(self):self.calculate('psi4',open_shell=True)
    def test_closed_shell_derivative_agreement(self):
        a=self.calculate('pyscf');b=self.calculate('psi4')
        sign=1 if np.vdot(a,b)>=0 else -1
        self.assertLess(np.linalg.norm(a-sign*b)/np.linalg.norm(a),1.e-3)

if __name__=='__main__':unittest.main()
