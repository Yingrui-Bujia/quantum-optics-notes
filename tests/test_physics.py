"""独立物理基准与容易被曲线重合掩盖的边界情况。"""
import unittest
import numpy as np
import qutip as qt
from week01_harmonic_oscillator import solve_stationary, analytic_state
from week01_matrix_exponential import calculate as matrix_calculate
from week02_qubit_rabi import rabi, zero_point, thermal, tensor_demo, OPTIONS
from week02_tutorials import exact_dephasing_z, tutorial001, tutorial002, tutorial003

class PhysicsTests(unittest.TestCase):
    def test_stationary_energies_parity_and_nodes(self):
        x,e,psi,ref=solve_stationary(800)
        np.testing.assert_allclose(e,np.arange(3)+0.5,atol=2e-4,rtol=0)
        for n in range(3):
            np.testing.assert_allclose(psi[n],(-1)**n*psi[n,::-1],atol=1e-10)
            central=psi[n,np.abs(x)<3]
            significant=central[np.abs(central)>1e-7]
            self.assertEqual(np.sum(significant[1:]*significant[:-1]<0),n)

    def test_second_order_grid(self):
        errors=[]
        for k in (200,400,800):
            _,e,_,_=solve_stationary(k)
            errors.append(np.max(abs(e-(np.arange(3)+0.5))))
        for a,b in zip(errors,errors[1:]):self.assertTrue(3.8<a/b<4.2)

    def test_expm_ode(self):
        _,_,_,metrics=matrix_calculate()
        self.assertLess(metrics['ode_vs_expm'],1e-8)

    def test_detuned_rabi_amplitude_and_pi_pulse(self):
        t,result,reference,metrics=rabi(detuning=2.,samples=1001)
        self.assertLess(metrics['population_error'],1e-7)
        self.assertLessEqual(np.max(result.expect[1]),1/5+1e-7)
        t,result,reference,metrics=rabi(samples=601)
        self.assertAlmostEqual(result.expect[1][150],1,places=7)

    def test_finite_commutator_boundary(self):
        data=zero_point(12,5)
        self.assertAlmostEqual(data['vacuum_energy'],0.5)
        self.assertAlmostEqual(data['vacuum_variance_product'],0.25)
        self.assertAlmostEqual(data['top_state_actual'],11)
        with self.assertRaises(ValueError):zero_point(6,5)

    def test_thermal_small_cutoff_is_not_infinite_distribution(self):
        _,p,infinite,d=thermal(10,3.)
        self.assertGreater(d['omitted_probability'],0.05)
        self.assertGreater(d['mean_number_error'],0.5)
        self.assertAlmostEqual(np.sum(p),1)
        self.assertLess(np.sum(infinite),1)

    def test_thermal_vacuum_and_convergence(self):
        _,p,_,d=thermal(10,0.)
        np.testing.assert_allclose(p,[1]+[0]*9)
        self.assertEqual(d['mean_number_error'],0)
        _,_,_,large=thermal(100,3.)
        self.assertLess(large['mean_number_error'],1e-8)

    def test_tensor_dimensions_and_subsystem(self):
        d=tensor_demo()
        self.assertEqual(d['dims'],[[6,2],[1]])
        self.assertAlmostEqual(d['mode_number'],2)
        self.assertEqual(d['qubit_sigma_z'],-1)
        self.assertLess(d['partial_trace_error'],1e-12)

    def test_tutorial001_002(self):
        tutorial001();tutorial002()

    def test_tutorial003_physicality(self):
        tutorial003()

    def test_critical_and_overdamped_dephasing(self):
        t=np.linspace(0,4,301)
        for delta,g in ((1.,1.),(1.,2.),(1.,0.),(0.,0.)):
            res=qt.mesolve(delta*qt.sigmax()/2,qt.basis(2,0).proj(),t,
                           [np.sqrt(g)*qt.sigmaz()],e_ops=[qt.sigmaz()],options=OPTIONS)
            np.testing.assert_allclose(res.expect[0],exact_dephasing_z(t,delta,g),atol=1e-8,rtol=0)

    def test_invalid_physical_parameters(self):
        for bad in (0,-1,np.nan,np.inf,True,'1',[1]):
            with self.subTest(omega=bad),self.assertRaises(ValueError):rabi(omega=bad)
        for bad in (-1,np.nan,np.inf,True,'1',1j):
            with self.subTest(nbar=bad),self.assertRaises(ValueError):thermal(nbar=bad)
        for bad in (0,-1,np.nan,True,'8',1j):
            with self.subTest(length=bad),self.assertRaises(ValueError):solve_stationary(length=bad)

if __name__=='__main__':unittest.main(verbosity=2)
