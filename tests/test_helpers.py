"""不依赖 pytest 的回归测试：python -m unittest -v test_helpers。

覆盖解析基准、复数、数组形状、错误输入、求解失败、绘图文件与风格恢复。
使用临时目录测试写文件，不污染真实笔记。
"""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from types import SimpleNamespace
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import helpers as h


class HelpersTests(unittest.TestCase):
    def setUp(self):
        self.sx = np.array([[0, 1], [1, 0]], complex)
        self.sy = np.array([[0, -1j], [1j, 0]], complex)
        self.sz = np.diag([1., -1.])

    def test_dagger(self):
        a = np.array([[1j, 2, 3], [4, 5j, 6]])
        np.testing.assert_allclose(h.dagger(a), a.conj().T)
        with self.assertRaises(ValueError): h.dagger([1, 2])

    def test_commutator(self):
        np.testing.assert_allclose(h.commutator(self.sx, self.sy), 2j*self.sz)
        with self.assertRaises(ValueError): h.commutator(self.sx, np.eye(3))

    def test_hermitian(self):
        self.assertTrue(h.is_hermitian(self.sy))
        self.assertFalse(h.is_hermitian([[0, 1], [0, 0]]))
        with self.assertRaises(ValueError): h.is_hermitian(self.sx, atol=-1)

    def test_discrete_normalization(self):
        original = np.array([1, 1j])
        before = original.copy()
        state = h.normalize_state(original)
        self.assertAlmostEqual(float(np.vdot(state, state).real), 1)
        np.testing.assert_array_equal(original, before)
        for scale in (1e-300, 1e300):
            np.testing.assert_allclose(h.normalize_state(original*scale), state)
        for bad in ([0, 0], [[1], [0]], [np.nan, 1], []):
            with self.subTest(bad=bad), self.assertRaises(ValueError): h.normalize_state(bad)

    def test_integral_and_broadcast(self):
        x = np.linspace(0, 1, 101)
        result = h.integrate_samples(np.stack([x**2, 1j*x**2]), x)
        np.testing.assert_allclose(result, [1/3, 1j/3], atol=1e-12)
        self.assertAlmostEqual(h.integrate_samples(x, x, method="trapezoid").real, 0.5)
        for xbad in ([0, 0, 1], [1, 0], [0, 1j], [0]):
            with self.assertRaises(ValueError): h.integrate_samples([1, 2], xbad)
        with self.assertRaises(ValueError): h.integrate_samples(x, x, method="bad")
        with self.assertRaises(ValueError): h.integrate_samples(x, x, axis=2)
        with self.assertRaises(ValueError): h.integrate_samples(x[:-1], x)

    def test_continuous_normalization(self):
        x = np.linspace(-8, 8, 1601)
        psi = h.normalize_wavefunction(np.exp(-x*x/2), x)
        self.assertAlmostEqual(h.integrate_samples(abs(psi)**2, x).real, 1)
        self.assertAlmostEqual(h.integrate_samples(x*x*abs(psi)**2, x).real, 0.5)
        # 不规则网格上的 Simpson 可给出负“概率”；明确拒绝，不静默改规则。
        with self.assertRaises(ValueError): h.normalize_wavefunction([1, 0, 0], [0, 1, 4])
        trap = h.normalize_wavefunction([1, 0, 0], [0, 1, 4], method="trapezoid")
        self.assertAlmostEqual(h.integrate_samples(abs(trap)**2, [0, 1, 4], method="trapezoid").real, 1)

    def test_expectation(self):
        self.assertAlmostEqual(h.expectation(self.sy, [3, 3j]), 1)
        self.assertAlmostEqual(h.expectation(1j*np.eye(2), [1, 0]), 1j)
        with self.assertRaises(ValueError): h.expectation(self.sx, [1, 2, 3])

    def test_eigensystem(self):
        values, vectors = h.diagonalize_hermitian(self.sy)
        np.testing.assert_allclose(values, [-1, 1])
        np.testing.assert_allclose(self.sy@vectors, vectors*values[None, :], atol=1e-14)
        np.testing.assert_allclose(h.dagger(vectors)@vectors, np.eye(2), atol=1e-14)
        with self.assertRaises(ValueError): h.diagonalize_hermitian([[0, 1], [0, 0]])

    def test_matrix_exponential(self):
        # 非厄米幂零矩阵也支持，exp(A)=I+A。
        a = np.array([[0, 2], [0, 0]])
        np.testing.assert_allclose(h.matrix_exp(a), np.eye(2)+a)
        np.testing.assert_allclose(h.matrix_exp(-0.3j*self.sx),
                                   np.cos(0.3)*np.eye(2)-1j*np.sin(0.3)*self.sx)
        for bad in ([1, 2], [[1, 2, 3]], [[np.inf]]):
            with self.assertRaises(ValueError): h.matrix_exp(bad)

    def test_time_origin_hbar_and_shape(self):
        t = np.linspace(5, 8, 51)
        states = h.evolve_constant(self.sx, [2, 0], t, hbar=2)
        analytical = 2*np.column_stack([np.cos((t-5)/2), -1j*np.sin((t-5)/2)])
        self.assertEqual(states.shape, (51, 2))
        np.testing.assert_allclose(states, analytical, atol=1e-14)
        np.testing.assert_allclose(h.evolve_constant(self.sx, [2, 0], [5]), [[2, 0]])
        numeric = h.solve_schrodinger(self.sx, [2, 0], t, hbar=2)
        np.testing.assert_allclose(numeric, analytical, atol=1e-8)

    def test_extreme_hbar_scaling(self):
        states = h.evolve_constant([[1e308]], [1], [0, 2], hbar=1e308)
        np.testing.assert_allclose(states[:, 0], [1, np.exp(-2j)], atol=1e-14)
        states = h.solve_schrodinger([[1e-320]], [1], [0, 1], hbar=1e-320)
        np.testing.assert_allclose(states[:, 0], [1, np.exp(-1j)], atol=1e-8)

    def test_time_dependent(self):
        t = np.linspace(0, 4, 81)
        numeric = h.solve_schrodinger(lambda s: (1+0.2*np.cos(s))*self.sx/2, [1, 0], t)
        theta = t+0.2*np.sin(t)
        np.testing.assert_allclose(numeric, np.column_stack([np.cos(theta/2), -1j*np.sin(theta/2)]), atol=1e-8)
        np.testing.assert_allclose(np.sum(abs(numeric)**2, axis=1), 1, atol=1e-8)

    def test_solver_errors(self):
        for kwargs in ({"hbar": 0}, {"rtol": -1}, {"method": "Radau"}):
            with self.assertRaises(ValueError): h.solve_schrodinger(self.sx, [1, 0], [0, 1], **kwargs)
        with self.assertRaises(ValueError): h.solve_schrodinger(self.sx, [1, 0], [0])
        with self.assertRaises(ValueError): h.solve_schrodinger(lambda t: [[0, 1], [0, 0]], [1, 0], [0, 1])
        with self.assertRaises(ValueError): h.evolve_constant(self.sx, [0, 0], [0, 1])
        with patch("helpers.solve_ivp", return_value=SimpleNamespace(success=False, message="test failure")):
            with self.assertRaises(RuntimeError): h.solve_schrodinger(self.sx, [1, 0], [0, 1])

    def test_subnormal_and_large_complex_states(self):
        for scale in (1e-320, 1.7e308):
            state = np.array([complex(scale, scale), 0])
            normalized = h.normalize_state(state)
            np.testing.assert_allclose(normalized, [(1+1j)/np.sqrt(2), 0])
            wave = h.normalize_wavefunction(state, [0, 1], method="trapezoid")
            self.assertAlmostEqual(h.integrate_samples(abs(wave)**2, [0, 1], method="trapezoid").real, 1)

    def test_scale_independent_hermitian_validation(self):
        # 公共检查保留用户绝对容差语义；物理求解不让微小单位掩盖错误。
        for scale in (1e-13, 1, 1e13):
            bad = np.array([[0, scale], [0, 0]])
            with self.assertRaises(ValueError): h.diagonalize_hermitian(bad)
            with self.assertRaises(ValueError): h.evolve_constant(bad, [1, 0], [0, 1])
            with self.assertRaises(ValueError): h.solve_schrodinger(bad, [1, 0], [0, 1])
        self.assertTrue(h.is_hermitian([[0, 1e-13], [0, 0]]))

    def test_narrow_pulse_and_max_step(self):
        # 高斯窄脉冲面积约为 angle；H(t)=pulse(t)*sigma_x/2。
        width, center, angle = 0.005, 0.537, 0.8
        def pulse(t):
            amplitude = 0.1 + angle/(np.sqrt(2*np.pi)*width)*np.exp(-0.5*((t-center)/width)**2)
            return amplitude*self.sx/2
        states = h.solve_schrodinger(pulse, [1, 0], [0, 1], max_step=width/5)
        total_angle = angle + 0.1  # 平滑小背景避免极小高斯尾部触发求解器下溢。
        np.testing.assert_allclose(states[-1], [np.cos(total_angle/2), -1j*np.sin(total_angle/2)], atol=1e-8)
        for bad in (0, -1, -np.inf, np.nan, True, "1", [1]):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                h.solve_schrodinger(self.sx, [1, 0], [0, 1], max_step=bad)

    def test_invalid_scalar_and_selector_types(self):
        for bad in (None, True, "1", 1j, [1], [[1], [2, 3]]):
            with self.subTest(bad=bad), self.assertRaises(ValueError): h.plot_style(dpi=bad)
        for bad in (None, np.array(1), [1], [[1], [2]]):
            with self.subTest(bad=bad), self.assertRaises(ValueError): h.new_figure(figsize=bad)
        for bad in ([], np.array(["simpson", "trapezoid"]), None):
            with self.assertRaises(ValueError): h.integrate_samples([1, 2], [0, 1], method=bad)
        with self.assertRaises(ValueError): h.integrate_samples([1, 2], [0, 1], axis=True)

    def test_floating_point_range_errors(self):
        with self.assertRaises(ValueError): h.matrix_exp([[1000]])
        with self.assertRaises(ValueError): h.evolve_constant([[1e308]], [1], [0, 1], hbar=1e-308)
        with self.assertRaises(ValueError): h.integrate_samples([1, 1], [-1e308, 1e308])
        # 仅一个时刻时不需要计算不可表示的 H/hbar。
        np.testing.assert_array_equal(h.evolve_constant([[1e308]], [2], [1], hbar=1e-308), [[2]])

    def test_atomic_save_failure_preserves_old_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"original.png"
            path.write_bytes(b"old-content")
            def failed_save(temp, **kwargs):
                Path(temp).write_bytes(b"incomplete")
                raise RuntimeError("simulated failure")
            fig = SimpleNamespace(savefig=failed_save)
            with self.assertRaises(RuntimeError): h.save_figure(fig, path, overwrite=True)
            self.assertEqual(path.read_bytes(), b"old-content")
            self.assertEqual(list(Path(tmp).iterdir()), [path])

    def test_atomic_save_race_does_not_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"race.png"
            def competing_save(temp, **kwargs):
                Path(temp).write_bytes(b"new-image")
                path.write_bytes(b"other-process")  # 模拟检查之后别的进程创建目标。
            fig = SimpleNamespace(savefig=competing_save)
            with self.assertRaises(FileExistsError): h.save_figure(fig, path)
            self.assertEqual(path.read_bytes(), b"other-process")
            self.assertEqual(list(Path(tmp).iterdir()), [path])

    def test_plots_and_style_restore(self):
        before = matplotlib.rcParams["font.size"]
        with tempfile.TemporaryDirectory() as tmp:
            with h.plot_style(font_size=17):
                self.assertEqual(matplotlib.rcParams["font.size"], 17)
                fig, ax = h.new_figure(xlabel="t", ylabel="p", title="test")
                ax.plot([0, 1], [0, 1])
                for suffix in ("png", "pdf", "svg"):
                    target = Path(tmp)/"nested"/f"test.{suffix}"
                    self.assertEqual(h.save_figure(fig, target), target)
                    self.assertGreater(target.stat().st_size, 0)
                    with self.assertRaises(FileExistsError): h.save_figure(fig, target)
                    h.save_figure(fig, target, overwrite=True)
                with self.assertRaises(ValueError): h.save_figure(fig, Path(tmp)/"bad.txt")
                plt.close(fig)
        self.assertEqual(matplotlib.rcParams["font.size"], before)


if __name__ == "__main__":
    unittest.main(verbosity=2)
