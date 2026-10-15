"""第一周第1–4次练习：语法、广播、线性代数、积分与经典谐振子。"""
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from helpers import diagonalize_hermitian, integrate_samples, plot_style
from workflow import run_directory, save_plot, save_json, require, max_error

def square(x):
    return x**2

def calculate():
    squares = [square(i) for i in range(1, 11)]
    array = np.arange(1, 7).reshape(2, 3)
    # (2,3) + (3,)：最后一维逐列相加；@ 才是矩阵乘法。
    broadcast = array + np.array([10, 20, 30])
    values, vectors = diagonalize_hermitian([[0, 1], [1, 0]])
    x = np.linspace(0, 1, 101)
    integral = integrate_samples(x*x, x).real
    t = np.linspace(0, 4*np.pi, 401)
    # y[0]=位置，y[1]=速度；m=omega=1，x''=-x。
    result = solve_ivp(lambda t, y: [y[1], -y[0]], (t[0], t[-1]), [1, 0],
                       t_eval=t, method='DOP853', rtol=1e-10, atol=1e-12)
    require(result.success, result.message)
    energy = (result.y[0]**2 + result.y[1]**2)/2
    metrics = {'squares_1_to_10': squares, 'broadcast': broadcast.tolist(),
               'sigma_x_eigenvalues': values.tolist(), 'integral_x2_0_1': float(integral),
               'classical_position_max_error': max_error(result.y[0], np.cos(t)),
               'classical_energy_max_drift': max_error(energy, 0.5)}
    require(metrics['classical_position_max_error'] < 1e-8, '经典谐振子位置误差过大')
    require(abs(integral-1/3) < 1e-12, '积分基准失败')
    return t, result.y, metrics

def main(output=None):
    folder = run_directory('week01-basics', output)
    t, y, metrics = calculate()
    with plot_style():
        fig, ax = plt.subplots(layout='constrained')
        ax.plot(t, np.sin(t), label='sin(t)')
        ax.set(xlabel='Time (dimensionless)', ylabel='Amplitude', title='A labelled sine curve')
        ax.legend(); save_plot(fig, folder, 'sine'); plt.close(fig)
        fig, ax = plt.subplots(layout='constrained')
        ax.plot(t, np.cos(t), label='Analytical cos(t)')
        ax.plot(t[::10], y[0, ::10], 'o', ms=3, label='solve_ivp')
        ax.set(xlabel='Time (1/omega)', ylabel='Position (x0)', title='Classical oscillator')
        ax.legend(); save_plot(fig, folder, 'classical_oscillator'); plt.close(fig)
    np.savetxt(folder/'classical.csv', np.column_stack([t, y.T]), delimiter=',', header='t,x,v', comments='')
    save_json(folder, 'metrics.json', metrics)
    return folder, metrics

if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('--output')
    args=parser.parse_args(); folder, metrics=main(args.output); print(metrics); print(folder)
