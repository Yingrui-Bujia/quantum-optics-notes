"""有限差分解一维谐振子定态，并量化前三态、网格与边界误差。"""
import argparse
from math import factorial
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.linalg import eigh_tridiagonal
from scipy.special import eval_hermite
from helpers import normalize_wavefunction, integrate_samples, plot_style
from workflow import run_directory, save_plot, save_json, require, max_error, real_parameter

def analytic_state(n, x):
    # 无量纲 x=x_phys/sqrt(hbar/(m*omega))；能量单位 hbar*omega。
    return eval_hermite(n, x)*np.exp(-x*x/2)/(np.pi**0.25*np.sqrt(2**n*factorial(n)))

def solve_stationary(intervals=800, length=8.0, count=3):
    if isinstance(intervals, bool) or not isinstance(intervals, (int, np.integer)) or intervals < 10:
        raise ValueError('intervals 必须是至少10的整数')
    length=real_parameter(length,'length',positive=True)
    if not isinstance(count, (int, np.integer)) or isinstance(count, bool) or not 1 <= count < intervals:
        raise ValueError('count 必须是1到intervals-1的整数')
    x = np.linspace(-length, length, intervals+1)
    dx = x[1]-x[0]
    # 端点psi=0；仅在内部格点上求解。二阶中心差分形成三对角矩阵。
    diagonal = 1/dx**2 + x[1:-1]**2/2
    off_diagonal = np.full(intervals-2, -1/(2*dx**2))
    energies, vectors = eigh_tridiagonal(diagonal, off_diagonal, select='i', select_range=(0, count-1))
    states, references = [], []
    for n in range(count):
        psi = normalize_wavefunction(np.pad(vectors[:, n], (1, 1)), x).real
        ref = analytic_state(n, x)
        # 本征矢允许整体负号；按重叠对齐，不能误判为物理差异。
        if integrate_samples(psi*ref, x).real < 0: psi = -psi
        states.append(psi); references.append(ref)
    return x, energies, np.asarray(states), np.asarray(references)

def calculate():
    x, energies, states, ref = solve_stationary()
    exact = np.arange(3)+0.5
    grid = []
    for intervals in (200, 400, 800, 1600):
        _, e, _, _ = solve_stationary(intervals)
        grid.append({'intervals': intervals, 'dx': 16/intervals, 'max_energy_error': max_error(e, exact)})
    # 只改盒长且固定 dx=0.02，分离有限边界与差分步长的影响。
    boundary = []
    for length in (3., 4., 6., 8.):
        _, e, _, _ = solve_stationary(round(2*length/0.02), length)
        boundary.append({'length': length, 'max_energy_error': max_error(e, exact)})
    metrics = {'energies': energies.tolist(), 'analytic_energies': exact.tolist(),
               'max_energy_error': max_error(energies, exact),
               'wavefunction_L2_errors': [float(np.sqrt(integrate_samples((p-r)**2, x).real)) for p,r in zip(states,ref)],
               'norm_errors': [abs(float(integrate_samples(p*p,x).real)-1) for p in states],
               'ground_x2': float(integrate_samples(x*x*states[0]**2, x).real),
               'grid_convergence': grid, 'boundary_convergence': boundary,
               'refinement_ratios': [grid[i]['max_energy_error']/grid[i+1]['max_energy_error'] for i in range(3)]}
    require(metrics['max_energy_error'] < 2e-4, '前三态能量误差超过2e-4')
    require(max(metrics['wavefunction_L2_errors']) < 2e-4, '波函数L2误差过大')
    require(all(3.8 < r < 4.2 for r in metrics['refinement_ratios']), '未观察到二阶收敛')
    return x, states, ref, metrics

def main(output=None):
    folder=run_directory('week01-oscillator', output)
    x, states, ref, metrics=calculate()
    with plot_style():
        fig, axes=plt.subplots(3,1,figsize=(8,9),layout='constrained',sharex=True)
        for n, ax in enumerate(axes):
            ax.plot(x,ref[n],label='Analytical')
            ax.plot(x[::15],states[n,::15],'o',ms=3,label='Finite difference')
            ax.set(ylabel=f'psi_{n}(x)',title=f'n={n}; E={metrics["energies"][n]:.7f} (hbar omega)')
            ax.legend()
        axes[-1].set_xlabel('x / sqrt(hbar / (m omega))')
        save_plot(fig,folder,'eigenfunctions');plt.close(fig)
        fig,ax=plt.subplots(layout='constrained')
        g=metrics['grid_convergence'];dx=np.array([r['dx'] for r in g]);err=np.array([r['max_energy_error'] for r in g])
        ax.loglog(dx,err,'o-',label='Maximum n=0,1,2 error')
        ax.loglog(dx,err[-1]*(dx/dx[-1])**2,'--',label='Second-order reference')
        ax.set(xlabel='Grid spacing (oscillator units)',ylabel='Energy error (hbar omega)',title='Grid convergence at L=8')
        ax.legend();save_plot(fig,folder,'convergence');plt.close(fig)
    np.savetxt(folder/'wavefunctions.csv',np.column_stack([x,states.T,ref.T]),delimiter=',',header='x,psi0,psi1,psi2,exact0,exact1,exact2',comments='')
    save_json(folder,'metrics.json',metrics)
    return folder,metrics

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output')
    args=parser.parse_args();folder,metrics=main(args.output);print(metrics);print(folder)
