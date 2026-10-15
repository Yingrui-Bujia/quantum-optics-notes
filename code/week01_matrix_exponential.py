"""一小时小测的参考解：矩阵指数、solve_ivp、解析解三方比较。"""
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from helpers import evolve_constant, solve_schrodinger, matrix_exp, plot_style
from workflow import run_directory,save_plot,save_json,require,max_error

def calculate():
    sx=np.array([[0,1],[1,0]],complex)
    h=sx/2;t=np.linspace(0,4*np.pi,501);initial=[1,0]
    exact=evolve_constant(h,initial,t)
    ode=solve_schrodinger(h,initial,t,rtol=1e-10,atol=1e-12)
    analytic=np.column_stack([np.cos(t/2),-1j*np.sin(t/2)])
    unitary=matrix_exp(-1j*h*t[-1])
    metrics={'expm_vs_analytic':max_error(exact,analytic), 'ode_vs_expm':max_error(ode,exact),
             'unitarity_error':max_error(unitary.conj().T@unitary,np.eye(2)),
             'ode_norm_drift':max_error(np.sum(abs(ode)**2,axis=1),1),
             'elementwise_exp_is_different':bool(not np.allclose(np.exp(-1j*h),matrix_exp(-1j*h)))}
    require(max(metrics[k] for k in ('expm_vs_analytic','ode_vs_expm','unitarity_error','ode_norm_drift'))<1e-8,'矩阵指数验收失败')
    return t,exact,ode,metrics

def main(output=None):
    folder=run_directory('week01-expm',output);t,exact,ode,metrics=calculate()
    with plot_style():
        fig,ax=plt.subplots(layout='constrained')
        ax.plot(t,np.sin(t/2)**2,label='sin^2(t/2)')
        ax.plot(t[::12],abs(exact[::12,1])**2,'o',ms=3,label='expm')
        ax.plot(t,abs(ode[:,1])**2,'--',label='solve_ivp')
        ax.set(xlabel='Time (1/Omega)',ylabel='P(1)',title='Matrix exponential examination reference')
        ax.legend();save_plot(fig,folder,'matrix_exponential');plt.close(fig)
    np.savetxt(folder/'populations.csv',np.column_stack([t,abs(exact[:,1])**2,abs(ode[:,1])**2]),delimiter=',',header='t,expm,solve_ivp',comments='')
    save_json(folder,'metrics.json',metrics);return folder,metrics

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output')
    args=parser.parse_args();folder,metrics=main(args.output);print(metrics);print(folder)
