"""官方Time Evolution 001–003的独立复现；版本与出处见docs/sources.json。

三种时间系数、Larmor进动、驱动退相干与倾斜场Bloch轨迹。
003使用精确阻尼解析解；原教程的exp(-g*t)*cos(delta*t)仅作近似。
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import qutip as qt
from helpers import plot_style
from workflow import run_directory, save_plot, save_json, require, max_error
from week02_qubit_rabi import OPTIONS

def tutorial001():
    a=qt.destroy(4);number=a.dag()*a;quadrature=a+a.dag()
    t=np.linspace(0,10,1001)
    def coefficient(t,frequency):return np.cos(frequency*t)
    function=qt.QobjEvo([number,[quadrature,coefficient]],args={'frequency':1.})
    # 显式用解释式字符串，无需额外安装Cython或编译器。
    from qutip.core.coefficient import CompilationOptions
    string_coeff=qt.coefficient('cos(frequency*t)',args={'frequency':1.},compile_opt=CompilationOptions(use_cython=False))
    string=number+quadrature*string_coeff
    array=qt.QobjEvo([number,[quadrature,np.cos(t)]],tlist=t)
    probes=np.linspace(0.015,9.985,137)
    function_errors=[];array_errors=[]
    for time in probes:
        expected=number+np.cos(time)*quadrature
        function_errors.append(max_error(function(time).full(),string(time).full()))
        # 先转数组再相减，防止Qobj稀疏相减时清理小元素，掩盖插值残差。
        array_errors.append(max_error(array(time).full(),expected.full()))
    override=(function(0.7,frequency=2.)-(number+np.cos(1.4)*quadrature)).norm()
    restored=(function(0.7)-(number+np.cos(0.7)*quadrature)).norm()
    # QobjEvo亦可做代数运算、dag、构造Liouvillian。
    algebra=((function.dag()*function)(0.7)-function(0.7).dag()*function(0.7)).norm()
    liouvillian=(qt.liouvillian(function,c_ops=[0.1*a])(0.7)-qt.liouvillian(function(0.7),[0.1*a])).norm()
    metrics={'function_vs_string':float(max(function_errors)), 'array_interpolation_error':float(max(array_errors)),
             'argument_override_error':float(override),'argument_restored_error':float(restored),
             'algebra_error':float(algebra),'liouvillian_error':float(liouvillian)}
    require(max(metrics.values())<1e-7,'001 QobjEvo复现失败')
    return metrics

def tutorial002():
    # 保留原教程初态(2|0>+|1>)/sqrt(5)和H=sigma_z。
    psi=(2*qt.basis(2,0)+qt.basis(2,1)).unit()
    t=np.linspace(0,10,501)
    def linear(t):return 0.3*t
    def periodic(t):return np.cos(0.5*t)
    models=[('constant',qt.sigmaz(),2*t),
            ('linear',qt.QobjEvo([[qt.sigmaz(),linear]]),0.3*t*t),
            ('periodic',qt.QobjEvo([[qt.sigmaz(),periodic]]),4*np.sin(t/2))]
    curves={};metrics={}
    for name,h,angle in models:
        result=qt.sesolve(h,psi,t,e_ops=[qt.sigmax(),qt.sigmay(),qt.sigmaz()],options=OPTIONS)
        xyz=np.asarray(result.expect)
        reference=np.array([0.8*np.cos(angle),0.8*np.sin(angle),np.full(t.size,0.6)])
        curves[name]=xyz;metrics[name+'_analytic_error']=max_error(xyz,reference)
        metrics[name+'_bloch_radius_error']=max_error(np.sum(xyz*xyz,axis=0),1)
    require(max(metrics.values())<1e-7,'002 Larmor解析基准失败')
    return t,curves,metrics

def exact_dephasing_z(t,delta,g):
    # z''+2g z'+delta² z=0，z(0)=1，z'(0)=0；包含欠/临界/过阻尼。
    discriminant=delta*delta-g*g
    if abs(discriminant)<1e-12:
        return np.exp(-g*t)*(1+g*t)
    if discriminant>0:
        rate=np.sqrt(discriminant)
        return np.exp(-g*t)*(np.cos(rate*t)+g/rate*np.sin(rate*t))
    rate=np.sqrt(-discriminant)
    # 指数形式比exp(-gt)*cosh(rate*t)更耐溢出。
    return 0.5*((1+g/rate)*np.exp((-g+rate)*t)+(1-g/rate)*np.exp((-g-rate)*t))

def density_diagnostics(states):
    return {'trace_error':max_error([s.tr() for s in states],1),
            'hermiticity_error':float(max((s-s.dag()).norm() for s in states)),
            'minimum_eigenvalue':float(min(np.min(s.eigenenergies()) for s in states))}

def tutorial003():
    delta=2*np.pi;g=0.25;t=np.linspace(0,5,601);psi=qt.basis(2,0)
    h=delta*qt.sigmax()/2
    result=qt.mesolve(h,psi.proj(),t,[np.sqrt(g)*qt.sigmaz()],e_ops=[qt.sigmaz()],options=OPTIONS)
    exact=exact_dephasing_z(t,delta,g);approx=np.exp(-g*t)*np.cos(delta*t)
    metrics={'exact_analytic_error':max_error(result.expect[0],exact),
             'upstream_approximation_error':max_error(result.expect[0],approx),
             'density':density_diagnostics(result.states)}
    theta=0.2*np.pi;tilted=delta*(np.cos(theta)*qt.sigmaz()+np.sin(theta)*qt.sigmax())
    trajectories={}
    for name,cops in [('unitary',[]),('dephasing',[np.sqrt(0.5)*qt.sigmaz()]),('relaxation',[np.sqrt(0.5)*qt.sigmam()])]:
        res=qt.mesolve(tilted,psi.proj(),t,cops,e_ops=[qt.sigmax(),qt.sigmay(),qt.sigmaz()],options=OPTIONS)
        trajectories[name]=np.asarray(res.expect)
        metrics[name+'_density']=density_diagnostics(res.states)
    # 另外单独检验无驱动的T1、T2：g=|0>, e=|1>, C=|g><e|。
    ground=qt.basis(2,0);excited=qt.basis(2,1);gamma=0.7;gamma_phi=0.4
    decay=qt.mesolve(0*qt.sigmaz(),excited.proj(),t,[np.sqrt(gamma)*ground*excited.dag()],e_ops=[excited.proj()],options=OPTIONS)
    plus=(ground+excited).unit()
    dephase=qt.mesolve(0*qt.sigmaz(),plus.proj(),t,[np.sqrt(gamma_phi/2)*qt.sigmaz()],e_ops=[qt.sigmax()],options=OPTIONS)
    metrics['T1_population_error']=max_error(decay.expect[0],np.exp(-gamma*t))
    metrics['T2_coherence_error']=max_error(dephase.expect[0],np.exp(-gamma_phi*t))
    require(metrics['exact_analytic_error']<1e-7 and metrics['T1_population_error']<1e-7 and metrics['T2_coherence_error']<1e-7,'003主方程解析基准失败')
    for key in ['density','unitary_density','dephasing_density','relaxation_density']:
        d=metrics[key];require(d['trace_error']<1e-7 and d['hermiticity_error']<1e-10 and d['minimum_eigenvalue']>-1e-7,'密度矩阵物理性失败')
    return t,np.asarray(result.expect[0]),exact,approx,trajectories,metrics

def main(output=None):
    folder=run_directory('week02-tutorials',output)
    first=tutorial001();t,curves,second=tutorial002();td,z,exact,approx,trajectories,third=tutorial003()
    with plot_style():
        fig,axes=plt.subplots(3,1,figsize=(8,9),layout='constrained',sharex=True)
        for ax,(name,xyz) in zip(axes,curves.items()):
            ax.plot(t,xyz[1],label='<sigma_y>');ax.set(ylabel='<sigma_y>',title=f'Larmor precession: {name} field');ax.legend()
        axes[-1].set_xlabel('Time (dimensionless)');save_plot(fig,folder,'tutorial002_larmor');plt.close(fig)
        fig,ax=plt.subplots(layout='constrained')
        ax.plot(td,exact,label='Exact damped oscillator')
        ax.plot(td[::10],z[::10],'o',ms=3,label='mesolve')
        ax.plot(td,approx,'--',label='Upstream approximation')
        ax.set(xlabel='Time (dimensionless)',ylabel='<sigma_z>',title='Driven qubit with dephasing')
        ax.legend();save_plot(fig,folder,'tutorial003_dephasing');plt.close(fig)
        fig=plt.figure(figsize=(12,4),layout='constrained')
        for i,(name,xyz) in enumerate(trajectories.items(),1):
            ax=fig.add_subplot(1,3,i,projection='3d')
            sphere=qt.Bloch(fig=fig,axes=ax);sphere.add_points(xyz,meth='l');sphere.render();ax.set_title(name)
        save_plot(fig,folder,'tutorial003_bloch');plt.close(fig)
    np.savetxt(folder/'dephasing.csv',np.column_stack([td,z,exact,approx]),delimiter=',',header='t,mesolve,exact,upstream_approx',comments='')
    metrics={'001':first,'002':second,'003':third};save_json(folder,'metrics.json',metrics);return folder,metrics

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output')
    args=parser.parse_args();folder,metrics=main(args.output);print(metrics);print(folder)
