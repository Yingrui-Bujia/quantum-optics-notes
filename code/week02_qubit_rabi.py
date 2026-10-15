"""第2周交付：Rabi振荡、零点涨落、热态分布、张量积与截断诊断。"""
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import qutip as qt
from helpers import plot_style
from workflow import run_directory, save_plot, save_json, require, max_error, real_parameter

OPTIONS = {'atol':1e-11,'rtol':1e-9,'method':'dop853','store_states':True,'normalize_output':False}

def rabi(omega=1., detuning=0., samples=601):
    omega=real_parameter(omega,'omega',positive=True)
    detuning=real_parameter(detuning,'detuning')
    if isinstance(samples,bool) or not isinstance(samples,(int,np.integer)) or samples<2:
        raise ValueError('samples 必须是至少2的整数')
    # basis(2,0)仅表示第0基矢。若H0=-omega0*sz/2，它才是能量基态。
    # 此处旋转系 H=Delta*sz/2+Omega*sx/2；P1=|1><1|。
    t=np.linspace(0,4*np.pi/omega,samples)
    initial=qt.basis(2,0);p1=qt.basis(2,1).proj()
    h=(detuning*qt.sigmaz()+omega*qt.sigmax())/2
    result=qt.sesolve(h,initial,t,e_ops=[qt.sigmaz(),p1],options=OPTIONS)
    rate=np.hypot(omega,detuning)
    reference=(omega/rate)**2*np.sin(rate*t/2)**2
    metrics={'population_error':max_error(result.expect[1],reference),
             'sigma_z_error':max_error(result.expect[0],1-2*reference),
             'norm_drift':max_error([s.norm()**2 for s in result.states],1)}
    require(max(metrics.values())<1e-7,'Rabi解析基准失败')
    return t,result,reference,metrics

def zero_point(cutoff=20, highest_n=5):
    if not isinstance(cutoff,int) or isinstance(cutoff,bool) or cutoff<2:
        raise ValueError('cutoff 必须是至少2的整数')
    if not isinstance(highest_n,int) or isinstance(highest_n,bool) or not 0 <= highest_n < cutoff-1:
        raise ValueError('最高验证态须小于cutoff-1，以免落在截断边界')
    a=qt.destroy(cutoff);ad=qt.create(cutoff)
    x=(a+ad)/np.sqrt(2);p=-1j*(a-ad)/np.sqrt(2)
    h=ad*a+0.5*qt.qeye(cutoff)
    rows=[]
    for n in range(highest_n+1):
        ket=qt.basis(cutoff,n)
        rows.append({'n':n,'quadrature_squared':float(qt.expect((a+ad)**2,ket)),
                     'expected':2*n+1,'energy':float(qt.expect(h,ket)),
                     'x_variance':float(qt.variance(x,ket)),'p_variance':float(qt.variance(p,ket))})
    correction=qt.qeye(cutoff)-cutoff*qt.basis(cutoff,cutoff-1).proj()
    defect=(qt.commutator(a,ad)-correction).norm()
    top=float(qt.expect((a+ad)**2,qt.basis(cutoff,cutoff-1)))
    metrics={'states':rows,'cutoff':cutoff,'finite_commutator_identity_error':float(defect),
             'top_state_actual':top,'top_state_infinite_formula':2*cutoff-1,
             'vacuum_energy':rows[0]['energy'],'vacuum_variance_product':rows[0]['x_variance']*rows[0]['p_variance']}
    require(max(abs(r['quadrature_squared']-r['expected']) for r in rows)<1e-12,'零点涨落基准失败')
    require(defect<1e-12,'截断对易关系验证失败')
    return metrics

def thermal(cutoff=100,nbar=3.):
    if not isinstance(cutoff,int) or isinstance(cutoff,bool) or cutoff<2:
        raise ValueError('cutoff 必须是至少2的整数')
    nbar=real_parameter(nbar,'nbar',nonnegative=True)
    n=np.arange(cutoff);q=nbar/(1+nbar);tail=q**cutoff
    if q>=1 or not np.isfinite(q):
        raise ValueError('nbar 过大，当前浮点精度无法分辨热态概率')
    infinite=(1-q)*q**n
    conditioned=infinite/(1-tail)
    rho=qt.thermal_dm(cutoff,nbar,method='operator')
    analytic=qt.thermal_dm(cutoff,nbar,method='analytic')
    p=rho.diag().real
    mean=float(qt.expect(qt.num(cutoff),rho))
    expected_mean=nbar-cutoff*tail/(1-tail)
    metrics={'cutoff':cutoff,'nbar_infinite':nbar,'omitted_probability':float(tail),
             'operator_trace':float(rho.tr()),'analytic_trace':float(analytic.tr()),
             'operator_vs_conditioned_max_error':max_error(p,conditioned),
             'analytic_vs_infinite_max_error':max_error(analytic.diag(),infinite),
             'mean_number':mean,'mean_number_error':abs(mean-nbar),
             'mean_vs_truncated_formula_error':abs(mean-expected_mean),
             'minimum_eigenvalue':float(np.min(rho.eigenenergies()))}
    require(metrics['operator_vs_conditioned_max_error']<1e-12,'热态归一化公式不符')
    require(metrics['mean_vs_truncated_formula_error']<1e-11,'热态平均数不符')
    return n,p,infinite,metrics

def tensor_demo():
    mode=qt.basis(6,2);qubit=qt.basis(2,1);joint=qt.tensor(mode,qubit)
    n=qt.tensor(qt.create(6)*qt.destroy(6),qt.qeye(2))
    sz=qt.tensor(qt.qeye(6),qt.sigmaz())
    # tensor 保留各子系统维数，ptrace(0)恢复模式的约化密度矩阵。
    return {'shape':list(joint.shape),'dims':joint.dims,'mode_number':float(qt.expect(n,joint)),
            'qubit_sigma_z':float(qt.expect(sz,joint)),
            'partial_trace_error':float((joint.ptrace(0)-mode.proj()).norm())}

def lab_frame():
    # 明确保留驱动幅度Omega；以多组omega0/Omega观察RWA误差的减小。
    curves=[]
    for ratio in (10.,20.,40.):
        omega=1.;w0=ratio*omega;t=np.linspace(0,4*np.pi,2001)
        def drive(t,w,amplitude):return amplitude*np.cos(w*t)
        h=qt.QobjEvo([-w0*qt.sigmaz()/2,[qt.sigmax(),drive]],args={'w':w0,'amplitude':omega})
        opts={**OPTIONS,'max_step':2*np.pi/w0/25}
        result=qt.sesolve(h,qt.basis(2,0),t,e_ops=[qt.basis(2,1).proj()],options=opts)
        reference=np.sin(t/2)**2
        curves.append((ratio,t,np.asarray(result.expect[0]),max_error(result.expect[0],reference)))
    require(curves[-1][3]<0.03,'高频极限未接近RWA基准')
    require(curves[-1][3]<curves[0][3],'增加频率未改善RWA误差')
    return curves

def main(output=None):
    folder=run_directory('week02-rabi',output)
    t,res,ref,rabi_metrics=rabi()
    td,detuned,dref,dm=rabi(detuning=1.2)
    zero=zero_point();n,p,pinf,th=thermal();tensor=tensor_demo();lab=lab_frame()
    convergence=[thermal(k)[3] for k in (20,40,80,100)]
    metrics={'rabi_resonant':rabi_metrics,'rabi_detuned':dm,'zero_point':zero,
             'thermal':th,'thermal_convergence':convergence,'tensor':tensor,
             'lab_frame':[{'omega0_over_Omega':r,'max_population_error_vs_RWA':err} for r,_,_,err in lab]}
    require(th['omitted_probability']<1e-10 and th['mean_number_error']<1e-8,'热态截断不够')
    with plot_style():
        fig,axes=plt.subplots(2,1,figsize=(8,7),sharex=True,layout='constrained')
        axes[0].plot(t,np.cos(t),label='Analytical cos(Omega t)')
        axes[0].plot(t[::12],np.asarray(res.expect[0])[::12],'o',ms=3,label='sesolve')
        axes[0].set(ylabel='<sigma_z>',title='Resonant Rabi oscillations');axes[0].legend()
        axes[1].plot(td,dref,label='Detuned analytical P1')
        axes[1].plot(td[::12],np.asarray(detuned.expect[1])[::12],'o',ms=3,label='sesolve')
        axes[1].set(xlabel='Time (1/Omega)',ylabel='P(1)',title='Detuning Delta / Omega = 1.2');axes[1].legend()
        save_plot(fig,folder,'rabi');plt.close(fig)
        fig,ax=plt.subplots(layout='constrained')
        ax.bar(n[:25],p[:25],label='thermal_dm (N=100)')
        ax.plot(n[:25],pinf[:25],'r.',label='Infinite geometric distribution')
        ax.set(xlabel='Photon number n',ylabel='Probability',title='Thermal distribution; nbar=3')
        ax.legend();save_plot(fig,folder,'thermal');plt.close(fig)
        fig,ax=plt.subplots(layout='constrained')
        ax.plot(t,np.sin(t/2)**2,'k--',label='Rotating wave approximation')
        for ratio,lt,lp,_ in lab:ax.plot(lt,lp,lw=1,label=f'Lab: omega0/Omega={ratio:g}')
        ax.set(xlabel='Time (1/Omega)',ylabel='P(1)',title='Lab frame vs. rotating wave approximation')
        ax.legend();save_plot(fig,folder,'rabi_lab');plt.close(fig)
    np.savetxt(folder/'rabi.csv',np.column_stack([t,res.expect[0],res.expect[1],ref]),delimiter=',',header='t,sigma_z,P1,P1_exact',comments='')
    np.savetxt(folder/'thermal.csv',np.column_stack([n,p,pinf]),delimiter=',',header='n,operator_probability,infinite_probability',comments='')
    save_json(folder,'metrics.json',metrics);return folder,metrics

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output')
    args=parser.parse_args();folder,metrics=main(args.output);print(metrics);print(folder)
