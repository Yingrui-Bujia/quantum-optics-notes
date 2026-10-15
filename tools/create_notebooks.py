"""创建中文跟做notebooks；已有笔记拒绝覆盖。以同一解释器逐格验证。"""
from pathlib import Path
from datetime import datetime,timezone
import base64
import contextlib
import io
import json
import os
import sys
import time
from unittest.mock import patch
import nbformat as nbf

ROOT=Path(__file__).resolve().parents[1]
BOOTSTRAP='''from pathlib import Path
import sys
import numpy as np
import qutip as qt
import matplotlib.pyplot as plt
from IPython.display import display
root = Path.cwd()
if not (root / 'code').is_dir():
    root = root.parent
if not (root / 'code').is_dir():
    raise RuntimeError('请在项目根目录或notebooks目录运行')
sys.path.insert(0, str(root / 'code'))
print('当前解释器：', sys.executable)
options = dict(atol=1e-11, rtol=1e-9, method='dop853', store_states=True, normalize_output=False)
'''

def markdown(text):return nbf.v4.new_markdown_cell(text)
def code(text):return nbf.v4.new_code_cell(text)

def notebooks():
    first=[markdown('''# 001 · 把时间放进算符

跟做官方QobjEvo教程的核心任务。先运行每格，再改frequency，预测变化。原文在third_party/qutip_tutorials/001_qobjevo.md；作者、许可证与固定版本见docs/sources.json。

Qobj是某个固定算符；QobjEvo像一个“给定时间就返回算符”的函数。'''),code(BOOTSTRAP),
    markdown('## 三种系数表示同一H(t)\nH(t)=a†a+cos(ωt)(a+a†)。数组系数在点之间需要插值；字符串使用解释执行，省去编译器。'),
    code('''from qutip.core.coefficient import CompilationOptions
a = qt.destroy(4)
number, quadrature = a.dag()*a, a+a.dag()
def drive(t, frequency):
    return np.cos(frequency*t)
function = qt.QobjEvo([number, [quadrature, drive]], args={'frequency':1.0})
string_coefficient = qt.coefficient('cos(frequency*t)', args={'frequency':1.0}, compile_opt=CompilationOptions(use_cython=False))
string = number + quadrature*string_coefficient
times = np.linspace(0,10,1001)
array = qt.QobjEvo([number,[quadrature,np.cos(times)]], tlist=times)
test_time = 0.137
print('函数和字符串的矩阵误差：', np.max(abs(function(test_time).full()-string(test_time).full())))
print('数组插值的矩阵误差：', np.max(abs(array(test_time).full()-function(test_time).full())))
np.testing.assert_allclose(function(test_time).full(), array(test_time).full(), atol=1e-7, rtol=0)
'''),markdown('## 改参数与量子对象\n临时传入frequency不应永久改变初始参数。再看一次shape与dims，练习dag、tensor、expect。'),
    code('''changed = function(0.7, frequency=2.0)
expected = number+np.cos(1.4)*quadrature
np.testing.assert_allclose(changed.full(),expected.full(),atol=1e-12)
np.testing.assert_allclose(function(0.7).full(),(number+np.cos(0.7)*quadrature).full(),atol=1e-12)
joint = qt.tensor(qt.basis(4,2),qt.basis(2,0))
observable = qt.tensor(number,qt.qeye(2))
print('shape:',joint.shape,'dims:',joint.dims,'模式光子数:',qt.expect(observable,joint))
liouvillian = qt.liouvillian(function,c_ops=[0.1*a])
print('Liouvillian shape:',liouvillian(0.7).shape)
'''),markdown('## 自己改一处\n把频率改成2，数组也必须改成cos(2*times)。然后尝试把数组采样降到21个点，检查非节点处的误差是否变大。不要只在已有采样点比较。')]
    second=[markdown('''# 002 · 看见Larmor进动

保留官方初态(2|0⟩+|1⟩)/sqrt(5)。原文在third_party/qutip_tutorials/002_larmor-precession.md。

对H=f(t)σz，Bloch向量绕z轴转动，角度是2∫f(t)dt。先预测恒定、线性、周期场的相位，再看图。'''),code(BOOTSTRAP),
    code('''psi=(2*qt.basis(2,0)+qt.basis(2,1)).unit()
times=np.linspace(0,10,501)
def linear(t):return 0.3*t
def periodic(t):return np.cos(0.5*t)
models=[('constant',qt.sigmaz(),2*times),
        ('linear',qt.QobjEvo([[qt.sigmaz(),linear]]),0.3*times**2),
        ('periodic',qt.QobjEvo([[qt.sigmaz(),periodic]]),4*np.sin(times/2))]
trajectories={}
for name,H,angle in models:
    result=qt.sesolve(H,psi,times,e_ops=[qt.sigmax(),qt.sigmay(),qt.sigmaz()],options=options)
    xyz=np.asarray(result.expect)
    reference=np.array([0.8*np.cos(angle),0.8*np.sin(angle),np.full(times.size,0.6)])
    np.testing.assert_allclose(xyz,reference,atol=1e-7,rtol=0)
    trajectories[name]=xyz
    print(name,'解析最大误差:',np.max(abs(xyz-reference)))
'''),markdown('## 同时比较相位与Bloch轨迹\n三个图的纵轴相同，便于比较；幺正演化的Bloch向量长度应保持1。'),
    code('''fig,axes=plt.subplots(3,1,figsize=(8,8),layout='constrained',sharex=True)
for ax,(name,xyz) in zip(axes,trajectories.items()):
    ax.plot(times,xyz[1],label='<sigma_y>')
    ax.set(ylabel='<sigma_y>',title=name)
    ax.legend()
axes[-1].set_xlabel('Time (dimensionless)')
display(fig);plt.close(fig)
fig=plt.figure(figsize=(6,5))
ax=fig.add_subplot(111,projection='3d')
sphere=qt.Bloch(fig=fig,axes=ax)
sphere.add_points(trajectories['constant'],meth='l');sphere.render()
display(fig);plt.close(fig)
'''),markdown('## 自己改一处\n把初态改成(|0⟩+i|1⟩)/sqrt(2)，先手算初始Bloch向量，再更改对应解析公式。恒定H若从σz改成σz/2，进动频率应减半。')]
    third=[markdown('''# 003 · 退相干与主方程

跟做官方qubit-dynamics，原文在third_party/qutip_tutorials/003_qubit-dynamics.md。本笔记补充精确阻尼解与密度矩阵检查。H=δσx/2，C=sqrt(g)σz。'''),code(BOOTSTRAP),
    markdown('## 精确解与近似解\nz满足z″+2g z′+δ²z=0，初始z=1、z′=0。下面计算欠阻尼解；教程的简单余弦作为近似显示。'),
    code('''delta=2*np.pi;g=0.25;times=np.linspace(0,5,601)
psi=qt.basis(2,0);H=delta*qt.sigmax()/2
result=qt.mesolve(H,psi.proj(),times,[np.sqrt(g)*qt.sigmaz()],e_ops=[qt.sigmaz()],options=options)
nu=np.sqrt(delta**2-g**2)
exact=np.exp(-g*times)*(np.cos(nu*times)+g/nu*np.sin(nu*times))
approx=np.exp(-g*times)*np.cos(delta*times)
z=np.asarray(result.expect[0])
np.testing.assert_allclose(z,exact,atol=1e-7,rtol=0)
print('精确解最大误差:',np.max(abs(z-exact)))
print('原文近似最大误差:',np.max(abs(z-approx)))
print('迹误差:',max(abs(rho.tr()-1) for rho in result.states))
print('最小本征值:',min(min(rho.eigenenergies()) for rho in result.states))
fig,ax=plt.subplots(figsize=(8,4),layout='constrained')
ax.plot(times,exact,label='Exact')
ax.plot(times[::10],z[::10],'o',ms=3,label='mesolve')
ax.plot(times,approx,'--',label='Upstream approximation')
ax.set(xlabel='Time (dimensionless)',ylabel='<sigma_z>');ax.legend()
display(fig);plt.close(fig)
'''),markdown('## 三种Bloch轨迹\n先设倾斜场；无耗散保持纯态，退相干和弛豫使轨迹发生不同变化。这里sigmam沿用官方约定，把第0基矢送到第1基矢。'),
    code('''theta=0.2*np.pi
H_tilt=delta*(np.cos(theta)*qt.sigmaz()+np.sin(theta)*qt.sigmax())
fig=plt.figure(figsize=(12,4),layout='constrained')
for index,(name,cops) in enumerate([('unitary',[]),('dephasing',[np.sqrt(0.5)*qt.sigmaz()]),('relaxation',[np.sqrt(0.5)*qt.sigmam()])],1):
    evolution=qt.mesolve(H_tilt,psi.proj(),times,cops,e_ops=[qt.sigmax(),qt.sigmay(),qt.sigmaz()],options=options)
    xyz=np.asarray(evolution.expect)
    for rho in evolution.states:
        np.testing.assert_allclose(rho.tr(),1,atol=1e-7)
        np.testing.assert_allclose(rho.full(),rho.full().conj().T,atol=1e-10)
        if min(rho.eigenenergies()) < -1e-7:raise RuntimeError('密度矩阵正性失败')
    ax=fig.add_subplot(1,3,index,projection='3d')
    sphere=qt.Bloch(fig=fig,axes=ax);sphere.add_points(xyz,meth='l');sphere.render();ax.set_title(name)
display(fig);plt.close(fig)
'''),markdown('## 独立核验衰减方向与因子2\n显式选g=|0⟩、e=|1⟩，直接用|g⟩⟨e|写collapse operator。'),
    code('''ground=qt.basis(2,0);excited=qt.basis(2,1);gamma=0.7
decay=qt.mesolve(0*qt.sigmaz(),excited.proj(),times,[np.sqrt(gamma)*ground*excited.dag()],e_ops=[excited.proj()],options=options)
np.testing.assert_allclose(decay.expect[0],np.exp(-gamma*times),atol=1e-7,rtol=0)
plus=(ground+excited).unit();gamma_phi=0.4
coherence=qt.mesolve(0*qt.sigmaz(),plus.proj(),times,[np.sqrt(gamma_phi/2)*qt.sigmaz()],e_ops=[qt.sigmax()],options=options)
np.testing.assert_allclose(coherence.expect[0],np.exp(-gamma_phi*times),atol=1e-7,rtol=0)
print('T1和T2解析基准通过')
'''),markdown('## 自己改一处\n把退相干算符从sqrt(gamma_phi/2)σz改成sqrt(gamma_phi)σz，预测衰减快多少，再改解析式验证。想增加g到大于delta时，改用讲义中的过阻尼公式。')]
    return {'001_qobjevo_中文跟做.ipynb':first,'002_larmor_中文跟做.ipynb':second,'003_qubit_dynamics_中文跟做.ipynb':third}

def execute(notebook):
    namespace={'__name__':'__main__'};count=0
    for cell in notebook.cells:
        if cell.cell_type!='code':continue
        count+=1;cell.execution_count=count;outputs=[];buffer=io.StringIO()
        def capture_display(obj):
            data={'text/plain':repr(obj)}
            if hasattr(obj,'savefig'):
                image=io.BytesIO();obj.savefig(image,format='png',dpi=110,bbox_inches='tight')
                data['image/png']=base64.b64encode(image.getvalue()).decode('ascii')
            outputs.append(nbf.v4.new_output('display_data',data=data,metadata={}))
        with contextlib.redirect_stdout(buffer),contextlib.redirect_stderr(buffer),patch('IPython.display.display',capture_display):
            exec(compile(cell.source,'<notebook-cell>','exec'),namespace)
        if buffer.getvalue():outputs.insert(0,nbf.v4.new_output('stream',name='stdout',text=buffer.getvalue()))
        cell.outputs=outputs
    return count

def main():
    folder=ROOT/'notebooks';folder.mkdir(exist_ok=True)
    definitions=notebooks()
    for name in definitions:
        if (folder/name).exists():raise FileExistsError(f'已有notebook不会被覆盖：{name}')
    os.chdir(ROOT);results=[]
    for name,cells in definitions.items():
        notebook=nbf.v4.new_notebook(cells=cells,metadata={'kernelspec':{'display_name':'Python 3 (main)','language':'python','name':'python3'},'language_info':{'name':'python','version':sys.version.split()[0]}})
        started=time.perf_counter();count=execute(notebook);nbf.validate(notebook)
        nbf.write(notebook,folder/name)
        results.append({'file':'notebooks/'+name,'code_cells':count,'status':'passed','seconds':time.perf_counter()-started})
    report={'method':'sequential Python execution in shared cell namespace; captured stdout and Figure display outputs; no Jupyter server or kernel registration',
            'executable':sys.executable,'utc':datetime.now(timezone.utc).isoformat(),'notebooks':results,'status':'passed'}
    (ROOT/'reports/notebook_execution.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
