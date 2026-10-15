# 第一周：自己的 helpers.py

## 1. 这次完成什么

根据桌面上的《量子光学与 QuTiP 预习计划》，第一周对应的是 **NumPy（数组、广播、线性代数）、SciPy（solve_ivp、积分）、Matplotlib（规范绘图）**。文中的代码片段按章节顺序阅读，后续片段会复用前面定义的 `np`、`sx`、`t` 等变量；需要独立可运行版本请执行 [example_usage.py](<example_usage.py>)。

这份工具模块只完成“封装自己的 helpers.py”这一项，并用小例子验证功能；不把第二周的 QuTiP、整套谐振子定态练习、闭卷小测混进交付范围。

交付文件：

- [helpers.py](<helpers.py>)：14 个公开工具函数，完整中文模块说明、参数/返回值注释及必要的数值检查。
- [example_usage.py](<example_usage.py>)：可直接运行的演示；包含解析解交叉验证与规范绘图。
- [test_helpers.py](<test_helpers.py>)：基于标准库 unittest 的回归测试，无须安装 pytest。
- [requirements.txt](<requirements.txt>)：环境依赖；不依赖 QuTiP。

设计原则：只封装重复出现的步骤，不替代理解 NumPy/SciPy；不在导入时计算或画图；不自动归一化演化结果来掩盖数值误差；遇到输入/形状错误明确报错。

## 2. 使用你现有的环境运行（Windows，不修改环境）

已使用你指定的 Python 3.12 环境完成验证，无须创建虚拟环境，无须执行 pip/conda 安装命令。在本任务目录打开 PowerShell：

```powershell
$python = 'D:\Master.s-work\Software-needed-for-master.s-work\python\main\python.exe'
# 禁用字节码写入；绘图缓存和测试临时文件放在当前任务目录。
$env:PYTHONDONTWRITEBYTECODE = '1'
$env:MPLCONFIGDIR = Join-Path $PWD '.mpl-cache'
$env:TEMP = Join-Path $PWD '.test-temp'
$env:TMP = $env:TEMP
New-Item -ItemType Directory -Path $env:TEMP -Force | Out-Null
& $python -B -X utf8 -m unittest -v test_helpers
& $python -B -X utf8 example_usage.py
```

这些变量仅影响当前 PowerShell 进程及其子进程，不修改系统环境变量，也不修改指定 Python 环境。测试会创建并清理临时文件；`.mpl-cache` 是项目侧绘图缓存，不是 Python 环境。目前 21 项回归测试已全部通过。新增 [验证.ps1](<验证.ps1>) 可一键完成测试、示例与临时缓存清理，并在结束后恢复进程环境变量：

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\验证.ps1
```

这里 Bypass 只作用于新启动的进程，不修改系统执行策略；脚本使用 UTF-8 BOM，兼容 Windows PowerShell 5.1 的中文编码。推荐日常使用此命令，不必手动管理缓存。

[requirements.txt](<requirements.txt>) 只作为依赖声明保留，不需要在现有环境执行安装。不要把自己写的脚本命名为 `numpy.py`、`scipy.py`、`matplotlib.py`，否则会遮蔽真正的库。

示例默认使用无窗口 Agg 后端，运行后把图和数值指标写入本目录下 `outputs/`。普通学习脚本可以不设置 Agg 并调用 `plt.show()`。再次运行示例会覆盖它自己的示例输出，工具函数本身默认拒绝覆盖。

## 3. 怎么把模块用在自己的周代码里

### 同目录脚本

把 [helpers.py](<helpers.py>) 放到你的脚本旁边即可：

```python
import numpy as np
from helpers import normalize_state, expectation

psi = normalize_state([1, 1j])
sigma_y = np.array([[0, -1j], [1j, 0]])
print(expectation(sigma_y, psi))  # 约为 (1+0j)
```

可以导入模块命名空间：`import helpers as hp`，然后使用 `hp.normalize_state(...)`。不建议 `from helpers import *`，不利于追踪函数来源。

### 放入计划要求的仓库

建议将此模块与各周脚本直接放在 `quantum-optics-notes/code/` 中；从仓库根目录运行 `python code/week01_your_script.py` 时，Python 会把脚本所在的 `code/` 加入导入搜索路径。若在别的文件夹或 Notebook 中工作，先确认工作目录，或将本模块所在文件夹加入导入路径；不建议长期硬编码桌面绝对路径。

这是单文件学习工具，不是完整 pip 包。以后功能多了再拆模块；当前不必引入复杂工程结构。

## 4. 数据形状与数学约定（最重要）

| 数据 | 形状 | 约定 |
|---|---|---|
| 离散态 `psi` | `(N,)` | 一维复数组，`psi[i]` 为第 i 个基态振幅 |
| 算符/哈密顿量 | `(N,N)` | 矩阵乘法使用 `@`；`*` 是逐元素乘法 |
| 多时刻态 `states` | `(Nt,N)` | **每一行是一个时刻**；`states[j]` 是 t[j] 的态 |
| 连续网格 `x` | `(Nx,)` | 有限实数、严格递增、至少 2 点 |
| 网格波函数 `psi_x` | `(Nx,)` | 每个点一个复振幅，使用积分范数 |
| 本征矢矩阵 `evecs` | `(N,N)` | **每一列是一个本征矢**；`evecs[:,j]` |

一维态不要写成 `[[1],[0]]`，应写 `np.array([1,0], dtype=complex)`；已有列向量可先 `.ravel()`，但不要对密度矩阵这么做。

离散基归一化：

```text
sum_i |psi_i|² = 1
<A> = psi† A psi / (psi† psi)
```

连续坐标归一化：

```text
integral |psi(x)|² dx = 1
<f(x)> = integral f(x) |psi(x)|² dx  （若已连续归一化）
```

两者不能混用。`normalize_state` 与 `normalize_wavefunction` 特意分成两个函数。如果在均匀网格上得到离散本征矢 v，常见近似为连续波函数 `v/sqrt(dx)`，也可明确采用本模块的积分归一化；但这仍需要和你的离散化边界条件、求积规则一致。

默认 ħ=1 是自然单位约定，**不是 ħ 的 SI 数值**。如果使用焦耳和秒，就应显式传实际 ħ 并保持 H、时间一致。为避免极端量级，初学阶段推荐先无量纲化。

## 5. API 一览

| 函数 | 用途 | 返回 |
|---|---|---|
| `dagger(a)` | 二维矩阵厄米共轭 | `(N,M)` 复数组 |
| `commutator(a,b)` | `a@b-b@a` | 同维复方阵 |
| `is_hermitian(a,atol=...,rtol=...)` | 近似厄米性检查 | bool |
| `normalize_state(psi)` | 离散态归一化 | `(N,)` 复数组 |
| `integrate_samples(y,x,axis=-1,method='simpson')` | 采样数据积分 | 标量或去掉积分轴的数组 |
| `normalize_wavefunction(psi,x,method='simpson')` | 积分范数归一化 | `(Nx,)` 复数组 |
| `expectation(operator,psi)` | 离散纯态期望值 | complex |
| `diagonalize_hermitian(a)` | 厄米矩阵对角化 | `(evals,evecs)` |
| `matrix_exp(a)` | 真正矩阵指数 | 同维复方阵 |
| `evolve_constant(h,psi0,times,hbar=1)` | 常量 H 的指数演化 | `(Nt,N)` |
| `solve_schrodinger(h,psi0,times,...)` | 常量或时间依赖 H 的 ODE 演化 | `(Nt,N)` |
| `plot_style(font_size=11,dpi=120)` | 局部绘图风格 | 上下文管理器 |
| `new_figure(xlabel=...,ylabel=...,title=...)` | 创建规范单面板 | `(fig,ax)` |
| `save_figure(fig,path,dpi=300,overwrite=False)` | 保存 PNG/PDF/SVG | Path |

表中省略的可选参数多为仅限关键字参数，完整签名与说明见 [helpers.py](<helpers.py>)。Python 中执行 `help(helpers.solve_schrodinger)` 可以直接查看中文文档。

## 6. 分功能详细示范

### 6.1 共轭、对易子、期望值

```python
import numpy as np
from helpers import dagger, commutator, is_hermitian, normalize_state, expectation

sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.diag([1., -1.])
print(is_hermitian(sy))                      # True
assert np.allclose(dagger(sy), sy)
assert np.allclose(commutator(sx, sy), 2j*sz)
psi = normalize_state([1, 1j])
assert np.isclose(expectation(sy, psi), 1)
```

`np.vdot(psi,phi)` 共轭第一个参数；`np.dot` 对复数态不自动共轭。`expectation` 也允许未归一化的非零态，并除以范数。返回复数是有意设计：厄米算符虚部应仅为数值残差，可用 `np.real_if_close(value)` 判断，而不要无条件 `.real` 掩盖非厄米错误。

### 6.2 本征分解与 NumPy 广播

```python
from helpers import diagonalize_hermitian
E, V = diagonalize_hermitian(sx)
# 每列分别乘以该列对应的本征值；E[None,:] 的形状为 (1,N)。
assert np.allclose(sx @ V, V * E[None, :])
ground_state = V[:, 0]
```

本征值按升序返回。整体相位任意，因此不要把数值本征矢与解析向量逐元素比较；非简并态比较重叠的绝对值更合理。简并情况下应比较子空间投影。非厄米矩阵应使用 `scipy.linalg.eig`，本函数故意拒绝。公共 `is_hermitian` 保留用户设定的绝对/相对容差语义；对角化和演化内部则先按矩阵分量尺度缩放再检查，避免 `1e-13` 能量单位下明显非厄米的矩阵被固定绝对容差误判。容许的舍入残差约为相对整体尺度 1e-10；更严格物理模型请另做专用验证。

### 6.3 采样积分与连续归一化

```python
from helpers import integrate_samples, normalize_wavefunction
x = np.linspace(-8, 8, 1601)
psi_x = normalize_wavefunction(np.exp(-x**2/2), x)
p = np.abs(psi_x)**2
norm = integrate_samples(p, x)               # 约 1
x2 = integrate_samples(x**2*p, x)            # 约 0.5
# 批量积分：每一行是一个被积函数，沿最后一轴求积。
moments = integrate_samples(np.stack([p, x*p, x**2*p]), x, axis=-1)
print(np.real_if_close(moments))             # 约 [1,0,0.5]
```

Simpson 更适合光滑且采样充分的数据，但不是所有输入都比梯形法准确。**不均匀/欠采样网格上的 Simpson 权重可能为负**，例如 x=[0,1,4]、p=[1,0,0] 可得到负的“概率积分”。归一化函数会拒绝非正积分，不会取绝对值掩盖问题。对不规则或尖峰数据优先试 `method='trapezoid'`，且归一化、后续求积必须使用同一规则。即使积分为正，也仍需检查网格收敛。

区间 [-8,8] 只是例子，不等于整个实轴。验证时应分别扩大区间和加密网格，避免把尾部截断与求积误差混在一起。仅有两个点虽然能调用积分，也不能保证描述真实波函数。

### 6.4 矩阵指数（第一周验收的常见易错点）

```python
from helpers import matrix_exp
alpha = 0.3
U = matrix_exp(-1j*alpha*sx)
assert np.allclose(U, np.cos(alpha)*np.eye(2)-1j*np.sin(alpha)*sx)
assert np.allclose(U.conj().T @ U, np.eye(2))
```

矩阵指数是 `I+A+A@A/2!+...`；`np.exp(A)` 只对每个元素取指数，不能代替 `matrix_exp(A)`。本函数也支持非厄米矩阵，但稠密存储适合小规模问题；大规模稀疏模型应另用 `expm_multiply`。有限输入并不保证 exp 的结果可表示（如 exp(1000) 会溢出）；现在会显式抛 ValueError 而不是返回 Inf/NaN。选择合适单位与数值尺度仍是你的建模责任。

### 6.5 指数演化与 solve_ivp

```python
from helpers import evolve_constant, solve_schrodinger
H = 0.5*sx
psi0 = np.array([1., 0.])
t = np.linspace(0, 10, 201)
exact = evolve_constant(H, psi0, t)
numeric = solve_schrodinger(H, psi0, t, rtol=1e-9, atol=1e-11)
P1 = np.abs(numeric[:, 1])**2                # 第二个基态的概率
norms = np.sum(np.abs(numeric)**2, axis=1)   # 每个时刻的范数
error = np.max(np.linalg.norm(numeric-exact, axis=1))
print('max state error:', error)
print('max norm drift:', np.max(np.abs(norms-1)))
```

方程约定：`i*hbar*dpsi/dt=H(t)@psi`。初态定义在 **t[0]**，不是一定在 0。非零起始时间时，指数演化采用 `t-t[0]`。两种函数都不替你归一化初态或输出。

时间依赖例子：

```python
def H_of_t(time):
    return (1+0.2*np.cos(time))*sx/2
states = solve_schrodinger(H_of_t, psi0, t)
```

本模块只支持厄米 H 的纯态薛定谔演化，不支持 Lindblad、密度矩阵或非厄米耗散。`solve_ivp` 的 `t_eval` 只是输出点，不是内部固定步长；加密输出点不等于提高积分精度。应把 rtol/atol 收紧后比较，必要时针对快振荡另写更细致的求解器配置。本封装提供 `max_step` 内部步长上限（默认 `np.inf`），但不暴露 events 等全部高级选项。窄脉冲尤其需要显式设置，例如 `solve_schrodinger(H_of_t, psi0, t, max_step=脉宽/10)`；更紧容差、更多输出点都不能保证内部步长采到脉冲。

默认 DOP853 支持复数；允许 RK23/RK45/BDF，不允许直接用不支持复数的 Radau/LSODA。ODE 方法不是严格幺正积分器，因此“算完后再归一化”不能当作误差控制。对于 H(t) 随时间变化，本函数在每次求导时验证矩阵，优先可读性而非性能。

### 6.6 规范绘图与保存

```python
import matplotlib.pyplot as plt
from helpers import plot_style, new_figure, save_figure
with plot_style(font_size=11):
    fig, ax = new_figure(xlabel='Time (1/Omega)', ylabel='Population', title='Evolution')
    ax.plot(t, P1, label='Numerical')
    ax.plot(t, np.sin(t/2)**2, '--', label='Analytical')
    ax.legend()
    ax.set_ylim(-0.05, 1.05)
    save_figure(fig, 'outputs/my_population.png')
    save_figure(fig, 'outputs/my_population.pdf')
    plt.close(fig)
```

风格仅在 `with` 内生效，不污染别人的绘图。坐标应写量与单位，图例区分数值与解析；显示概率不要无故截断误差。PNG 用于预览，PDF/SVG 适合矢量插图。

`save_figure` 默认禁止覆盖；确实要替换文件时显式传 `overwrite=True`。先写同目录临时文件再原子提交，保存中途失败不破坏原图，另一进程同时创建目标也不会被误覆盖。默认防覆盖提交需要文件系统支持硬链接（如本机 NTFS）；不支持时保留 OSError，不会静默退回不安全覆盖。原子提交不等同于断电持久性保证。父目录自动创建；相对路径相对于**运行进程的工作目录**。`dpi` 对 PNG 有意义，PDF 中线条一般仍为矢量。函数不会自动关闭图，批量计算请自行 `plt.close(fig)`。

本工具不硬编码中文字体。需要中文时，在 `with plot_style():` 内设置已安装的中文字体，例如 `plt.rcParams['font.sans-serif']=['Microsoft YaHei']`；务必确认该字体实际存在，不能把缺字警告当作正常排版。

## 7. 错误与自查

| 现象 | 优先检查 |
|---|---|
| 导入提示缺少 scipy/matplotlib | 是否调用了指定的现有 Python；不要为了排错直接修改该环境 |
| 矩阵/态维数不匹配 | H 是 `(N,N)`，态是 `(N,)`，不是列矩阵 |
| 零态归一化失败 | 零向量不能表示量子态，不应除以零 |
| 网格报错 | x/t 是否含重复、倒序、NaN/Inf 或虚部 |
| eigh 被拒绝 | H 是否真的厄米，是否把矩阵乘法写成了 `*` |
| 连续积分不为 1 | 是否混用离散范数/连续范数，积分规则是否一致 |
| solve_ivp 失败或漂移 | 单位、容差、模型尺度与快振荡；比较解析或指数结果 |
| 结果相位与书上不同 | 是否只是整体相位；简并时是否比较了任意基 |
| 保存 FileExistsError | 防覆盖机制；确认后显式 overwrite=True |

非法数值/形状一般抛 ValueError，ODE 失败抛 RuntimeError，已有图文件抛 FileExistsError；文件系统权限问题保留原生异常。回调函数自身的异常也会向上传递。模块不是通用数据清洗器，不负责自动纠正用户任意输入类型。

## 8. 验收标准与边界

测试包括：Pauli 对易关系、复共轭、极大/极小尺度态归一化、连续高斯范数与二阶矩、复数批量积分、厄米本征方程、矩阵指数解析式、非零时间起点与 ħ、常量/时间依赖 ODE 与解析解、范数漂移、非法输入、求解失败、图像保存与风格恢复。

实际执行结果见 [验证记录.md](<验证记录.md>)；完整演示逻辑见 [example_usage.py](<example_usage.py>)。测试只验证这些有限基准，不意味着任意物理模型都正确。后续学习仍要针对自己的边界条件、Hilbert 空间截断、网格、步长/容差做独立收敛测试。

刻意没有加入：QuTiP 包装器、谐振子专用求解器、密度矩阵演化、控制优化/EP 工具。这些超出本周这一项任务。你现在可以先读懂并运行本模块，再在自己的练习中复用。
