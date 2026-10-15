# 第1周：从数组到物理计算

目标：看懂并亲自改写本项目的 `helpers.py` 与谐振子脚本。建议5次，每次约3小时。运行结果已经由电脑验证；你的理解程度用最后的闭卷题确认。

## 第1次：Python最小语法集

先打开 VS Code → 文件 → 打开文件夹 → 选择 `D:\Master.s-work\quantum-optics-notes`。按 `Ctrl+Shift+P`，选择 `Python: Select Interpreter`，选现有 main 环境的 Python 3.12。若列表没有它，输入解释器路径：

```text
D:\Master.s-work\Software-needed-for-master.s-work\python\main\python.exe
```

项目右下角显示的解释器应与上面一致。终端输入 `python` 可能命中 Windows 的商店占位程序；运行项目时优先使用根目录 `run_all.ps1` 或下面的完整路径，不需要改全局 PATH。

在 `code/my_practice.py` 中亲自输入：

```python
def square(x):             # def定义函数；冒号后缩进4个空格
    return x**2           # **是乘方

values = []               # 列表：按顺序装东西
for i in range(1, 11):     # 含1，不含11
    values.append(square(i))
print(values)             # 应为[1,4,9,...,100]
```

知识点：`=`是赋值，`==`是比较；`#`后是注释；`import numpy as np`表示导入工具库并取短名。函数先接收输入，再返回结果。只写函数定义不会运行它，需要调用 `square(3)`。

动手：把上面的平方改成立方；再写 `energy(n)` 返回 `n+0.5`。不看例子写出1到10的平方列表，再解释为什么 `range(1,10)` 少一个数。

## 第2次：NumPy数组与广播

```python
import numpy as np
x = np.linspace(0, 10, 101)  # 101个点，两端均包含，间距0.1
a = np.array([1., 2., 3.])
print(a.shape)               # (3,)，一维数组
print(a[1:3])                # [2.,3.]，切片右端不含
print(a*2)                   # [2.,4.,6.]，逐元素运算
b = np.arange(6).reshape(2, 3)
print(b + np.array([10,20,30]))
print(b @ b.T)               # 矩阵乘法，结果shape=(2,2)
```

广播不是复制物理空间，而是数组运算的对齐规则：从最右一维比较，长度相同或其中一方为1就可相容。`(2,3)+(3,)`成立，`(2,3)+(2,)`失败。若想给两行分别加1、2，使用 `np.array([1,2])[:,None]`，形状变成 `(2,1)`。

动手：打印 `(3,)`、`(3,1)`、`(1,3)` 的 `shape`；预测它们相加后的形状再运行。比较 `[1,2]*2` 与 `np.array([1,2])*2`，前者重复列表，后者乘数值。

## 第3次：复数、共轭、本征值、矩阵指数

```python
sx = np.array([[0,1],[1,0]], dtype=complex)
energies, vectors = np.linalg.eigh(sx)
print(energies)                    # [-1.,1.]
print(sx @ vectors[:,0])           # 等于energies[0]*vectors[:,0]
psi = np.array([1,1j])/np.sqrt(2)
print(np.vdot(psi,psi))            # vdot自动共轭第一个参数，结果1
print(sx.conj().T)                 # 厄米共轭
```

本征方程 `Hv=Ev` 的意思是：H作用后方向不变，只乘一个数E。厄米矩阵的本征值是实数，使用 `eigh`；普通非厄米问题再考虑 `eig`。返回的本征矢在列上；整体正负号或复相位可以变化，并不代表求错。

矩阵指数是 `I+A+A²/2!+…`。SciPy 的 `expm(A)`计算这个级数代表的矩阵函数；`np.exp(A)`只对每个元素取指数，不能替代它。第1周小测专门验证这一点。

动手：手算σx的本征值±1与本征矢；用本征方程残差验证代码。把 `dtype=complex` 去掉，观察为什么初始实矩阵仍可用于乘复数，但向已有实数组直接写入虚数需要小心。

## 第4次：积分、微分方程、规范画图

积分：采样的函数值与真实坐标要一起传入。`integrate_samples(x**2,x)`在[0,1]应得到1/3。简单地 `sum(y)`没有乘步长，也没有处理积分规则，不能当积分。

经典谐振子 `x''=-x` 改成两个一阶方程：`x'=v, v'=-x`。`solve_ivp` 的输入 `y=[x,v]`，输出形状是 `(变量数,时间点数)`。`t_eval`决定你在哪里取结果，内部积分器仍会自己选步长。

运行 `week01_python_basics.py`，对照数值 `x(t)` 与 `cos(t)`，并检查能量 `(x²+v²)/2`。图必须有横轴、纵轴、单位或无量纲说明、图例。保存 PNG 便于查看，PDF 便于放报告。

```powershell
& 'D:\Master.s-work\Software-needed-for-master.s-work\python\main\python.exe' -B -X utf8 'D:\Master.s-work\quantum-optics-notes\code\week01_python_basics.py'
```

动手：把初速度改成0.5，先写解析解 `cos(t)+0.5*sin(t)`，再对照图；不要继续用原来 `cos(t)` 作为答案。更改物理参数后也要更新验收基准。

## 第5次：谐振子与工具封装

先读下一篇推导，再打开 `code/week01_harmonic_oscillator.py`。按“坐标网格→三对角矩阵→求前三态→连续归一化→相位对齐→对照与误差”六步找代码。

工具函数放在 `code/helpers.py`，物理模型放在周脚本。现有工具来自你的桌面文件，保持原API：

| 工具 | 输入和用途 | 常见误区 |
|---|---|---|
| `normalize_state` | 离散复向量；使Σ|ψᵢ|²=1 | 不能拿它代替连续积分归一化 |
| `normalize_wavefunction` | 连续网格波函数与x；使∫|ψ|²dx=1 | 不能忘记坐标间距 |
| `integrate_samples` | 积分；支持复数与批量轴 | 坐标必须严格递增 |
| `diagonalize_hermitian` | 厄米本征问题 | 本征矢按列读取 |
| `matrix_exp` | 小型稠密矩阵指数 | `np.exp`不是矩阵指数 |
| `evolve_constant` | 恒定H的expm演化 | 初态对应times[0]，不一定是0 |
| `solve_schrodinger` | solve_ivp量子演化 | 不是保幺正算法，要检查范数 |
| `plot_style` | `with`包住绘图，退出恢复风格 | 导入工具不应改全局样式 |
| `save_figure` | PNG/PDF/SVG，默认拒绝覆盖 | 接受完整文件路径，非旧对话的save_fig名称 |

完整工具手册在 `docs/helpers_guide.md`。该手册原来的相邻文件位置描述针对桌面原包；本仓库工具在code/，测试在tests/，使用根目录一键验证即可。

## 卡住时按这个顺序

1. `ModuleNotFoundError`：核对完整解释器路径与main，不要直接重装库。
2. 数组形状错误：在报错上一行打印各对象的 `shape`。
3. 图不重合：先查单位、网格边界、归一化和整体相位。
4. 数值只显示1.00：提高打印精度，读取JSON误差，不凭四舍五入下结论。
5. 仍卡住：记“输入、期望、实际、完整报错、试过什么”；当天自查、次日查资料、第三日求助。

资料：NumPy用户指南 `https://numpy.org/doc/stable/user/`；SciPy积分 `https://docs.scipy.org/doc/scipy/reference/integrate.html`；矩阵指数 `https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.expm.html`；Matplotlib入门 `https://matplotlib.org/stable/tutorials/pyplot.html`。本项目已用当前安装版本实测这些练习，不要求额外刷完课程。
