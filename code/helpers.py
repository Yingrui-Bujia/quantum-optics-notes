"""第一周科学计算工具：NumPy + SciPy + Matplotlib（不依赖 QuTiP）。

约定
----
* 离散基中的态是 shape=(N,) 的复数组，范数为 sum(|psi|**2)。
* 连续坐标网格上的波函数使用积分范数 ∫|psi(x)|² dx；不要混用两种范数。
* 演化公式为 i*hbar*dpsi/dt = H*psi；时间、能量、hbar 必须单位自洽。
* 本模块不在 import 时画图、创建目录、修改全局绘图风格或运行计算。

快速例子：
    import numpy as np
    from helpers import normalize_state, expectation, evolve_constant
    psi = normalize_state([1, 1j])
    H = np.diag([0., 1.])
    states = evolve_constant(H, psi, np.linspace(0, 2, 101))
    print(expectation(H, states[-1]))

详细教程见同目录 README.md；可运行示例见 example_usage.py。
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable
import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import simpson, solve_ivp
from scipy.linalg import eigh, expm

__all__ = [
    "dagger", "commutator", "is_hermitian", "normalize_state", "integrate_samples",
    "normalize_wavefunction", "expectation", "diagonalize_hermitian", "matrix_exp",
    "evolve_constant", "solve_schrodinger", "plot_style", "new_figure", "save_figure",
]


def _array(value: ArrayLike, name: str) -> NDArray:
    """统一转换数值数组，并拒绝 NaN/Inf，避免错误悄悄传入求解器。"""
    try:
        arr = np.asarray(value, dtype=complex)
    except (ValueError, TypeError) as exc:
        raise ValueError(f"{name} 必须是数值数组") from exc
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} 不能含 NaN 或 Inf")
    return arr


def _matrix(value: ArrayLike, name: str = "matrix") -> NDArray:
    arr = _array(value, name)
    if arr.ndim != 2 or arr.shape[0] != arr.shape[1] or arr.shape[0] == 0:
        raise ValueError(f"{name} 必须是非空方阵，收到 shape={arr.shape}")
    return arr


def _state(value: ArrayLike) -> NDArray:
    arr = _array(value, "psi")
    if arr.ndim != 1 or arr.size == 0:
        raise ValueError("psi 必须是非空一维态矢；列向量请先用 ravel() 展平")
    return arr


def _real_scalar(value: float, name: str, *, allow_zero: bool = False) -> float:
    """验证有限实标量；拒绝 bool、字符串、复数和数组，统一抛 ValueError。"""
    try:
        arr = np.asarray(value)
        if arr.ndim != 0 or arr.dtype.kind not in "iuf":
            raise ValueError("不是实数标量")
        number = float(arr)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} 必须是有限实数标量") from exc
    if not np.isfinite(number) or number < 0 or (number == 0 and not allow_zero):
        condition = "非负" if allow_zero else "正"
        raise ValueError(f"{name} 必须是有限{condition}实数")
    return number


def _positive(value: float, name: str) -> float:
    return _real_scalar(value, name)


def _finite_result(value: ArrayLike, name: str):
    """输入有限不保证运算结果有限；让溢出成为可定位的显式错误。"""
    if not np.all(np.isfinite(value)):
        raise ValueError(f"{name} 超出浮点数范围；请使用合适的无量纲单位")
    return value


def _scaled_state(arr: NDArray) -> NDArray:
    """按最大实/虚分量缩放，连 |z| 自身溢出的有限复数也可处理。"""
    scale = max(np.max(np.abs(arr.real)), np.max(np.abs(arr.imag)))
    if scale == 0:
        raise ValueError("不能归一化零态或零波函数")
    return arr.real / scale + 1j * (arr.imag / scale)


def _scaled_hamiltonian(matrix: NDArray, hbar: float) -> NDArray:
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        result = matrix.real / hbar + 1j * (matrix.imag / hbar)
    return _finite_result(result, "H/hbar")


def _grid(value: ArrayLike, name: str, minimum: int = 2) -> NDArray:
    arr = _array(value, name)
    if arr.ndim != 1 or arr.size < minimum or np.any(arr.imag != 0):
        raise ValueError(f"{name} 必须是至少含 {minimum} 个点的一维实数组")
    real = arr.real
    if np.any(real[1:] <= real[:-1]):
        raise ValueError(f"{name} 必须严格递增，不能重复或倒序")
    with np.errstate(over="ignore", invalid="ignore"):
        span = real[-1] - real[0]
    if not np.isfinite(span):
        raise ValueError(f"{name} 的区间跨度超出浮点数范围")
    return real


def _require_hermitian(matrix: NDArray) -> None:
    """物理求解用尺度无关检查，防止微小能量单位被绝对容差掩盖。"""
    scale = max(np.max(np.abs(matrix.real)), np.max(np.abs(matrix.imag)))
    if scale == 0:
        return
    scaled = matrix.real / scale + 1j * (matrix.imag / scale)
    if not np.allclose(scaled, scaled.conj().T, atol=1e-12, rtol=1e-10):
        raise ValueError("矩阵必须厄米；请检查共轭、单位与相对误差")


def dagger(a: ArrayLike) -> NDArray:
    """返回二维矩阵的厄米共轭 A† = (A*)ᵀ，允许非方阵。

    参数 a：shape=(M,N) 的数值矩阵。
    返回：shape=(N,M) 的复数组。不是普通转置；复矩阵必须先共轭。
    态矢的内积请使用 np.vdot(psi, phi)，它自动共轭第一个参数。
    """
    arr = _array(a, "a")
    if arr.ndim != 2:
        raise ValueError("dagger 接收二维矩阵；态矢内积请用 np.vdot")
    return arr.conj().T


def commutator(a: ArrayLike, b: ArrayLike) -> NDArray:
    """计算对易子 [A,B]=AB-BA；a、b 必须是同维方阵。

    用 @ 做矩阵乘法；* 是逐元素乘法，不能用来代替 AB。
    返回同维复数组，可用于验证 Pauli 矩阵对易关系。
    """
    a, b = _matrix(a, "a"), _matrix(b, "b")
    if a.shape != b.shape:
        raise ValueError("a 和 b 必须有相同形状")
    return a @ b - b @ a


def is_hermitian(a: ArrayLike, *, atol: float = 1e-12, rtol: float = 1e-10) -> bool:
    """按 np.allclose(A,A†,atol,rtol) 检查厄米性。

    非方阵或非有限数值会抛 ValueError；合法非厄米方阵返回 False。
    atol 是绝对容差，rtol 是相对容差，应按问题尺度选择。
    """
    arr = _matrix(a)
    absolute = _real_scalar(atol, "atol", allow_zero=True)
    relative = _real_scalar(rtol, "rtol", allow_zero=True)
    return bool(np.allclose(arr, arr.conj().T, atol=absolute, rtol=relative))


def normalize_state(psi: ArrayLike) -> NDArray:
    """归一化离散基态矢：返回 psi/sqrt(Σ|psi_i|²)，不修改输入。

    psi：shape=(N,)，不接受密度矩阵、批量态或连续波函数积分范数。
    零态没有可定义的归一化，抛 ValueError。采用尺度缩放减少溢出。
    """
    arr = _state(psi)
    scaled = _scaled_state(arr)
    return scaled / np.linalg.norm(scaled)


def integrate_samples(y: ArrayLike, x: ArrayLike, *, axis: int = -1, method: str = "simpson"):
    """积分采样数据 ∫y(x)dx（不是无限区间符号积分）。

    x：严格递增实网格，至少 2 点；允许非均匀网格。
    y：实/复数组；沿 axis 的长度必须等于 len(x)。例如 y.shape=(3,N)
       且 axis=-1 时，返回三个积分；复数实部/虚部都会保留。
    method：'simpson'（默认，光滑数据通常更准确）或 'trapezoid'。
    返回：标量或去掉积分轴后的数组。精度仍需靠加密网格检查。
    """
    grid = _grid(x, "x")
    values = _array(y, "y")
    if isinstance(axis, (bool, np.bool_)) or not isinstance(axis, (int, np.integer)) or values.ndim == 0:
        raise ValueError("axis 必须是有效整数轴，y 不能是标量")
    if not -values.ndim <= axis < values.ndim:
        raise ValueError("axis 超出 y 的维度")
    if values.shape[axis] != grid.size:
        raise ValueError("y 沿积分轴的长度必须等于 len(x)")
    if not isinstance(method, str) or method not in {"simpson", "trapezoid"}:
        raise ValueError("method 只能是 'simpson' 或 'trapezoid'")
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        result = (simpson(values, x=grid, axis=axis) if method == "simpson"
                  else np.trapezoid(values, x=grid, axis=axis))
    return _finite_result(result, "积分结果")


def normalize_wavefunction(psi: ArrayLike, x: ArrayLike, *, method: str = "simpson") -> NDArray:
    """按连续范数归一化采样波函数：psi / sqrt(∫|psi|²dx)。

    psi 与 x 均为一维等长数组；x 是实际坐标，不是数组下标。
    输出使指定积分规则得到的概率约为 1；有限区间遗漏的尾部仍需检查。
    不规则/欠采样网格的 Simpson 权重可能为负，不能保证概率积分为正；
    此时会拒绝非正范数，请加密网格或显式选 method='trapezoid'，不要取绝对值。
    注意：离散对角化器输出的本征矢通常满足 sum(|v|²)=1，作连续
    波函数图或积分期望值之前可用本函数转换；不能再用离散 expectation。
    """
    arr = _state(psi)
    scaled = _scaled_state(arr)
    norm2 = float(np.real(integrate_samples(np.abs(scaled)**2, x, method=method)))
    if not np.isfinite(norm2) or norm2 <= 0:
        raise ValueError("积分范数必须为有限正数；请检查采样网格")
    return scaled / np.sqrt(norm2)


def expectation(operator: ArrayLike, psi: ArrayLike) -> complex:
    """计算离散态期望值 <psi|A|psi>/<psi|psi>，接受未归一化非零态。

    operator：shape=(N,N)；psi：shape=(N,)；返回 complex。
    即使 A 厄米也不静默删除虚部，便于发现错误/舍入残差。
    显示实物理量时可用 np.real_if_close(result)，不要盲目取 .real。
    本函数不接收密度矩阵，也不替代连续波函数积分。
    """
    op = _matrix(operator, "operator")
    state = normalize_state(psi)
    if op.shape[0] != state.size:
        raise ValueError("operator 与 psi 维数不匹配")
    return complex(np.vdot(state, op @ state))


def diagonalize_hermitian(a: ArrayLike) -> tuple[NDArray, NDArray]:
    """厄米矩阵对角化，返回 (evals, evecs)，本征值按升序排列。

    evals.shape=(N,)，为实数；evecs.shape=(N,N)，第 j 列 evecs[:,j]
    对应 evals[j]（不是第 j 行）。使用 scipy.linalg.eigh。
    非厄米矩阵拒绝计算；该情形应另用 scipy.linalg.eig。
    本征矢整体相位任意；简并子空间内基也任意，不能逐元素硬比较。
    """
    arr = _matrix(a)
    _require_hermitian(arr)
    return eigh(arr)


def matrix_exp(a: ArrayLike) -> NDArray:
    """返回矩阵指数 exp(A)=I+A+A²/2!+…，使用 scipy.linalg.expm。

    a：任意非空有限方阵（也可非厄米）。返回同维复数组。
    np.exp(A) 是逐元素指数，不是矩阵指数；二者一般完全不同。
    大型稀疏矩阵宜另用 scipy.sparse.linalg.expm_multiply，本函数是稠密版。
    """
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        result = expm(_matrix(a))
    return _finite_result(result, "矩阵指数")


def evolve_constant(h: ArrayLike, psi0: ArrayLike, times: ArrayLike, *, hbar: float = 1.0) -> NDArray:
    """常量厄米 H 的精确矩阵指数演化，输出 shape=(Nt,N)。

    psi0 是 times[0] 时刻的态；第 j 行等于
    exp[-i*H*(times[j]-times[0])/hbar] @ psi0。
    times 严格递增，可只含 1 点；hbar 为正实数，默认自然单位 1。
    不自动归一化输入/输出：这样不会隐藏范数漂移或初态错误。
    每个时刻独立调用 expm，适合小矩阵学习，不适合大规模生产计算。
    """
    matrix, initial = _matrix(h, "h"), _state(psi0)
    grid, hb = _grid(times, "times", minimum=1), _positive(hbar, "hbar")
    if matrix.shape[0] != initial.size:
        raise ValueError("h 必须与 psi0 同维")
    _require_hermitian(matrix)
    if not np.any(initial != 0):
        raise ValueError("psi0 不能是零态")
    # 只有初始时刻时无需形成 H/hbar，直接返回独立副本。
    if grid.size == 1:
        return initial[np.newaxis, :].copy()
    scaled_h = _scaled_hamiltonian(matrix, hb)
    with np.errstate(over="ignore", invalid="ignore"):
        states = np.stack([
            matrix_exp(-1j * scaled_h * (t-grid[0])) @ initial for t in grid
        ])
    return _finite_result(states, "演化结果")


def solve_schrodinger(
    h: ArrayLike | Callable[[float], ArrayLike], psi0: ArrayLike, times: ArrayLike,
    *, hbar: float = 1.0, rtol: float = 1e-9, atol: float = 1e-11,
    method: str = "DOP853", max_step: float = np.inf,
) -> NDArray:
    """用 solve_ivp 解 dpsi/dt=-i H(t)psi/hbar，返回 shape=(Nt,N)。

    h：常量厄米矩阵，或只接收时间 t 并返回矩阵的函数 h(t)。
    psi0：times[0] 时刻的一维非零态；times：至少两个严格递增时间。
    hbar：正实数；rtol/atol：正实求解容差，不是采样时间步长。
    method：'RK23'、'RK45'、'DOP853' 或 'BDF'，均支持复数初态。
    这些方法不是保结构幺正积分器，务必检查范数并做容差收敛测试。
    时间采样 times 只指定输出点，求解器自行选择内部步长。
    max_step：内部步长上限，正实数或 np.inf（默认）；窄脉冲需按脉宽
    设置，例如脉宽/10，否则再密的输出点或更紧容差也可能漏采脉冲。
    每次回调验证 H 的维数/有限性/厄米性；不用于非厄米耗散模型。
    求解失败抛 RuntimeError；输入不合法抛 ValueError；不重新归一化。
    """
    grid, initial = _grid(times, "times"), _state(psi0)
    hb = _positive(hbar, "hbar")
    rel, absolute = _positive(rtol, "rtol"), _positive(atol, "atol")
    if not isinstance(method, str) or method not in {"RK23", "RK45", "DOP853", "BDF"}:
        raise ValueError("请选择支持复数的 RK23/RK45/DOP853/BDF")
    try:
        raw_step = np.asarray(max_step)
        unlimited = raw_step.ndim == 0 and raw_step.dtype.kind in "iuf" and float(raw_step) == np.inf
    except (TypeError, ValueError, OverflowError):
        unlimited = False
    step = np.inf if unlimited else _positive(max_step, "max_step")
    if not np.any(initial != 0):
        raise ValueError("psi0 不能是零态")

    def checked_matrix(t: float) -> NDArray:
        matrix = _matrix(h(t) if callable(h) else h, "h(t)")
        if matrix.shape != (initial.size, initial.size):
            raise ValueError("每个 h(t) 都必须与 psi0 同维")
        _require_hermitian(matrix)
        return _scaled_hamiltonian(matrix, hb)

    # 常量矩阵只验证一次；时间依赖矩阵必须在每次求导时验证。
    constant = None if callable(h) else checked_matrix(float(grid[0]))

    def rhs(t: float, psi: NDArray) -> NDArray:
        matrix = checked_matrix(t) if constant is None else constant
        with np.errstate(over="ignore", invalid="ignore"):
            derivative = -1j * (matrix @ psi)
        if not np.all(np.isfinite(derivative)):
            raise ValueError("求导超出浮点数范围；请改用合适的无量纲单位")
        return derivative

    result = solve_ivp(rhs, (grid[0], grid[-1]), initial, t_eval=grid,
                       method=method, rtol=rel, atol=absolute, max_step=step)
    if not result.success:
        raise RuntimeError(f"solve_ivp 失败：{result.message}")
    if not np.all(np.isfinite(result.y)):
        raise RuntimeError("solve_ivp 产生非有限数值")
    # SciPy 是 (N,Nt)，转置后统一为每一行一个时刻。
    return result.y.T


def plot_style(*, font_size: float = 11, dpi: float = 120):
    """返回 Matplotlib 局部 rc_context；用 with plot_style(): 包住绘图。

    font_size、dpi 为正实数；不永久改全局 rcParams，不依赖 LaTeX。
    默认不强设中文字体（各机器字体不同）；示例使用英文标签避免缺字。
    中文图可在此上下文内设置 plt.rcParams['font.sans-serif'] 为已安装字体。
    """
    import matplotlib as mpl  # 延迟导入：纯数值调用不要求加载绘图库。
    size, resolution = _positive(font_size, "font_size"), _positive(dpi, "dpi")
    return mpl.rc_context({
        "font.size": size, "figure.dpi": resolution,
        "lines.linewidth": 1.8, "axes.grid": True,
        "grid.alpha": 0.25, "axes.spines.top": False,
        "axes.spines.right": False, "savefig.bbox": "tight",
    })


def new_figure(*, xlabel: str = "", ylabel: str = "", title: str = "",
               figsize: tuple[float, float] = (6.4, 4.0)):
    """创建单面板 (fig,ax)，设置标签/标题与 constrained_layout。

    figsize 是英寸宽高；图中坐标请注明物理量与单位，例如 'Time (s)'。
    不调用 show()；由调用者决定显示、保存或 plt.close(fig)。
    """
    import matplotlib.pyplot as plt
    try:
        valid_size = isinstance(figsize, (tuple, list, np.ndarray)) and len(figsize) == 2
    except TypeError:
        valid_size = False
    if not valid_size:
        raise ValueError("figsize 必须是 (宽,高)")
    width, height = (_positive(v, "figsize") for v in figsize)
    fig, ax = plt.subplots(figsize=(width, height), layout="constrained")
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
    return fig, ax


def save_figure(fig, path: str | Path, *, dpi: float = 300, overwrite: bool = False) -> Path:
    """保存指定 Figure，返回最终 Path；支持 .png/.pdf/.svg。

    path 相对当前进程工作目录，不是 helpers.py 所在目录。
    自动建立父目录；默认拒绝覆盖已有文件（overwrite=True 才允许）。
    PNG 的 dpi 影响像素数；PDF/SVG 中线条/文字通常为矢量。
    先写同目录临时文件再原子提交，保存失败不破坏已有图；禁止覆盖时
    Windows 使用原生 rename 的不覆盖语义；其他平台使用硬链接原子创建。
    非 Windows 且不支持硬链接时保留原生 OSError。
    本函数不关闭图、不调用 show；批量画图后调用 plt.close(fig) 释放内存。
    """
    target = Path(path)
    resolution = _positive(dpi, "dpi")
    if target.suffix.lower() not in {".png", ".pdf", ".svg"}:
        raise ValueError("输出扩展名必须是 .png/.pdf/.svg")
    if target.exists() and not overwrite:
        raise FileExistsError(f"文件已存在，拒绝覆盖：{target}")
    if not isinstance(overwrite, bool):
        raise ValueError("overwrite 必须是 bool")
    target.parent.mkdir(parents=True, exist_ok=True)
    # 同目录临时文件保证提交在同一文件系统；先完整写成，失败不碰原图。
    import os
    import tempfile
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{target.stem}-", suffix=target.suffix.lower(), dir=target.parent,
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        fig.savefig(temporary, dpi=resolution, bbox_inches="tight",
                    format=target.suffix[1:].lower())
        if overwrite:
            os.replace(temporary, target)
        else:
            # Windows rename 原子移动且拒绝已存在目标；工作区可能禁止硬链接。
            # POSIX rename 会覆盖目标，所以其他平台仍需 link 的不覆盖语义。
            if os.name == 'nt':
                os.rename(temporary, target)
            else:
                os.link(temporary, target)
    finally:
        temporary.unlink(missing_ok=True)
    return target
