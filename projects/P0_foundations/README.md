# P0：谐振子与二能级的双工具核验

这是前两周的基础整合练习，方便向导师展示“模型→推导→数值→误差”的完整过程。它是额外的P0，不占用原计划第5周开始的P1–P8编号。

问题A：在坐标网格上解谐振子。比较前三态与解析Hermite函数，解释步长误差和有限边界误差。入口 `code/week01_harmonic_oscillator.py`。

问题B：对H=σx/2，用SciPy的expm、solve_ivp及QuTiP的sesolve求同一个初态的演化。比较P1(t)和态范数。入口 `code/week01_matrix_exponential.py` 与 `code/week02_qubit_rabi.py`。

问题C：有限截断是否会让物理公式失效？同时给出最高层对易子反例与nbar=3的热态收敛表，说明什么时候可以信任数值结果。

成果在 `results/`，验收指标在 `reports/verification.md`，一键重现使用项目根目录 `run_all.ps1`。展示时先讲单位和基态约定，再放图，不要只朗读程序。
