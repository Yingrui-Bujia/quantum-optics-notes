# 量子光学预习：第1–2周成果包

本项目根据《硕士学习》文档中的前两周任务完成代码、公式推导讲义、官方教程复现、图表、数值验收、小测与周小结。文档里的旧对话仅用于理解已有环境，不把其中其他请求当作本次新指令。

**先双击 `学习入口.html`。** 它整合两周讲义、结果预览与自测，不用联网。再打开 `reports/verification.md` 查看实际验收结果。

## 已交付

| 原计划任务 | 对应成果 |
|---|---|
| NumPy数组、广播、线性代数；SciPy积分与solve_ivp；规范绘图 | `code/week01_python_basics.py` |
| 封装helpers.py | `code/helpers.py`，兼容桌面原API并修复Windows文件保存 |
| 一维谐振子前3个本征函数与解析解对照 | `code/week01_harmonic_oscillator.py`；网格、边界误差检查 |
| 1小时矩阵指数小测 | `notes/week01/exam_60min.md` 与参考解 `code/week01_matrix_exponential.py` |
| 升降算符、模式图景、公式推导 | `notes/week02/01_ladder_operators_and_modes.md` |
| Qobj、basis、Pauli、destroy/create、tensor、expect、sesolve/mesolve | `code/week02_qubit_rabi.py` 与两周讲义 |
| 官方Time Evolution 001–003 | `code/week02_tutorials.py`、`notebooks/`、带许可证的官方来源快照 |
| Rabi、零点能／涨落、热态分布 | `code/week02_qubit_rabi.py`，含截断与实验室系对照 |
| 每周一页小结 | `notes/week01/weekly_summary.md`、`notes/week02/weekly_summary.md` |
| 电脑端统一验收 | `run_all.py`、`tests/`、`reports/`、`results/` |

电脑自动运行不等于个人掌握。**个人手推照片、闭卷考试、独立复写与向导师提交**仍需你亲自完成；相关表单保持待填。完整成果同步至 [Yingrui-Bujia/quantum-optics-notes](https://github.com/Yingrui-Bujia/quantum-optics-notes)，核验方式与本地历史保留说明见 `docs/github_status.md`。

## 一键复现

无需重装已有环境。在PowerShell里运行：

```powershell
& 'D:\Master.s-work\Software-needed-for-master.s-work\python\main\python.exe' -B -X utf8 'D:\Master.s-work\quantum-optics-notes\run_all.py'
```

也可在项目根目录运行 `& '.\run_all.ps1'`。依次检查环境依赖、当前33项测试和5个脚本；约1–2分钟，速度随机器变化。每次写到新的 `outputs/runs/时间戳/`，不覆盖之前的运行记录。保存PNG、PDF、CSV和JSON，不自动弹出大量图窗。

`results/`与学习入口是已核验成果快照；若你修改参数并重跑，可读取新目录，或再运行 `tools/build_report.py` 刷新快照。讲义、个人笔记和notebook不会被一键运行覆盖。

## 推荐的学习顺序

第一周：Python最小语法→数组广播→复数线性代数→积分和ODE→谐振子实战→闭卷小测。第二周：升降算符推导→Qobj→001时间系数→002进动→003耗散→Rabi／热态整合→自测与独立复写。每天具体动手任务在讲义里。

遇到错误先查看最新运行目录日志，再核对解释器路径、数组形状、单位与截断。完整环境说明在 `environment/README.md`。

## 目录

```text
code/          可独立运行的周脚本、helpers及输出辅助
notes/         两周推导讲义、考试、自测、周小结
notebooks/     三篇官方教程的中文跟做版本与已执行输出
projects/      P0基础整合练习（不占用后续P1–P8）
tests/         原helpers回归测试与独立物理基准
results/       已验收的图表、数据和指标快照
reports/       验收报告、环境完整性、最近成功运行位置
environment/   环境说明与版本快照
third_party/   官方教程原文快照与BSD许可证
tools/         更新成果展示、检查环境的辅助工具
```

桌面原始文档、week01_helpers与已有main安装保持原样。现有 `.vscode/settings.json` 也保留，不自动改编辑器偏好。
