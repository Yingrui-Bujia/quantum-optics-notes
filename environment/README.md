# 环境使用与保护

本机解释器：`D:\Master.s-work\Software-needed-for-master.s-work\python\main\python.exe`。

本次直接使用已有main，没有安装、升级、降级或删除包，没有改系统PATH、镜像、Git全局配置及原VS Code设置。开始前与结束后对比完整已安装包的名称和版本；同时核验桌面原始helpers、任务文档、环境history及已有配置的SHA-256。结果见 `reports/environment_integrity.json`。

`installed_versions.json`是当前机器的版本快照；`requirements.txt`给出项目核心依赖范围，`requirements-reference.txt`记录实测核心版本。这些文件用于了解环境或将来在新环境恢复，当前机器无需执行安装命令。Python科学计算工具的完整Conda构建不能仅靠pip版本清单逐字节复原。

## 常用启动

在项目根目录右键 → 在终端打开，运行：

```powershell
& '.\run_all.ps1'
```

若终端执行策略限制脚本，直接调用Python，不修改执行策略：

```powershell
& 'D:\Master.s-work\Software-needed-for-master.s-work\python\main\python.exe' -B -X utf8 '.\run_all.py'
```

换机器时，把 `run_all.ps1` 的 `-PythonPath` 参数指向自己的现有解释器。脚本不负责装包、不调用网络、不改变用户安装。

## 输出与临时文件

每次生成独立 `outputs/runs/时间戳/`，成功才更新 `reports/latest_run.json`。失败日志保留。测试临时目录放在这次结果内部，避免默认临时目录在受限工作区出现清理错误。图片由指定Figure保存，Windows无覆盖保存采用原生rename，仍保留文件存在和竞态保护。

QuTiP字符串系数采用 `CompilationOptions(use_cython=False)`；原版本会发出一次回退到解释执行的提示，这不是数值失败。无需为了去掉提示安装编译器。

`notebooks/`提供已逐格执行的学习笔记。VS Code打开后仍应选main内核；若换了内核，先运行第一格确认 `sys.executable`。本次未注册或修改全局Jupyter内核。
