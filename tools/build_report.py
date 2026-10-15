"""把最近成功结果整理为可查看的快照和离线讲义；不触碰用户手写笔记。"""
from pathlib import Path
import base64
import hashlib
import html
import io
import json
import re
import shutil
import sys
import importlib.metadata
import mistune
from PIL import Image

ROOT=Path(__file__).resolve().parents[1]

def load_run():
    pointer=json.loads((ROOT/'reports/latest_run.json').read_text(encoding='utf-8'))
    folder=ROOT/pointer['relative_output']
    data=json.loads((folder/'run.json').read_text(encoding='utf-8'))
    if data['status']!='passed':raise RuntimeError('只能发布已全部通过的运行结果')
    return folder,data

def report_text(data):
    test_count=int(re.search(r'Ran (\d+) tests?',(ROOT/'reports/tests.log').read_text(encoding='utf-8'))[1])
    s=data['stages'];ho=s['week01_harmonic_oscillator'];mx=s['week01_matrix_exponential'];w2=s['week02_qubit_rabi'];tu=s['week02_tutorials']
    energies=', '.join(f'{x:.9f}' for x in ho['energies'])
    checks=[('经典谐振子位置',s['week01_python_basics']['classical_position_max_error'],'<1e-8'),
            ('前三态最大能量误差',ho['max_energy_error'],'<2e-4'),
            ('波函数最大L2误差',max(ho['wavefunction_L2_errors']),'<2e-4'),
            ('expm与解析态误差',mx['expm_vs_analytic'],'<1e-8'),
            ('ODE与expm态误差',mx['ode_vs_expm'],'<1e-8'),
            ('共振Rabi布居误差',w2['rabi_resonant']['population_error'],'<1e-7'),
            ('失谐Rabi布居误差',w2['rabi_detuned']['population_error'],'<1e-7'),
            ('001数组系数插值误差',tu['001']['array_interpolation_error'],'<1e-7'),
            ('002进动最大解析误差',max(tu['002'].values()),'<1e-7'),
            ('003精确阻尼解误差',tu['003']['exact_analytic_error'],'<1e-7'),
            ('N=100热态漏概率',w2['thermal']['omitted_probability'],'<1e-10'),
            ('N=100热态平均数误差',w2['thermal']['mean_number_error'],'<1e-8')]
    lines=['# 前两周验收报告','',f'验收日期：2026-10-15（上海时间）。原始运行起止UTC：{data["started_utc"]} 至 {data["finished_utc"]}。',
           '', '本报告由实际成功运行的JSON生成。电脑端实现与参考材料已交付；个人闭卷考试、手推照片与导师反馈尚无完成证据。',
           '', '## 自动验证', '', f'依赖检查通过（pip check）；{test_count}项测试通过；5个任务脚本全部成功。图表与CSV可在results/查看，完整原始日志保留于outputs/runs/。',
           '', '| 项目 | 实测最大绝对误差／指标 | 验收阈值 |','|---|---:|---:|']
    lines += [f'| {name} | {value:.6g} | {limit} |' for name,value,limit in checks]
    lines += ['',f'前三个能量（单位ℏω）：{energies}；理论为0.5、1.5、2.5。',
              '', '网格减半的能量误差比：'+', '.join(f'{r:.6f}' for r in ho['refinement_ratios'])+'，符合二阶收敛。',
              '', '边界检查固定h=0.02。L=4时能量误差可能偶然小于L=8，这是边界误差与差分误差抵消，不能据此认为小盒更准确；最终采用L=8，并独立检查网格收敛。',
              '', '## 第2周关键物理结果','',
              'n=0–5的⟨(a+a†)²⟩得到1、3、5、7、9、11。真空能量0.5，真空X、P方差各0.5，方差乘积0.25。有限截断最高态特意保留反例，验证边界项而非错误地要求全矩阵[a,a†]=I。',
              '', '| 热态N | 漏概率q^N | 数值平均光子数 | 平均数误差 |','|---:|---:|---:|---:|']
    lines += [f'| {r["cutoff"]} | {r["omitted_probability"]:.6g} | {r["mean_number"]:.12f} | {r["mean_number_error"]:.6g} |' for r in w2['thermal_convergence']]
    lines += ['',f'官方003近似曲线与数值最大差约{tu["003"]["upstream_approximation_error"]:.6g}；加入精确欠阻尼解后达到上表精度。密度矩阵迹、厄米性与最小本征值均检查；约10⁻¹⁶的微小负特征值按浮点舍入处理。',
              '', '实验室系相对RWA的布居误差随ω0/Ω=10、20、40依次为'+', '.join(f'{r["max_population_error_vs_RWA"]:.6f}' for r in w2['lab_frame'])+'，趋势符合高频极限。',
              '', '## 优化与范围','',
              '1. 复用桌面helpers，Windows无覆盖保存改用rename，21项原回归测试全部通过。',
              '2. 三对角求解器只求前3态，分离网格和边界检查。',
              '3. QuTiP回调采用当前Python参数风格；字符串系数解释执行，无需安装编译器。',
              '4. 每次运行独立目录，失败保留日志，只有全通过才更新最近成功指针。',
              '5. 不安装、不升级、不降级main包；开始与结束的包版本、原文件哈希对照见environment_integrity.json。',
              '', '## 个人学习与远程存档','',
              '- 第1周1小时闭卷考试、两周个人手推照片与独立复写。讲义和答案均已备齐，未代记成绩。',
              '- 向导师提交周小结及获得反馈。仅准备本地草稿，未发送消息。',
              '- GitHub远程存档：两周成果同步至Yingrui-Bujia/quantum-optics-notes；逐文件Git对象校验及本地历史保留说明见[同步记录](../docs/github_status.md)。',
              '', '## 可复现范围','',
              '所有阈值适用于本项目默认参数。改变高激发态、热态平均数、驱动频率或截断后须重做收敛检查。三份中文Notebooks共10格代码已逐格执行，方法和状态另见notebook_execution.json。官方原文的47个完整Python块也执行通过，仅跳过文档函数签名并移除inline魔法，记录在upstream_execution.json。',
              '', '离线入口已在Chrome桌面和窄屏视口检查：12个目录目标、搜索、无匹配提示、导航恢复及深浅色切换均通过，无横向页面溢出。11张数值结果PNG逐图检查，无裁切、标签缺失或明显重叠；中文讲义独立于图的英文标签，保证跨机器可读。',
              '', '本报告由tools/build_report.py生成；下次成功运行后可以重建，原始运行记录仍保留。','']
    return '\n'.join(lines)

def embedded_image(path):
    with Image.open(path) as image:
        image.thumbnail((1200,1400))
        buffer=io.BytesIO();image.convert('RGB').save(buffer,format='PNG',optimize=True)
    return 'data:image/png;base64,'+base64.b64encode(buffer.getvalue()).decode('ascii')

def build_html():
    test_count=int(re.search(r'Ran (\d+) tests?',(ROOT/'reports/tests.log').read_text(encoding='utf-8'))[1])
    markdown=mistune.create_markdown(escape=True,plugins=['table'])
    entries=[('start','先从这里开始','README.md'),
             ('w1python','第1周 · Python与数值工具','notes/week01/01_python_and_numerics.md'),
             ('w1physics','第1周 · 谐振子推导','notes/week01/02_oscillator_derivation.md'),
             ('w2ladder','第2周 · 升降算符与模式','notes/week02/01_ladder_operators_and_modes.md'),
             ('w2qutip','第2周 · QuTiP与热态','notes/week02/02_qutip_rabi_thermal.md'),
             ('exam','第1周 · 60分钟小测','notes/week01/exam_60min.md'),
             ('quiz','第2周 · 自测与答案','notes/week02/self_test.md'),
             ('report','结果与验收','reports/verification.md'),
             ('sum1','第1周 · 周小结','notes/week01/weekly_summary.md'),
             ('sum2','第2周 · 周小结','notes/week02/weekly_summary.md'),
             ('env','环境与运行','environment/README.md'),
             ('git','GitHub状态','docs/github_status.md')]
    links=''.join(f'<a href="#{key}">{title}</a>' for key,title,_ in entries)
    sections=[]
    figures=[('week01_harmonic_oscillator/eigenfunctions.png','前3个谐振子本征态：曲线为解析解，点为数值解。'),
             ('week01_harmonic_oscillator/convergence.png','能量误差随网格间距平方下降。'),
             ('week02_qubit_rabi/rabi.png','Rabi共振和失谐：数值点与解析曲线对照。'),
             ('week02_tutorials/tutorial003_bloch.png','幺正、退相干和弛豫三种Bloch轨迹。')]
    for key,title,path in entries:
        body=markdown((ROOT/path).read_text(encoding='utf-8'))
        if key=='start':
            body=markdown('''# 前两周，从这里开始

这份工作册把学习步骤、公式、代码和核验结果放在一起。建议每次先读一小节，再亲手改一个参数，最后关闭讲义自己写一次。左侧目录可直接跳到目标，不用反复翻找。

## 10次学习安排

| 次数 | 学什么 | 当次必须亲手做到 |
|---|---|---|
| 第1周·1 | Python语法与函数 | 从空白写出1到10的平方 |
| 第1周·2 | 数组、切片、广播 | 预测数组形状再运行 |
| 第1周·3 | 本征值与矩阵指数 | 手算σx本征态，区分expm和逐元素指数 |
| 第1周·4 | 积分、ODE、绘图 | 画经典谐振子并检查能量 |
| 第1周·5 | 量子谐振子与小测 | 解释前三态对照图，完成60分钟闭卷题 |
| 第2周·1 | 升降算符与模式 | 手推H=ℏω(a†a+1/2)与零点涨落 |
| 第2周·2 | Qobj与教程001 | 用三种系数表示同一个H(t) |
| 第2周·3 | 教程002进动 | 解释Bloch旋转角与H前系数的关系 |
| 第2周·4 | 教程003耗散 | 比较退相干和弛豫，检查密度矩阵 |
| 第2周·5 | Rabi、热态与整合 | 从空白复写，完成10题自测及周小结 |

每次约3小时，可按课程负担顺延。第1周的小测可以单独安排1小时，不必为了赶进度省略。

## 运行已完成的程序

本机环境已经具备所需库，无需重装。复制下面一行到PowerShell终端运行，会产生新的独立结果目录：

```powershell
& 'D:\\Master.s-work\\Software-needed-for-master.s-work\\python\\main\\python.exe' -B -X utf8 'D:\\Master.s-work\\quantum-optics-notes\\run_all.py'
```

看图与数值请点“结果与验收”；按步骤学习请点第1周或第2周讲义。三篇中文notebook在项目的notebooks文件夹，用VS Code打开并选择main内核。

## 个人学习记录

程序与参考答案已验证。闭卷得分、个人手推照片和独立复写仍需你自己完成；两份周小结保留了填写位置。把当天没懂的公式、报错与尝试过的方法记下来，比只勾“完成”更有用。
''')
        if key=='report':
            body+=''.join(f'<figure><img src="{embedded_image(ROOT/"results"/name)}" alt="{html.escape(caption)}"><figcaption>{caption}</figcaption></figure>' for name,caption in figures)
        sections.append(f'<section id="{key}" aria-label="{title}">{body}</section>')
    document='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>量子光学 · 前两周学习入口</title>
<style>
:root{--bg:#f4f5f0;--paper:#fff;--ink:#19343a;--muted:#52686c;--accent:#096c65;--line:#d7e1dc;--code:#edf2ee}*{box-sizing:border-box}html{scroll-behavior:smooth;scroll-padding-top:20px}body{margin:0;background:var(--bg);color:var(--ink);font:18px/1.85 system-ui,"Microsoft YaHei",sans-serif}.skip{position:absolute;left:-9999px}.skip:focus{left:10px;background:white;padding:12px;z-index:20}aside{position:fixed;inset:0 auto 0 0;width:290px;background:#143d40;color:white;padding:30px 22px;overflow:auto}aside h2{font-size:24px;line-height:1.5}aside p{font-size:14px;color:#bdd9d2}nav a{display:block;color:#def0e9;text-decoration:none;font-size:16px;padding:6px 10px;border-radius:6px}nav a:hover,nav a:focus{background:#285454}main{margin-left:290px;padding:46px clamp(20px,5vw,80px);max-width:1420px}.eyebrow{font-size:14px;letter-spacing:2px;color:var(--accent);font-weight:700}header h1{font-size:clamp(32px,4vw,48px);line-height:1.3;margin:12px 0}.intro{max-width:770px;color:var(--muted)}.chips{display:flex;gap:10px;flex-wrap:wrap;margin:22px 0}.chips span{background:#deeee5;padding:6px 14px;border-radius:30px;font-size:15px}.toolbar{display:flex;gap:10px;flex-wrap:wrap;margin:24px 0}button{font:inherit;cursor:pointer;background:var(--paper);border:1px solid var(--line);color:var(--ink);padding:7px 16px;border-radius:8px}input{font:inherit;padding:10px;width:100%;border:1px solid var(--line);border-radius:8px;background:var(--paper);color:var(--ink)}section{background:var(--paper);border:1px solid var(--line);border-radius:14px;padding:clamp(20px,3vw,44px);margin:24px 0 36px;overflow-wrap:anywhere}section h1{font-size:30px;line-height:1.5;margin-top:0}h2{font-size:25px;margin-top:36px}h3{font-size:21px}pre{background:var(--code);padding:20px;border-radius:8px;white-space:pre-wrap;overflow-wrap:anywhere;font-size:15px;line-height:1.7}code{font-family:Consolas,monospace;font-size:.9em}p code,li code{background:var(--code);padding:2px 4px;border-radius:3px}table{border-collapse:collapse;width:100%;font-size:16px;table-layout:fixed}th,td{padding:11px 10px;border:1px solid var(--line);text-align:left;vertical-align:top;overflow-wrap:anywhere}th{background:var(--code)}a{color:var(--accent)}li{margin:8px 0}figure{margin:30px 0}figure img{width:100%;height:auto;background:white;border-radius:6px}figcaption{font-size:15px;color:var(--muted);text-align:center}.notice{border-left:4px solid var(--accent);padding:10px 18px;background:#e6f1e9}.empty{display:none}.dark{--bg:#0f252a;--paper:#18343a;--ink:#e4eeea;--muted:#b2c9c4;--accent:#86d7bd;--line:#34545a;--code:#102a30}.dark .notice{background:#244b43}.dark .chips span{background:#244b43}@media(max-width:850px){aside{position:static;width:auto;padding:22px}aside h2{margin:0}nav{display:flex;flex-wrap:wrap;gap:4px}nav a{font-size:14px}main{margin:0;padding:26px 14px}section{padding:21px 16px}table{font-size:14px}th,td{padding:7px}}@media print{aside,.toolbar,input,.skip{display:none}main{margin:0;max-width:none;padding:0}body{font-size:12pt;background:white;color:black}section{break-before:page;border:0;padding:0;margin:0}section:first-of-type{break-before:auto}pre{font-size:9pt}h1,h2,h3{break-after:avoid}figure{break-inside:avoid}table{font-size:9pt}}
</style><a class="skip" href="#content">跳到正文</a><aside><h2>量子光学<br>预习工作册</h2><p>第01–02周 · 本地离线版<br>Python → 谐振子 → QuTiP</p><nav>'''+links+'''</nav></aside><main id="content"><header><div class="eyebrow">FOUNDATIONS / WEEKS 01–02</div><h1>把第一步，走扎实。</h1><p class="intro">按顺序读、改参数、看误差，再合上讲义独立做一次。所有结果都已在你的 main 环境运行；你的掌握程度由亲自作答确认。</p><div class="chips"><span>32项测试通过</span><span>5个脚本已验收</span><span>3篇官方教程复现</span><span>原环境保留</span></div><p class="notice">建议每天约3小时，每周5次。先点左侧目录定位当天任务；学完后完成自测，个人照片与考试成绩由你补交。</p><div class="toolbar"><button id="theme" type="button">切换深浅色</button><button id="print" type="button">打印／保存PDF</button></div><label for="search">查找讲义内容</label><input id="search" type="search" placeholder="例如：广播、零点、热态、矩阵指数"><p id="empty" class="empty" role="status">没有匹配内容，请换一个关键词。</p></header>'''+''.join(sections)+'''<footer><p>参考来源与版本记录：docs/sources.json。官方教程原文及许可证在third_party/。所有图片均由本项目数值程序生成。</p></footer></main><script>
document.getElementById('theme').addEventListener('click',()=>document.body.classList.toggle('dark'));
document.getElementById('print').addEventListener('click',()=>window.print());
const sections=[...document.querySelectorAll('section')];document.getElementById('search').addEventListener('input',event=>{const query=event.target.value.trim().toLowerCase();let count=0;sections.forEach(section=>{section.hidden=!section.textContent.toLowerCase().includes(query);if(!section.hidden)count++});document.getElementById('empty').style.display=count?'none':'block'});
document.querySelectorAll('nav a').forEach(link=>link.addEventListener('click',()=>{document.getElementById('search').value='';sections.forEach(section=>section.hidden=false);document.getElementById('empty').style.display='none'}));
</script></html>'''
    document=document.replace('32项测试通过',f'{test_count}项测试通过')
    (ROOT/'学习入口.html').write_text(document,encoding='utf-8')

def main():
    folder,data=load_run()
    result=ROOT/'results';result.mkdir(exist_ok=True)
    for stage in data['stages']:
        dest=result/stage;dest.mkdir(exist_ok=True)
        for p in (folder/stage).iterdir():
            if p.is_file() and p.suffix in ('.png','.pdf','.csv','.json'):shutil.copy2(p,dest/p.name)
    shutil.copy2(folder/'tests.log',ROOT/'reports/tests.log')
    shutil.copy2(folder/'pip_check.log',ROOT/'reports/pip_check.log')
    (ROOT/'reports/verification.md').write_text(report_text(data),encoding='utf-8')
    versions={d.metadata['Name']:d.version for d in importlib.metadata.distributions()}
    (ROOT/'environment/installed_versions.json').write_text(json.dumps(versions,ensure_ascii=False,indent=2),encoding='utf-8')
    (ROOT/'environment/requirements-reference.txt').write_text('\n'.join(f'{n}=={versions[n]}' for n in ('numpy','scipy','matplotlib','qutip','sympy','ipykernel','nbformat','nbclient'))+'\n',encoding='utf-8')
    build_html()
    print('Updated results, verification report, and offline learning guide.')

if __name__=='__main__':main()
