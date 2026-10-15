"""一键复现前两周；仅使用当前Python，不安装、升级或删除环境内容。"""
from pathlib import Path
from datetime import datetime, timezone
import os
import sys
import json
import subprocess
import importlib.metadata

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'code'))
# 子进程继承配置；绘图库缓存与测试临时文件留在本次输出目录。

def main():
    output=ROOT/'outputs'/'runs'/datetime.now().strftime('%Y%m%d-%H%M%S-%f')
    output.mkdir(parents=True,exist_ok=False)
    temporary=output/'cache'/'tmp'
    temporary.mkdir(parents=True)
    for key in ('TMP','TEMP','TMPDIR'):
        os.environ[key]=str(temporary)
    os.environ['MPLCONFIGDIR']=str(output/'cache'/'matplotlib')
    os.environ['PYTHONDONTWRITEBYTECODE']='1'
    os.environ['MPLBACKEND']='Agg'
    os.environ['PYTHONIOENCODING']='utf-8'
    env=os.environ.copy();env['PYTHONPATH']=str(ROOT/'code')
    packages={d.metadata['Name']:d.version for d in importlib.metadata.distributions()}
    result={'started_utc':datetime.now(timezone.utc).isoformat(),'executable':sys.executable,
            'python':sys.version,'packages':packages,'output':str(output),'stages':{},'status':'running'}
    report=output/'run.json'
    try:
        pip_check=subprocess.run([sys.executable,'-m','pip','check'],cwd=ROOT,env=env,capture_output=True,text=True,encoding='utf-8')
        (output/'pip_check.log').write_text(pip_check.stdout+pip_check.stderr,encoding='utf-8')
        if pip_check.returncode:raise RuntimeError('现有环境依赖检查未通过，详见pip_check.log')
        with (output/'tests.log').open('w',encoding='utf-8') as log:
            test=subprocess.run([sys.executable,'-B','-m','unittest','discover','-s','tests','-v'],cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=180)
        result['tests']={'returncode':test.returncode,'log':'tests.log'}
        if test.returncode:raise RuntimeError('测试未通过，详见tests.log')
        stages=['week01_python_basics','week01_harmonic_oscillator','week01_matrix_exponential','week02_tutorials','week02_qubit_rabi']
        for name in stages:
            print('Running:',name,flush=True)
            folder=output/name
            execution=subprocess.run([sys.executable,'-B',str(ROOT/'code'/f'{name}.py'),'--output',str(folder)],cwd=ROOT,env=env,capture_output=True,text=True,encoding='utf-8')
            (output/f'{name}.log').write_text(execution.stdout+execution.stderr,encoding='utf-8')
            if execution.returncode:raise RuntimeError(f'{name}运行失败，详见日志')
            result['stages'][name]=json.loads((folder/'metrics.json').read_text(encoding='utf-8'))
        result['status']='passed'
        # 最近成功结果只在所有阶段通过后更新；失败运行不取代成功结果。
        (ROOT/'reports'/'latest_run.json').write_text(json.dumps({'output':str(output),'relative_output':output.relative_to(ROOT).as_posix()},ensure_ascii=False,indent=2),encoding='utf-8')
        print('All checks passed:',output,flush=True)
    except Exception as exc:
        result['status']='failed';result['error']=str(exc)
        print(str(exc),file=sys.stderr)
    finally:
        result['finished_utc']=datetime.now(timezone.utc).isoformat()
        report.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    return 0 if result['status']=='passed' else 1

if __name__=='__main__':sys.exit(main())
