"""执行官方三篇原始教程的完整Python块；保留原文，只移除inline魔法。"""
from pathlib import Path
import contextlib
import io
import json
import re
import sys
import time
import warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]

def main():
    folder=ROOT/'reports';results=[]
    for source in sorted((ROOT/'third_party/qutip_tutorials').glob('*.md')):
        text=source.read_text(encoding='utf-8')
        blocks=re.findall(r'^```python\s*\n(.*?)^```\s*$',text,re.M|re.S)
        namespace={'__name__':'__main__'};log=io.StringIO();executed=0;skipped=[];started=time.perf_counter()
        with contextlib.redirect_stdout(log),contextlib.redirect_stderr(log),warnings.catch_warnings(record=True) as recorded:
            warnings.simplefilter('always')
            for i,block in enumerate(blocks,1):
                block='\n'.join(line for line in block.splitlines() if not line.strip().startswith('%matplotlib'))
                if not block.strip():continue
                if block.strip()=='(t: float, ...) -> complex':
                    skipped.append({'block':i,'reason':'documentation-only function signature'})
                    continue
                exec(compile(block,f'{source.name}:block{i}','exec'),namespace)
                executed+=1
                plt.close('all')
            messages=sorted({str(w.message) for w in recorded})
        (folder/(source.stem+'_upstream.log')).write_text(log.getvalue(),encoding='utf-8')
        results.append({'source':source.relative_to(ROOT).as_posix(),'python_blocks_executed':executed,
                        'documentation_blocks_skipped':skipped,'warnings':messages,'status':'passed','seconds':time.perf_counter()-started})
    report={'method':'complete Python code fences executed in source order; notebook inline magic removed; documentation-only type signature skipped; Agg backend',
            'executable':sys.executable,'status':'passed','tutorials':results}
    (folder/'upstream_execution.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
