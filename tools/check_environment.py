"""对照开工前快照，核查包版本与受保护原文件；不修改任何安装。"""
from pathlib import Path
import argparse
import hashlib
import importlib.metadata
import json
import sys

ROOT=Path(__file__).resolve().parents[1]

def check(baseline):
    before=json.loads(Path(baseline).read_text(encoding='utf-8'))
    after={d.metadata['Name']:d.version for d in importlib.metadata.distributions()}
    package_changes={k:{'before':before['packages'].get(k),'after':after.get(k)}
                     for k in set(before['packages'])|set(after) if before['packages'].get(k)!=after.get(k)}
    files=[]
    for name,digest in before['protected_files'].items():
        p=Path(name)
        actual=hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else None
        files.append({'path':name,'unchanged':actual==digest,'sha256_before':digest,'sha256_after':actual})
    result={'status':'passed' if not package_changes and all(f['unchanged'] for f in files) else 'failed',
            'python_executable_unchanged':sys.executable==before['executable'],
            'installed_package_count':len(after),'package_changes':package_changes,
            'protected_file_count':len(files),'protected_files':files}
    if not result['python_executable_unchanged']:result['status']='failed'
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline',type=Path)
    args=parser.parse_args();result=check(args.baseline)
    (ROOT/'reports/environment_integrity.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(result['status'],f"packages={result['installed_package_count']}",f"protected_files={result['protected_file_count']}")
    sys.exit(0 if result['status']=='passed' else 1)
