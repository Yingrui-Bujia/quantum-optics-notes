"""仅管理本项目生成的结果；每次运行使用独立目录，避免覆盖学习记录。"""
from pathlib import Path
from datetime import datetime
import json
import numpy as np
from helpers import save_figure

ROOT = Path(__file__).resolve().parents[1]

def run_directory(name, output=None):
    folder = Path(output) if output else ROOT / 'outputs' / 'runs' / (datetime.now().strftime('%Y%m%d-%H%M%S-%f') + '-' + name)
    folder.mkdir(parents=True, exist_ok=False)
    return folder

def save_plot(fig, folder, name):
    for ext in ('png', 'pdf'):
        save_figure(fig, folder / f'{name}.{ext}')

def save_json(folder, name, data):
    (folder / name).write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

def max_error(a, b):
    return float(np.max(np.abs(np.asarray(a) - np.asarray(b))))

def real_parameter(value, name, *, positive=False, nonnegative=False):
    raw=np.asarray(value)
    if raw.ndim!=0 or raw.dtype.kind not in 'iuf':
        raise ValueError(f'{name} 必须是有限实数标量')
    number=float(raw)
    if not np.isfinite(number) or (positive and number<=0) or (nonnegative and number<0):
        raise ValueError(f'{name} 超出允许范围')
    return number
