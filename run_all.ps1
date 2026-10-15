param([string]$PythonPath = 'D:\Master.s-work\Software-needed-for-master.s-work\python\main\python.exe')
$ErrorActionPreference = 'Stop'
if (-not (Test-Path -LiteralPath $PythonPath -PathType Leaf)) {
    throw "Python not found: $PythonPath. Pass -PythonPath with your existing interpreter."
}
& $PythonPath -B -X utf8 (Join-Path $PSScriptRoot 'run_all.py')
if ($LASTEXITCODE -ne 0) { throw 'Verification failed. Read the newest outputs/runs log.' }
