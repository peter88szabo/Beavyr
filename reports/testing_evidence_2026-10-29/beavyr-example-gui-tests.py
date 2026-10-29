import json
import os
from pathlib import Path
import re
import signal
import subprocess
import time

repo = Path('/home/peter/DevCode/Beavyr')
root = Path('/tmp/beavyr-example-gui-tests-20261029')
root.mkdir(exist_ok=True)
binary = '/tmp/beavyr-release-package-check-20261029-2135/usr/bin/beavyr'
examples = repo / 'examples'
trajectories = sorted(p for p in examples.rglob('*.xyz') if 'trajectory' in p.name)
structures = sorted(p for p in examples.rglob('*.xyz') if p not in trajectories)
cases = [('structures', structures)]
cases += [('trajectory_' + p.stem, [p]) for p in trajectories]
cases += [('hessian_' + p.stem, [p]) for p in sorted(examples.rglob('*.hess'))]
cases += [('molden_' + str(i), [p]) for i, p in enumerate(sorted(examples.rglob('*.molden')))]
cases += [('output_' + str(i), [p]) for i, p in enumerate(sorted(examples.rglob('*.out')))]
cases += [('gaussian_log', [examples / 'c164-ts1-2-1001.log'])]
cases += [('hessian_with_gradient', [examples / 'Exc_root_OnlyFreq_from_optimized.hess', examples / 'Exc_root_OnlyFreq_from_optimized.engrad'])]
results = []

def capture_own_window(pid, destination):
    try:
        tree = subprocess.run(['xwininfo', '-root', '-tree'], capture_output=True, text=True, timeout=3).stdout
        for line in tree.splitlines():
            if 'Beavyr' not in line:
                continue
            match = re.search(r'0x[0-9a-f]+', line)
            if not match:
                continue
            wid = match.group()
            owner = subprocess.run(['xprop', '-id', wid, '_NET_WM_PID'], capture_output=True, text=True, timeout=3).stdout
            if re.search(r'=\s*' + str(pid) + r'\s*$', owner):
                image = subprocess.run(['import', '-window', wid, str(destination)], capture_output=True, timeout=5)
                return image.returncode == 0
    except (OSError, subprocess.TimeoutExpired):
        pass
    return False

for name, files in cases:
    directory = root / name
    directory.mkdir(exist_ok=True)
    env = os.environ.copy()
    env.update(XDG_CONFIG_HOME=str(directory / 'config'), XDG_CACHE_HOME=str(directory / 'cache'),
               XDG_DATA_HOME=str(directory / 'data'), OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='2', RUST_LOG='info')
    log_path = directory / 'application.log'
    with log_path.open('w') as log:
        proc = subprocess.Popen([binary, *map(str, files)], cwd=repo, env=env,
                                stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        started = time.monotonic()
        screenshot = False
        try:
            while time.monotonic() - started < 14 and proc.poll() is None:
                time.sleep(0.2)
            alive = proc.poll() is None
            if alive:
                screenshot = capture_own_window(proc.pid, directory / 'window.png')
        finally:
            if proc.poll() is None:
                os.killpg(proc.pid, signal.SIGTERM)
                try:
                    proc.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL)
                    proc.wait(timeout=5)
    text = log_path.read_text(errors='replace')
    failures = [line for line in text.splitlines() if 'beavyr:' in line or 'panicked at' in line or 'Failed to initialize' in line]
    # Import messages come from the real production command-line loaders.
    loaded = [line for line in text.splitlines() if re.search(r': (\d+ (?:frames|atoms)|excited states loaded|loaded)', line)
              or ('Loaded' in line and ('Hessian' in line or 'gradient' in line or 'Gaussian' in line))]
    entry = dict(case=name, files=[str(p.relative_to(repo)) for p in files],
                 alive_after_14_seconds=alive, exit_code=proc.returncode,
                 own_window_screenshot=screenshot, load_messages=loaded,
                 reported_errors=failures, log=str(log_path))
    results.append(entry)
    (root / 'results.json').write_text(json.dumps(results, indent=2))
    print(json.dumps(entry), flush=True)
