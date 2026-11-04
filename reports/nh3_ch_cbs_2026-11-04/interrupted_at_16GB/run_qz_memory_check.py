"""One new aVQZ point; never repeat an attempted calculation."""
import json
import os
from pathlib import Path
import signal
import subprocess
import time

root = Path(__file__).resolve().parent
record = root / 'status.json'
if record.exists():
    raise SystemExit('Already attempted: no retry permitted.')
exe = '/home/peter/DevCode/Behemoth/target/release/behemoth'
env = dict(os.environ)
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'RAYON_NUM_THREADS'):
    env[key] = '8'
state = dict(status='STARTING', basis='aug-cc-pVQZ', distance_angstrom=10.0,
             reference='UHF', cores=8, memory_limit_bytes=16000000000,
             reuse_aVTZ=True, retries=False)
record.write_text(json.dumps(state, indent=2))
start = time.monotonic()
peak = 0
reason = None
with (root/'output.out').open('w') as out, (root/'error.err').open('w') as err:
    child = subprocess.Popen(['taskset', '-c', '0-7', exe, 'input.inp'], cwd=root,
                             env=env, stdout=out, stderr=err, start_new_session=True)
    state.update(status='RUNNING', pid=child.pid)
    record.write_text(json.dumps(state, indent=2))
    while child.poll() is None:
        try:
            fields = dict(line.split(':', 1) for line in Path(f'/proc/{child.pid}/status').read_text().splitlines() if ':' in line)
            rss = int(fields.get('VmRSS', '0 kB').split()[0])*1024
            peak = max(peak, rss)
            mem = dict(line.split(':',1) for line in Path('/proc/meminfo').read_text().splitlines())
            available = int(mem['MemAvailable'].split()[0])*1024
            if rss > 16000000000 or available < 4*1024**3:
                reason = 'MEMORY_LIMIT: stop CBS series; no retry'
                os.killpg(child.pid, signal.SIGTERM)
                try:
                    child.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    os.killpg(child.pid, signal.SIGKILL)
                break
        except FileNotFoundError:
            pass
        time.sleep(0.25)
    code = child.wait()
log = (root/'output.out').read_text(errors='replace')
success = code == 0 and 'Total UCCSD(T) energy:' in log and 'SCF successfully converged' in log
state.update(status='STOPPED_MEMORY' if reason else ('FINISHED_REQUIRES_VALIDATION' if success else 'FAILED_DISCARDED'),
             exit_code=code, wall_seconds=time.monotonic()-start,
             peak_rss_bytes=peak, reason=reason)
record.write_text(json.dumps(state, indent=2))
print(json.dumps(state), flush=True)
