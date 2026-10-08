"""One-shot standard-library capture for a separately authorized worker run."""
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
REVIEW = Path(__file__).resolve().parent
EXPERIMENT = ROOT / 'pgfrag/experiments/mixed93_moving_front_reanalysis_20261007'
ENTRY = EXPERIMENT / 'mixed93_moving_front_reanalysis.py'
CASE = EXPERIMENT / 'case.json'
PINS = REVIEW / 'SOURCE_PINS.json'
RUN = REVIEW / 'run01'
CAPTURE = REVIEW / 'capture_once01'
LOCK = REVIEW / 'run01.capture.lock'
PYTHON = ROOT / 'third_party/vesselforge-runtime/venv/Scripts/python.exe'
MAX_CHILD_SECONDS = 60
MAX_CAPTURE_SECONDS = 64
MAX_TOTAL_BYTES = 12 * 1024 * 1024


def _record(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return {'sha256': digest.hexdigest(), 'bytes': path.stat().st_size}


def _input_paths(case, pins):
    paths = {str(ENTRY.relative_to(ROOT)): ENTRY, str(CASE.relative_to(ROOT)): CASE,
        str(PINS.relative_to(ROOT)): PINS, str(Path(__file__).resolve().relative_to(ROOT)): Path(__file__).resolve()}
    for item in pins['references'].values():
        paths[item['path']] = ROOT / item['path']
    for item in pins.get('artifacts', {}).values():
        paths[item['path']] = ROOT / item['path']
    for item in case['dependency_pins'].values():
        paths[item['path']] = ROOT / item['path']
    return paths


def _snapshot(paths):
    return {name: _record(path) for name, path in sorted(paths.items())}


def _verified_inputs():
    case = json.loads(CASE.read_text(encoding='utf-8'))
    pins = json.loads(PINS.read_text(encoding='utf-8'))
    contract = case['capture_contract']
    if contract != {'script': str(Path(__file__).resolve().relative_to(ROOT)).replace('\\', '/'), 'exact_child_invocations': 1, 'max_child_seconds': MAX_CHILD_SECONDS, 'max_capture_seconds': MAX_CAPTURE_SECONDS, 'max_total_output_bytes': MAX_TOTAL_BYTES, 'pre_post_all_input_hashes_must_match': True, 'raw_stdout_stderr': True, 'automatic_retry': False}:
        raise RuntimeError('capture script settings differ from the reviewed case contract')
    paths = _input_paths(case, pins)
    if case['execution_authorization']['authorized'] is not True or pins['execution_authorized'] is not True:
        raise RuntimeError('capture blocked: exact case and source pins are not authorized')
    if case.get('source_entry') != pins.get('entry_path') or pins.get('entry_path') != ENTRY.relative_to(ROOT).as_posix():
        raise RuntimeError('capture blocked: entry path differs between case and pins')
    if _record(ENTRY)['sha256'] != pins['entry_sha256'] or _record(CASE)['sha256'] != pins['current_case_sha256']:
        raise RuntimeError('capture blocked: current entry/case does not match source pins')
    capture_pin = pins['artifacts']['capture_wrapper']
    if capture_pin['path'] != str(Path(__file__).resolve().relative_to(ROOT)).replace('\\', '/') or _record(Path(__file__).resolve())['sha256'] != capture_pin['sha256']:
        raise RuntimeError('capture blocked: wrapper differs from its artifact pin')
    for key, item in pins['references'].items():
        if _record(ROOT / item['path'])['sha256'] != item['sha256']:
            raise RuntimeError(f'capture blocked: pinned reference changed: {key}')
    for key, item in pins.get('artifacts', {}).items():
        if _record(ROOT / item['path'])['sha256'] != item['sha256']:
            raise RuntimeError(f'capture blocked: pinned artifact changed: {key}')
    for key, item in case['dependency_pins'].items():
        if _record(ROOT / item['path'])['sha256'] != item['sha256']:
            raise RuntimeError(f'capture blocked: case dependency changed: {key}')
    return (case, pins, paths)


def _tree_records(path):
    if not path.exists():
        return {}
    return {p.relative_to(path).as_posix(): _record(p) for p in sorted(path.rglob('*')) if p.is_file()}


def _write_receipt(receipt, receipt_path, index_path, started):
    total_bytes = 0
    for _ in range(6):
        receipt['combined_output_bytes_including_receipt_index'] = total_bytes
        receipt['capture_elapsed_seconds_after_posthash_tally'] = time.monotonic() - started
        receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n', encoding='utf-8')
        index_path.write_text(json.dumps({'receipt_path': receipt_path.name, **_record(receipt_path)}, indent=2, sort_keys=True) + '\n', encoding='utf-8')
        actual = sum(p.stat().st_size for folder in (RUN, CAPTURE) if folder.exists() for p in folder.rglob('*') if p.is_file())
        if actual == total_bytes:
            return (actual, time.monotonic() - started)
        total_bytes = actual
    raise RuntimeError('capture receipt/index byte tally did not stabilize')


def main():
    started = time.monotonic()
    if LOCK.exists() or RUN.exists() or CAPTURE.exists():
        raise RuntimeError('capture blocked: lock or exclusive output already exists')
    LOCK.mkdir()
    try:
        case, pins, paths = _verified_inputs()
        before = _snapshot(paths)
        child_timeout = min(MAX_CHILD_SECONDS, MAX_CAPTURE_SECONDS - (time.monotonic() - started) - 4.0)
        if child_timeout <= 0.0:
            raise RuntimeError('capture preflight exhausted the reserved four-second receipt budget')
        command = [str(PYTHON), '-B', str(ENTRY), '--execute', '--output', str(RUN / 'result.json')]
        timed_out = False; child_return_code = None
        try:
            child = subprocess.run(command, cwd=ROOT, capture_output=True, timeout=child_timeout, check=False)
            child_return_code = child.returncode; stdout = child.stdout; stderr = child.stderr
        except subprocess.TimeoutExpired as exc:
            timed_out = True; stdout = exc.stdout or b''; stderr = exc.stderr or b''
        after = _snapshot(paths)
        inputs_stable = before == after
        CAPTURE.mkdir()
        stdout_path = CAPTURE / 'stdout.bin'; stdout_path.write_bytes(stdout)
        stderr_path = CAPTURE / 'stderr.bin'; stderr_path.write_bytes(stderr)
        run_files = _tree_records(RUN)
        receipt_path = CAPTURE / 'capture_receipt.json'
        status = 'CAPTURED_EXACTLY_ONCE' if inputs_stable and not timed_out and child_return_code == 0 else 'CAPTURE_NEEDS_REVIEW'
        receipt = {'schema': 'mixed93_moving_front_capture.v1', 'status': status,
            'case_id': case['case_id'], 'child_invocations': 1, 'command': command,
            'actual_child_return_code': child_return_code, 'child_timed_out': timed_out, 'child_timeout_seconds': child_timeout,
            'whole_capture_timeout_seconds': MAX_CAPTURE_SECONDS,
            'inputs_stable_pre_post': inputs_stable, 'inputs_pre': before, 'inputs_post': after,
            'run01_files': run_files, 'stdout': {'path': stdout_path.name, **_record(stdout_path)},
            'stderr': {'path': stderr_path.name, **_record(stderr_path)},
            'capture_source': {'path': Path(__file__).name, **_record(Path(__file__).resolve())},
            'primary_result_present': (RUN / 'result.json').is_file(), 'automatic_retry': False,
            'combined_output_limit_bytes': MAX_TOTAL_BYTES,
            'combined_output_bytes_including_receipt_index': 0, 'capture_elapsed_seconds_after_posthash_tally': 0.0}
        index_path = CAPTURE / 'receipt_index.json'
        total_bytes, elapsed = _write_receipt(receipt, receipt_path, index_path, started)
        return_code = child_return_code if child_return_code is not None else 124
        if child_return_code == 0 and not inputs_stable: return_code = 3
        if total_bytes > MAX_TOTAL_BYTES:
            receipt['status'] = 'CAPTURE_OUTPUT_LIMIT_EXCEEDED'
            return_code = 3
        if elapsed > MAX_CAPTURE_SECONDS:
            receipt['status'] = 'CAPTURE_TIME_LIMIT_EXCEEDED'
            return_code = 124
        receipt['wrapper_return_code'] = return_code
        if receipt['status'] != 'CAPTURED_EXACTLY_ONCE' or return_code != 0:
            total_bytes, elapsed = _write_receipt(receipt, receipt_path, index_path, started)
        if elapsed > MAX_CAPTURE_SECONDS and return_code != 124:
            receipt['status'] = 'CAPTURE_TIME_LIMIT_EXCEEDED'; return_code = 124; receipt['wrapper_return_code'] = return_code
            _write_receipt(receipt, receipt_path, index_path, started)
        return return_code
    finally:
        LOCK.rmdir()


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f'capture failed before/around the single child run: {type(exc).__name__}: {exc}', file=sys.stderr)
        raise SystemExit(2)
