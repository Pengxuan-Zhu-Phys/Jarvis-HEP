"""Regression contracts for the 2.0.15 reliability fixes."""
import contextlib
import io
import json
import tempfile
import threading
import time
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from jarvishep2.cli.client import main, _dispatch_monitor_snapshot, dispatch_plot
from jarvishep2.dashboard import attach_reader
from jarvishep2.process_cleanup import JarvisProcess, JarvisScan
from jarvishep2.redis_queue import make_fakeredis_queue
from jarvishep2.task_validation import validate_task_config
from test_task_validation import _minimal_dynesty_config


def test_missing_operas_names_report_each_expression_path():
    cfg = _minimal_dynesty_config()
    cfg['Sampling']['LogLikelihood'] = [{'name': 'LogL', 'expression': 'nosuch.f(x)'}]
    cfg['Sampling']['selection'] = 'helper.nosuch(x) > 0'
    report = validate_task_config(cfg, strict=True)
    errors = [item for item in report.errors() if item.code == 'JV2-EXP-001']
    assert {item.path for item in errors} == {'Sampling.LogLikelihood[0].expression', 'Sampling.selection'}
    assert all(item.suggestion and item.hint for item in errors)


def test_validate_compiles_plugin_without_numeric_execution():
    import jarvis_operas
    from jarvishep2.operas_functions import reset_operas_discovery_for_tests
    registry = jarvis_operas.OperasRegistry()
    calls = []
    @jarvis_operas.oper('probe', namespace='review', registry=registry)
    def probe(x):
        calls.append(x)
        raise AssertionError('validation executed numeric physics')
    cfg = _minimal_dynesty_config()
    cfg['Sampling']['LogLikelihood'] = [{'name': 'LogL', 'expression': 'review.probe(x)'}]
    reset_operas_discovery_for_tests()
    with patch.object(jarvis_operas, 'get_global_operas_registry', return_value=registry):
        report = validate_task_config(cfg)
    assert report.ok, report.issues
    assert not calls


def test_mapper_keeps_its_existing_diagnostic():
    cfg = _minimal_dynesty_config()
    cfg['Sampling']['Mapper'] = [{'name': 'z', 'expression': 'nosuch.f(x)'}]
    errors = validate_task_config(cfg).errors()
    assert any(item.code == 'JV2-MAP-007' for item in errors)
    assert not any(item.code == 'JV2-EXP-001' for item in errors)


@pytest.mark.parametrize('argv', [
    ['--monitor', 'validate', 'task.yaml'],
    ['validate', 'task.yaml', '--monitor'],
    ['task.yaml', '--validate', '--monitor'],
    ['task.yaml', '--plot', '--validate'],
])
def test_conflicts_before_and_after_command(argv):
    with contextlib.redirect_stderr(io.StringIO()) as err:
        assert main(argv) == 2
    assert 'conflicting CLI intents' in err.getvalue()


def test_dependency_cli_tail_is_opaque():
    with patch('jarvishep2.cli.client.dispatch_plot', return_value=0) as forward:
        assert main(['plot', 'scene.yaml', '--monitor']) == 0
    forward.assert_called_once_with(['scene.yaml', '--monitor'])


def test_plot_task_card_explains_generation(tmp_path):
    path = tmp_path / 'task.yaml'
    path.write_text('Scan: {name: review}\nSampling: {Method: Random}\n')
    with contextlib.redirect_stderr(io.StringIO()) as err:
        assert dispatch_plot([str(path)]) == 2
    assert 'DataSet' in err.getvalue()
    assert f'Jarvis gen-plot-yaml {path}' in err.getvalue()


def _scan():
    return JarvisScan('R1', 'review', (
        JarvisProcess(900, 'Jarvis-Core:review'),
        JarvisProcess(901, 'Jarvis-Worker-0:review'),
        JarvisProcess(902, 'Jarvis-Worker-1:review'),
    ))


def test_monitor_json_ref_and_name_use_os_inventory():
    queue = make_fakeredis_queue()
    queue.connect()
    for wid, pid in [('0', 901), ('1', 902)]:
        queue.heartbeat(wid, status='idle', pid=pid, ts=time.time())
    scan = _scan()
    metadata = {'control_pid': 900, 'redis': {}}
    with patch('jarvishep2.process_cleanup.list_active_scans', return_value=[scan]), \
         patch('jarvishep2.process_cleanup.runtime_metadata_for_scan', return_value=metadata), \
         patch('jarvishep2.redis_queue.RedisQueue', return_value=queue), patch.object(queue, 'close'):
        for ref in ('R1', 'review'):
            with contextlib.redirect_stdout(io.StringIO()) as out:
                assert _dispatch_monitor_snapshot(ref, as_json=True) == 0
            payload = json.loads(out.getvalue())
            assert payload['ok'] is True
            assert payload['factory'] == {'workers_alive': 2, 'workers_total': 2}
            assert {row['pid'] for row in payload['workers']} == {901, 902}
            assert all(row['status'] == 'idle' for row in payload['workers'])
            assert abs(payload['timestamp'] - time.time()) < 2


def test_monitor_errors_remain_json():
    with patch('jarvishep2.process_cleanup.list_active_scans', return_value=[]):
        with contextlib.redirect_stdout(io.StringIO()) as out:
            assert _dispatch_monitor_snapshot('R99', as_json=True) == 2
        assert json.loads(out.getvalue())['ok'] is False


@pytest.mark.parametrize('mode', ['missing', 'stale', 'zero', 'unknown'])
def test_monitor_missing_boards_and_unknown_counts(mode):
    queue = make_fakeredis_queue()
    queue.connect()
    if mode == 'stale':
        queue.heartbeat('0', pid=12, status='busy', ts=1)
    workers = [] if mode == 'zero' else [{'worker_id': '0', 'pid': 22, 'alive': True}]
    view = attach_reader(redis=queue, owner_ids=['0'], workers=None if mode == 'unknown' else workers).read()
    if mode == 'unknown':
        assert view.factory['workers_alive'] is None
    elif mode == 'zero':
        assert view.factory['workers_alive'] == 0
        assert not view.workers
    else:
        assert view.factory['workers_alive'] == 1
        assert view.workers[0]['pid'] == 22
        assert view.workers[0]['status'] is None


@pytest.mark.slow
def test_spawn_initialization_failure_is_fatal_and_cleans_children(tmp_path):
    from jarvishep2.factory import TaskFactory
    from test_worker_calculator import _start_tcp_fakeredis, _worker_config
    server, redis_config = _start_tcp_fakeredis()
    factory = TaskFactory(redis_config)
    factory.init_redis()
    config = _worker_config(str(tmp_path))
    config['likelihood_expressions'] = [{'name': 'LogL', 'expression': 'nosuch.f(z)'}]
    started = time.monotonic()
    try:
        with pytest.raises(RuntimeError, match="initialization failed: AttributeError.*no attribute 'f'"):
            factory.start_workers(2, **config)
            factory.wait_for_worker_startup(timeout=20)
        assert time.monotonic() - started < 30
        assert factory._respawn_count == 0
        assert all(not worker.is_alive() for worker in factory.workers)
        assert int(factory.redis.fetch_sample_stats().get('failed', 0)) == 0
        assert not factory.redis.read_proc_board('worker', owner_id='0')
        assert not factory.redis.read_proc_board('worker', owner_id='1')
    finally:
        factory.shutdown()
        server.shutdown()
        server.server_close()


@pytest.mark.parametrize('method', ['MultiNest', 'Dynesty'])
def test_nested_budget_return_does_not_claim_convergence(tmp_path, method):
    import numpy as np
    from jarvishep2.sampling.dynesty_sampler import DynestySampler
    from jarvishep2.sampling.multinest_sampler import MultiNestSampler
    queue = make_fakeredis_queue()
    stop = threading.Event()
    def feedback_worker():
        while not stop.is_set():
            task = queue.pull_task(timeout=1)
            if task:
                u = np.asarray(task['u_coords'], dtype=float)
                queue.publish_feedback({'uuid': task['uuid'], 'logL': -float(np.sum((u - 0.5)**2))*1000})
    thread = threading.Thread(target=feedback_worker, daemon=True)
    thread.start()
    sampler = MultiNestSampler() if method == 'MultiNest' else DynestySampler()
    sampler.set_config({
        'Sampling': {'Method': method, 'Variables': [{'name': 'x', 'distribution': {'type': 'Flat', 'parameters': {'min': 0, 'max': 1}}}],
                     'Bounds': {'seed': 7, 'nlive': 10, 'dlogz': 1e-8, 'run_nested': {'maxiter': 1, 'print_progress': False}}},
        'Runtime': {'workers': 1, 'batch_size': 1}, 'EnvReqs': {'V2': {'checkpoint': {'enabled': False}}},
        'task_root': str(tmp_path), 'task_result_dir': str(tmp_path), 'scan_name': 'budget-test',
    })
    sampler.set_redis(queue)
    try:
        sampler.run_adaptive(generation_timeout=30)
        reason = sampler.log_stop_reason()
        assert sampler._finished
        assert 'engine returned' in reason and 'maxiter' in reason
        assert 'converged' not in reason
    finally:
        stop.set()
        thread.join(timeout=2)


def test_monitor_no_ref_keeps_json_array():
    with patch('jarvishep2.process_cleanup.list_active_scans', return_value=[_scan()]):
        with contextlib.redirect_stdout(io.StringIO()) as out:
            assert _dispatch_monitor_snapshot('', as_json=True) == 0
        payload = json.loads(out.getvalue())
        assert isinstance(payload, list) and payload[0]['reference'] == 'R1'


@pytest.mark.parametrize('failure', ['metadata', 'redis'])
def test_monitor_attach_failures_are_json(failure):
    scan = _scan()
    metadata = None if failure == 'metadata' else {'control_pid': 900, 'redis': {}}
    broken = SimpleNamespace(connect=lambda: (_ for _ in ()).throw(ConnectionError('offline')), close=lambda: None)
    with patch('jarvishep2.process_cleanup.list_active_scans', return_value=[scan]), \
         patch('jarvishep2.process_cleanup.runtime_metadata_for_scan', return_value=metadata), \
         patch('jarvishep2.redis_queue.RedisQueue', return_value=broken):
        with contextlib.redirect_stdout(io.StringIO()) as out:
            assert _dispatch_monitor_snapshot('R1', as_json=True) == 1
        assert json.loads(out.getvalue())['ok'] is False


def test_partial_spawn_failure_stops_already_started_workers():
    from jarvishep2.runtime.factory import TaskFactory
    factory = TaskFactory()
    factory.redis = make_fakeredis_queue()
    started = []
    class FakeWorker:
        def __init__(self, worker_id, *args):
            self.worker_id = worker_id
        def start(self):
            if self.worker_id == 1:
                raise OSError('spawn failed')
            started.append(self)
        def is_alive(self):
            return True
    with patch('jarvishep2.runtime.factory.Worker', FakeWorker), patch.object(factory, 'stop_all_workers') as stop:
        with pytest.raises(OSError, match='spawn failed'):
            factory.start_workers(2)
        assert factory.workers == started
        stop.assert_called_once()
