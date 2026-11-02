"""Prepare and run a resumable local SOL144 dataset without an online agent.

python -m core.batch_run prepare --output-dir data_2026-10-17 --count 100 --seed 20261017
python -m core.batch_run run --output-dir data_2026-10-17
"""
import argparse
import contextlib
import csv
import ctypes
from datetime import datetime, timedelta
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
import zipfile

from .data_generation import generate_data_sets
from .file_io import atomic_text, replace_file
from .nodal_info import get_box_nodes, get_floating_nodes
from .paths import BASE_DIR, CORE_DIR, DATA_DIR, PROJECT_ROOT, WORKSPACE_DIR, require_working_output
from .plate_strain import export_plate_strains


def timestamp():
    return datetime.now().astimezone().isoformat(timespec='seconds')


def write_json(path, data):
    atomic_text(path,json.dumps(data,indent=2,ensure_ascii=True))


def digest(path):
    with Path(path).open('rb') as source:
        return hashlib.file_digest(source,'sha256').hexdigest()


def input_hashes(folder):
    return {p.name:digest(p) for p in Path(folder).iterdir() if p.suffix.lower() in {'.bdf','.inc'}}


def batch_work(output):
    output = require_working_output(output)
    if output.parent != PROJECT_ROOT or not output.name.startswith('data_'):
        raise ValueError('Dataset must be a data_<date> directory directly under the project root')
    return output,WORKSPACE_DIR/f'batch_{output.name}'


def make_fiber_base(target):
    target.mkdir(parents=True,exist_ok=True)
    for p in BASE_DIR.iterdir():
        if p.suffix.lower() in {'.bdf','.inc'}:
            shutil.copy2(p,target/p.name)
    control = target/'subcases.inc'
    text = re.sub(r'^\s*(?:VECTOR|SPCFORCES|STRESS|APRES).*\n','',control.read_text(),flags=re.M)
    text = re.sub(r'STRAIN\([^\n]+\)=ALL','STRAIN(PRINT,SORT1,REAL,FIBER,CENTER)=ALL',text)
    if text.count('STRAIN(PRINT,SORT1,REAL,FIBER,CENTER)=ALL')!=10:
        raise ValueError('Expected the current ten-subcase template')
    control.write_text(text)


def prepare(output, count, seed):
    output,work = batch_work(output)
    if (output/'batch.json').exists():
        plan = json.loads((output/'batch.json').read_text())
        if plan['count']!=count or plan['master_seed']!=seed:
            raise ValueError('Existing batch has a different count/seed; use a new output directory')
        return plan
    output.mkdir(exist_ok=True)
    work.mkdir(parents=True,exist_ok=True)
    base = work/'base'
    make_fiber_base(base)
    shared = output/'shared'
    shared.mkdir(exist_ok=True)
    shutil.copytree(base,shared/'solver_base',dirs_exist_ok=True)
    shutil.copytree(BASE_DIR,shared/'original_base',dirs_exist_ok=True)
    for name in ('nodes.csv','elements.csv','fixed_nodes.csv'):
        shutil.copy2(DATA_DIR/name,shared/name)
    with (DATA_DIR/'nodes.csv').open() as f:
        nodes = {int(row['node_id']):[float(row[a]) for a in ('x','y','z')] for row in csv.DictReader(f)}
    excluded = get_floating_nodes(DATA_DIR) | get_box_nodes(DATA_DIR)
    eligible = sorted(set(nodes)-excluded)
    cases = []
    with (work/'generation.log').open('a',encoding='utf-8') as log,contextlib.redirect_stdout(log):
        for index in range(1,count+1):
            case_id = f'case_{index:03d}'
            case = work/case_id/'sample_1'
            accepted = None
            for rejection in range(100):
                case_seed = seed+index-1+count*rejection
                if not (case/'case_manifest.json').exists():
                    try:
                        generate_data_sets(n=1,base_dir=base,temp_dir=case.parent,info_dir=DATA_DIR,seed=case_seed,d_n=1)
                    except ValueError as exc:
                        if 'No elements changed' not in str(exc):
                            raise
                        continue
                manifest = json.loads((case/'case_manifest.json').read_text())
                if manifest['seed']!=case_seed or manifest['base_input_sha256']!=input_hashes(base):
                    raise ValueError(f'Existing inputs differ from this sampling plan: {case_id}')
                accepted = manifest
                break
            if accepted is None:
                raise ValueError(f'Could not sample a nonempty damage case: {case_id}')
            region = accepted['regions'][0]
            case_output = output/case_id
            case_output.mkdir(exist_ok=True)
            record = {'case_id':case_id,'seed':case_seed,'empty_selection_rejections':rejection,
                      'damage_regions':1,'center_node':region['center_node'],'center_xyz':nodes[region['center_node']],
                      'semi_axes':region['axes'],'rotation_matrix':region['orientation'],
                      'stiffness_scale':accepted['stiffness_scale'],'stiffness_loss_fraction':1-accepted['stiffness_scale'],
                      'modified_element_count':accepted['modified_element_count'],'modified_element_ids':accepted['modified_element_ids'],
                      'selected_node_count':len(region['selected_node_ids']),
                      'case_input_sha256':accepted['case_input_sha256']}
            write_json(case_output/'case.json',record)
            cases.append({'case_id':case_id,'seed':case_seed,'state':'pending','attempts':0})
    plan = {'schema_version':1,'created_at':timestamp(),'count':count,'master_seed':seed,'subcases_per_case':10,
            'output_dir':str(output),'work_dir':str(work),'solver':r'D:\Siemens\Femap 2406\nastran\bin\nastran.exe',
            'license_server':'29000@license-siemens.eng.uci.edu','memory':'1GB','threads':2,'workers':1,
            'original_base_sha256':input_hashes(BASE_DIR),'solver_base_sha256':input_hashes(base),'cases':cases}
    sampling = {'design':'100 independent random damage configurations; one ellipsoid per case; same 10 flight conditions.',
                'rng':'NumPy default_rng / PCG64. Case i seed = master_seed + i - 1.',
                'master_seed':seed,'numpy_version':__import__('numpy').__version__,
                'centers':f'Uniform over {len(eligible)} eligible mesh nodes, independently across cases; node-uniform, not volume-uniform.',
                'exclusions':'Unused nodes and nodes inside the fixed-node bounding box expanded by 30 model length units.',
                'semi_axes':'Independent Normal(30,5) draws; redraw nonpositive values. Model length units.',
                'orientation':'QR of a 3x3 standard-normal matrix; correct diagonal signs and determinant to +1. q=(xyz-center)@R.',
                'stiffness_scale':'clip(Normal(0.5,0.1),0.1,0.7); remaining stiffness fraction, not damage fraction.',
                'selection':'Modify an element only when all its nodes are inside the ellipsoid and outside the excluded node set.',
                'empty_selection':'Reject empty cases; retry with seed + count*rejection_index. Record accepted seed and rejection count.',
                'material_update':'Scale MAT1 E/G in all active shell membrane/bending/shear/coupling references and PBAR beams. Keep density and geometry.',
                'eligible_node_ids_sha256':hashlib.sha256(','.join(map(str,eligible)).encode()).hexdigest(),
                'coordinates':'Center and semi-axes use the mesh coordinate system and length units.',
                'strain_definition':'Centroid element-axis exx, eyy, engineering gamma_xy. Top=+element-z; bottom=-element-z; mid=(top+bottom)/2.',
                'solver_output_change':'Request FIBER strains; omit unused displacement/stress/SPC/APRES output. Loads and solution settings are unchanged.'}
    sampling['design']=sampling['design'].replace('100 ',f'{count} ',1)
    sampling['source_reference']='Installed Siemens Simcenter Nastran User Guide (Femap 2406), section 4-79: STRAIN(FIBER), property Z1/Z2.'
    write_json(output/'sampling.json',sampling)
    write_json(output/'batch.json',plan)
    (output/'README.txt').write_text(
        'SOL144 DAMAGE DATASET\n'
        f'{count} damage cases x 10 flight conditions. Sampling: sampling.json.\n'
        'Each case: top_plate_strain.csv, bottom_plate_strain.csv, mid_plate_strain.csv.\n'
        'Rows: 12,475 shell elements. Columns: ID + exx, eyy, gamma_xy for each subcase.\n'
        'Strains are dimensionless, in element axes at the centroid. gamma_xy is engineering shear.\n'
        'Top/bottom follow the element normal; mid is the component-wise arithmetic mean.\n'
        'case.json records position, ellipsoid, stiffness scale and modified element IDs.\n'
        'Inputs are archived per case. Shared model/mesh snapshots are in shared/.\n'
        'metadata.json records validation, file hashes and solver diagnostics.\n'
        'A previously traced empty-output DGEMM diagnostic is retained; other runtime errors fail validation.\n'
        'After verified export, redundant working BDF/F04/F06/XDB and scratch files are removed.\n'
        'The first raw F06 is retained compressed for audit; failed cases retain diagnostic files.\n'
        'status.json gives progress and estimated finish time. Resume with the same batch run command.\n',encoding='utf-8')
    publish(output,plan,'prepared')
    return plan


def publish(output, plan, state, current=None, worker_pid=None, solver_pid=None, message=None):
    done = [c for c in plan['cases'] if c['state']=='completed']
    failed = [c['case_id'] for c in plan['cases'] if c['state']=='failed']
    durations = plan.get('run_cycle_seconds') or [c['cycle_seconds'] for c in sorted(done,key=lambda c:c.get('completed_at','')) if c.get('cycle_seconds')]
    seconds = statistics.median(durations[-10:]) if durations else 60.
    remaining = max(0,plan['count']-len(done))*seconds
    now = datetime.now().astimezone()
    status = {'state':state,'updated_at':now.isoformat(timespec='seconds'),'completed':len(done),'target':plan['count'],
              'failed_cases':failed,'current_case':current,'worker_pid':worker_pid,'solver_pid':solver_pid,
              'seconds_per_case_estimate':round(seconds,1),'remaining_seconds_estimate':round(remaining),
              'estimated_finish_local':(now+timedelta(seconds=remaining)).isoformat(timespec='seconds') if state=='running' else None,
              'message':message}
    # The checkpoint is authoritative. A reader locking a progress view must
    # not abort a solver or turn a valid case into a failed calculation.
    write_json(output/'batch.json',plan)
    try:
        write_json(output/'status.json',status)
    except PermissionError as exc:
        print(f'{timestamp()} Progress view locked; checkpoint saved: {exc}',flush=True)
    temporary = output/'cases.csv.tmp'
    with temporary.open('w',newline='') as stream:
        fields = ['case_id','state','seed','center_node','x','y','z','a','b','c','stiffness_scale','modified_elements','cycle_seconds']
        writer = csv.DictWriter(stream,fieldnames=fields)
        writer.writeheader()
        for item in plan['cases']:
            info = json.loads((output/item['case_id']/'case.json').read_text())
            writer.writerow({'case_id':item['case_id'],'state':item['state'],'seed':item['seed'],
                             'center_node':info['center_node'],**dict(zip(('x','y','z'),info['center_xyz'])),
                             **dict(zip(('a','b','c'),info['semi_axes'])),'stiffness_scale':info['stiffness_scale'],
                             'modified_elements':info['modified_element_count'],'cycle_seconds':item.get('cycle_seconds','')})
    try:
        replace_file(temporary,output/'cases.csv')
    except PermissionError as exc:
        print(f'{timestamp()} Case index locked; checkpoint saved: {exc}',flush=True)
    return status


def cleanup_case(case_parent, work):
    """Delete only one verified batch work directory; never the dataset or base."""
    work,case_parent = Path(work).resolve(),Path(case_parent).resolve()
    if not work.is_relative_to(WORKSPACE_DIR.resolve()) or case_parent.parent!=work or not re.fullmatch(r'case_\d{3,}',case_parent.name):
        raise ValueError(f'Unsafe cleanup target: {case_parent}')
    if case_parent.exists():
        for path in [case_parent,*case_parent.rglob('*')]:
            if path.is_symlink() or (hasattr(path,'is_junction') and path.is_junction()):
                raise ValueError(f'Cleanup will not traverse links: {path}')
        shutil.rmtree(case_parent)


def archive_case(case, destination, keep_raw=False):
    temporary = destination/'inputs.zip.part'
    with zipfile.ZipFile(temporary,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
        for path in case.iterdir():
            if path.suffix.lower() in {'.bdf','.inc'} or path.name in {'params.txt','case_manifest.json'}:
                archive.write(path,path.name)
    with zipfile.ZipFile(temporary) as archive:
        if archive.testzip() is not None:
            raise ValueError('Input archive integrity check failed')
    replace_file(temporary,destination/'inputs.zip')
    for name in ('sol144.log','run-summary.json','run-invocation.json','launcher.stdout.txt','launcher.stderr.txt'):
        if (case/name).exists():
            shutil.copy2(case/name,destination/name)
    if keep_raw:
        temporary = destination/'sol144.f06.gz.part'
        with (case/'sol144.f06').open('rb') as source,gzip.open(temporary,'wb',compresslevel=1) as target:
            shutil.copyfileobj(source,target,1024*1024)
        replace_file(temporary,destination/'sol144.f06.gz')


@contextlib.contextmanager
def batch_lock(output):
    with (output/'batch.lock').open('a+b') as lock:
        if os.name=='nt':
            import msvcrt
            if lock.tell()==0:
                lock.write(b'0'); lock.flush()
            lock.seek(0)
            msvcrt.locking(lock.fileno(),msvcrt.LK_NBLCK,1)
        else:
            import fcntl
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        try:
            yield
        finally:
            if os.name=='nt':
                lock.seek(0); msvcrt.locking(lock.fileno(),msvcrt.LK_UNLCK,1)


def stop_solver(process):
    """Reap our own solver tree before any retry or fatal worker exit."""
    if process.poll() is None:
        if os.name=='nt':
            subprocess.run(['taskkill','/PID',str(process.pid),'/T','/F'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,
                           creationflags=subprocess.CREATE_NO_WINDOW)
        else:
            process.kill()
        process.wait(timeout=30)


def saved_solver_case(case, item):
    candidate = Path(item.get('solver_directory',case)).resolve()
    case = case.resolve()
    if candidate!=case and (candidate.parent!=case or not candidate.name.startswith('solve_')):
        raise ValueError(f'Invalid solver working directory: {candidate}')
    return candidate


def available_memory():
    if os.name!='nt':
        return None
    class MemoryStatus(ctypes.Structure):
        _fields_ = [('length',ctypes.c_uint),('load',ctypes.c_uint)] + [
            (name,ctypes.c_ulonglong) for name in ('total_phys','avail_phys','total_commit','avail_commit',
                                                 'total_virtual','avail_virtual','avail_extended')]
    memory = MemoryStatus()
    memory.length = ctypes.sizeof(memory)
    if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(memory)):
        raise OSError('Cannot read Windows memory availability')
    return {'available_commit':memory.avail_commit,'available_physical':memory.avail_phys}


def wait_for_memory(output, plan, item):
    size = re.fullmatch(r'(\d+(?:\.\d+)?)([KMGT])B',plan['memory'].upper())
    if not size:
        raise ValueError('Batch solver memory must use KB, MB, GB or TB')
    requested = float(size[1])*1024**('KMGT'.index(size[2])+1)
    required = requested + 1024**3
    notified = False
    while True:
        memory = available_memory()
        if memory is None or (memory['available_commit']>=required and memory['available_physical']>=512*1024**2):
            if notified:
                print(f'{timestamp()} Memory available; resuming {item["case_id"]}',flush=True)
            return
        message = (f'Waiting for memory: {memory["available_commit"]/1024**3:.2f} GiB commit available; '
                   f'{required/1024**3:.2f} GiB required including headroom. Automatic retry in 30 seconds.')
        publish(output,plan,'waiting_resources',item['case_id'],os.getpid(),message=message)
        if not notified:
            print(f'{timestamp()} {message}',flush=True)
            notified = True
        time.sleep(30)


def memory_allocation_failed(case):
    log = case/'sol144.log'
    if not log.exists():
        return False
    with log.open(errors='replace') as stream:
        return any('Unable to allocate requested memory' in line or 'OPEN CORE Allocation Failed' in line for line in stream)


def retain_attempt_diagnostics(case, destination, attempt):
    target = destination/'attempts'/f'attempt_{attempt:03d}'
    target.mkdir(parents=True,exist_ok=True)
    for name in ('sol144.log','run-summary.json','run-invocation.json','launcher.stdout.txt','launcher.stderr.txt','batch-launcher.log'):
        if (case/name).is_file():
            shutil.copy2(case/name,target/name)


def run_solver(case, output, work, plan, item):
    wait_for_memory(output,plan,item)
    # Never reuse an XDB or scratch database, including after an interrupted run.
    active = Path(tempfile.mkdtemp(prefix='solve_',dir=case))
    for path in case.iterdir():
        if path.is_file() and (path.suffix.lower() in {'.bdf','.inc'} or path.name in {'params.txt','case_manifest.json'}):
            shutil.copy2(path,active/path.name)
    item['solver_directory'] = str(active)
    publish(output,plan,'running',item['case_id'],os.getpid())
    command = [sys.executable,'-m','core.run_local_nastran',str(active),'--solver',plan['solver'],
               '--license-server',plan['license_server'],'--memory',plan['memory'],'--threads',str(plan['threads'])]
    env = os.environ.copy()
    env.update(PYTHONIOENCODING='utf-8',PYTHONDONTWRITEBYTECODE='1',MPLBACKEND='Agg')
    env.pop('MKL_VERBOSE',None)
    started = time.monotonic()
    with (active/'batch-launcher.log').open('wb') as stream:
        process = subprocess.Popen(command,cwd=PROJECT_ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,
                                   creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
        try:
            while process.poll() is None:
                publish(output,plan,'running',item['case_id'],os.getpid(),process.pid)
                if time.monotonic()-started>1800:
                    raise TimeoutError('Solver exceeded the 30-minute case timeout')
                time.sleep(5)
        finally:
            stop_solver(process)
    return active


def run(output):
    output,work = batch_work(output)
    plan = json.loads((output/'batch.json').read_text())
    if input_hashes(output/'shared/solver_base')!=plan['solver_base_sha256']:
        raise ValueError('Shared solver template changed')
    with batch_lock(output):
        source_hashes={p.name:digest(p) for p in CORE_DIR.glob('*.py')}
        if plan.get('pipeline_source_sha256') and plan['pipeline_source_sha256']!=source_hashes:
            raise ValueError('Pipeline code differs from the saved batch; review before resuming')
        plan['pipeline_source_sha256']=source_hashes
        plan['python_version']=sys.version.split()[0]
        plan['run_cycle_seconds']=[]
        plan['run_started_at']=timestamp()
        snapshot=output/'shared/pipeline_source.zip'
        if not snapshot.exists():
            with zipfile.ZipFile(snapshot,'w',compression=zipfile.ZIP_DEFLATED) as archive:
                for p in CORE_DIR.glob('*.py'):
                    archive.write(p,'core/'+p.name)
        awake = False
        if os.name=='nt':
            ctypes.windll.kernel32.SetThreadExecutionState.argtypes=[ctypes.c_uint]
            ctypes.windll.kernel32.SetThreadExecutionState.restype=ctypes.c_uint
            awake = bool(ctypes.windll.kernel32.SetThreadExecutionState(0x80000001))
        print(f'{timestamp()} Started worker {os.getpid()}; idle-sleep prevention={awake}',flush=True)
        try:
            for item in plan['cases']:
                if item['state']=='completed':
                    metadata=json.loads((output/item['case_id']/'metadata.json').read_text())
                    if any(digest(output/item['case_id']/name)!=expected for name,expected in metadata['csv_sha256'].items()):
                        raise ValueError(f'Completed result changed: {item["case_id"]}')
                    if (work/item['case_id']).exists():
                        if digest(output/item['case_id']/'inputs.zip')!=metadata['inputs_archive_sha256']:
                            raise ValueError('Input archive changed; refusing cleanup')
                        cleanup_case(work/item['case_id'],work)
                    continue
                if {p.name:digest(p) for p in CORE_DIR.glob('*.py')}!=source_hashes:
                    raise ValueError('Pipeline code changed during the batch; resume after reviewing the change')
                case_id = item['case_id']
                case = work/case_id/'sample_1'
                destination = output/case_id
                if not case.exists():
                    raise FileNotFoundError(f'Missing pending case inputs: {case}')
                reuse_existing = item['state']!='failed'
                item['state'] = 'running'
                publish(output,plan,'running',case_id,os.getpid())
                last_error = None
                attempt = 1
                while attempt<=3:
                    item['attempts'] += 1
                    try:
                        active = saved_solver_case(case,item)
                        saved = active/'run-summary.json'
                        reusable = reuse_existing and attempt==1 and saved.exists() and json.loads(saved.read_text()).get('solver_completed',False)
                        if not reusable:
                            active = run_solver(case,output,work,plan,item)
                        export_started = time.monotonic()
                        metadata = export_plate_strains(active,destination)
                        archive_case(active,destination,keep_raw=(case_id=='case_001'))
                        if case_id=='case_001':
                            with gzip.open(destination/'sol144.f06.gz','rb') as raw:
                                if hashlib.file_digest(raw,'sha256').hexdigest()!=metadata['source_f06_sha256']:
                                    raise ValueError('Raw pilot archive integrity check failed')
                        metadata['export_seconds'] = round(time.monotonic()-export_started,3)
                        metadata['completed_at'] = timestamp()
                        metadata['inputs_archive_sha256'] = digest(destination/'inputs.zip')
                        metadata['raw_f06_retained'] = case_id=='case_001'
                        metadata['pipeline_source_sha256'] = source_hashes
                        write_json(destination/'metadata.json',metadata)
                        if (destination/'failure.json').exists():
                            history = destination/'failure-history.jsonl'
                            with history.open('a',encoding='utf-8') as stream:
                                stream.write(json.dumps(json.loads((destination/'failure.json').read_text()))+'\n')
                            (destination/'failure.json').unlink()
                        item.pop('last_error',None)
                        item.update(state='completed',cycle_seconds=round(metadata['solver_seconds']+metadata['export_seconds'],3),completed_at=timestamp())
                        plan['run_cycle_seconds'].append(item['cycle_seconds']+5)
                        publish(output,plan,'running',case_id,os.getpid())
                        try:
                            cleanup_case(case.parent,work)
                            cleanup_status='temporary files removed'
                        except OSError as cleanup_error:
                            item['cleanup_error']=str(cleanup_error)
                            cleanup_status='temporary files retained after cleanup error'
                        print(f'{timestamp()} {case_id} completed in {item["cycle_seconds"]:.1f}s; {cleanup_status}',flush=True)
                        last_error = None
                        break
                    except Exception as exc:
                        last_error = f'{type(exc).__name__}: {exc}'
                        item['last_error'] = last_error
                        active = saved_solver_case(case,item)
                        retain_attempt_diagnostics(active,destination,item['attempts'])
                        if memory_allocation_failed(active):
                            item['resource_retries']=item.get('resource_retries',0)+1
                            reuse_existing = False
                            publish(output,plan,'waiting_resources',case_id,os.getpid(),message='Solver memory allocation failed; waiting and retrying automatically.')
                            print(f'{timestamp()} {case_id}: memory allocation failed; retrying after resource wait',flush=True)
                            time.sleep(30)
                            continue
                        write_json(destination/'failure.json',{'at':timestamp(),'attempt':attempt,'error':last_error})
                        publish(output,plan,'running',case_id,os.getpid(),message=last_error)
                        print(f'{timestamp()} {case_id} attempt {attempt} failed: {last_error}',flush=True)
                        if attempt<3:
                            time.sleep(60)
                        attempt += 1
                if last_error:
                    item['state']='failed'
                    # Stop when authorization/connectivity is unavailable; do not burn through all 100 cases.
                    summary_path=saved_solver_case(case,item)/'run-summary.json'
                    failed_summary=json.loads(summary_path.read_text()) if summary_path.exists() else {}
                    diagnostic=' '.join(failed_summary.get('fatal_messages',[])+failed_summary.get('runtime_errors',[]))+last_error
                    if re.search(r'authoriz|licens|connection|network|FATAL MESSAGE 1259',diagnostic,re.I):
                        publish(output,plan,'paused',case_id,os.getpid(),message='Restore VPN/license access, then resume the same command.')
                        return 2
            state='completed' if all(c['state']=='completed' for c in plan['cases']) else 'incomplete'
            publish(output,plan,state,worker_pid=os.getpid())
            print(f'{timestamp()} Batch {state}',flush=True)
            return 0 if state=='completed' else 2
        except BaseException as exc:
            publish(output,plan,'interrupted',worker_pid=os.getpid(),message=f'{type(exc).__name__}: {exc}')
            raise
        finally:
            if awake:
                ctypes.windll.kernel32.SetThreadExecutionState(0x80000000)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['prepare','run'])
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--count',type=int,default=100)
    parser.add_argument('--seed',type=int,default=20261017)
    args=parser.parse_args()
    if args.count<1:
        parser.error('--count must be positive')
    if args.action=='prepare':
        plan=prepare(args.output_dir,args.count,args.seed)
        print(f'Prepared {plan["count"]} cases in {plan["output_dir"]}')
        return 0
    return run(args.output_dir)


if __name__=='__main__':
    raise SystemExit(main())
