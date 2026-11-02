"""Export shell top/bottom fiber strains and their arithmetic mean.

Input must request STRAIN(PRINT,SORT1,REAL,FIBER,CENTER)=ALL.
Top is +element-z, bottom is -element-z; gamma_xy is engineering shear strain.
"""
import csv
import hashlib
import json
import math
from pathlib import Path
import re

from .bdf_damage import read_cards, nastran_float
from .file_io import replace_file
from .run_local_nastran import inspect_results
from .paths import require_working_output

KNOWN_MKL = 'Intel MKL ERROR: Parameter 13 was incorrect on entry to DGEMM .'


def read_fiber_strains(path):
    rows, pending = {}, None
    section, fiber, subcase, kind, ended = False, False, None, None, False
    with Path(path).open(errors='replace') as stream:
        for line in stream:
            ended |= '* * * END OF JOB * * *' in line
            match = re.search(r'^0\s+SUBCASE\s+(\d+)\s*$',line)
            if match:
                subcase = int(match[1])
            compact = ''.join(line.split())
            if 'STRAINSINQUADRILATERALELEMENTS(QUAD4)' in compact:
                section, fiber, kind = True, False, 'CQUAD4'
            elif 'STRAINSINTRIANGULARELEMENTS(TRIA3)' in compact:
                section, fiber, kind = True, False, 'CTRIA3'
            elif any(title in compact for title in ('STRESSESIN','STRAINSINBAR','DISPLACEMENTVECTOR','ANALYSISSUMMARY','FORCESOFSINGLE')):
                section = False
            if not section:
                continue
            if 'DISTANCE' in line:
                fiber = True
            if 'CURVATURE' in line:
                raise ValueError('Expected fiber strains; received membrane/curvature output')
            fields = line.split()
            if len(fields)==10 and fields[0]=='0' and fields[1].isdigit():
                if pending is not None or not fiber or subcase is None:
                    raise ValueError('Incomplete/unrecognized fiber-strain pair')
                pending = (subcase,int(fields[1]),kind,float(fields[2]),tuple(map(float,fields[3:6])))
            elif pending is not None and len(fields)==8:
                try:
                    z2 = float(fields[0])
                    values2 = tuple(map(float,fields[1:4]))
                except ValueError:
                    continue
                s,eid,element_type,z1,values1 = pending
                if z1==z2 or not all(math.isfinite(v) for v in (z1,z2,*values1,*values2)):
                    raise ValueError(f'Invalid fiber results: {s}/{eid}')
                bottom,top = sorted([(z1,values1),(z2,values2)])
                if not bottom[0]<0<top[0]:
                    raise ValueError(f'Fibers do not straddle the midsurface: {s}/{eid}')
                key = (s,eid)
                if key in rows:
                    raise ValueError(f'Duplicate fiber strain: {key}')
                rows[key] = {'type':element_type,'z_bottom':bottom[0],'z_top':top[0],
                             'bottom':bottom[1],'top':top[1],
                             'mid':tuple((a+b)/2 for a,b in zip(bottom[1],top[1]))}
                pending = None
    if not ended or not rows or pending is not None:
        raise ValueError('Incomplete fiber-strain output')
    return rows


def export_plate_strains(case_dir, output_dir):
    case_dir, output_dir = Path(case_dir).resolve(), require_working_output(output_dir)
    summary = inspect_results(case_dir)
    if not summary['solver_completed']:
        raise ValueError('Solver output is incomplete or contains FATAL messages')
    recorded = json.loads((case_dir/'run-summary.json').read_text())
    hashes = {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in case_dir.iterdir() if p.suffix.lower() in {'.bdf','.inc'}}
    if recorded.get('fresh_f06') is not True or recorded.get('input_sha256')!=hashes:
        raise ValueError('Output freshness/input hash check failed')
    errors = set(recorded.get('runtime_errors',[])) | set(summary['runtime_errors'])
    if errors - {KNOWN_MKL,'Launcher returned exit code 1'}:
        raise ValueError(f'Unreviewed runtime diagnostics: {sorted(errors)}')
    # An unexplained launcher failure is not accepted on its own.
    if 'Launcher returned exit code 1' in errors and KNOWN_MKL not in errors:
        raise ValueError('Unexplained nonzero launcher return code')
    for warning in summary['warning_messages']:
        if not re.search(r'WARNING MESSAGE (324|4699)\b',warning):
            raise ValueError(f'Unreviewed solver warning: {warning}')
    rows = read_fiber_strains(case_dir/'sol144.f06')
    subcases = sorted({int(s) for s in re.findall(r'^\s*SUBCASE\s+(\d+)',(case_dir/'subcases.inc').read_text(),re.M|re.I)})
    cards = read_cards((case_dir/'FEM_only.bdf').read_text().splitlines(True))
    props = {int(c.fields[1]):c for c in cards if c.name=='PSHELL'}
    shell_props = {int(c.fields[1]):int(c.fields[2]) for c in cards if c.name in {'CQUAD4','CTRIA3'}}
    if set(rows)!={(s,eid) for s in subcases for eid in shell_props}:
        raise ValueError('Missing/extra shell IDs or subcases')
    for (s,eid),row in rows.items():
        prop = props[shell_props[eid]]
        thickness = nastran_float(prop.fields[3])
        z1 = nastran_float(prop.fields[9]) if len(prop.fields)>9 and prop.fields[9] else -thickness/2
        z2 = nastran_float(prop.fields[10]) if len(prop.fields)>10 and prop.fields[10] else thickness/2
        for actual,expected in zip((row['z_bottom'],row['z_top']),sorted((z1,z2))):
            if not math.isclose(actual,expected,rel_tol=2e-6,abs_tol=1e-10):
                raise ValueError(f'Fiber distance differs from PSHELL: {s}/{eid}')
    output_dir.mkdir(parents=True,exist_ok=True)
    ids = sorted(shell_props)
    file_hashes = {}
    for surface in ('top','bottom','mid'):
        name = f'{surface}_plate_strain.csv'
        destination = output_dir/name
        temporary = destination.with_suffix('.csv.part')
        with temporary.open('w',newline='') as stream:
            writer = csv.writer(stream)
            writer.writerow(['ID']+[f'subcase_{s}_{component}' for s in subcases for component in ('exx','eyy','gamma_xy')])
            writer.writerows([eid]+[value for s in subcases for value in rows[(s,eid)][surface]] for eid in ids)
        replace_file(temporary,destination)
        file_hashes[name] = hashlib.sha256(destination.read_bytes()).hexdigest()
    with (case_dir/'sol144.f06').open('rb') as stream:
        f06_hash = hashlib.file_digest(stream,'sha256').hexdigest()
    return {'shell_elements':len(ids),'subcases':subcases,'rows_per_file':len(ids),'columns_per_file':1+3*len(subcases),
            'components':['exx','eyy','gamma_xy'],'coordinates':'Element axes; centroid values; no nodal averaging.',
            'surfaces':{'top':'+element-z fiber','bottom':'-element-z fiber','mid':'(top + bottom) / 2, component-wise'},
            'units':'Dimensionless strain; gamma_xy is engineering shear strain.',
            'precision':'Nastran F06 printed precision.','csv_sha256':file_hashes,'source_f06_sha256':f06_hash,
            'input_sha256':hashes,'solver_seconds':summary['elapsed_seconds'],'runtime_diagnostics':sorted(errors),
            'warning_messages':summary['warning_messages'],'solver_completed':True}
