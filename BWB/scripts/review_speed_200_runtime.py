"""Record the actual at-speed MKL trace and compare all printed result values."""
from __future__ import annotations

import json
from pathlib import Path
import re

from BWB.paths import REPOSITORY_ROOT as ROOT
from .export_sol144_strain import export_case
from .run_aoa_sweep import sha, EXPECTED_DIAGNOSTICS
from .build_bwb_aero_panels import provenance


def review():
    normal = ROOT / 'workspace/BWB/speed_200/run_01/aoa_06/part_1'
    trace = ROOT / 'workspace/BWB/speed_200/aoa_06_part_1_trace'
    normal_summary = json.loads((normal / 'run-summary.json').read_text())
    trace_summary = json.loads((trace / 'run-summary.json').read_text())
    for summary in (normal_summary, trace_summary):
        assert summary['solver_completed'] and summary['fresh_f06']
        assert not summary['fatal_messages'] and not summary['warning_messages']
        assert set(summary['runtime_errors']) == EXPECTED_DIAGNOSTICS
    assert normal_summary['input_sha256'] == trace_summary['input_sha256']
    context, waiting = [], False
    with (trace / 'bwb_sol144.log').open(errors='replace') as stream:
        for number, line in enumerate(stream, 1):
            if 'Intel MKL ERROR:' in line:
                assert 'Parameter 13 was incorrect on entry to DGEMM' in line
                context.append(dict(line=number, text=line.strip()))
                waiting = True
            elif waiting and 'MKL_VERBOSE DGEMM' in line:
                context.append(dict(line=number, text=line.strip()))
                args = re.search(r'DGEMM\(([^)]+)\)', line).group(1).split(',')
                assert args[:5] == ['T', 'N', '0', '1', '3'] and args[-1] == '0'
                waiting = False
    assert len(context) == 2 and not waiting
    trace_export = trace / 'export_verification'
    export_case(trace, trace_export, make_plot=False)
    normal_export = ROOT / 'BWB/results/speed_200/aoa_06/part_1'
    compared = {}
    for name in ('shell_strain_long.csv', 'grid_displacements.csv', 'solver_generated_aero_displacements.csv'):
        compared[name] = dict(normal_sha256=sha(normal_export / name), trace_sha256=sha(trace_export / name))
        assert compared[name]['normal_sha256'] == compared[name]['trace_sha256']
    result = {
        'status': 'at_speed_empty_output_call_reviewed; raw runtime diagnostic retained',
        'velocity_m_s': 200., 'angle_of_attack_deg': 6, 'part': 1,
        'log_path': (trace / 'bwb_sol144.log').relative_to(ROOT).as_posix(),
        'log_sha256': sha(trace / 'bwb_sol144.log'), 'context': context,
        'decoded_DGEMM': dict(m=0, n=1, k=3, lda=3, ldb=3, ldc=0),
        'interpretation': 'The diagnosed DGEMM invocation has zero output rows and invalid ldc=0. This is the same empty-output diagnostic observed at 20 m/s.',
        'same_inputs_between_normal_and_trace': True, 'all_exported_result_csvs_byte_identical': True,
        'result_comparison': compared,
        'f06_sha256': dict(normal=sha(normal / 'bwb_sol144.f06'), trace=sha(trace / 'bwb_sol144.f06')),
        'trace_run_summary': trace_summary,
        'scope': 'Directly traced at 200 m/s, AOA 6, part 1. Other seven cases share solver, Mach, q, trim-variable count and error signature but were not individually verbose-traced.',
        'limitations': 'Numerical repeatability does not validate large deformation, yielding, flow separation, buckling or divergence stability.',
        'review_script': provenance(Path(__file__)),
    }
    output = ROOT / 'BWB/results/speed_200/runtime_diagnostic_review.json'
    output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    for review_path in output.parent.glob('aoa_*/part_*/case_review.json'):
        case_review = json.loads(review_path.read_text())
        assert set(case_review['runtime_errors_retained']) == EXPECTED_DIAGNOSTICS
        case_review['runtime_diagnostic_review'] = provenance(output)
        case_review['runtime_review_scope'] = result['scope']
        # This update documents the runtime call; it never overrides physical failures.
        review_path.write_text(json.dumps(case_review, indent=2) + '\n', encoding='utf-8')
    print(output)


if __name__ == '__main__':
    review()
