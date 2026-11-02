"""Create separate Femap models containing actual SOL144 AOA 3--6 results.

Each OP2 may use SUBCASE 1. Independent Femap result sets are required and are
checked after every import. No analysis input is exported or overwritten.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import shutil

from femap_connection import ROOT, check, connect, ids
from view_femap_aero_panels import entity_ids
from view_femap_sol144_results import (
    catalog_results, choose_strain, require_completed_run, sha256,
    strain_view, summarize_vector,
)

ANGLES = (3, 4, 5, 6)


def result_ids(app) -> list[int]:
    results = app.feResults
    if not results.NumberOfSets():
        return []
    selected = app.feSet
    check(results.Sets(selected.ID, True), 'list result sets')
    return ids(selected)


def sample_preserved(app, record: dict) -> None:
    output = app.feOutputSet
    check(output.Get(record['output_set_id']), 'read retained sweep output set')
    if output.title != record['output_set_title']:
        raise RuntimeError('A later OP2 import replaced an earlier result-set label')
    sample = record['strain_statistics']['maximum_absolute']
    rc, value = app.feResults.EntityValueV2(record['output_set_id'], record['vector_id'], sample['element_id'])
    check(rc, 'read retained sweep strain sample')
    if not math.isclose(value, sample['strain'], rel_tol=1e-10, abs_tol=1e-14):
        raise RuntimeError('A later OP2 import replaced an earlier result value')


def check_inputs(source_model: Path, run_root: Path, part: int, allow_diagnostic: bool) -> list[dict]:
    inputs = []
    for angle in ANGLES:
        op2 = run_root / f'aoa_{angle:02d}/part_{part}/bwb_sol144.op2'
        run = require_completed_run(op2, source_model, allow_diagnostic)
        if (not run.get('fresh_f06') or run.get('fatal_messages')
                or not run.get('output_sections', {}).get('strains')):
            raise RuntimeError(f'AOA {angle}: incomplete or invalid solver evidence')
        if run['status'] == 'completed' and run.get('runtime_errors'):
            raise RuntimeError(f'AOA {angle}: clean status conflicts with recorded errors')
        inputs.append({'aoa_deg': angle, 'op2': op2, 'run': run,
                       'op2_sha256': sha256(op2),
                       'run_summary_sha256': sha256(op2.parent / 'run-summary.json')})
    return inputs


def build(app, part: int, run_root: Path, output: Path, allow_diagnostic: bool = False) -> dict:
    source = ROOT / f'BWB/models/aero_panels/part_{part}/part_{part}_aero_structure.modfem'
    inputs = check_inputs(source, run_root, part, allow_diagnostic)
    source_hash = sha256(source)
    if output.exists():
        raise FileExistsError(f'Use a new output path to preserve previous sweep results: {output}')
    output.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, output)
    check(app.feFileOpen(False, str(output)), 'open separate sweep model copy')
    if result_ids(app):
        raise RuntimeError('Reference model contains old results')
    counts_before = {name: len(entity_ids(app, kind)) for name, kind in
                     [('nodes', 7), ('elements', 8), ('aero_panels', 95), ('aero_splines', 98)]}
    records = []
    for case in inputs:
        prior_ids = set(result_ids(app))
        check(app.feFileReadNastranResults(8, str(case['op2'])), f"import AOA {case['aoa_deg']} OP2")
        current_ids = set(result_ids(app))
        new_ids = current_ids - prior_ids
        if len(new_ids) != 1 or not prior_ids.issubset(current_ids):
            raise RuntimeError('Each OP2 must append exactly one independent result set')
        set_id = new_ids.pop()
        selected, vector = choose_strain(catalog_results(app), set_id, None)
        statistics = summarize_vector(app, set_id, vector['id'])
        if statistics['element_count'] != counts_before['elements']:
            raise RuntimeError('Strain output does not cover all structural shell elements')
        title = f"AOA {case['aoa_deg']} deg | Part {part} | SOL144"
        output_set = app.feOutputSet
        check(output_set.Get(set_id), 'read imported sweep result metadata')
        output_set.title = title
        output_set.value = float(case['aoa_deg'])
        output_set.notes = (
            f"Angle of attack = {case['aoa_deg']} degrees. Output-set value records AOA in degrees.\n"
            f"Imported actual OP2: {case['op2'].relative_to(ROOT).as_posix()}\n"
            f"Original solver status: {case['run']['status']}.\n"
            + '\n'.join(case['run'].get('runtime_errors', []))
        )
        check(output_set.Put(set_id), 'name independent sweep output set')
        record = {
            'aoa_deg': case['aoa_deg'], 'output_set_id': set_id, 'output_set_title': title,
            'vector_id': vector['id'], 'vector_title': vector['title'],
            'imported_vector_count': selected['vector_count'],
            'source_op2': case['op2'].relative_to(ROOT).as_posix(),
            'source_op2_sha256': case['op2_sha256'],
            'source_run_summary_sha256': case['run_summary_sha256'],
            'source_case': case['run']['source_case'],
            'original_solver_status': case['run']['status'],
            'runtime_diagnostics': case['run'].get('runtime_errors', []),
            'solver_completed': case['run']['solver_completed'],
            'fresh_f06': case['run']['fresh_f06'],
            'strain_statistics': statistics,
        }
        records.append(record)
        for prior in records:
            sample_preserved(app, prior)
        print(f"Part {part}, AOA {case['aoa_deg']}: native set {set_id}, {statistics['element_count']} shell strains", flush=True)

    final_case = records[-1]
    view = strain_view(app, final_case['output_set_id'], final_case['vector_id'], final_case['vector_title'])
    native_view = app.feView
    check(native_view.Get(view['id']), 'read final sweep view')
    native_view.title = f'BWB Part {part} | AOA 6 deg | Plate Top MajorPrn Strain [m/m]'
    check(native_view.Put(view['id']), 'name final sweep view')
    app.feViewRegenerate(view['id'])
    counts_after = {name: len(entity_ids(app, kind)) for name, kind in
                    [('nodes', 7), ('elements', 8), ('aero_panels', 95), ('aero_splines', 98)]}
    if counts_before != counts_after or len(result_ids(app)) != len(ANGLES):
        raise RuntimeError('Sweep import changed the model or did not preserve four result sets')
    check(app.feFileSaveAs(False, str(output)), 'save four-case sweep result model')
    check(app.feFileClose(False), 'close saved sweep result model')
    check(app.feFileOpen(False, str(output)), 'reopen saved sweep result model')
    if set(result_ids(app)) != {record['output_set_id'] for record in records}:
        raise RuntimeError('Sweep result-set IDs changed after reopening')
    for record in records:
        rc, location = app.feResults.SetLocation(record['output_set_id'])
        check(rc, 'check persisted sweep result location')
        if location != 1:
            raise RuntimeError('Sweep results depend on an external result attachment')
        sample_preserved(app, record)
    rc, active_view_id = app.feAppGetActiveView()
    check(rc, 'get persisted sweep view')
    native_view = app.feView
    check(native_view.Get(active_view_id), 'read persisted sweep view')
    if (native_view.OutputSet != final_case['output_set_id'] or native_view.Contour != 1
            or native_view.ContourDataV2 != final_case['vector_id']):
        raise RuntimeError('Default AOA 6 strain contour did not persist')
    if sha256(source) != source_hash or any(sha256(case['op2']) != case['op2_sha256'] for case in inputs):
        raise RuntimeError('A source model or OP2 changed during sweep import')
    report = {
        'part': part, 'result_model': output.relative_to(ROOT).as_posix(),
        'source_model': source.relative_to(ROOT).as_posix(), 'source_model_sha256': source_hash,
        'angles_deg': list(ANGLES), 'result_sets': records, 'entity_counts': counts_after,
        'independent_result_sets': 4, 'all_prior_result_samples_preserved': True,
        'reviewed_empty_matrix_diagnostic_allowed': allow_diagnostic,
        'diagnostic_policy': 'Original statuses retained; only the explicitly reviewed DGEMM parameter 13 and launcher 1 pair can bypass the clean-status gate',
        'default_display': {'aoa_deg': 6, **view},
        'persistence_check': {'reopened_from_disk': True, 'all_four_sets_in_femap_database': True,
                              'all_four_labels_and_strain_samples_retained': True,
                              'default_aoa_6_strain_contour_retained': True},
        'source_model_unchanged': True, 'source_op2_files_unchanged': True,
        'analysis_deck_exported': False,
    }
    report_path = output.with_name(output.stem + '_femap_report.json')
    report_path.write_text(json.dumps(report, indent=2), encoding='utf-8')
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-root', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'BWB/results/aoa_sweep')
    parser.add_argument('--parts', nargs='+', type=int, choices=(1, 2), default=[1, 2])
    parser.add_argument('--allow-reviewed-empty-matrix-diagnostic', action='store_true')
    args = parser.parse_args()
    run_root, output_dir = args.run_root.resolve(strict=True), args.output_dir.resolve()
    if not output_dir.is_relative_to(ROOT / 'BWB/results'):
        parser.error('Use a dedicated BWB/results output directory')
    for part in args.parts:
        check_inputs(ROOT / f'BWB/models/aero_panels/part_{part}/part_{part}_aero_structure.modfem',
                     run_root, part, args.allow_reviewed_empty_matrix_diagnostic)
    app = connect(True)
    previous = (app.DialogAutoAnswer, app.DialogAutoSkip, app.DialogAutoSkipMsg,
                app.Pref_UseSubcaseIdForOutput)
    try:
        app.DialogAutoAnswer = app.DialogAutoSkip = app.DialogAutoSkipMsg = 1
        app.Pref_UseSubcaseIdForOutput = False
        for part in args.parts:
            report = build(app, part, run_root, output_dir / f'part_{part}_aoa_3to6.modfem',
                           args.allow_reviewed_empty_matrix_diagnostic)
            print(json.dumps({key: value for key, value in report.items() if key != 'result_sets'}, indent=2))
    finally:
        (app.DialogAutoAnswer, app.DialogAutoSkip, app.DialogAutoSkipMsg,
         app.Pref_UseSubcaseIdForOutput) = previous


if __name__ == '__main__':
    main()
