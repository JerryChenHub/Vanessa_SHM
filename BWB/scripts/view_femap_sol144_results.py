"""Import a completed BWB SOL144 OP2 into a separate native Femap result model.

The source model and solver inputs remain unchanged. This script never exports
an analysis deck or launches a solver. The full external BDF, including W2GJ,
remains the authoritative analysis input. Femap contour values are strain in
m/m; use the reported vector title to distinguish shell side and component.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import shutil

from femap_connection import ROOT, check, connect, ids
from view_femap_aero_panels import entity_ids


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require_completed_run(op2: Path, source_model: Path, allow_reviewed_diagnostic: bool = False) -> dict:
    summary_path = op2.parent / 'run-summary.json'
    summary = json.loads(summary_path.read_text(encoding='utf-8'))
    reviewed_messages = {
        'Intel MKL ERROR: Parameter 13 was incorrect on entry to DGEMM .',
        'Launcher returned exit code 1',
    }
    diagnostic_allowed = (
        allow_reviewed_diagnostic
        and summary.get('status') == 'completed_with_runtime_errors'
        and len(summary.get('runtime_errors', [])) == 2
        and set(summary['runtime_errors']) == reviewed_messages
        and summary.get('fresh_f06') is True
        and not summary.get('fatal_messages')
        and summary.get('output_sections', {}).get('strains') is True
    )
    if not summary.get('solver_completed') or not (summary.get('status') == 'completed' or diagnostic_allowed):
        raise RuntimeError('Only a verified completed solver run can be imported as a result')
    if not op2.is_file() or not op2.stat().st_size:
        raise FileNotFoundError(f'Nonempty OP2 missing: {op2}')
    if op2.stat().st_mtime < float(summary['started_unix']) - 2:
        raise RuntimeError('OP2 timestamp predates the recorded solver invocation')
    model_part = re.search(r'part_([12])', str(source_model))
    source_case = str(summary.get('source_case', '')).replace('\\', '/')
    if model_part and f'/part_{model_part.group(1)}' not in source_case:
        raise RuntimeError('Solver case part does not match the Femap reference model')
    return summary


def catalog_results(app) -> list[dict]:
    result = []
    output = app.feOutputSet
    vectors = app.feResults
    selected_sets = app.feSet
    check(vectors.Sets(selected_sets.ID, True), 'list imported result sets')
    for set_id in ids(selected_sets):
        check(output.Get(set_id), 'read imported result set')
        rc, count, vector_ids, titles = vectors.VectorTitlesV2(set_id, False, 0, 0)
        check(rc, 'read imported result vector titles')
        inventory = [{'id': int(vector_id), 'title': title}
                     for vector_id, title in zip(vector_ids or (), titles or ())]
        if count != len(inventory):
            raise RuntimeError('Result vector title inventory is incomplete')
        result.append({'id': set_id, 'title': output.title,
                       'nastran_subcase': output.nas_case,
                       'analysis_type': output.analysis,
                       'vector_count': count, 'vectors': inventory})
    return result


def choose_strain(catalog: list[dict], set_id: int | None, vector_id: int | None) -> tuple[dict, dict]:
    candidates = []
    for output in catalog:
        if set_id is not None and output['id'] != set_id:
            continue
        for vector in output['vectors']:
            title = vector['title'].lower()
            if 'strain' not in title:
                continue
            if vector_id is not None:
                if vector['id'] == vector_id:
                    candidates.append((output, vector))
            elif 'top' in title and 'major' in title and 'corner' not in title:
                candidates.append((output, vector))
    if len(candidates) != 1:
        choices = [(output['id'], vector['id'], vector['title'])
                   for output in catalog for vector in output['vectors']
                   if 'strain' in vector['title'].lower()]
        raise RuntimeError(f'Choose an unambiguous strain set/vector with CLI arguments. Available: {choices}')
    return candidates[0]


def summarize_vector(app, set_id: int, vector_id: int) -> dict:
    results = app.feResults
    vector_info = results.VectorInfoV2(set_id, vector_id)
    check(vector_info[0], 'read strain vector data location')
    if vector_info[1] != 8:
        raise RuntimeError('Expected element strain output, not a nodal approximation')
    valid = app.feSet
    check(results.VectorEntitiesV2(set_id, vector_id, valid.ID, True), 'query strain result locations')
    valid_ids = ids(valid)
    if not valid_ids:
        raise RuntimeError('Selected strain vector has no actual result locations')
    check(results.clear(), 'clear result extraction workspace')
    check(results.DataNeeded(8, valid.ID), 'select element strain data')
    rc, count, columns = results.AddColumnV2(set_id, vector_id, False)
    check(rc, 'select strain result column')
    if count != 1:
        raise RuntimeError('Expected one scalar strain column')
    check(results.Populate(), 'read actual strain values')
    rc, result_ids, values = results.GetColumn(int(columns[0]))
    check(rc, 'extract imported strain values')
    valid_id_set = set(valid_ids)
    pairs = [(int(eid), float(value)) for eid, value in zip(result_ids, values)
             if int(eid) in valid_id_set]
    if len(pairs) != len(valid_ids) or any(not math.isfinite(value) for _, value in pairs):
        raise RuntimeError('Strain data are incomplete or nonfinite')
    low, high = min(pairs, key=lambda pair: pair[1]), max(pairs, key=lambda pair: pair[1])
    absolute = max(pairs, key=lambda pair: abs(pair[1]))
    return {'element_count': len(pairs), 'units': 'm/m',
            'minimum': {'element_id': low[0], 'strain': low[1], 'microstrain': low[1] * 1e6},
            'maximum': {'element_id': high[0], 'strain': high[1], 'microstrain': high[1] * 1e6},
            'maximum_absolute': {'element_id': absolute[0], 'strain': absolute[1],
                                 'microstrain': absolute[1] * 1e6}}


def strain_view(app, set_id: int, vector_id: int, title: str) -> dict:
    rc, view_id = app.feAppGetActiveView()
    check(rc, 'get result view')
    view = app.feView
    check(view.Get(view_id), 'read result view')
    view.title = ('BWB SOL144 | ' + title + ' | strain m/m')[:79]
    view.OutputSet, view.ContourDataV2 = set_id, vector_id
    view.Contour, view.Deformed, view.Mode = 1, 0, 3
    view.ContourElemUnaveraged = True
    view.ContourCornerData = False
    view.group, view.BC, view.load = 0, 0, 0
    for option in (9, 40, 48):
        view.SetDraw(option, True)
    for option in (2, 3, 4, 5, 7, 8, 102, 103, 104, 105, 108):
        view.SetDraw(option, False)
    view.SetLabel(9, 0)
    check(view.Put(view_id), 'save strain contour view')
    orientation = app.feViewOrient
    check(orientation.Get(view_id), 'get result view orientation')
    orientation.Perspective = False
    orientation.SetRotationAngles(-65, 20, -25)
    check(orientation.Put(view_id), 'set result view orientation')
    app.feAppVisible(True)
    app.feViewAutoscaleAll(view_id, True)
    app.feViewRegenerate(view_id)
    check(view.Get(view_id), 'verify saved strain view')
    if view.OutputSet != set_id or view.ContourDataV2 != vector_id or view.Contour != 1:
        raise RuntimeError('Femap did not retain the selected strain contour')
    return {'id': view_id, 'result_set': set_id, 'vector_id': vector_id, 'vector_title': title,
            'contour_on': True, 'unaveraged_element_contour': True,
            'deformation_display': False, 'units': 'm/m',
            'shell_side_note': 'Top is the positive shell-normal face, not the aircraft upper skin'}


def verify_saved_result(app, output: Path, set_id: int, vector_id: int, sample: dict) -> dict:
    check(app.feFileClose(False), 'close saved result model for persistence check')
    check(app.feFileOpen(False, str(output)), 'reopen result model from disk')
    results = app.feResults
    rc, location = results.SetLocation(set_id)
    check(rc, 'read persisted result storage location')
    if location != 1:
        raise RuntimeError('Result model depends on an external result attachment')
    rc, value = results.EntityValueV2(set_id, vector_id, sample['element_id'])
    check(rc, 'read persisted actual strain sample')
    if not math.isclose(value, sample['strain'], rel_tol=1e-10, abs_tol=1e-14):
        raise RuntimeError('Strain value changed after model save/reopen')
    rc, view_id = app.feAppGetActiveView()
    check(rc, 'read persisted active view')
    view = app.feView
    check(view.Get(view_id), 'read persisted contour options')
    if view.OutputSet != set_id or view.ContourDataV2 != vector_id or view.Contour != 1:
        raise RuntimeError('Strain contour selection was not retained on disk')
    return {'reopened_from_disk': True, 'results_in_femap_database': True,
            'sample_strain_retained': True, 'strain_contour_retained': True}


def build(app, source_model: Path, op2: Path, output: Path,
          set_id: int | None = None, vector_id: int | None = None,
          allow_reviewed_diagnostic: bool = False) -> dict:
    run = require_completed_run(op2, source_model, allow_reviewed_diagnostic)
    if output.exists():
        raise FileExistsError(f'Use a new result model path to preserve previous results: {output}')
    source_hash, op2_hash = sha256(source_model), sha256(op2)
    output.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source_model, output)
    check(app.feFileOpen(False, str(output)), 'open separate copy of corresponding Femap model')
    before = {name: len(entity_ids(app, kind)) for name, kind in
              [('nodes', 7), ('elements', 8), ('aero_panels', 95), ('aero_splines', 98)]}
    if app.feResults.NumberOfSets():
        raise RuntimeError('Reference model already contains results; use a clean matching model')
    check(app.feFileReadNastranResults(8, str(op2)), 'import completed SOL144 OP2 results')
    catalog = catalog_results(app)
    (output.parent / 'femap_result_report.json').write_text(json.dumps(
        {'state': 'imported_awaiting_view_validation', 'catalog': catalog,
         'original_solver_status': run['status'], 'runtime_diagnostics': run.get('runtime_errors', [])},
        indent=2), encoding='utf-8')
    selected_set, selected_vector = choose_strain(catalog, set_id, vector_id)
    summary = summarize_vector(app, selected_set['id'], selected_vector['id'])
    view = strain_view(app, selected_set['id'], selected_vector['id'], selected_vector['title'])
    after = {name: len(entity_ids(app, kind)) for name, kind in
             [('nodes', 7), ('elements', 8), ('aero_panels', 95), ('aero_splines', 98)]}
    if before != after:
        raise RuntimeError('Result import changed geometry/entity counts')
    check(app.feFileSaveAs(False, str(output)), 'save independent Femap result model')
    persistence = verify_saved_result(app, output, selected_set['id'], selected_vector['id'],
                                      summary['maximum_absolute'])
    if sha256(source_model) != source_hash or sha256(op2) != op2_hash:
        raise RuntimeError('A source artifact changed during results import')
    report = {'source_model': source_model.relative_to(ROOT).as_posix(),
              'source_model_sha256': source_hash, 'source_op2': op2.relative_to(ROOT).as_posix(),
              'source_op2_sha256': op2_hash, 'source_case': run['source_case'],
              'solver_completed': True, 'result_model': output.relative_to(ROOT).as_posix(),
              'original_solver_status': run['status'], 'runtime_diagnostics': run.get('runtime_errors', []),
              'reviewed_empty_matrix_diagnostic_allowed': allow_reviewed_diagnostic,
              'diagnostic_review_note': ('Explicit review identified DGEMM(T,N,0,1,2,...,ldc0) '
                                         'with an empty output matrix; original solver diagnostic is retained')
                                        if allow_reviewed_diagnostic else None,
              'entity_counts': after, 'imported_result_sets': len(catalog),
              'catalog': catalog,
              'selected_vector_statistics': summary, 'view': view,
              'persistence_check': persistence,
              'source_model_unchanged': True, 'source_op2_unchanged': True,
              'analysis_deck_exported': False,
              'analysis_authority': 'External SOL144 BDF and W2GJ input used by the completed solver run'}
    (output.parent / 'femap_result_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', type=Path, required=True)
    parser.add_argument('--op2', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--output-set', type=int)
    parser.add_argument('--vector-id', type=int)
    parser.add_argument('--allow-reviewed-empty-matrix-diagnostic', action='store_true',
                        help='Explicitly accept only the reviewed DGEMM parameter 13 + launcher 1 diagnostic pair')
    args = parser.parse_args()
    source, op2, output = args.model.resolve(strict=True), args.op2.resolve(strict=True), args.output.resolve()
    if output.suffix.lower() != '.modfem' or not output.is_relative_to(ROOT / 'BWB'):
        parser.error('Use a new BWB .modfem output path')
    if output == source:
        parser.error('Source and result models must be different files')
    require_completed_run(op2, source, args.allow_reviewed_empty_matrix_diagnostic)
    app = connect(True)
    previous = (app.DialogAutoAnswer, app.DialogAutoSkip, app.DialogAutoSkipMsg)
    try:
        app.DialogAutoAnswer = app.DialogAutoSkip = app.DialogAutoSkipMsg = 1
        report = build(app, source, op2, output, args.output_set, args.vector_id,
                       args.allow_reviewed_empty_matrix_diagnostic)
        print(json.dumps({key: value for key, value in report.items() if key != 'catalog'}, indent=2))
    finally:
        app.DialogAutoAnswer, app.DialogAutoSkip, app.DialogAutoSkipMsg = previous


if __name__ == '__main__':
    main()
