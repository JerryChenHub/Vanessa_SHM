"""View a new speed sweep without changing previous Femap result models.

Reuses the verified four-OP2 importer. Adds explicit speed/angle labels and
stress contours while retaining every imported OP2 vector and the original
solver diagnostics. Strength and small-deformation screens are interpretation
checks, not conditions for discarding an available linear solution.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from femap_connection import ROOT, check, connect, ids
from view_femap_aoa_sweep import ANGLES, build as build_angle_sweep, check_inputs, result_ids, sample_preserved
from view_femap_sol144_results import catalog_results, sha256, strain_view

LINEAR_NOTE = (
    'LINEAR SOL144 RESULT ONLY. Linear elastic material and small-displacement structural theory; '
    'no plasticity, geometric nonlinearity, or failure prediction. '
    'A strength or deformation screen may invalidate physical interpretation without deleting this computed field.'
)


def read_review_evidence(review_dir: Path, part: int, run_root: Path, speed: float) -> tuple[dict, dict[int, dict]]:
    runtime_path = review_dir / 'runtime_diagnostic_review.json'
    runtime = json.loads(runtime_path.read_text(encoding='utf-8'))
    if (runtime.get('status') != 'at_speed_empty_output_call_reviewed; raw runtime diagnostic retained'
            or not runtime.get('same_inputs_between_normal_and_trace')
            or not runtime.get('all_exported_result_csvs_byte_identical')
            or runtime.get('decoded_DGEMM', {}).get('m') != 0
            or runtime.get('decoded_DGEMM', {}).get('ldc') != 0):
        raise RuntimeError('Expected the explicit at-speed empty-matrix diagnostic review')
    if not math.isclose(runtime.get('velocity_m_s', math.nan), speed, rel_tol=1e-12):
        raise RuntimeError('Runtime diagnostic review velocity disagrees with the requested speed')
    traced_f06 = run_root / f"aoa_{runtime['angle_of_attack_deg']:02d}/part_{runtime['part']}/bwb_sol144.f06"
    if runtime.get('f06_sha256', {}).get('normal') != sha256(traced_f06):
        raise RuntimeError('Runtime diagnostic review does not match the actual run F06')
    review = {'path': runtime_path.relative_to(ROOT).as_posix(), 'sha256': sha256(runtime_path),
              'status': runtime['status'], 'scope': runtime['scope'], 'limitations': runtime['limitations']}
    cases = {}
    for angle in ANGLES:
        folder = review_dir / f'aoa_{angle:02d}/part_{part}'
        case_path, reason_path = folder / 'case_review.json', folder / 'reasonableness.json'
        case = json.loads(case_path.read_text(encoding='utf-8'))
        reason = json.loads(reason_path.read_text(encoding='utf-8'))
        if not case.get('result_integrity_checks_pass'):
            raise RuntimeError(f'Result integrity review did not pass: {case_path}')
        if (case.get('part') != part or case.get('aoa_deg') != angle
                or not math.isclose(case.get('velocity_m_s', math.nan), speed, rel_tol=1e-12)
                or reason.get('part') != part):
            raise RuntimeError(f'Incorrect case review identity: {case_path}')
        actual_f06_hash = sha256(run_root / f'aoa_{angle:02d}/part_{part}/bwb_sol144.f06')
        if (case.get('source_f06_sha256') != actual_f06_hash
                or reason.get('source_f06_sha256') != actual_f06_hash):
            raise RuntimeError(f'Review source F06 hash disagrees with the actual run: {case_path}')
        if case.get('reasonableness_sha256') != sha256(reason_path):
            raise RuntimeError(f'Reasonableness file changed since case review: {reason_path}')
        physical = reason['subcases'][0]
        cases[angle] = {
            'case_review': case_path.relative_to(ROOT).as_posix(), 'case_review_sha256': sha256(case_path),
            'reasonableness': reason_path.relative_to(ROOT).as_posix(), 'reasonableness_sha256': sha256(reason_path),
            'status': case['status'], 'result_integrity_checks_pass': True,
            'source_f06_sha256': actual_f06_hash,
            'linear_physics_screens_pass': case['linear_physics_screens_pass'],
            'failed_physical_screens': case['failed_physical_screens'],
            'max_translation_m': physical['max_translation_m'],
            'max_translation_over_reference_span': physical['max_translation_over_reference_span'],
            'max_rotation_deg': physical['max_rotation_deg'],
            'assumed_yield_MPa': physical['assumed_yield_MPa'],
            'peak_von_mises_MPa_from_F06_strains': physical['peak_von_mises_MPa'],
            'limitations': case['limitations'],
        }
    return review, cases


def verify_flight_labels(run_root: Path, part: int, speed: float) -> list[dict]:
    records = []
    for angle in ANGLES:
        folder = run_root / f'aoa_{angle:02d}/part_{part}'
        manifest_path = folder / 'case_manifest.json'
        data = json.loads(manifest_path.read_text(encoding='utf-8'))
        flight = data['flight']
        if (data['part'] != part or not math.isclose(flight['velocity_m_s'], speed, rel_tol=1e-12)
                or not math.isclose(flight['angle_of_attack_deg'], angle, abs_tol=1e-12)):
            raise RuntimeError(f'Speed/angle/part label disagrees with actual case manifest: {manifest_path}')
        records.append({'flight': flight, 'manifest': manifest_path,
                        'manifest_sha256': sha256(manifest_path)})
    return records


def find_scalar(catalog: list[dict], set_id: int, result_name: str) -> dict:
    needles = {
        'top-vonmises-stress': ('plate', 'top', 'vonmises', 'stress'),
        'bottom-vonmises-stress': ('plate', 'bot', 'vonmises', 'stress'),
        'top-major-strain': ('plate', 'top', 'major', 'strain'),
    }[result_name]
    matches = [vector for output in catalog if output['id'] == set_id for vector in output['vectors']
               if all(word in vector['title'].lower() for word in needles)
               and 'corner' not in vector['title'].lower()]
    if len(matches) != 1:
        raise RuntimeError(f'Expected one {result_name} vector in output set {set_id}: {matches}')
    return matches[0]


def scalar_statistics(app, set_id: int, vector_id: int, is_stress: bool) -> dict:
    results = app.feResults
    vector_info = results.VectorInfoV2(set_id, vector_id)
    check(vector_info[0], 'read scalar result metadata')
    if vector_info[1] != 8:
        raise RuntimeError('Expected an element-based shell scalar')
    selected = app.feSet
    check(results.VectorEntitiesV2(set_id, vector_id, selected.ID, True), 'find actual scalar result locations')
    valid = set(ids(selected))
    check(results.clear(), 'clear result extraction workspace')
    check(results.DataNeeded(8, selected.ID), 'select shell scalar data')
    rc, count, columns = results.AddColumnV2(set_id, vector_id, False)
    check(rc, 'select scalar result column')
    if count != 1:
        raise RuntimeError('Expected one scalar result column')
    check(results.Populate(), 'read actual scalar result values')
    rc, element_ids, values = results.GetColumn(int(columns[0]))
    check(rc, 'extract scalar values')
    pairs = [(int(eid), float(value)) for eid, value in zip(element_ids, values) if int(eid) in valid]
    if not valid or len(pairs) != len(valid) or any(not math.isfinite(value) for _, value in pairs):
        raise RuntimeError('Scalar field is missing or nonfinite')

    def item(pair):
        return {'element_id': pair[0], 'value': pair[1],
                'MPa' if is_stress else 'microstrain': pair[1] / 1e6 if is_stress else pair[1] * 1e6}

    return {'element_count': len(pairs), 'units': 'Pa' if is_stress else 'm/m',
            'minimum': item(min(pairs, key=lambda pair: pair[1])),
            'maximum': item(max(pairs, key=lambda pair: pair[1])),
            'maximum_absolute': item(max(pairs, key=lambda pair: abs(pair[1])))}


def check_display_sample(app, record: dict) -> None:
    sample = record['display_statistics']['maximum_absolute']
    rc, value = app.feResults.EntityValueV2(record['output_set_id'], record['display_vector_id'], sample['element_id'])
    check(rc, 'read retained display scalar')
    if not math.isclose(value, sample['value'], rel_tol=1e-10, abs_tol=1e-14):
        raise RuntimeError('Display scalar value changed after save/reopen')


def build(app, part: int, run_root: Path, output: Path, speed: float,
          default_result: str, allow_diagnostic: bool = False,
          review_dir: Path = ROOT / 'BWB/results/speed_200') -> dict:
    flights = verify_flight_labels(run_root, part, speed)
    runtime_review, physical_reviews = read_review_evidence(review_dir, part, run_root, speed)
    report = build_angle_sweep(app, part, run_root, output, allow_diagnostic)
    catalog = catalog_results(app)
    for record, manifest in zip(report['result_sets'], flights):
        set_id = record['output_set_id']
        title = f"{speed:g} m/s | AOA {record['aoa_deg']} deg | Part {part} | LINEAR SOL144"
        native = app.feOutputSet
        check(native.Get(set_id), 'read speed sweep output set')
        native.title = title
        physical = physical_reviews[record['aoa_deg']]
        native.notes += (f"\nVelocity = {speed:g} m/s. Mach = {manifest['flight']['mach']:.9g}. "
                         f"Dynamic pressure = {manifest['flight']['dynamic_pressure_Pa']:.9g} Pa.\n"
                         f"Failed physical screens: {', '.join(physical['failed_physical_screens'])}.\n"
                         f"Maximum translation = {physical['max_translation_m']:.6g} m; "
                         f"maximum rotation = {physical['max_rotation_deg']:.6g} deg.\n"
                         f"Assumed yield = {physical['assumed_yield_MPa']:.6g} MPa. "
                         "Stress vectors are in Pa; MPa = Pa / 1,000,000.\n"
                         "Top/Bottom identify shell-normal fibers, not aircraft upper/lower skins.\n" + LINEAR_NOTE)
        record['output_set_notes'] = native.notes
        check(native.Put(set_id), 'label speed and applicability in output set')
        record['output_set_title'] = title
        record['physical_review'] = physical
        record['flight'] = manifest['flight']
        record['case_manifest_sha256'] = manifest['manifest_sha256']
        display_vector = find_scalar(catalog, set_id, default_result)
        record['display_vector_id'] = display_vector['id']
        record['display_vector_title'] = display_vector['title']
        record['display_statistics'] = scalar_statistics(app, set_id, display_vector['id'], 'stress' in default_result)
        record['von_mises_stress_fibers'] = {}
        for side, key in (('top', 'top-vonmises-stress'), ('bottom', 'bottom-vonmises-stress')):
            vector = find_scalar(catalog, set_id, key)
            statistics = scalar_statistics(app, set_id, vector['id'], True)
            if statistics['element_count'] != report['entity_counts']['elements']:
                raise RuntimeError('Stress output does not cover every shell element')
            record['von_mises_stress_fibers'][side] = {'vector_id': vector['id'], 'vector_title': vector['title'], **statistics}
        sample_preserved(app, record)
        check_display_sample(app, record)

    final = report['result_sets'][-1]
    # The existing helper configures an unaveraged scalar contour; correct its
    # strain-specific caption/units before saving this stress-capable model.
    view = strain_view(app, final['output_set_id'], final['display_vector_id'], final['display_vector_title'])
    view['units'] = final['display_statistics']['units']
    view['vector_title'] = final['display_vector_title']
    view['shell_side_note'] = 'Top/Bottom follow positive/negative shell-normal fibers, not aircraft upper/lower skins.'
    native_view = app.feView
    check(native_view.Get(view['id']), 'read scalar contour view')
    native_view.title = (f'{speed:g} m/s | AOA 6 deg | Part {part} | '
                         f"{final['display_vector_title']} [{view['units']}]")[:79]
    check(native_view.Put(view['id']), 'save speed-aware scalar contour caption')
    app.feViewRegenerate(view['id'])
    check(app.feFileSaveAs(False, str(output)), 'save speed sweep result model')
    check(app.feFileClose(False), 'close saved speed sweep model')
    check(app.feFileOpen(False, str(output)), 'reopen saved speed sweep model')
    if set(result_ids(app)) != {record['output_set_id'] for record in report['result_sets']}:
        raise RuntimeError('Result sets changed after reopening')
    persisted_catalog = {item['id']: item for item in catalog_results(app)}
    for record in report['result_sets']:
        rc, location = app.feResults.SetLocation(record['output_set_id'])
        check(rc, 'read persisted result storage')
        if location != 1:
            raise RuntimeError('Speed sweep results are externally attached, not embedded')
        if persisted_catalog[record['output_set_id']]['vector_count'] != record['imported_vector_count']:
            raise RuntimeError('An imported OP2 result vector was lost')
        saved_set = app.feOutputSet
        check(saved_set.Get(record['output_set_id']), 'read persisted applicability notes')
        if saved_set.notes != record['output_set_notes']:
            raise RuntimeError('Physical applicability notes did not persist')
        sample_preserved(app, record)
        check_display_sample(app, record)
    rc, active = app.feAppGetActiveView()
    check(rc, 'read persisted speed view')
    native_view = app.feView
    check(native_view.Get(active), 'read persisted scalar view selection')
    if (native_view.OutputSet != final['output_set_id'] or native_view.Contour != 1
            or native_view.ContourDataV2 != final['display_vector_id']):
        raise RuntimeError('Default AOA 6 scalar contour did not persist')
    report.update({
        'velocity_m_s': speed, 'linear_solution_only': True,
        'applicability_note': LINEAR_NOTE, 'physical_screen_is_not_an_import_gate': True,
        'runtime_diagnostic_review': runtime_review,
        'stress_unit_conversion': 'Native stress values are Pa; MPa = Pa / 1000000.',
        'shell_fiber_convention': 'Top/Bottom follow the shell normal, not aircraft upper/lower skins.',
        'default_result_kind': default_result,
        'default_display': {'aoa_deg': 6, 'velocity_m_s': speed, **view},
        'all_imported_op2_vectors_retained': True,
        'result_model_sha256': sha256(output),
        'persistence_check': {'reopened_from_disk': True, 'all_four_sets_in_femap_database': True,
                              'all_four_labels_and_strain_samples_retained': True,
                              'all_four_display_samples_retained': True,
                              'all_four_applicability_notes_retained': True,
                              'default_aoa_6_scalar_contour_retained': True},
    })
    output.with_name(output.stem + '_femap_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-root', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'BWB/results/speed_200')
    parser.add_argument('--review-dir', type=Path, default=ROOT / 'BWB/results/speed_200',
                        help='Existing reviewed evidence directory, independent of the new model output directory')
    parser.add_argument('--speed-mps', type=float, default=200.0)
    parser.add_argument('--parts', nargs='+', type=int, choices=(1, 2), default=[1, 2])
    parser.add_argument('--default-result', choices=('top-vonmises-stress', 'bottom-vonmises-stress', 'top-major-strain'),
                        default='bottom-vonmises-stress')
    parser.add_argument('--allow-reviewed-empty-matrix-diagnostic', action='store_true')
    args = parser.parse_args()
    if not math.isfinite(args.speed_mps) or args.speed_mps <= 0:
        parser.error('--speed-mps must be positive and finite')
    run_root, output_dir = args.run_root.resolve(strict=True), args.output_dir.resolve()
    review_dir = args.review_dir.resolve(strict=True)
    if not output_dir.is_relative_to(ROOT / 'BWB/results'):
        parser.error('Use a dedicated BWB/results output directory')
    for part in args.parts:
        verify_flight_labels(run_root, part, args.speed_mps)
        read_review_evidence(review_dir, part, run_root, args.speed_mps)
        check_inputs(ROOT / f'BWB/models/aero_panels/part_{part}/part_{part}_aero_structure.modfem',
                     run_root, part, args.allow_reviewed_empty_matrix_diagnostic)
    app = connect(True)
    previous = (app.DialogAutoAnswer, app.DialogAutoSkip, app.DialogAutoSkipMsg, app.Pref_UseSubcaseIdForOutput)
    try:
        app.DialogAutoAnswer = app.DialogAutoSkip = app.DialogAutoSkipMsg = 1
        app.Pref_UseSubcaseIdForOutput = False
        speed_tag = f'{args.speed_mps:g}'.replace('.', 'p')
        for part in args.parts:
            report = build(app, part, run_root, output_dir / f'part_{part}_{speed_tag}mps_aoa_3to6.modfem',
                           args.speed_mps, args.default_result, args.allow_reviewed_empty_matrix_diagnostic,
                           review_dir=review_dir)
            print(json.dumps({key: value for key, value in report.items() if key != 'result_sets'}, indent=2))
    finally:
        (app.DialogAutoAnswer, app.DialogAutoSkip, app.DialogAutoSkipMsg,
         app.Pref_UseSubcaseIdForOutput) = previous


if __name__ == '__main__':
    main()
