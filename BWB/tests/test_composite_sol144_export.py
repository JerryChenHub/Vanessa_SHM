"""Synthetic PCOMP result fixtures; no trial data or Nastran solver required."""
import csv
import json
from pathlib import Path
import tempfile
import unittest

from BWB.scripts.export_composite_sol144 import export_case, read_composite_model, validate_case


STRUCTURE = """GRID,1,,0.,0.,0.
GRID,2,,1.,0.,0.
GRID,3,,1.,0.,1.
GRID,4,,0.,0.,1.
GRID,5,,.5,0.,2.
CQUAD4,10,11,1,2,3,4,30.
CTRIA3,20,12,4,3,5,30.
MAT8,101,5.4E10,5.2E10,.05,4.E9,3.E9,3.E9,1550.
PCOMP,11,,,,,,,
+,101,.0002,0.,YES,101,.0002,45.,YES
PCOMP,12,,,,,,,
+,101,.0003,-45.,YES
SPC1,1,123456,1,2
"""


def f06_text():
    text = """0 COMPOSITE_TEST SUBCASE 1
 D I S P L A C E M E N T   V E C T O R
 1 G 0.0 0.0 0.0 0.0 0.0 0.0
 2 G 0.0 0.0 0.0 0.0 0.0 0.0
 3 G 0.0 .001 0.0 0.0 0.0 0.0
 4 G 0.0 .001 0.0 0.0 0.0 0.0
 5 G 0.0 .002 0.0 0.0 0.0 0.0
"""
    for quantity in ('STRESSES', 'STRAINS'):
        for kind, pairs in (('QUAD4', ((10, 1), (10, 2))), ('TRIA3', ((20, 1),))):
            text += f' {quantity} IN LAYERED COMPOSITE ELEMENTS ( {kind} )\n'
            text += ' ELEMENT PLY NORMAL-1 NORMAL-2 SHEAR-12 SHEAR-XZ SHEAR-YZ ANGLE MAJOR MINOR MAX SHEAR\n'
            for eid, ply in pairs:
                # Deliberately independent printed stress and strain: the exporter
                # must preserve solver stress, never synthesize it using MAT1.
                values = '7000000. -2000000. 300000. 0. 0. 0. 7010000. -2010000. 4510000.' if quantity == 'STRESSES' else '0. 0. .002 0. 0. 45. .001 -.001 .002'
                text += f'0 {eid} {ply} {values}\n'
    return text + '* * * END OF JOB * * *\n'


def write_case(folder, text=None, structure=STRUCTURE):
    folder.mkdir(parents=True, exist_ok=True)
    (folder / 'structural_cards.bdf').write_text(structure)
    (folder / 'bwb_sol144.bdf').write_text("SOL 144\nCEND\nSUBCASE 1\nBEGIN BULK\nINCLUDE 'structural_cards.bdf'\nENDDATA\n")
    (folder / 'bwb_sol144.f06').write_text(f06_text() if text is None else text)


class CompositeExportTests(unittest.TestCase):
    def test_ply_midpoints_coordinates_and_direct_stress_export(self):
        with tempfile.TemporaryDirectory() as temporary:
            case, output = Path(temporary) / 'case', Path(temporary) / 'output'
            write_case(case)
            summary = export_case(case, output)
            self.assertEqual(summary['ply_rows'], 3)
            self.assertEqual(summary['shell_elements'], 2)
            self.assertEqual(summary['root_fixed_nodes'], 2)
            self.assertEqual(summary['max_abs_root_dof'], 0.)
            with (output / 'ply_stress_strain.csv').open(newline='') as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(float(rows[0]['sigma_11_Pa']), 7e6)
            self.assertEqual(float(rows[0]['principal_max_microstrain']), 1000.)
            self.assertEqual(float(rows[0]['principal_min_microstrain']), -1000.)
            self.assertEqual(rows[0]['element_theta_or_mcid'], '30.')
            self.assertAlmostEqual(float(rows[0]['ply_mid_z_m']), -.0001)
            self.assertAlmostEqual(float(rows[1]['ply_mid_z_m']), .0001)
            self.assertEqual(float(rows[1]['ply_angle_deg']), 45.)

    def test_missing_duplicate_nonfinite_and_incomplete_results_rejected(self):
        original = f06_text()
        stress_row = next(line for line in original.splitlines(True) if line.startswith('0 10 2 7000000.'))
        variants = (original.replace(stress_row, ''),
                    original.replace(stress_row, stress_row + stress_row),
                    original.replace('7000000.', 'NaN', 1),
                    original.replace('* * * END OF JOB * * *', ''))
        with tempfile.TemporaryDirectory() as temporary:
            case, output = Path(temporary) / 'case', Path(temporary) / 'output'
            for text in variants:
                with self.subTest(text=text[:60]):
                    write_case(case, text)
                    with self.assertRaises(ValueError):
                        export_case(case, output)
                    self.assertFalse(output.exists())

    def test_moving_root_and_all_zero_strain_are_rejected(self):
        variants = [f06_text().replace('1 G 0.0 0.0', '1 G 0.0 .001'),
                    f06_text().replace('0. 0. .002 0. 0. 45. .001 -.001 .002', '0. 0. 0. 0. 0. 0. 0. 0. 0.')]
        with tempfile.TemporaryDirectory() as temporary:
            case = Path(temporary)
            for text in variants:
                write_case(case, text)
                with self.assertRaises(ValueError):
                    validate_case(case)

    def test_ply_output_must_be_requested_for_full_explicit_stack(self):
        with tempfile.TemporaryDirectory() as temporary:
            case = Path(temporary)
            for structure in (STRUCTURE.replace('0.,YES', '0.,NO'),
                              STRUCTURE.replace('PCOMP,11,,,,,,,', 'PCOMP,11,,,,,,,SYM')):
                write_case(case, structure=structure)
                with self.assertRaises(ValueError):
                    read_composite_model(case / 'bwb_sol144.bdf')

    def test_runtime_diagnostics_are_retained(self):
        with tempfile.TemporaryDirectory() as temporary:
            case, output = Path(temporary) / 'case', Path(temporary) / 'output'
            write_case(case)
            (case / 'bwb_sol144.log').write_text('Intel MKL ERROR: Parameter 13 was incorrect on entry to DGEMM .\n')
            (case / 'run-summary.json').write_text(json.dumps({'launcher_returncode': 1, 'runtime_errors': ['Launcher returned exit code 1']}))
            summary = export_case(case, output)
            self.assertTrue(summary['diagnostic_only'])
            self.assertEqual(summary['recorded_launcher_returncode'], 1)
            self.assertEqual(len(summary['diagnostics']['runtime_errors']), 2)

    def test_real_f06_output_block_trailer_is_not_a_ply_row(self):
        # Simcenter inserts this block between the last stress page and strains.
        trailer = """ *** USER INFORMATION MESSAGE 4114 (OUTPBN2)
 DATA BLOCK OES1C WRITTEN ON FORTRAN UNIT 12, TRL =
 101 2 1 0 0 0 0
 (TOTAL DATA WRITTEN FOR DATA BLOCK = 145601 WORDS.)
0 COMPOSITE_TEST SUBCASE 1
"""
        text = f06_text().replace(' STRAINS IN LAYERED COMPOSITE ELEMENTS ( QUAD4 )',
                                  trailer + ' STRAINS IN LAYERED COMPOSITE ELEMENTS ( QUAD4 )', 1)
        # Real composite rows have no separate carriage-control zero token.
        text = text.replace('\n0 10 ', '\n 10 ').replace('\n0 20 ', '\n 20 ')
        with tempfile.TemporaryDirectory() as temporary:
            case = Path(temporary)
            write_case(case, text)
            result = validate_case(case)
            self.assertEqual(len(result['results']['stress']), 3)
            self.assertEqual(len(result['results']['strain']), 3)


if __name__ == '__main__':
    unittest.main()
