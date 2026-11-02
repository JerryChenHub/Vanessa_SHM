import tempfile
from pathlib import Path
import unittest
import json

from core.export_results import read_tables, export_case
from core.run_local_nastran import inspect_results


F06 = '''0                          SUBCASE 1
 D I S P L A C E M E N T   V E C T O R
 42 G 1.0 2.0 3.0 0.0 0.0 0.0
 S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )
 ID. CURVATURE NORMAL-X NORMAL-Y SHEAR-XY ANGLE MAJOR MINOR VON MISES
0 60 0.0 0.001 -0.0002 0.0006 12.0 0.0011 -0.0003 0.0009
 -1.0 0.1 0.2 0.3 12.0 0.3 0.1 0.2
 * * * END OF JOB * * *
'''


class ResultTablesTest(unittest.TestCase):
    def test_export_rejects_missing_requested_subcase_and_stale_input(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            (path/'sol144.f06').write_text(F06)
            (path/'sol144.log').write_text('Analysis complete  0\n')
            (path/'subcases.inc').write_text('SUBCASE 1\nSUBCASE 2\n')
            with self.assertRaisesRegex(ValueError,'subcases'):
                export_case(path)
            (path/'run-summary.json').write_text(json.dumps({'fresh_f06':False}))
            with self.assertRaisesRegex(ValueError,'freshness'):
                export_case(path)
            self.assertFalse((path/'results').exists())

    def test_membrane_row_is_not_averaged_with_curvature(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'sol144.f06'
            path.write_text(F06)
            displacement, strain = read_tables(path)
            self.assertEqual(displacement[(1,42)], (1.,2.,3.,0.,0.,0.))
            self.assertEqual(strain[(1,60)], ('CQUAD4',0.001,-0.0002,0.0006))

    def test_incomplete_and_unknown_fiber_output_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'sol144.f06'
            for text in (F06.replace('* * * END OF JOB * * *',''), F06.replace('CURVATURE','FIBER DISTANCE')):
                path.write_text(text)
                with self.assertRaises(ValueError):
                    read_tables(path)

    def test_zero_exit_message_does_not_hide_fatal_or_mkl_error(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            (path/'sol144.f06').write_text(F06)
            (path/'sol144.log').write_text('Analysis complete  0\nIntel MKL ERROR: Parameter 13 was incorrect on entry to DGEMM .\n')
            self.assertEqual(inspect_results(path)['status'], 'completed_with_runtime_errors')
            (path/'sol144.f06').write_text(F06+'*** USER FATAL MESSAGE 315\n')
            self.assertFalse(inspect_results(path)['solver_completed'])


if __name__ == '__main__':
    unittest.main()
