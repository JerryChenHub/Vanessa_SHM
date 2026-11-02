import tempfile
from pathlib import Path
import unittest
from core.plate_strain import read_fiber_strains
from core.batch_run import cleanup_case
from core.paths import WORKSPACE_DIR

FIBER = '''0 SUBCASE 1
 S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )
 ID. DISTANCE NORMAL-X NORMAL-Y SHEAR-XY ANGLE MAJOR MINOR VON MISES
0 60 -0.1 0.001 -0.002 0.003 0.0 0.0 0.0 0.0
 0.1 0.005 0.006 -0.007 0.0 0.0 0.0 0.0
 * * * END OF JOB * * *
'''


class PlateStrainTest(unittest.TestCase):
    def test_surface_order_and_componentwise_mean(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)/'test.f06'
            p.write_text(FIBER)
            row=read_fiber_strains(p)[(1,60)]
            self.assertEqual(row['bottom'],(.001,-.002,.003))
            self.assertEqual(row['top'],(.005,.006,-.007))
            for a,b in zip(row['mid'],(.003,.002,-.002)):
                self.assertAlmostEqual(a,b)
            # Reverse output fiber order; top still follows positive element-z.
            lines=FIBER.splitlines(True)
            lines[3]='0 60 0.1 0.005 0.006 -0.007 0.0 0.0 0.0 0.0\n'
            lines[4]=' -0.1 0.001 -0.002 0.003 0.0 0.0 0.0 0.0\n'
            p.write_text(''.join(lines))
            self.assertEqual(read_fiber_strains(p)[(1,60)],row)

    def test_incomplete_nonfinite_and_curvature_rows_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)/'test.f06'
            for content in (FIBER.replace('* * * END OF JOB * * *',''),FIBER.replace('0.005','nan'),
                            FIBER.replace('DISTANCE','CURVATURE'),FIBER.replace(' 0.1 0.005 0.006 -0.007 0.0 0.0 0.0 0.0\n','')):
                p.write_text(content)
                with self.assertRaises(ValueError): read_fiber_strains(p)

    def test_cleanup_rejects_shared_dataset_and_workspace_root(self):
        work=WORKSPACE_DIR/'batch_test'
        for path in (work,work/'base',WORKSPACE_DIR.parent/'data_test/case_001'):
            with self.assertRaises(ValueError): cleanup_case(path,work)


if __name__=='__main__':
    unittest.main()
