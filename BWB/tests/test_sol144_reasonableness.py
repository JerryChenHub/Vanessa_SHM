"""Small numerical guards for independent SOL144 coarse checks."""
import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

from BWB.scripts.check_sol144_reasonableness import (
    plane_stress_von_mises, read_load_tables, shift_moment,
)


class ReasonablenessTests(unittest.TestCase):
    def test_plane_stress_uniaxial_and_engineering_shear(self):
        young, poisson, stress = 69e9, .33, 12e6
        self.assertAlmostEqual(plane_stress_von_mises(stress / young, -poisson * stress / young, 0), stress)
        gamma = .002
        expected = math.sqrt(3) * young / (2 * (1 + poisson)) * gamma
        self.assertAlmostEqual(plane_stress_von_mises(0, 0, gamma), expected)

    def test_part2_root_origin_moment_shift(self):
        force = np.array([0., 7., 0.])
        point, root = np.array([0., 0., 2.56]), np.array([0., 0., 1.56])
        moment_zero = np.cross(point, force)
        np.testing.assert_allclose(shift_moment(moment_zero, force, root), [-7., 0., 0.])

    def test_labeled_apres_and_spc_tables(self):
        text = """0 ASSUMED_LOW_SPEED_ELASTIC SUBCASE 1
AERODYNAMIC PRESSURES ON THE AERODYNAMIC ELEMENTS
1 LS .2 49.0 1000001
0 SPCFORCE RESULTANT
TOTALS 0.0 -7.0 0.0 7.0 0.0 0.0
F O R C E S O F S I N G L E - P O I N T C O N S T R A I N T
1 G 0.0 -7.0 0.0 7.0 0.0 0.0
* * * END OF JOB * * *
"""
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "fixture.f06"
            path.write_text(text)
            pressures, resultants, reactions = read_load_tables(path)
            self.assertEqual(pressures[(1, 1000001)], (.2, 49.))
            np.testing.assert_allclose(resultants[1], [0, -7, 0, 7, 0, 0])
            np.testing.assert_allclose(reactions[(1, 1)], resultants[1])


if __name__ == "__main__":
    unittest.main()
