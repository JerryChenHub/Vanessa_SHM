"""Regressions for the flight-condition bugs found by actual SOL144 runs."""
from pathlib import Path
import tempfile
import unittest

from BWB.paths import WORKSPACE_DIR
from BWB.scripts.build_sol144_demo import build
from BWB.scripts.inspect_shell_bdf import cards, nastran_float


class Sol144CaseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.scratch_root = WORKSPACE_DIR.resolve()
        cls.scratch_root.mkdir(parents=True, exist_ok=True)
        cls.temporary = tempfile.TemporaryDirectory(prefix='test_sol144_case_', dir=cls.scratch_root)
        cls.fixture = Path(cls.temporary.name).resolve()
        cls.addClassCleanup(cls.cleanup_fixture)
        if cls.fixture == cls.scratch_root or not cls.fixture.is_relative_to(cls.scratch_root):
            raise RuntimeError('Test fixture must stay strictly inside workspace/BWB')
        # Use the retained model sources; no saved trial inputs or results are needed.
        build(cls.fixture, velocity=20., density=1.225, sound_speed=340.3, alpha_deg=2.)

    @classmethod
    def cleanup_fixture(cls):
        target = Path(cls.temporary.name).resolve()
        if target == cls.scratch_root or not target.is_relative_to(cls.scratch_root):
            raise RuntimeError('Refusing to clean a test fixture outside workspace/BWB')
        cls.temporary.cleanup()

    def test_explicit_zero_sideslip_preserves_elastic_flag(self):
        for part in (1, 2):
            path = self.fixture / f'part_{part}/flight_condition.bdf'
            parsed = list(cards(path))
            trim = next(fields for name, fields, _ in parsed if name == 'TRIM')
            self.assertEqual(trim[4], 'ANGLEA')
            self.assertEqual(trim[6], 'SIDES')
            self.assertEqual(nastran_float(trim[7]), 0.)
            self.assertEqual(nastran_float(trim[8]), 1.)
            self.assertAlmostEqual(nastran_float(trim[2]), 20 / 340.3, places=10)
            self.assertAlmostEqual(nastran_float(trim[3]), .5 * 1.225 * 20**2, places=10)
            self.assertEqual({f[2] for n, f, _ in parsed if n == 'AESTAT'}, {'ANGLEA', 'SIDES'})

    def test_si_case_does_not_inherit_nasa_weight_conversion(self):
        parsed = list(cards(self.fixture / 'part_1/bwb_sol144.bdf'))
        names = [name for name, _, _ in parsed]
        self.assertEqual(names.count('GRID'), 2940)
        self.assertEqual(names.count('CQUAD4') + names.count('CTRIA3'), 2915)
        self.assertFalse(any(name == 'PARAM' and fields[1] in ('WTMASS', 'AUNITS')
                             for name, fields, _ in parsed))
        self.assertEqual(names.count('DMI'), 2)  # header plus one continued column


if __name__ == '__main__':
    unittest.main()
