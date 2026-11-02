"""Regressions for the flight-condition bugs found by actual SOL144 runs."""
from pathlib import Path
import unittest

from BWB.scripts.inspect_shell_bdf import cards, nastran_float


ROOT = Path(__file__).resolve().parents[2]


class Sol144CaseTests(unittest.TestCase):
    def test_explicit_zero_sideslip_preserves_elastic_flag(self):
        for part in (1, 2):
            path = ROOT / f'BWB/cases/sol144_demo/part_{part}/flight_condition.bdf'
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
        parsed = list(cards(ROOT / 'BWB/cases/sol144_demo/part_1/bwb_sol144.bdf'))
        names = [name for name, _, _ in parsed]
        self.assertEqual(names.count('GRID'), 2940)
        self.assertEqual(names.count('CQUAD4') + names.count('CTRIA3'), 2915)
        self.assertFalse(any(name == 'PARAM' and fields[1] in ('WTMASS', 'AUNITS')
                             for name, fields, _ in parsed))
        self.assertEqual(names.count('DMI'), 2)  # header plus one continued column


if __name__ == '__main__':
    unittest.main()
