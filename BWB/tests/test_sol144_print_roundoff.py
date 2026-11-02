"""Printed F06 precision must distinguish rounding from inconsistent high loads."""
import unittest

import numpy as np

from BWB.scripts.check_sol144_reasonableness import (
    pressure_consistency, printed_roundoff, reaction_consistency,
)


def printed(values):
    tokens = [f"{value:.6E}" for value in values]
    return np.array([float(token) for token in tokens]), np.array([printed_roundoff(token) for token in tokens])


class PrintedRoundoffTests(unittest.TestCase):
    def test_exponent_and_printed_precision_are_preserved(self):
        self.assertAlmostEqual(printed_roundoff("3.024690E+04"), .005)
        self.assertAlmostEqual(printed_roundoff("-1.234567D-03"), 5e-10)
        self.assertAlmostEqual(printed_roundoff("1.234567+03"), .0005)
        self.assertAlmostEqual(printed_roundoff("49.000"), .0005)

    def test_high_q_rounding_passes_but_changed_pressure_fails(self):
        q, exact_cp = 24500., 1.2345674499
        values, errors = printed([exact_cp, q * exact_cp])
        residual, tolerance = pressure_consistency(*values, q, *errors)
        self.assertGreater(residual, 1e-3)  # The old fixed tolerance rejected this valid pair.
        self.assertLessEqual(residual, tolerance)
        changed_residual, changed_tolerance = pressure_consistency(values[0], values[1] + .1, q, *errors)
        self.assertGreater(changed_residual, changed_tolerance)

    def test_small_cp_does_not_inherit_large_pressure_tolerance(self):
        q, exact_cp = 24500., 1.2345674499e-8
        values, errors = printed([exact_cp, q * exact_cp])
        residual, tolerance = pressure_consistency(*values, q, *errors)
        self.assertLessEqual(residual, tolerance)
        changed_residual, changed_tolerance = pressure_consistency(values[0], values[1] + 1e-6, q, *errors)
        self.assertGreater(changed_residual, changed_tolerance)

    def test_cancelling_nodal_forces_use_operand_rounding_not_small_resultant(self):
        points = np.array([[.4, .02, 1.56], [-.3, -.01, 1.56]])
        exact = np.array([[123456.74499, -45678.12349, 1000.123449, 102.123449, -89.123449, 65.123449],
                          [-123455.25501, 45676.63351, -999.033551, -101.033551, 88.033551, -64.033551]])
        exact_resultant = np.r_[exact[:, :3].sum(axis=0),
                                (exact[:, 3:] + np.cross(points, exact[:, :3])).sum(axis=0)]
        values, budgets = zip(*(printed(row) for row in exact))
        listed, listed_budget = printed(exact_resultant)
        difference, tolerance = reaction_consistency(points, values, budgets, listed, listed_budget)
        self.assertGreater(max(abs(difference)), 1e-4)
        self.assertTrue(np.all(abs(difference) <= tolerance))
        altered = listed.copy()
        altered[0] += 1.
        altered_difference, altered_tolerance = reaction_consistency(points, values, budgets, altered, listed_budget)
        self.assertGreater(abs(altered_difference[0]), altered_tolerance[0])

    def test_moment_uses_printed_force_uncertainty_times_lever_arm(self):
        points = np.array([[0., 0., 100.]])
        values, budgets = printed([123456.74499, 0., 0., 0., 0., 0.])
        exact_total = np.array([123456.74499, 0., 0., 0., 12345674.499, 0.])
        listed, listed_budget = printed(exact_total)
        difference, tolerance = reaction_consistency(points, [values], [budgets], listed, listed_budget)
        self.assertTrue(np.all(abs(difference) <= tolerance))
        altered = listed.copy()
        altered[4] += 100.
        difference, tolerance = reaction_consistency(points, [values], [budgets], altered, listed_budget)
        self.assertGreater(abs(difference[4]), tolerance[4])


if __name__ == "__main__":
    unittest.main()
