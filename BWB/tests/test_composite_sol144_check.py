"""Independent material/diagnostic fixtures; no saved trial results are required."""
from pathlib import Path
import tempfile
import unittest

from BWB.scripts.check_composite_sol144 import constitutive_residual, coupling_evidence


MATERIALS = {'7': {'E1_Pa': 100e9, 'E2_Pa': 10e9, 'nu12': .25, 'G12_Pa': 5e9}}


def row(stress, strain):
    return dict(material_id='7', sigma_11_Pa=stress[0], sigma_22_Pa=stress[1], tau_12_Pa=stress[2],
                epsilon_11=strain[0], epsilon_22=strain[1], gamma_12=strain[2])


def write_coupling_fixture(folder, density='5.86273E-04', increment='2.500000E-03', include_monitor=True):
    (folder / 'fixture.f04').write_text(
        ' *8** Module DMAP Matrix Cols Rows F T NzWds Density\n'
        f' MPYAD 247 KALL 12 12 1 1 3 {density} 1 1 1 1 1 0 *8**\n')
    text = ('STRUCTURAL MONITOR POINT TOTAL VEHICLE COEFFICIENTS:\n'
            ' WIND CL 9.000000E-01 0.000000E+00 0. 0. 0. 9.000000E-01\n')
    if include_monitor:
        text += ('AERODYNAMIC MONITOR POINT TOTAL VEHICLE COEFFICIENTS:\n'
                 ' BODY CZ 1.500000E-01 9.000000E-04 0. 0. 0. 1.509000E-01\n'
                 f' WIND CL 1.000000E-01 {increment} 0. 0. 0. 1.025000E-01\n')
    (folder / 'fixture.f06').write_text(text)


class CompositeConstitutiveChecks(unittest.TestCase):
    def test_unequal_moduli_obey_poisson_reciprocity_in_both_loading_directions(self):
        # nu21 = .025, not .25. Equal-modulus woven material alone would
        # not detect using E2 in the off-diagonal compliance by mistake.
        rows = [row((100e6, 0., 0.), (.001, -.00025, 0.)),
                row((0., 8e6, 0.), (-.00002, .0008, 0.))]
        result = constitutive_residual(rows, MATERIALS)
        self.assertTrue(result['passed'])
        self.assertLess(result['relative_error'], 1e-12)
        wrong = row((100e6, 0., 0.), (.001, -.0025, 0.))
        self.assertFalse(constitutive_residual([wrong], MATERIALS)['passed'])

    def test_engineering_shear_passes_and_tensor_shear_substitution_fails(self):
        correct = row((0., 0., 4e6), (0., 0., .0008))
        incorrect = row((0., 0., 4e6), (0., 0., .0004))
        self.assertTrue(constitutive_residual([correct], MATERIALS)['passed'])
        self.assertFalse(constitutive_residual([incorrect], MATERIALS)['passed'])

    def test_six_significant_digit_printing_passes_but_stress_unit_error_fails(self):
        # Independent stiffness calculation, then separate F06-like rounding
        # of strains and stresses before checking the compliance relationship.
        e1, e2, nu, g = 100e9, 10e9, .25, 5e9
        exact_strain = (1.23456789e-4, -3.45678912e-5, 8.7654321e-5)
        determinant = 1 - nu * nu * e2 / e1
        exact_stress = ((e1 * exact_strain[0] + nu * e2 * exact_strain[1]) / determinant,
                        (nu * e2 * exact_strain[0] + e2 * exact_strain[1]) / determinant,
                        g * exact_strain[2])
        rounding = lambda values: tuple(float(f'{v:.5E}') for v in values)
        valid = row(rounding(exact_stress), rounding(exact_strain))
        self.assertTrue(constitutive_residual([valid], MATERIALS)['passed'])
        wrong_units = dict(valid, sigma_11_Pa=valid['sigma_11_Pa'] / 1e6)
        self.assertFalse(constitutive_residual([wrong_units], MATERIALS)['passed'])

    def test_nonfinite_or_zero_response_is_not_considered_valid(self):
        for invalid in (row((float('nan'), 0., 0.), (.001, 0., 0.)),
                        row((0., 0., 0.), (0., 0., 0.))):
            with self.subTest(invalid=invalid):
                with self.assertRaises(ValueError):
                    constitutive_residual([invalid], MATERIALS)


class CompositeCouplingChecks(unittest.TestCase):
    def test_selects_aerodynamic_wind_monitor_and_nonzero_matrix_density(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            write_coupling_fixture(folder)
            result = coupling_evidence(folder, 'fixture')
            self.assertTrue(result['nonzero_KALL'])
            self.assertTrue(result['nonzero_elastic_lift_increment'])
            self.assertAlmostEqual(result['rigid_CL'], .1)
            self.assertAlmostEqual(result['elastic_CL_increment'], .0025)
            self.assertAlmostEqual(result['total_CL'], .1025)

    def test_zero_kall_density_cannot_pass_the_coupling_screen(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            write_coupling_fixture(folder, density='0.00000E+00')
            result = coupling_evidence(folder, 'fixture')
            self.assertFalse(result['nonzero_KALL'])

    def test_zero_lift_increment_cannot_pass_the_feedback_screen(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            write_coupling_fixture(folder, increment='0.000000E+00')
            result = coupling_evidence(folder, 'fixture')
            self.assertTrue(result['nonzero_KALL'])
            self.assertFalse(result['nonzero_elastic_lift_increment'])

    def test_missing_aerodynamic_monitor_cannot_use_a_structural_table(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            write_coupling_fixture(folder, include_monitor=False)
            with self.assertRaisesRegex(ValueError, 'Missing aeroelastic monitor'):
                coupling_evidence(folder, 'fixture')


if __name__ == '__main__':
    unittest.main()
