"""Independent constitutive, thickness and structural-coverage checks."""
from __future__ import annotations

from collections import Counter
import copy
import json
import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

from BWB.paths import REPOSITORY_ROOT as ROOT, WORKSPACE_DIR
from BWB.scripts.build_composite_skin import (
    DEFAULT_CONFIG, REGIONS, build, element_assignment, laminate, material_values,
    report_band_boundaries, transformed_stiffness,
)
from BWB.scripts.inspect_shell_bdf import cards, load_model, nastran_float


def independent_rotated_stiffness(material, angle):
    """Rotate physical strain/stress tensors, not closed-form Qbar entries."""
    e1, e2, nu, shear = (material[key] for key in ('E1_Pa', 'E2_Pa', 'nu12', 'G12_Pa'))
    compliance = np.array([[1/e1, -nu/e1, 0.], [-nu/e1, 1/e2, 0.], [0., 0., 1/shear]])
    local_q = np.linalg.inv(compliance)
    c, s = math.cos(math.radians(angle)), math.sin(math.radians(angle))
    rotation = np.array([[c, -s], [s, c]])
    columns = []
    for unit in np.eye(3):
        global_strain = np.array([[unit[0], unit[2]/2], [unit[2]/2, unit[1]]])
        strain = rotation.T @ global_strain @ rotation
        stress = local_q @ np.array([strain[0, 0], strain[1, 1], 2*strain[0, 1]])
        tensor = rotation @ np.array([[stress[0], stress[2]], [stress[2], stress[1]]]) @ rotation.T
        columns.append([tensor[0, 0], tensor[1, 1], tensor[0, 1]])
    return np.array(columns).T


class LaminateTests(unittest.TestCase):
    def setUp(self):
        self.config = json.loads(DEFAULT_CONFIG.read_text(encoding='utf-8-sig'))
        self.material = material_values(self.config)

    def test_rotated_stiffness_matches_independent_tensor_operation(self):
        material = dict(self.material, E1_Pa=131e9, E2_Pa=11e9, nu12=.27, G12_Pa=4.7e9)
        for angle in (0., 17., -31., 45., 90., 137.):
            with self.subTest(angle=angle):
                actual = transformed_stiffness(material, angle)
                np.testing.assert_allclose(actual, independent_rotated_stiffness(material, angle), rtol=2e-14, atol=1e-4)
                self.assertGreater(np.linalg.eigvalsh(actual).min(), 0.)

    def test_uniform_plate_abd_matches_closed_form_and_mass(self):
        stack = laminate(dict(ply_count=4, angle_deg=23., ply_thickness_m=.0002), self.material, 1)
        thickness = .0008
        q = independent_rotated_stiffness(self.material, 23.)
        np.testing.assert_allclose(stack['A_N_per_m'], q*thickness, rtol=1e-14)
        np.testing.assert_allclose(stack['B_N'], 0., atol=1e-11)
        np.testing.assert_allclose(stack['D_Nm'], q*thickness**3/12, rtol=1e-14)
        self.assertEqual(stack['ply_thicknesses_m'], [.0002]*4)
        self.assertEqual(stack['ply_thickness_m'], .0002)
        self.assertAlmostEqual(stack['areal_mass_kg_m2'], 1350*thickness)

    def test_nonuniform_abd_scales_with_thickness_and_reversing_stack(self):
        region = dict(ply_angles_deg=[0., 32., -18.], ply_thicknesses_m=[.0001, .0003, .0002])
        first = laminate(region, self.material, 1)
        factor = 2.5
        scaled = laminate(dict(region, ply_thicknesses_m=[factor*t for t in region['ply_thicknesses_m']]), self.material, 1)
        self.assertGreater(np.linalg.norm(first['B_N']), 1.)
        for name, exponent in (('A_N_per_m', 1), ('B_N', 2), ('D_Nm', 3)):
            np.testing.assert_allclose(scaled[name], np.array(first[name])*factor**exponent, rtol=1e-13, atol=1e-10)
        reverse = laminate(dict(ply_angles_deg=region['ply_angles_deg'][::-1],
                                ply_thicknesses_m=region['ply_thicknesses_m'][::-1]), self.material, 1)
        np.testing.assert_allclose(reverse['A_N_per_m'], first['A_N_per_m'], rtol=1e-14)
        np.testing.assert_allclose(reverse['B_N'], -np.array(first['B_N']), rtol=1e-13, atol=1e-11)
        np.testing.assert_allclose(reverse['D_Nm'], first['D_Nm'], rtol=1e-13, atol=1e-12)
        self.assertNotIn('ply_thickness_m', first)

    def test_symmetry_accounts_for_thickness_and_angles(self):
        region = dict(ply_angles_deg=[0., 45., 45., 0.], ply_thicknesses_m=[.0001, .0002, .0002, .0001])
        symmetric = laminate(region, self.material, 1)
        self.assertTrue(symmetric['symmetric'])
        np.testing.assert_allclose(symmetric['B_N'], 0., atol=1e-11)
        asymmetric = laminate(dict(region, ply_thicknesses_m=[.0001, .0002, .0002, .0003]), self.material, 1)
        self.assertFalse(asymmetric['symmetric'])
        self.assertGreater(np.linalg.norm(asymmetric['B_N']), 1.)
        self.assertFalse(laminate(dict(region, ply_angles_deg=[0., 45., 20., 0.]), self.material, 1)['symmetric'])

    def test_counts_totals_redundant_consistent_inputs_and_mirror_sign(self):
        region = dict(ply_count=3, total_thickness_m=.0006, angle_deg=13.)
        stack = laminate(region, self.material, -1)
        self.assertEqual(stack['angles_deg'], [-13.]*3)
        np.testing.assert_allclose(stack['ply_thicknesses_m'], [.0002]*3)
        explicit = dict(region, ply_angles_deg=[0., 2., 0.], ply_thickness_m=.0002,
                        ply_thicknesses_m=[.0002]*3)
        self.assertEqual(laminate(explicit, self.material, -1)['angles_deg'], [-13., -15., -13.])
        inferred = laminate(dict(ply_thicknesses_m=[.0001, .0002], angle_deg=5.), self.material, 1)
        self.assertEqual(inferred['angles_deg'], [5., 5.])
        self.assertNotIn('ply_thickness_m', inferred)

    def test_conflicts_nonfinite_and_invalid_counts_are_rejected(self):
        base = dict(ply_count=2, ply_thickness_m=.0002)
        invalid = [dict(base, total_thickness_m=.0006),
                   dict(base, ply_thicknesses_m=[.0002, .0003]),
                   dict(ply_count=2, ply_thicknesses_m=[.0001, .0002], total_thickness_m=.0004),
                   dict(base, ply_angles_deg=[0.]), dict(base, ply_angles_deg=[]),
                   dict(base, ply_thicknesses_m=[.0002]), dict(base, ply_thicknesses_m=[]),
                   dict(base, ply_angles_deg=[0., float('nan')]),
                   dict(base, angle_deg=float('inf')), dict(ply_count=2),
                   dict(ply_thickness_m=.0002), dict(base, ply_thickness_m=True)]
        invalid += [dict(base, ply_count=count) for count in (0, -1, 2.5, True, float('nan'))]
        invalid += [dict(base, ply_thickness_m=value) for value in (0, -.01, float('inf'), float('nan'))]
        for region in invalid:
            with self.subTest(region=region), self.assertRaises(ValueError):
                laminate(region, self.material, 1)
        with self.assertRaises(ValueError):
            laminate(base, self.material, 0)

    def test_nonpositive_orthotropic_compliance_is_rejected(self):
        config = copy.deepcopy(self.config)
        config['material'].update(E1_Pa=1e9, E2_Pa=100e9, nu12=.5)
        with self.assertRaisesRegex(ValueError, 'positive definite'):
            material_values(config)


class CompositeGeometryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = json.loads(DEFAULT_CONFIG.read_text(encoding='utf-8-sig'))
        cls.models, cls.assignments = {}, {}
        for part in (1, 2):
            model = load_model(ROOT/f'BWB/models/aero_panels/part_{part}/structural_cards.bdf')
            layout = json.loads((ROOT/f'BWB/models/skin_structured/part_{part}/mesh_layout.json').read_text())
            cls.models[part] = model
            cls.assignments[part] = element_assignment(model, layout, cls.config)

    def test_all_shells_have_unique_region_and_tangent_material_frame(self):
        expected = {'leading_edge': 238, 'lower': 1062, 'trailing_edge': 535, 'upper': 1080}
        for part, assignment in self.assignments.items():
            self.assertEqual(set(assignment), set(self.models[part]['elements']))
            self.assertEqual(len(assignment), 2915)
            self.assertEqual(Counter(item['region'] for item in assignment.values()), expected)
            self.assertEqual(len({item['mcid'] for item in assignment.values()}), len(assignment))
            for item in assignment.values():
                self.assertEqual(item['pid'], REGIONS[item['region']])
                self.assertGreater(item['area_m2'], 0.)
                self.assertAlmostEqual(np.linalg.norm(item['normal']), 1.)
                self.assertAlmostEqual(np.linalg.norm(item['reference_axis']), 1.)
                self.assertAlmostEqual(np.dot(item['normal'], item['reference_axis']), 0.)
                self.assertEqual(item['reference_kind'], 'tip_chord' if item['patch']=='tip_cap' else 'outward_span')

    def test_mirror_geometry_reference_axes_and_region_labels(self):
        for eid, first in self.assignments[1].items():
            second = self.assignments[2][eid]
            expected = np.array(first['centroid_m'])*[1, 1, -1] + [0, 0, 1.56]
            np.testing.assert_allclose(second['centroid_m'], expected, atol=1e-6, rtol=0)
            np.testing.assert_allclose(second['reference_axis'], np.array(first['reference_axis'])*[1, 1, -1], atol=1e-4, rtol=0)
            self.assertEqual(first['region'], second['region'])

    def test_band_report_uses_outer_boundary_not_internal_leading_edge_columns(self):
        assignment = self.assignments[1]
        report = report_band_boundaries(assignment, self.config)
        for patch, record in report.items():
            for region, expected_fraction in (('leading_edge', 2/24), ('trailing_edge', 4/24)):
                band = record[region]
                self.assertFalse(band['missing_span_rows'], patch)
                self.assertEqual(band['represented_span_rows'], record['source_span_rows'])
                self.assertAlmostEqual(band['achieved_fraction_min'], expected_fraction, delta=3e-4)
                self.assertAlmostEqual(band['achieved_fraction_max'], expected_fraction, delta=3e-4)
                self.assertAlmostEqual(band['fraction_error_min'], band['achieved_fraction_min']-band['target_fraction'])
                self.assertAlmostEqual(band['fraction_error_max'], band['achieved_fraction_max']-band['target_fraction'])
        # Every skin node uses its own source row's endpoints, including the taper.
        for item in assignment.values():
            if 'chord_node_fractions' in item:
                self.assertTrue(all(-1e-12 <= xi <= 1+1e-12 for xi in item['chord_node_fractions']))

    def test_pcomp_records_preserve_nonuniform_thicknesses(self):
        WORKSPACE_DIR.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='composite-test-', dir=WORKSPACE_DIR) as temporary:
            folder = Path(temporary).resolve()
            self.assertTrue(folder.is_relative_to(WORKSPACE_DIR.resolve()))
            config = copy.deepcopy(self.config)
            config['regions']['upper'] = dict(ply_count=4, angle_deg=17.,
                                             ply_thicknesses_m=[.0001, .0002, .0003, .0004])
            config_path = folder/'config.json'
            config_path.write_text(json.dumps(config), encoding='utf-8')
            output = build(config_path, folder/'generated')
            for part, angle in ((1, 17.), (2, -17.)):
                pcomp = next(fields for name, fields, _ in cards(output/f'part_{part}/structural_cards.bdf')
                             if name == 'PCOMP' and int(fields[1]) == REGIONS['upper'])
                self.assertAlmostEqual(nastran_float(pcomp[2]), -.0005)
                np.testing.assert_allclose([nastran_float(pcomp[10+4*i]) for i in range(4)], [.0001, .0002, .0003, .0004])
                self.assertEqual([nastran_float(pcomp[11+4*i]) for i in range(4)], [angle]*4)
                manifest = json.loads((output/f'part_{part}/composite_assignment.json').read_text())
                self.assertNotIn('ply_thickness_m', manifest['laminates']['upper'])
                self.assertIn('ply_thickness_m', manifest['laminates']['lower'])
                self.assertIn('band_boundaries', manifest)
            # Invalid input must be rejected before creating any trial deck.
            config['regions']['upper']['total_thickness_m'] = .002
            config_path.write_text(json.dumps(config), encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'conflicts'):
                build(config_path, folder/'rejected')
            self.assertFalse((folder/'rejected').exists())


if __name__ == '__main__':
    unittest.main()
