"""Synthetic parser fixtures only; these tests do not represent a solved wing."""
from __future__ import annotations

import csv
import json
from pathlib import Path
import tempfile
import unittest

from BWB.scripts.export_sol144_strain import export_case, principal_strains, validate_case


STRUCTURE = """GRID,1,,0.,0.,0.
GRID,2,,1.,0.,0.
GRID,3,,1.,.005,1.
GRID,4,,0.,.005,1.
GRID,5,,.5,.01,2.
MAT1,1,6.9E10,,.33,2700.
PSHELL,1,1,.002
CQUAD4,10,1,1,2,3,4
CTRIA3,11,1,4,3,5
SPC1,1,123456,1,2
"""

F06 = """0 SUBCASE 1
 D I S P L A C E M E N T   V E C T O R
 POINT ID. TYPE T1 T2 T3 R1 R2 R3
 1 G 0.0 0.0 0.0 0.0 0.0 0.0
 2 G 0.0 0.0 0.0 0.0 0.0 0.0
 3 G 0.0 0.001 0.0 0.0 0.0 0.0
 4 G 0.0 0.001 0.0 0.0 0.0 0.0
 5 G 0.0 0.002 0.0 0.0 0.0 0.0
 S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )
 ID. DISTANCE NORMAL-X NORMAL-Y SHEAR-XY ANGLE MAJOR MINOR VON MISES
0 10 -0.001 0.0 0.001 0.0 0.0 0.0 0.0 0.0
 0.001 0.001 0.0 0.0 0.0 0.0 0.0 0.0
 S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )
 ID. DISTANCE NORMAL-X NORMAL-Y SHEAR-XY ANGLE MAJOR MINOR VON MISES
0 11 -0.001 0.0 0.0 -0.002 0.0 0.0 0.0 0.0
 0.001 0.0 0.0 0.002 0.0 0.0 0.0 0.0
 * * * END OF JOB * * *
"""


def write_case(folder, f06=F06):
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "structural_cards.bdf").write_text(STRUCTURE, encoding="ascii")
    (folder / "bwb_sol144.bdf").write_text("SOL 144\nCEND\nSUBCASE 1\nBEGIN BULK\nINCLUDE 'structural_cards.bdf'\nENDDATA\n", encoding="ascii")
    (folder / "bwb_sol144.f06").write_text(f06, encoding="ascii")


class Sol144StrainTests(unittest.TestCase):
    def test_engineering_shear_uses_half_shear_in_tensor(self):
        first, second, angle = principal_strains(0.0, 0.0, 0.002)
        self.assertAlmostEqual(first, .001)
        self.assertAlmostEqual(second, -.001)
        self.assertAlmostEqual(angle, 45.)

    def test_long_export_mid_components_and_provenance(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary) / "case"
            write_case(folder)
            output = Path(temporary) / "result"
            summary = export_case(folder, output, make_plot=False)
            self.assertEqual(summary["strain_csv_rows"], 6)
            self.assertEqual(summary["displacement_csv_rows"], 5)
            self.assertEqual(summary["root_fixed_nodes"], 2)
            self.assertEqual(summary["max_abs_root_dof"], 0.)
            self.assertEqual(summary["nasa_reuse"]["symbol"], "read_fiber_strains")
            with (output / "shell_strain_long.csv").open(newline="") as stream:
                rows = {(int(r["element_id"]), r["surface"]): r for r in csv.DictReader(stream)}
            # Principal strains must be recomputed after the component average.
            self.assertAlmostEqual(float(rows[(10, "mid")]["principal_max_microstrain"]), 500.)
            self.assertAlmostEqual(float(rows[(10, "mid")]["principal_min_microstrain"]), 500.)
            self.assertAlmostEqual(float(rows[(11, "top")]["principal_max_microstrain"]), 1000.)
            self.assertEqual(float(rows[(11, "mid")]["dominant_principal_microstrain"]), 0.)

    def test_fatal_incomplete_and_nonfinite_output_are_rejected_before_export(self):
        variants = (
            F06 + "*** USER FATAL MESSAGE 123\n",
            F06.replace("* * * END OF JOB * * *", ""),
            F06.replace(" 5 G 0.0 0.002", " 5 G 0.0 nan"),
            F06.replace("0.001 0.001 0.0 0.0", "0.001 nan 0.0 0.0"),
        )
        with tempfile.TemporaryDirectory() as temporary:
            folder, output = Path(temporary) / "case", Path(temporary) / "output"
            for content in variants:
                with self.subTest(content=content[-80:]):
                    write_case(folder, content)
                    with self.assertRaises(ValueError):
                        export_case(folder, output, make_plot=False)
                    self.assertFalse(output.exists())

    def test_missing_shell_and_fixed_root_motion_are_rejected(self):
        missing_shell = F06.replace("0 11 -0.001 0.0 0.0 -0.002 0.0 0.0 0.0 0.0\n", "").replace(
            " 0.001 0.0 0.0 0.002 0.0 0.0 0.0 0.0\n", "")
        moving_root = F06.replace(" 1 G 0.0 0.0", " 1 G 0.0 0.001")
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            write_case(folder, missing_shell)
            with self.assertRaisesRegex(ValueError, "Missing/extra shell"):
                validate_case(folder)
            write_case(folder, moving_root)
            with self.assertRaisesRegex(ValueError, "Fixed-root displacement"):
                validate_case(folder)

    def test_stale_input_hash_or_fiber_thickness_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            write_case(folder, F06.replace("-0.001", "-0.002"))
            with self.assertRaisesRegex(ValueError, "Fiber distance differs"):
                validate_case(folder)
            write_case(folder)
            (folder / "run-summary.json").write_text(json.dumps({"input_sha256": {"structural_cards.bdf": "0" * 64}}))
            with self.assertRaisesRegex(ValueError, "Solved input hashes differ"):
                validate_case(folder)

    def test_labeled_subcase_header_is_adapted_without_modifying_source(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder, output = Path(temporary) / "case", Path(temporary) / "output"
            write_case(folder, F06.replace("0 SUBCASE 1", "0     ASSUMED_LOW_SPEED_ELASTIC                SUBCASE 1"))
            before = (folder / "bwb_sol144.f06").read_bytes()
            (folder / "run-summary.json").write_text(json.dumps({"runtime_errors": ["Diagnostic fixture error"]}))
            summary = export_case(folder, output, make_plot=False)
            self.assertEqual(before, (folder / "bwb_sol144.f06").read_bytes())
            self.assertEqual(summary["f06_parser_adaptation"]["normalized_subcase_headers"], 1)
            self.assertTrue(summary["diagnostic_only"])
            self.assertEqual(summary["status"], "exported_with_runtime_diagnostics")

    def test_generated_aero_points_are_checked_and_preserved_separately(self):
        extra = "".join(f" {nid} G 0.0 0.0 0.001 0.0 0.002 0.0\n" for nid in range(6, 11))
        marker = " S T R A I N S   I N   Q U A D R I L A T E R A L"
        with tempfile.TemporaryDirectory() as temporary:
            folder, output = Path(temporary) / "case", Path(temporary) / "output"
            write_case(folder, F06.replace(marker, extra + marker))
            # One box center plus four drawing corners: five generated points.
            (folder / "aero_cards.bdf").write_text("CAERO1,100,1,0,1,1,,,1\n+,0.,0.,1.,1.,0.,0.,0.,1.\n")
            deck = folder / "bwb_sol144.bdf"
            deck.write_text(deck.read_text().replace("ENDDATA", "INCLUDE 'aero_cards.bdf'\nENDDATA"))
            summary = export_case(folder, output, make_plot=False)
            self.assertEqual(summary["solver_generated_aero_points_per_subcase"], 5)
            self.assertEqual(summary["displacement_csv_rows"], 5)
            self.assertTrue((output / "solver_generated_aero_displacements.csv").is_file())
            altered = (folder / "bwb_sol144.f06").read_text().replace(" 10 G 0.0 0.0 0.001 0.0 0.002 0.0\n", "")
            (folder / "bwb_sol144.f06").write_text(altered)
            with self.assertRaisesRegex(ValueError, "Unexpected nonstructural GRID"):
                validate_case(folder)


if __name__ == "__main__":
    unittest.main()
