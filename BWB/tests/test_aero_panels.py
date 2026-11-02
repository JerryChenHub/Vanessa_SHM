"""Meaningful input-validation regressions for the BWB aerodynamic add-on."""
from __future__ import annotations

from contextlib import contextmanager
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import unittest

import numpy as np

from BWB.scripts.check_aero_panels import (
    ValidationError, check_aero, check_physical_records, grouped_cards, read_aero_geometry, read_w2gj,
    source_path,
)
from BWB.scripts.build_bwb_aero_panels import normalwash_from_mean_section


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "BWB/models/aero_panels"


def change_field(path, card_name, field_index, replacement, predicate=lambda fields: True):
    """Mutate an actual card field without depending on its storage format."""
    fields, location = next((f, loc) for f, loc in grouped_cards(path)[card_name] if predicate(f))
    del fields
    lines = path.read_text(encoding="ascii").splitlines()
    line_index = int(location.rsplit(":", 1)[1]) - 1
    first = lines[line_index]
    if "," in first[:8]:
        line_index += (field_index - 1) // 8
        tokens = lines[line_index].split(",")
        tokens[1 + (field_index - 1) % 8] = replacement
        lines[line_index] = ",".join(tokens)
    else:
        width = 16 if first[:8].strip().endswith("*") else 8
        per_line = 64 // width
        line_index += (field_index - 1) // per_line
        start = 8 + width * ((field_index - 1) % per_line)
        line = lines[line_index].ljust(start + width)
        if len(replacement) > width:
            raise ValueError("Mutation value does not fit field")
        lines[line_index] = line[:start] + replacement.rjust(width) + line[start + width:]
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


class AeroCardTests(unittest.TestCase):
    def parse(self, text):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "cards.bdf"
            path.write_text(text, encoding="ascii")
            return grouped_cards(path)

    def geometry_text(self):
        return (
            "CORD2R,1,0,0.,0.,0.,0.,1.,0.\n"
            "+,1.,0.,0.\n"
            "AEROS,1,0,.75,1.,.75,0,0\n"
            "PAERO1,10\n"
            "CAERO1,100,10,0,2,2,0,0,1\n"
            "+,0.,0.,1.,1.,0.,0.,0.,.5\n"
        )

    def test_trapezoid_area_box_order_and_right_handed_aero_frame(self):
        result = read_aero_geometry(self.parse(self.geometry_text()))
        self.assertEqual(sorted(result["boxes"]), [100, 101, 102, 103])
        self.assertEqual(result["boxes"][101]["chord_index"], 1)
        self.assertEqual(result["boxes"][102]["span_index"], 1)
        self.assertAlmostEqual(sum(b["area"] for b in result["boxes"].values()), 0.75)
        np.testing.assert_allclose(result["up"], [0, 1, 0])
        np.testing.assert_allclose(result["macros"][100]["normal"], [0, 1, 0])

    def test_overlapping_box_id_ranges_are_rejected(self):
        text = self.geometry_text() + "CAERO1,102,10,0,2,2,0,0,1\n+,0.,0.,2.,1.,0.,0.,1.,.5\n"
        with self.assertRaisesRegex(ValidationError, "Overlapping aerodynamic box"):
            read_aero_geometry(self.parse(text))

    def test_w2gj_reserved_field_and_missing_row_are_rejected(self):
        with self.assertRaisesRegex(ValidationError, "reserved field"):
            read_w2gj(self.parse("DMI,W2GJ,0,2,2,0,2,1\nDMI,W2GJ,1,1,.1\n"), 2)
        with self.assertRaisesRegex(ValidationError, "Missing W2GJ rows"):
            read_w2gj(self.parse("DMI,W2GJ,0,2,2,0,,2,1\nDMI,W2GJ,1,1,.1\n"), 2)

    def test_initial_normalwash_sign_scale_and_mirror_invariance(self):
        # Nose-up h = -0.02 * X, despite a two-metre chord and translation.
        x = np.linspace(3.0, 5.0, 9)
        points = np.column_stack((x, -0.02 * x + 1, np.full_like(x, 4.0)))
        for xi in (0.075, 0.475, 0.975):
            self.assertAlmostEqual(normalwash_from_mean_section(points, [0, 1, 0], xi), 0.02)
        points[:, 2] = 0.03 * x + 4
        normal = np.array([0.0, 0.6, 0.8])
        original = normalwash_from_mean_section(points, normal, 0.75)
        points[:, 2] *= -1
        normal[2] *= -1
        self.assertAlmostEqual(normalwash_from_mean_section(points, normal, 0.75), original)
        points[:, 1:] = 0
        self.assertEqual(normalwash_from_mean_section(points, [0, 1, 0], 0.75), 0.0)

    def test_excessive_free_field_width_and_record_length_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "cards.bdf"
            path.write_text("DMI,W2GJ,1,1,1.234567890123E-02\n", encoding="ascii")
            with self.assertRaisesRegex(ValidationError, "eight characters"):
                check_physical_records(path)
            path.write_text("GRID," + "1," * 40 + "1\n", encoding="ascii")
            with self.assertRaisesRegex(ValidationError, "80 columns"):
                check_physical_records(path)

    def test_provenance_normalizes_only_explicit_crlf_mode(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / "source.txt"
            lf = b"alpha\nbeta\n"
            record = dict(path="source.txt", sha256=hashlib.sha256(lf).hexdigest(),
                          hash_mode="crlf_to_lf")
            for content in (lf, b"alpha\r\nbeta\r\n"):
                path.write_bytes(content)
                self.assertEqual(source_path(root, record), path.resolve())
            raw_record = {key: value for key, value in record.items() if key != "hash_mode"}
            with self.assertRaisesRegex(ValidationError, "Source SHA256 mismatch"):
                source_path(root, raw_record)
            path.write_bytes(b"alpha\r\ngamma\r\n")
            with self.assertRaisesRegex(ValidationError, "Source SHA256 mismatch"):
                source_path(root, record)
            record["hash_mode"] = "unknown"
            with self.assertRaisesRegex(ValidationError, "Unknown provenance hash mode"):
                source_path(root, record)

    def test_invalid_mean_section_order_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "increase strictly"):
            normalwash_from_mean_section([[0, 0, 0], [0, 0.1, 0], [1, 0, 0]], [0, 1, 0], 0.75)


class GeneratedAeroTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not (GENERATED / "part_1/aero_manifest.json").is_file():
            raise unittest.SkipTest("Generate BWB aerodynamic inputs before artifact regressions")

    @contextmanager
    def copied_part(self):
        with tempfile.TemporaryDirectory() as temporary:
            target = Path(temporary) / "part_1"
            shutil.copytree(GENERATED / "part_1", target)
            yield target

    def test_both_exported_parts_pass_independent_checks(self):
        for part in ("1", "2"):
            with self.subTest(part=part):
                result = check_aero(GENERATED / f"part_{part}", ROOT)
                self.assertEqual(result["status"], "PASS")
                self.assertEqual(result["caero1_macros"], 11)
                self.assertEqual(result["aerodynamic_boxes"], 590)
                self.assertEqual(result["w2gj_rows"], 590)

    def test_stale_source_hash_is_rejected(self):
        with self.copied_part() as folder:
            path = folder / "aero_manifest.json"
            manifest = json.loads(path.read_text(encoding="utf-8"))
            manifest["source_files"][0]["sha256"] = "0" * 64
            path.write_text(json.dumps(manifest), encoding="utf-8")
            with self.assertRaisesRegex(ValidationError, "Source SHA256 mismatch"):
                check_aero(folder, ROOT)

    def test_unknown_spline_support_node_is_rejected(self):
        with self.copied_part() as folder:
            manifest = json.loads((folder / "aero_manifest.json").read_text(encoding="utf-8"))
            path = folder / manifest["files"]["spline"]
            change_field(path, "SET1", 2, "9999999")
            with self.assertRaisesRegex(ValidationError, "unknown structural nodes"):
                check_aero(folder, ROOT)

    def test_changed_camber_value_is_rejected(self):
        with self.copied_part() as folder:
            manifest = json.loads((folder / "aero_manifest.json").read_text(encoding="utf-8"))
            path = folder / manifest["files"]["camber"]
            change_field(path, "DMI", 4, "0.56789", lambda f: f[1] == "W2GJ" and f[2] == "1")
            with self.assertRaisesRegex(ValidationError, "W2GJ"):
                check_aero(folder, ROOT)


if __name__ == "__main__":
    unittest.main()
