import contextlib
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from core.data_generation import generate_data_sets
from core.paths import BASE_DIR, CORE_DIR, DATA_DIR, WORKSPACE_DIR, require_working_output
from core.run_local_nastran import main as solve


class ProjectLayoutTest(unittest.TestCase):
    def test_generation_from_another_directory_matches_the_verified_case(self):
        WORKSPACE_DIR.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=WORKSPACE_DIR) as directory:
            target = Path(directory).resolve()
            self.assertTrue(target.is_relative_to(WORKSPACE_DIR))
            original_cwd = Path.cwd()
            try:
                os.chdir(target)
                with contextlib.redirect_stdout(io.StringIO()):
                    generate_data_sets(n=1, temp_dir=target/'cases', seed=20261017)
            finally:
                os.chdir(original_cwd)
            case = target/'cases/sample_1'
            manifest = json.loads((case/'case_manifest.json').read_text())
            self.assertEqual(manifest['modified_element_count'],110)
            # Same verified seed, model and damage logic as before the directory move.
            self.assertEqual(hashlib.sha256((case/'FEM_only.bdf').read_bytes()).hexdigest(),
                             'be327b3cb26f49246d0de527dd660c7175fde1b71666ca8d7f28e745d281e484')
            self.assertEqual(manifest['case_input_sha256']['sol144.bdf'],
                             hashlib.sha256((BASE_DIR/'sol144.bdf').read_bytes()).hexdigest())
            self.assertTrue((DATA_DIR/'nodes.csv').exists())

    def test_core_is_not_a_generated_output_directory(self):
        for target in (CORE_DIR, BASE_DIR, DATA_DIR/'sample_1'):
            with self.assertRaises(ValueError):
                require_working_output(target)
        with patch.object(sys,'argv',['core.run_local_nastran',str(BASE_DIR),'--inspect-only']), \
             contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as raised:
            solve()
        self.assertEqual(raised.exception.code,2)

    def test_solver_dispatch_uses_the_new_case_location(self):
        with tempfile.TemporaryDirectory(dir=WORKSPACE_DIR) as directory:
            case = Path(directory).resolve()
            self.assertTrue(case.is_relative_to(WORKSPACE_DIR))
            (case/'sol144.bdf').write_text('SOL 144\nCEND\nBEGIN BULK\nENDDATA\n')
            solver = case/'solver.exe'
            solver.touch()

            def fake_solver(command, **kwargs):
                (case/'sol144.f06').write_text('* * * END OF JOB * * *\n')
                (case/'sol144.log').write_text('Analysis complete  0\n')
                return subprocess.CompletedProcess(command,0)

            with patch.object(sys,'argv',['core.run_local_nastran',str(case),'--solver',str(solver)]), \
                 patch('core.run_local_nastran.subprocess.run',side_effect=fake_solver) as launch, \
                 contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(solve(),0)
            command = launch.call_args.args[0]
            self.assertEqual(command[:2],[str(solver),'sol144.bdf'])
            self.assertIn(f'sdirectory={case / "scratch"}',command)
            self.assertEqual(launch.call_args.kwargs['cwd'],case)
            self.assertTrue(json.loads((case/'run-summary.json').read_text())['solver_completed'])


if __name__ == '__main__':
    unittest.main()
