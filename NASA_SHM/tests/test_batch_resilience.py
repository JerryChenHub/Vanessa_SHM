import contextlib
import ctypes
import io
import json
import os
from pathlib import Path
import tempfile
import threading
import time
import unittest
from unittest.mock import Mock, patch

from core import batch_run
from core.file_io import replace_file


class BatchResilienceTest(unittest.TestCase):
    @unittest.skipUnless(os.name=='nt','Windows sharing-lock regression')
    def test_atomic_replace_waits_for_a_real_windows_reader(self):
        kernel = ctypes.WinDLL('kernel32',use_last_error=True)
        kernel.CreateFileW.argtypes = [ctypes.c_wchar_p,ctypes.c_uint,ctypes.c_uint,ctypes.c_void_p,
                                       ctypes.c_uint,ctypes.c_uint,ctypes.c_void_p]
        kernel.CreateFileW.restype = ctypes.c_void_p
        kernel.CloseHandle.argtypes = [ctypes.c_void_p]
        with tempfile.TemporaryDirectory() as directory:
            target,source = Path(directory)/'status.json',Path(directory)/'status.tmp'
            target.write_text('old'); source.write_text('new')
            # Read/write sharing allowed; deletion/rename sharing withheld.
            handle = kernel.CreateFileW(str(target),0x80000000,3,None,3,0,None)
            self.assertNotEqual(handle,ctypes.c_void_p(-1).value)
            try:
                with self.assertRaises(PermissionError):
                    replace_file(source,target,timeout=0)
                self.assertEqual(target.read_text(),'old')
            finally:
                kernel.CloseHandle(handle)
            handle = kernel.CreateFileW(str(target),0x80000000,3,None,3,0,None)
            release = threading.Timer(0.15,kernel.CloseHandle,args=(handle,))
            release.start()
            try:
                replace_file(source,target,timeout=2)
            finally:
                release.join()
            self.assertEqual(target.read_text(),'new')
            self.assertFalse(source.exists())

    def test_locked_progress_does_not_fail_a_checkpoint(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            plan = {'count':1,'cases':[]}
            original = batch_run.write_json

            def locked_status(path,data):
                if path.name=='status.json':
                    raise PermissionError('Reader temporarily denies replacement')
                original(path,data)

            with patch.object(batch_run,'write_json',side_effect=locked_status), \
                 contextlib.redirect_stdout(io.StringIO()):
                status = batch_run.publish(output,plan,'running')
            self.assertEqual(status['state'],'running')
            self.assertEqual(json.loads((output/'batch.json').read_text()),plan)
            self.assertTrue((output/'cases.csv').exists())

    def test_solver_exception_reaps_child_and_uses_clean_attempt(self):
        with tempfile.TemporaryDirectory() as directory:
            case = Path(directory)
            (case/'sol144.bdf').write_text('SOL 144\n')
            (case/'subcases.inc').write_text('SUBCASE 1\n')
            (case/'sol144.xdb').write_text('old database')
            (case/'sol144.f06').write_text('old output')
            process = Mock(pid=123456)
            process.poll.return_value = None
            plan = {'solver':'solver.exe','license_server':'test','memory':'1GB','threads':2}
            item = {'case_id':'case_001'}
            with patch.object(batch_run,'publish',side_effect=[None,RuntimeError('checkpoint failed')]), \
                 patch.object(batch_run,'wait_for_memory'), \
                 patch.object(batch_run.subprocess,'Popen',return_value=process) as start, \
                 patch.object(batch_run.subprocess,'run') as terminate:
                with self.assertRaisesRegex(RuntimeError,'checkpoint failed'):
                    batch_run.run_solver(case,case,case,plan,item)
            active = Path(start.call_args.args[0][3])
            self.assertEqual(active.parent,case)
            self.assertTrue((active/'sol144.bdf').exists())
            self.assertTrue((active/'subcases.inc').exists())
            self.assertFalse((active/'sol144.xdb').exists())
            self.assertFalse((active/'sol144.f06').exists())
            self.assertEqual((case/'sol144.xdb').read_text(),'old database')
            process.wait.assert_called_once_with(timeout=30)
            if os.name=='nt':
                self.assertEqual(terminate.call_args.args[0],['taskkill','/PID','123456','/T','/F'])
            else:
                process.kill.assert_called_once()

    def test_resume_rejects_a_solver_directory_outside_its_case(self):
        with tempfile.TemporaryDirectory() as directory:
            case = Path(directory)/'case_001'
            self.assertEqual(batch_run.saved_solver_case(case,{}),case.resolve())
            with self.assertRaises(ValueError):
                batch_run.saved_solver_case(case,{'solver_directory':str(case.parent/'solve_other')})

    def test_low_commit_waits_without_failing_the_case(self):
        low = {'available_commit':600*1024**2,'available_physical':4*1024**3}
        ready = {'available_commit':3*1024**3,'available_physical':2*1024**3}
        plan = {'memory':'1GB'}
        item = {'case_id':'case_001','state':'running','attempts':1}
        with patch.object(batch_run,'available_memory',side_effect=[low,ready]), \
             patch.object(batch_run,'publish') as publish, \
             patch.object(batch_run.time,'sleep') as sleep, \
             contextlib.redirect_stdout(io.StringIO()):
            batch_run.wait_for_memory(Path('.'),plan,item)
        self.assertEqual(publish.call_args.args[2],'waiting_resources')
        self.assertEqual(item,{'case_id':'case_001','state':'running','attempts':1})
        sleep.assert_called_once_with(30)

    def test_memory_failure_is_recognized_even_without_a_fatal_keyword(self):
        with tempfile.TemporaryDirectory() as directory:
            case=Path(directory)
            self.assertFalse(batch_run.memory_allocation_failed(case))
            (case/'sol144.log').write_text('Unable to allocate requested memory.\nRequested size = 1 Gbytes.\nOperating system error code = 12.\n')
            self.assertTrue(batch_run.memory_allocation_failed(case))
            (case/'sol144.log').write_text('Analysis complete  0\n')
            self.assertFalse(batch_run.memory_allocation_failed(case))

    def test_eta_uses_current_run_speed_after_a_restart(self):
        with tempfile.TemporaryDirectory() as directory:
            plan={'count':10,'cases':[],'run_cycle_seconds':[80,60]}
            status=batch_run.publish(Path(directory),plan,'running')
            self.assertEqual(status['seconds_per_case_estimate'],70)
            self.assertEqual(status['remaining_seconds_estimate'],700)
            status=batch_run.publish(Path(directory),plan,'waiting_resources')
            self.assertIsNone(status['estimated_finish_local'])


if __name__=='__main__':
    unittest.main()
