"""Atomic file publication with bounded retries for Windows sharing locks."""
import os
from pathlib import Path
import time


def replace_file(source, destination, timeout=5.0):
    deadline = time.monotonic() + timeout
    delay = 0.05
    while True:
        try:
            os.replace(source, destination)
            return
        except PermissionError:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise
            time.sleep(min(delay, remaining))
            delay = min(delay * 2, 0.5)


def atomic_text(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(text, encoding='utf-8')
    replace_file(temporary, path)
