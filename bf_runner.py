"""Background workers and filesystem helpers for the simulator run.

``BuildRunWorker`` builds (cmake/make when needed) and runs the C++ ``Brain_growth``
solver in a ``QThread`` so the GUI stays responsive and the Stop button can
terminate it. ``ParaviewWorker`` regenerates a result video/image under
ParaView's ``pvpython`` when it is missing. Both replace the subprocess+psutil
handoff of the old ``make_run.py``/``progress.py``.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from os.path import exists
from pathlib import Path

from PySide6.QtCore import QThread, Signal

OUTPUT_DIR = "Folder_Output"


def find_paraview() -> str | None:
    """Return the path to a ParaView.app under /Applications, or None."""
    apps = Path("/Applications")
    if not apps.is_dir():
        return None
    for f in apps.iterdir():
        if re.match("ParaView", f.name):
            return str(f)
    return None


def pvpython_path() -> str | None:
    pv = find_paraview()
    return f"{pv}/Contents/bin/pvpython" if pv else None


def collect_outputs() -> str:
    """Recreate Folder_Output and move solver artifacts into it (mirrors make_run.py)."""
    parent = os.getcwd()
    path_folder = os.path.join(parent, OUTPUT_DIR)
    if os.path.exists(path_folder):
        shutil.rmtree(path_folder)
    os.mkdir(path_folder)

    for dirs, _subdirs, files in os.walk(parent):
        if os.path.abspath(dirs) == os.path.abspath(path_folder):
            continue
        for file in files:
            if file.endswith(".vtk"):
                src = os.path.join(dirs, file)
                dst = os.path.join(path_folder, file)
                if not os.path.exists(dst):
                    shutil.move(src, path_folder)

    for extra in ("timeing.csv", "Parameters.prm"):
        src = os.path.join(parent, extra)
        if os.path.exists(src):
            shutil.copy(src, path_folder)
    return path_folder


def save_outputs_to(target_dir: str) -> None:
    """Move Folder_Output into the user-chosen directory, then clear the working copy."""
    parent = os.getcwd()
    path_folder = os.path.join(parent, OUTPUT_DIR)
    if target_dir:
        target_path = os.path.join(target_dir, OUTPUT_DIR)
        os.makedirs(target_path, exist_ok=True)
        for dirs, _subdirs, files in os.walk(path_folder):
            for file in files:
                shutil.move(os.path.join(dirs, file), target_path)
    if os.path.exists(path_folder):
        shutil.rmtree(path_folder)


def discard_outputs() -> None:
    path_folder = os.path.join(os.getcwd(), OUTPUT_DIR)
    if os.path.exists(path_folder):
        shutil.rmtree(path_folder)


class BuildRunWorker(QThread):
    """cmake/make (if needed) → run ./Brain_growth Parameters.prm <case> → collect outputs."""

    finished_ok = Signal()
    failed = Signal(str)
    output = Signal(str)            # one line of terminal output at a time
    progress = Signal(float, float)  # (current simulation time, total time)

    # Brain_growth prints e.g. "Timestep 5 @ 12.5s  @  delta t 0.1"
    _TIME_RE = re.compile(r"Timestep\s+\d+\s+@\s+([-+0-9.eE]+)s")

    def __init__(self, case: str, total_time: str = "", parent=None):
        super().__init__(parent)
        self.case = case
        try:
            self.total = float(total_time)
        except (TypeError, ValueError):
            self.total = 0.0
        self._proc: subprocess.Popen | None = None
        self._stopped = False

    def _stream(self, cmd: list[str]) -> int | None:
        """Run a command, forwarding its combined stdout/stderr line-by-line."""
        self.output.emit("$ " + " ".join(cmd))
        self._proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, bufsize=1)
        for line in self._proc.stdout:
            if self._stopped:
                break
            line = line.rstrip("\n")
            self.output.emit(line)
            m = self._TIME_RE.search(line)
            if m:
                try:
                    self.progress.emit(float(m.group(1)), self.total)
                except ValueError:
                    pass
        self._proc.wait()
        return self._proc.returncode

    def run(self) -> None:
        try:
            if not exists("Makefile"):
                self._stream(["cmake", "CMakeLists.txt"])
            if self._stopped:
                return
            if not exists("Brain_growth"):
                self._stream(["make"])
            if self._stopped:
                return
            if not exists("Brain_growth"):
                self.failed.emit("Build failed: 'Brain_growth' executable not found.\n"
                                 "Make sure deal.II and the build toolchain are installed.")
                return

            self.output.emit("")
            rc = self._stream(["./Brain_growth", "Parameters.prm", self.case])

            if self._stopped:
                return
            if rc not in (0, None):
                # The solver frequently exits non-zero when it auto-stops at the
                # mechanical-instability point, yet still writes valid output.
                # Mirror the original behaviour: proceed to the results instead
                # of treating a non-zero exit code as a failure.
                self.output.emit(f"\nSolver exited with code {rc}; collecting any output produced …")

            self.output.emit("\nGathering output files …")
            collect_outputs()
            self.finished_ok.emit()
        except Exception as exc:  # noqa: BLE001 - surface any failure to the UI
            if not self._stopped:
                self.failed.emit(str(exc))

    def stop(self) -> None:
        """Terminate the running solver (Stop button)."""
        self._stopped = True
        if self._proc and self._proc.poll() is None:
            self._proc.terminate()
            try:
                self._proc.wait(timeout=3)
            except subprocess.TimeoutExpired:
                self._proc.kill()


class ParaviewWorker(QThread):
    """Run the ParaView pvpython script for a result if its output file is missing."""

    ready = Signal(str)     # emits the output path when available
    failed = Signal(str)

    def __init__(self, output_path: str, script: str, parent=None):
        super().__init__(parent)
        self.output_path = output_path
        self.script = script

    def run(self) -> None:
        try:
            if not exists(self.output_path):
                pv = pvpython_path()
                if not pv:
                    self.failed.emit("ParaView not found in /Applications - cannot generate results.")
                    return
                subprocess.run([pv, self.script], check=False)
            if exists(self.output_path):
                self.ready.emit(self.output_path)
            else:
                self.failed.emit(f"Result file was not generated: {self.output_path}")
        except Exception as exc:  # noqa: BLE001
            self.failed.emit(str(exc))
