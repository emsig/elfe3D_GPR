"""
runner.py
=========
  1. Knows where the executable lives and handles WSL  →  ProjectPaths
  2. Calls TetGen on the .poly file from survey.io     →  run_tetgen()
  3. Calls elfe3d_gpr                                  →  run_solver()

Usage in a notebook
-------------------
    from runner  import ProjectPaths, run_tetgen, run_solver
    from survey  import GPRSurvey

    paths  = ProjectPaths(master_dir=MASTER)
    survey = GPRSurvey.build(experiment_name='air_only', base_dir=..., ...)
    survey.generate()

    run_tetgen(paths, survey.io.poly_file)
    run_solver(paths)

    # results live at:
    # survey.io.output_dir / 'electric_fields_receiver_line.txt'
"""

from __future__ import annotations

import subprocess
import time
from dataclasses import dataclass
from pathlib import Path


@dataclass
class ProjectPaths:
    """
    Machine-level path configuration.  Set once at the top of every notebook.

    Parameters
    ----------
    master_dir : absolute path to the master\\ folder
    exec_rel   : path to elfe3d_gpr relative to master_dir
                 default: elfe3D_GPR\\elfe3d_gpr
    use_wsl    : True  → prefix subprocess calls with "wsl" and convert
                         Windows paths to /mnt/<drive>/... format
                 False → run natively (already inside WSL)
    """
    master_dir: str
    exec_rel:   str  = r"elfe3D_GPR\elfe3d_gpr"
    use_wsl:    bool = True

    def exec_path(self) -> Path:
        """Absolute path to the elfe3d_gpr executable."""
        return Path(self.master_dir) / self.exec_rel

    def to_wsl(self, p) -> str:
        """Convert a Windows path to its WSL equivalent.
        C:\\Users\\foo\\bar  →  /mnt/c/Users/foo/bar
        """
        s = str(p).replace("\\", "/")
        if len(s) >= 2 and s[1] == ":":
            s = f"/mnt/{s[0].lower()}/{s[3:]}"
        return s

    def _prefix(self, *args) -> list[str]:
        return (["wsl"] + list(args)) if self.use_wsl else list(args)


def _stream(cmd: list[str], label: str) -> None:
    """
    Run *cmd* and stream every output line to the notebook cell in real time.
    Raises RuntimeError on non-zero exit.
    """
    import sys
    print(f"── {label} {'─' * (50 - len(label))}")
    print(f"   {' '.join(cmd)}\n")

    t0   = time.time()
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,   # merge stderr into stdout
        text=True,
        bufsize=1,                  # line-buffered
    )

    for line in proc.stdout:
        print(line, end="")
        sys.stdout.flush()          # forces Jupyter to render immediately

    proc.wait()
    elapsed = time.time() - t0

    if proc.returncode != 0:
        raise RuntimeError(f"{label} failed (exit {proc.returncode})")

    print(f"\n✓  done in {elapsed:.1f}s")


def run_tetgen(
    paths,
    poly_file,
    flags:  str = "-pq1.2kAaen",
    binary: str = "/usr/bin/tetgen",
) -> None:
    """
    Mesh the geometry with TetGen.

    tetgen15 is an alias in ~/.bashrc pointing to /usr/bin/tetgen — aliases
    are invisible to bash -c, so we use the real binary path directly.

        run_tetgen(paths, survey.io.poly_file)

    Output is streamed line-by-line.  Raises RuntimeError on failure.
    """
    import sys
    poly_path = Path(poly_file)
    wsl_dir   = paths.to_wsl(poly_path.parent)
    wsl_poly  = poly_path.name
    bash_cmd  = f"cd {wsl_dir} && {binary} {flags} {wsl_poly}"
    cmd       = ["wsl", "bash", "-c", bash_cmd]

    print(f"\u2500\u2500 TetGen {chr(9472) * 44}")
    print(f"   {bash_cmd}\n")

    t0   = time.time()
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    for line in proc.stdout:
        print(line, end="")
        sys.stdout.flush()
    proc.wait()

    if proc.returncode != 0:
        raise RuntimeError(f"TetGen failed (exit {proc.returncode})")
    print(f"\n\u2713  done in {time.time() - t0:.1f}s")


def run_solver(paths: ProjectPaths, survey) -> None:
    """
    Run the elfe3d_gpr FEM solver.

    Runs from survey.io.base_dir so the exe finds elfe3D_input.txt in its
    working directory — no file copying needed.  The exe is called by its
    absolute WSL path.

        run_solver(paths, survey)

    Output is streamed line-by-line.  Raises RuntimeError on failure.
    """
    input_dir = paths.to_wsl(survey.io.base_dir)
    exe_path  = paths.to_wsl(paths.exec_path())
    bash_cmd  = f"cd {input_dir} && {exe_path}"
    cmd       = ["wsl", "bash", "-c", bash_cmd]

    _stream(cmd, "elfe3d_gpr")
