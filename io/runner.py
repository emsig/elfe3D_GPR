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

import shlex
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence, Union


@dataclass
class ProjectPaths:
    """
    Machine-level path configuration.  Set once at the top of every notebook.

    Parameters
    ----------
    master_dir : absolute path to the master\\ folder
    exec_rel   : path to elfe3d_gpr relative to master_dir
                 default: elfe3D_GPR\\elfe3d_gpr
    use_wsl    : True  : prefix subprocess calls with "wsl" and convert
                         Windows paths to /mnt/<drive>/... format
                 False : run natively on Linux
    """
    master_dir: Union[str, Path]
    exec_rel:   str = ""
    use_wsl:    bool = True

    def exec_path(self) -> Path:
        """Absolute path to the elfe3d_gpr executable."""
        exec_rel = self.exec_rel if self.exec_rel != "" else "elfe3d_gpr"
        return Path(self.master_dir) / exec_rel

    def to_target_path(self, p: Union[str, Path]) -> str:
        """Convert a path to the runtime target format."""
        p = Path(p)
        if self.use_wsl:
            s = str(p).replace("\\", "/")
            if len(s) >= 2 and s[1] == ":":
                s = f"/mnt/{s[0].lower()}/{s[3:]}"
            return s
        return str(p)

    def _prefix(self, *args) -> list[str]:
        return (["wsl"] + list(args)) if self.use_wsl else list(args)


def _stream(
    cmd: Union[list[str], str],
    label: str,
    cwd: Union[str, Path, None] = None,
    shell: bool = False,
) -> None:
    """
    Run *cmd* and stream every output line to the notebook cell in real time.
    Raises RuntimeError on non-zero exit.
    """
    import sys
    label_line = f"── {label} {'─' * max(1, 50 - len(label))}"
    # print(label_line)
    # if isinstance(cmd, list):
    #     print(f"   {' '.join(cmd)}\n")
    # else:
    #     print(f"   {cmd}\n")

    t0 = time.time()
    proc = subprocess.Popen(
        cmd,
        cwd=str(cwd) if cwd is not None else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        shell=shell,
    )

    for line in proc.stdout:
        print(line, end="")
        sys.stdout.flush()

    proc.wait()
    elapsed = time.time() - t0

    if proc.returncode != 0:
        raise RuntimeError(f"{label} failed (exit {proc.returncode})")

    print(f"\nDone in {elapsed:.1f}s")


def _command_tokens(binary: str, flags: str, extra: Sequence[str] | None = None) -> list[str]:
    tokens = [binary] + shlex.split(flags)
    if extra:
        tokens += list(extra)
    return tokens


def run_tetgen(
    paths: ProjectPaths,
    poly_file,
    flags: str = "-pq1.2kAaen",
    binary: str = "/usr/bin/tetgen",
) -> None:
    """
    Mesh the geometry with TetGen.

    On Windows with WSL enabled the command is run inside WSL.
    On native Linux the command is run directly from the mesh folder.

        run_tetgen(paths, survey.io.poly_file)
    """
    poly_path = Path(poly_file)
    cmd = _command_tokens(binary, flags, [poly_path.name])

    if paths.use_wsl:
        wsl_dir = paths.to_target_path(poly_path.parent)
        wsl_cmd = shlex.join(cmd)
        bash_cmd = f"cd {wsl_dir} && {wsl_cmd}"
        cmd = ["wsl", "bash", "-c", bash_cmd]
        cwd = None
    else:
        cwd = poly_path.parent

    _stream(cmd, "TetGen", cwd=cwd)


def run_solver(paths: ProjectPaths, survey) -> None:
    """
    Run the elfe3d_gpr FEM solver.

    On Windows with WSL enabled the executable is run inside WSL.
    On native Linux it is run directly from the survey input folder.
    """
    input_dir = Path(survey.io.base_dir)
    exe_path = paths.exec_path()

    if paths.use_wsl:
        wsl_dir = paths.to_target_path(input_dir)
        wsl_exe = paths.to_target_path(exe_path)
        bash_cmd = f"cd {wsl_dir} && {wsl_exe}"
        cmd = ["wsl", "bash", "-c", bash_cmd]
        cwd = None
    else:
        cmd = [str(exe_path)]
        cwd = input_dir

    _stream(cmd, "elfe3d_gpr", cwd=cwd)


def run_custom_command(
    paths: ProjectPaths,
    command,
    cwd=None,
    use_wsl: bool | None = None,
    shell: bool = False,
) -> None:
    """
    Run a custom command with the same WSL/native routing rules.

    Parameters
    ----------
    command : list[str] or str
        If list, it is executed directly. If str and shell=False, it is split.
    cwd : path-like or None
        Working directory for the command.
    use_wsl : bool | None
        Use WSL when True, native Linux when False, otherwise use paths.use_wsl.
    shell : bool
        Pass to subprocess.Popen when running natively.
    """
    if use_wsl is None:
        use_wsl = paths.use_wsl

    if isinstance(command, str) and not shell:
        command = shlex.split(command)

    if use_wsl:
        cwd_path = Path(cwd) if cwd is not None else Path.cwd()
        wsl_cwd = paths.to_target_path(cwd_path)
        command_text = shlex.join(command) if isinstance(command, list) else command
        bash_cmd = f"cd {wsl_cwd} && {command_text}"
        cmd = ["wsl", "bash", "-c", bash_cmd]
        cwd = None
        shell = False
    else:
        cmd = command
        cwd = Path(cwd) if cwd is not None else None

    _stream(cmd, "custom command", cwd=cwd, shell=shell)
