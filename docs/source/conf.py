import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

project = 'elfe3D_GPR'
author = 'elfe3D_GPR developers'
release = '0.9.0-beta'

extensions = [
    'myst_parser',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

myst_enable_extensions = [
    'deflist',
    'colon_fence',
]

import subprocess
from pathlib import Path


def run_doxygen():
    repo_root = Path(__file__).resolve().parents[2]
    doxyfile = repo_root / 'Doxyfile'
    if not doxyfile.exists():
        print('Doxygen config not found at', doxyfile)
        return
    try:
        subprocess.run(['doxygen', str(doxyfile)], cwd=str(repo_root), check=True)
    except FileNotFoundError:
        print('Doxygen executable not found; skipping Fortran docs generation.')
    except subprocess.CalledProcessError:
        print('Doxygen failed; Fortran docs may be outdated.')


def setup(app):
    run_doxygen()
