#!/usr/bin/python3

import glob
import os
import sys
from shutil import which

def setup_gromacs_env():
    home = os.path.expanduser("~")
    gromacs_root = os.path.join(home, "Software", "gromacs-2024.2")
    gromacs_bin = os.path.join(gromacs_root, "bin")

    if gromacs_bin not in os.environ["PATH"]:
        os.environ["PATH"] += os.pathsep + gromacs_bin

    if which("gmx") is None:
        raise EnvironmentError(f"'gmx' not found in PATH. Expected in: {gromacs_bin}")

    py_version = f"python{sys.version_info.major}.{sys.version_info.minor}"
    gmxapi_dir = os.path.join(gromacs_root, "lib", py_version, "site-packages")
    if os.path.isdir(gmxapi_dir) and gmxapi_dir not in sys.path:
        sys.path.insert(0, gmxapi_dir)

setup_gromacs_env()

try:
    import gromacs
except ImportError:
    try:
        import gmx as gromacs
    except ImportError:
        try:
            import gmxapi as gromacs
        except ImportError:
            raise ImportError("Could not import 'gromacs', 'gmx', or 'gmxapi'. Check your PYTHONPATH and installation.")

import pandas as pd
import re
from subprocess import check_call

##############################################################################################
class Bash_steps(object):
    def __init__(self):
        pass
##############################################################################################
    def idrive_molecule_dir(self,system_title):
        text = 'cd ../\n'
        text += f'sudo cp -r {system_title} ~/idrive/MaterialsModelling/Zoe/SCEE_Results/{system_title}\n'
        text += f'cd {system_title}\n'
        text += f'ln -s ~/idrive/MaterialsModelling/Zoe/SCEE_Results/{system_title}/Simulations/ Results\n'
        
        f=open(f'idrive_{system_title}_dir.sh','w')
        f.write(text)
        f.close()

        cmd = ['bash', f'create_{system_title}_dir.sh']
        check_call(cmd)
##############################################################################################
