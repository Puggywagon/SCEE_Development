#!/bin/bash
# launch.sh: run from the folder containing Run_SCEE.py and Settings.yml
# Edit the two paths below for your machine.
source ~/venvs/scee/bin/activate
source /path/to/gromacs/bin/GMXRC
nohup python3 -u Run_SCEE.py > run_$(date +%Y%m%d_%H%M).log 2>&1 &
echo "Started PID $!"
