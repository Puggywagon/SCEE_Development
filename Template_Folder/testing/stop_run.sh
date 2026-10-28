#!/bin/bash
# stop_run.sh: cleanly stop a running SCEE run on THIS machine.
# Stops the launcher first so no new Gaussian jobs start, then running jobs.
# Finished gas-phase outputs are kept; a restart only reruns unfinished ones.
pgrep -af "Run_SCEE|run_.sh|g09|gmx mdrun" || { echo "Nothing running."; exit 0; }
read -p "Stop all of the above? [y/N] " ans
[ "$ans" = "y" ] || exit 0
pkill -f Run_SCEE.py
pkill -f "run_.sh"
pkill g09
pkill -f 'l[0-9]+\.exe'
pkill -f "gmx mdrun"
sleep 2
pgrep -af "g09|l[0-9]+\.exe|gmx mdrun" && echo "Some processes survived; rerun or use kill -9." || echo "All stopped."
