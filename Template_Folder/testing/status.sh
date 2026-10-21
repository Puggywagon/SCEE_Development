#!/bin/bash
# status.sh: one-screen summary of all SCEE runs (MD + Gaussian)
# Usage: watch -n 60 ./status.sh
BASE=~/Research/SCEE_Project/Aromatics     # folder containing the molecule run folders

md_progress() {   # $1 = a GROMACS .log file
    local log=$1 nsteps step
    nsteps=$(grep -m1 -E '^\s*nsteps\s*=' "$log" | awk '{print $3}')
    step=$(grep -A1 -E '^\s*Step\s+Time' "$log" | tail -1 | awk '{print $1}')
    if [[ -n $nsteps && -n $step && $nsteps -gt 0 ]]; then
        echo "$step / $nsteps steps ($(( 100 * step / nsteps ))%)"
    else
        echo "starting"
    fi
}

tick() { [ -e "$1" ] && echo "✓" || echo "·"; }

echo "=== $(date '+%a %H:%M:%S') ==="
echo "Running scripts:"
pgrep -af Run_SCEE.py | sed 's/^/  /' || echo "  none"
echo "Active jobs: g09 = $(pgrep -c g09)   mdrun = $(pgrep -fc 'gmx mdrun')"
echo

for res in "$BASE"/*/*/Results; do
    sys=$(dirname "$res"); run=$(dirname "$sys"); name=$(basename "$sys")
    echo "--- $name ---"

    echo "  Vacuum:    em $(tick "$sys/em.gro")  NVT $(tick "$sys/nvt_vacuum2.gro")  gas configs $(ls "$sys"/Gas_Configs/gas_conf_*.gro 2>/dev/null | wc -l)"
    echo "  Box:       em $(tick "$sys/em1.gro")  NVT $(tick "$sys/nvt_eq.gro")  NPT $(tick "$sys/npt.gro")"

    live=$(find "$sys" -name '*.log' -mmin -5 -printf '%T@ %p\n' 2>/dev/null | sort -n | tail -1 | cut -d' ' -f2-)
    [ -n "$live" ] && echo "  MD now:    ${live#$sys/} — $(md_progress "$live")"

    line="  MD prod:  "
    for leaf in "$sys"/Simulations/replica_*/*/*/; do
        [ -d "$leaf" ] || continue
        r=$(basename "$(dirname "$(dirname "$leaf")")" | sed 's/replica_/r/')
        nconf=$(ls "$leaf"/conf_*.gro 2>/dev/null | wc -l)
        if ls "$leaf"/*_QMMM_md3.gro >/dev/null 2>&1; then line+=" $r ✓($nconf cfg)"
        elif [ -e "$leaf/md.log" ];                   then line+=" $r running"
        else                                               line+=" $r ·"; fi
    done
    echo "$line"

    v0=$(ls "$sys"/Simulations/replica_*/*/*/SCEE/*_q?.out 2>/dev/null | wc -l)
    v1=$(ls "$sys"/Simulations/replica_*/*/*/SCEE/*_v1.out 2>/dev/null | wc -l)
    gas=$(ls "$sys"/{Vacuum,PCM1,PCM2}/*.out 2>/dev/null | wc -l)
    echo "  Gaussian:  gas $gas/150   SCEE v0 $v0/3000   SCEE v1 $v1/3000"

    log=$(ls -t "$run"/*.log 2>/dev/null | head -1)
    [ -n "$log" ] && echo "  last log:  $(grep -v '^\s*$' "$log" | tail -1 | cut -c1-90)"
    echo
done
