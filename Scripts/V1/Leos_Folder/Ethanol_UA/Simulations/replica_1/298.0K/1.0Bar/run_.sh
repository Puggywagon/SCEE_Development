#!/bin/bash
cd opt-test-SCEE/
export g09root=/home/zoe/Software/Gaussian/g09_pgi/
source $g09root/g09/bsd/g09.profile
SCRATCH_DIRS=(/home/zoe/Research/Gaussian/scratch)
MAXJOBS=2
i=0
for datfile in $(ls *_c*_q?.dat); do
  outfile=${datfile%.dat}.out
  idx=$(( i % ${#SCRATCH_DIRS[@]} ))
  GAUSS_SCRDIR=${SCRATCH_DIRS[$idx]}
  echo "Preparing $datfile -> $outfile on $GAUSS_SCRDIR"
  while true; do
      g09c=$(pgrep -c g09 2>/dev/null)
      mdc=$(pgrep -c mdrun 2>/dev/null)
      g09c=${g09c:-0}
      mdc=${mdc:-0}
      if [ $((g09c + mdc)) -lt "$MAXJOBS" ]; then
          break
      fi
      sleep 5
  done
  GAUSS_SCRDIR=$GAUSS_SCRDIR g09 < "${datfile}" > "${outfile}" &
  i=$((i+1))
done
wait
