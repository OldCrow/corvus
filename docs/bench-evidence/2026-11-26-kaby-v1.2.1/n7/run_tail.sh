#!/bin/bash
# bench_tail over several libraries, alternating, REPS fresh processes each.
# The user's background processes are stopped (SIGSTOP) for the whole pass
# and resumed (SIGCONT) at exit; ambient CPU is sampled before and after.
# usage: run_tail.sh <outdir> <bindir> <reps> <name>...
set -u
O=$1; B=$2; REPS=$3; shift 3
HERE=$(cd "$(dirname "$0")" && pwd)
SUSPEND=$(sed -n 's/^suspended during each run: //p' "$HERE/../environment.txt")
hold() { for d in $SUSPEND; do pkill "-$1" -u "$(id -u)" -x "$d"; done; }
amb() { /usr/bin/top -l 2 -s 10 -n 0 | grep 'CPU usage' | tail -1 | sed -E 's/.* ([0-9.]+)% idle.*/\1/' | awk '{printf "%.2f", 100-$1}'; }
trap 'hold CONT' EXIT
trap 'exit 1' INT TERM
mkdir -p "$O"
caffeinate -is -w $$ &
hold STOP
echo "$(date -u '+%FT%TZ') ambient before: $(amb)%" >> "$O/run.log"
for rep in $(seq 1 "$REPS"); do
  for name in "$@"; do
    "$B/bt_$name" > "$O/${name}_$rep.txt" 2>&1 || echo "$name rep $rep: rc=$?" >> "$O/run.log"
  done
  echo "$(date -u '+%FT%TZ') rep $rep done, ambient: $(amb)%" >> "$O/run.log"
done
hold CONT
