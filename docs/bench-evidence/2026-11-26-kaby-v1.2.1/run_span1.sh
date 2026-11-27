#!/bin/bash
# Kaby Lake, v1.2.1 fleet leg 3 (PLAN.md on main): bench_span1, `main` against
# `dev/v1.2.1`, two separate build trees, alternating, three fresh processes
# each. Each run is tools/quiet_bench.sh (gate: ambient < 5% over two
# consecutive 10 s windows), from the branch's copy of the script for both.
# The user's own background processes in SUSPEND are stopped (SIGSTOP) before
# each gate and resumed (SIGCONT) after each run [user, 2026-11-26; the
# technique of libstats' 2026-11-24 M1 leg]. Only processes owned by the
# invoking user are signalled. A gate that does not pass resumes them, waits
# two minutes and tries again, up to TRIES times per run.
set -u
BR=/Users/wolfman/Development/corvus/.claude/worktrees/corvus-v1-2-1-setup-0a6c67
MN=/Users/wolfman/Development/corvus-main
O=$BR/docs/bench-evidence/2026-11-26-kaby-v1.2.1
LOG=$O/run.log
TRIES=6
SUSPEND=${SUSPEND:-'mediaanalysisd corespotlightd spotlightknowledged mdworker_shared mdbulkimport storeuid CGPDFService VTEncoderXPCService VTDecoderXPCService photolibraryd photoanalysisd OneDrive contactsd AddressBookSourceSync cloudd bird suggestd ScreenSaverEngine'}
caffeinate -is -w $$ &

log() { echo "$(date -u '+%FT%TZ') load=$(sysctl -n vm.loadavg | awk '{print $2}') $*" >> "$LOG"; }
hold() { for d in $SUSPEND; do pkill "-$1" -u "$(id -u)" -x "$d"; done; }
trap 'hold CONT' EXIT
trap 'exit 1' INT TERM

{
  echo "branch:   $(git -C $BR rev-parse --short HEAD) dirty=$(git -C $BR status --porcelain --untracked-files=no | wc -l | tr -d ' ')"
  echo "main:     $(git -C $MN rev-parse --short HEAD) (+ tests/bench_span1.cpp from the branch and its two CMake lines)"
  echo "machine:  $(sysctl -n hw.model) $(sysctl -n machdep.cpu.brand_string), macOS $(sw_vers -productVersion) ($(sw_vers -buildVersion))"
  echo "compiler: $(clang++ --version | head -1)"
  echo "power:    $(pmset -g batt | head -1)"
  echo "suspended during each run: $SUSPEND"
  echo "libcorvus.a bytes: branch $(stat -f %z $BR/build/libcorvus.a) main $(stat -f %z $MN/build/libcorvus.a)"
  echo "bench_span1 sha1: branch $(shasum $BR/build/tests/bench_span1 | cut -c1-12) main $(shasum $MN/build/tests/bench_span1 | cut -c1-12)"
} > "$O/environment.txt"

one() {  # one <name> <tree> <rep>
  local name=$1 tree=$2 rep=$3 try rc
  for try in $(seq 1 $TRIES); do
    hold STOP
    ( cd "$tree" && "$BR/tools/quiet_bench.sh" -b build -t AVX2 -m 5 -s 10 -w 5 bench_span1 ) > /dev/null 2>&1; rc=$?
    hold CONT
    mkdir -p "$O/$name$rep"
    if [ $rc -eq 3 ]; then
      cp "$tree/build/quiet_bench/quiet_bench.log" "$O/$name$rep/gate_failed_try$try.log"
      rm -f "$tree/build/quiet_bench/quiet_bench.log"
      log "$name$rep try $try: no quiet window in 5 min; resumed, retry in 2 min"; sleep 120; continue
    fi
    cp "$tree/build/quiet_bench/quiet_bench.log" "$tree/build/quiet_bench/quiet_bench_bench_span1.txt" "$O/$name$rep/"
    rm -f "$tree/build/quiet_bench/quiet_bench.log"
    log "$name$rep try $try: rc=$rc"; return $rc
  done
  log "$name$rep: SKIPPED, no quiet window in $TRIES tries"; return 3
}

log "start"
for rep in 1 2 3; do
  one main "$MN" $rep
  one branch "$BR" $rep
done
log "DONE"
