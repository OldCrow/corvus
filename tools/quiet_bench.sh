#!/bin/zsh
# Self-gating benchmark runner (macOS/Linux) — the Unix counterpart of
# tools/quiet_bench.ps1, same protocol: refuse to measure a busy machine,
# sample ambient CPU between targets, write a noise-annotated log beside the
# raw per-target output. See the .ps1 header for the rationale; nothing here
# stops services for you — when the gate fails it prints the top consumers.
#
# usage: tools/quiet_bench.sh [-b build] [-t AVX2] [-m 5] [-s 10] [-w 10] [bench_lgamma ...]
#   -b  build directory holding tests/bench_* (must be a Release build)
#   -t  CORVUS_EXPECT_TARGET value; unset means report-only
#   -m  max ambient CPU percent to start (default 5)
#   -s  ambient window seconds (default 10)
#   -w  give up after this many minutes (default 10)
set -u
BUILD=build; EXPECT=""; MAXAMB=5; SAMPLE=10; MAXWAIT=10
while getopts "b:t:m:s:w:" o; do case $o in b) BUILD=$OPTARG;; t) EXPECT=$OPTARG;; m) MAXAMB=$OPTARG;; s) SAMPLE=$OPTARG;; w) MAXWAIT=$OPTARG;; esac; done
shift $((OPTIND-1))
TARGETS=("$@"); [ ${#TARGETS[@]} -eq 0 ] && TARGETS=(bench_erf bench_erfc bench_lgamma bench_erfinv bench_gamma bench_beta bench_digamma bench_trigamma bench_gammainv bench_betainv bench_bessel bench_lbeta bench_trig bench_exp bench_log)
OUT=$BUILD/quiet_bench; mkdir -p $OUT; LOG=$OUT/quiet_bench.log
log() { echo "$(date '+%H:%M:%S') $*" | tee -a $LOG; }
# ambient busy % over one window: 100 - idle, from top's CPU-usage line
# (two samples; the first is since boot and is discarded)
ambient() {
  if [ "$(uname)" = Darwin ]; then
    top -l 2 -s $SAMPLE -n 0 2>/dev/null | grep 'CPU usage' | tail -1 | sed -E 's/.* ([0-9.]+)% idle.*/\1/' | awk '{printf "%.2f", 100-$1}'
  else
    vmstat $SAMPLE 2 | tail -1 | awk '{printf "%.2f", 100-$15}'
  fi
}
consumers() { ps -Ao pcpu,comm -r | head -6 | tail -5 | sed 's/^/    /'; }
grep -q 'CMAKE_BUILD_TYPE:STRING=Release' $BUILD/CMakeCache.txt || { log "ABORT: $BUILD is not a Release build"; exit 2; }
log "quiet_bench.sh  build=$BUILD  expect=${EXPECT:-<unset>}  targets=${#TARGETS[@]}"
log "gate          : ambient < $MAXAMB% over two consecutive ${SAMPLE}s windows"
log "machine       : $(sysctl -n machdep.cpu.brand_string 2>/dev/null || uname -m); $(sw_vers -productName 2>/dev/null) $(sw_vers -productVersion 2>/dev/null)"
deadline=$(( $(date +%s) + MAXWAIT*60 )); quiet=0; gate=""
while [ $quiet -lt 2 ]; do
  a=$(ambient); log "ambient sample : ${a}%"
  if awk -v a=$a -v m=$MAXAMB 'BEGIN{exit !(a<m)}'; then quiet=$((quiet+1)); gate="$gate $a"; else quiet=0; gate=""; log "  top consumers:"; consumers | tee -a $LOG; fi
  [ $(date +%s) -gt $deadline ] && { log "ABORT: machine did not go quiet within $MAXWAIT min"; exit 3; }
done
log "gate passed at:$gate%"
[ -n "$EXPECT" ] && export CORVUS_EXPECT_TARGET=$EXPECT
fail=0
for t in "${TARGETS[@]}"; do
  exe=$BUILD/tests/$t; [ -x $exe ] || { log "$t: MISSING ($exe)"; fail=1; continue; }
  n=$(ambient); log "$t: noise before ${n}%  top: $(ps -Ao pcpu,comm -r | sed -n '2,4p' | awk '{printf "%s%%%s ", $1, substr($2, match($2, /[^\/]+$/))}')"
  start=$(date +%s); $exe > $OUT/quiet_bench_$t.txt 2>&1; rc=$?; dur=$(( $(date +%s) - start ))
  tgt=$(grep -oE 'SIMD target: *[A-Z0-9_]+' $OUT/quiet_bench_$t.txt | head -1 | grep -oE '[A-Z][A-Z0-9_]+$')
  log "$t: rc=$rc ${dur}s target=${tgt:-?}"; [ $rc -eq 0 ] || fail=1
done
log "done fail=$fail. Raw per-target output: $OUT/quiet_bench_<target>.txt"
exit $fail
