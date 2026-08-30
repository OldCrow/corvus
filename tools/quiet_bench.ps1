<#
.SYNOPSIS
    Self-gating benchmark runner (Windows): refuses to measure a busy machine,
    and records the ambient noise alongside every number it produces.

.DESCRIPTION
    Benchmark numbers from a loaded machine are not wrong so much as
    meaningless -- they measure the machine's other work as much as the
    kernel's. This runner will not start until the box is demonstrably quiet,
    samples the noise again between every target, and writes both the raw
    per-target output and a noise-annotated log. The log is the evidence chain:
    a number without its noise context cannot be published later, because
    nobody can tell afterwards what the machine was doing.

    Gate: ambient CPU must be below -MaxAmbient percent across TWO consecutive
    windows before the first target runs. Two windows rather than one because a
    single quiet sample is regularly just the gap between a background task's
    bursts.

    Reconstructed 2026-08-15 from the protocol recorded in PLAN.md for the
    2026-08-12 pass, after the original was lost with its build directory --
    it had never been checked in. That is why this lives in tools/ and not in
    a build tree.

    QUIETING THE MACHINE. The 2026-08-12 pass needed APSDaemon killed (~1 core,
    constant) and WSearch (0.83 core) plus LightingService (0.13) temporarily
    stopped, on the Performance power plan. This script does NOT stop services
    for you: on your own machine that is your call, not a script's. When the
    gate fails it prints the top consumers so you know what to go after.

    Tier attribution: pass -ExpectTarget to set CORVUS_EXPECT_TARGET, and each
    bench aborts rather than reporting numbers under the wrong tier's name.
    Strongly recommended -- an uncapped Ryzen build silently reports AVX3_ZEN4
    under MSVC's blocklist as AVX2, and the numbers look perfectly plausible.

.PARAMETER BuildDir
    Build directory holding the bench executables. Must be a Release build:
    anything else measures the wrong program.

.PARAMETER Targets
    Bench targets to run. Default: all twelve. Narrow it to keep the quiet
    window short -- e.g. just bench_lgamma when only its per-region
    zone/recurrence/Stirling breakdown is wanted.

.PARAMETER ExpectTarget
    Value for CORVUS_EXPECT_TARGET, e.g. AVX3_ZEN4. Unset means report-only.

.PARAMETER MaxAmbient
    Ambient CPU percent the machine must be under to start. Default 5.

.PARAMETER SampleSeconds
    Length of each ambient window, in seconds. Default 10.

.PARAMETER MaxWaitMinutes
    Give up if the machine has not gone quiet within this long. Default 10.

.PARAMETER OutDir
    Where the logs land. Default: the build directory.

.EXAMPLE
    tools\quiet_bench.ps1 -Targets bench_lgamma -ExpectTarget AVX3_ZEN4
    Per-region lgamma numbers against libm, tier asserted.

.EXAMPLE
    tools\quiet_bench.ps1 -BuildDir build-clangcl -ExpectTarget AVX3_ZEN4
    Full twelve-target sweep.
#>
[CmdletBinding()]
param(
  [string] $BuildDir = "build-clangcl",
  [string[]] $Targets = @(
    "bench_erf", "bench_erfc", "bench_lgamma", "bench_erfinv",
    "bench_gamma", "bench_beta", "bench_digamma", "bench_trigamma",
    "bench_gammainv", "bench_betainv", "bench_bessel", "bench_lbeta",
    "bench_trig"
  ),
  [string] $ExpectTarget = "",
  [double] $MaxAmbient = 5.0,
  [int] $SampleSeconds = 10,
  [int] $MaxWaitMinutes = 10,
  [string] $OutDir = ""
)

$ErrorActionPreference = "Stop"

if ([string]::IsNullOrWhiteSpace($OutDir)) { $OutDir = $BuildDir }
if (-not (Test-Path $BuildDir)) { throw "Build directory not found: $BuildDir" }
if (-not (Test-Path $OutDir)) { New-Item -ItemType Directory -Path $OutDir | Out-Null }

$logPath = Join-Path $OutDir "quiet_bench.log"

function Write-Log {
  param([string] $Message)
  $stamp = (Get-Date).ToString("HH:mm:ss")
  $line = "[$stamp] $Message"
  Write-Host $line
  Add-Content -Path $logPath -Value $line
}

# Average total CPU over a window. One-second samples averaged, rather than a
# single long sample, so a short spike is visible in the mean instead of being
# swallowed by the counter's own smoothing.
function Get-AmbientLoad {
  param([int] $Seconds)
  $samples = (Get-Counter -Counter '\Processor(_Total)\% Processor Time' `
                          -SampleInterval 1 -MaxSamples $Seconds).CounterSamples
  $values = $samples | ForEach-Object { $_.CookedValue }
  return ($values | Measure-Object -Average).Average
}

# Current RATE, not Get-Process's CPU property -- that is cumulative
# processor-seconds since the process started, so it ranks whatever has been
# open longest rather than whatever is burning cycles now. Percentages are of a
# single core, so a process pinning two cores reads ~200.
function Show-TopConsumers {
  Write-Log "  top CPU consumers right now (% of one core):"
  $samples = (Get-Counter -Counter '\Process(*)\% Processor Time' `
                          -SampleInterval 1 -MaxSamples 3 -ErrorAction SilentlyContinue).CounterSamples
  if (-not $samples) { Write-Log "    (per-process counters unavailable)"; return }
  $samples |
    Where-Object { $_.InstanceName -notin @('_total', 'idle') } |
    Group-Object InstanceName |
    ForEach-Object {
      [PSCustomObject]@{
        Name = $_.Name
        Pct  = ($_.Group | Measure-Object CookedValue -Average).Average
      }
    } |
    Where-Object { $_.Pct -ge 1 } |
    Sort-Object Pct -Descending | Select-Object -First 8 |
    ForEach-Object { Write-Log ("    {0,-24} {1,6:N1}%" -f $_.Name, $_.Pct) }
}

# --- Gate -------------------------------------------------------------------

Write-Log "=== quiet_bench start ==="
Write-Log "build dir     : $BuildDir"
Write-Log "targets       : $($Targets -join ', ')"
Write-Log "expect target : $(if ($ExpectTarget) { $ExpectTarget } else { '(report only)' })"
Write-Log "gate          : ambient < $MaxAmbient% over two consecutive ${SampleSeconds}s windows"

if (-not $ExpectTarget) {
  Write-Log "WARNING: no -ExpectTarget. Numbers will not be tier-asserted."
}

$deadline = (Get-Date).AddMinutes($MaxWaitMinutes)
$gateValue = $null
while ($true) {
  if ((Get-Date) -gt $deadline) {
    Write-Log "ABORT: machine did not go quiet within $MaxWaitMinutes minutes."
    Show-TopConsumers
    exit 1
  }

  $first = Get-AmbientLoad -Seconds $SampleSeconds
  Write-Log ("gate window 1: {0:N2}%" -f $first)
  if ($first -ge $MaxAmbient) {
    Show-TopConsumers
    continue
  }

  $second = Get-AmbientLoad -Seconds $SampleSeconds
  Write-Log ("gate window 2: {0:N2}%" -f $second)
  if ($second -ge $MaxAmbient) {
    Show-TopConsumers
    continue
  }

  $gateValue = [Math]::Max($first, $second)
  Write-Log ("GATE PASSED at {0:N2}%" -f $gateValue)
  break
}

# --- Run --------------------------------------------------------------------

if ($ExpectTarget) { $env:CORVUS_EXPECT_TARGET = $ExpectTarget }

$noise = New-Object System.Collections.Generic.List[double]
$failed = @()

foreach ($target in $Targets) {
  # The executable may sit at the build root (Ninja) or under a config
  # subdirectory (multi-config generators).
  $exe = Get-ChildItem -Path $BuildDir -Filter "$target.exe" -Recurse -File -ErrorAction SilentlyContinue |
         Select-Object -First 1
  if (-not $exe) {
    Write-Log "SKIP $target - executable not found under $BuildDir"
    $failed += "$target (missing)"
    continue
  }

  $before = Get-AmbientLoad -Seconds 3
  $noise.Add($before)
  Write-Log ("--- {0} --- noise before: {1:N2}%" -f $target, $before)

  $outFile = Join-Path $OutDir "quiet_bench_$target.txt"
  & $exe.FullName *> $outFile
  $code = $LASTEXITCODE

  $after = Get-AmbientLoad -Seconds 3
  $noise.Add($after)

  if ($code -ne 0) {
    # Exit 2 is the tier assertion in expect_target.h. Either way the numbers
    # in the file are not usable, and saying so here beats discovering it when
    # the log is read months later.
    Write-Log ("FAIL {0} exited {1} - output NOT publishable (tier mismatch is exit 2)" -f $target, $code)
    $failed += "$target (exit $code)"
  } else {
    Write-Log ("ok   {0} -> {1}   noise after: {2:N2}%" -f $target, (Split-Path $outFile -Leaf), $after)
  }
}

# --- Summary ----------------------------------------------------------------

$stats = $noise | Measure-Object -Average -Minimum -Maximum
Write-Log "=== quiet_bench summary ==="
Write-Log ("gated at      : {0:N2}%" -f $gateValue)
Write-Log ("noise held    : {0:N2}% - {1:N2}% (avg {2:N2}%)" -f $stats.Minimum, $stats.Maximum, $stats.Average)
Write-Log ("targets run   : {0} of {1}" -f ($Targets.Count - $failed.Count), $Targets.Count)

if ($failed.Count -gt 0) {
  Write-Log "FAILURES      : $($failed -join ', ')"
  Write-Log "Results are incomplete. Do not publish this pass."
  exit 1
}

Write-Log "Raw per-target output: $OutDir\quiet_bench_<target>.txt"
Write-Log "This log is the evidence chain - keep it with any published number."
exit 0
