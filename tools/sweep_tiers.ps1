<#
.SYNOPSIS
    Per-tier native validation sweep (Windows), the PowerShell counterpart of
    the shell recipe in AGENTS.md.

.DESCRIPTION
    Caps compilation to one SIMD tier at a time with CORVUS_DISABLED_TARGETS,
    rebuilds, and runs every ULP gate with CORVUS_EXPECT_TARGET set, so a cap
    that fails to bite is a hard failure instead of a green suite measuring
    the wrong tier.

    Caps only ever REMOVE targets, so the sweep starts below the machine's
    native ceiling. On the Ryzen box, run the uncapped build separately for
    the AVX3* tiers -- and with GCC or clang-cl, never MSVC, whose Highway
    blocklist silently removes every AVX-512 target (see AGENTS.md).

    Aborts on the first configure, build, or gate failure. That matters: a
    build failure leaves the PREVIOUS tier's binaries in place, so continuing
    would re-measure the last tier under the new tier's name.

.PARAMETER Tier
    Restrict the sweep to one tier (AVX2, SSE4, SSSE3, SSE2). Default: all.

.PARAMETER BuildDir
    Build directory, reused across iterations so Highway is built once.

.EXAMPLE
    tools\sweep_tiers.ps1
    tools\sweep_tiers.ps1 -Tier SSE2
#>
[CmdletBinding()]
param(
  [ValidateSet("AVX2", "SSE4", "SSSE3", "SSE2")]
  [string[]] $Tier,
  [string] $BuildDir = "build-cap",
  [string] $CxxCompiler = "g++",
  # "cc" is a Unix convention that mingw-w64 does not ship; Highway's own
  # CMakeLists enables C, so an unresolvable name fails at configure.
  [string] $CCompiler = "gcc"
)

$ErrorActionPreference = "Stop"
$repo = Resolve-Path (Join-Path $PSScriptRoot "..")
$bd = if ([System.IO.Path]::IsPathRooted($BuildDir)) { $BuildDir }
      else { Join-Path $repo $BuildDir }

# Cap the WHOLE AVX-512 family including AVX10_2: leaving it out works only
# while HWY_BROKEN_AVX10_2's compiler-version gate holds, and that one expires.
$base = "HWY_AVX10_2|HWY_AVX3_SPR|HWY_AVX3_ZEN4|HWY_AVX3_DL|HWY_AVX3"
$caps = [ordered]@{
  "AVX2"  = $base
  "SSE4"  = "$base|HWY_AVX2"
  "SSSE3" = "$base|HWY_AVX2|HWY_SSE4"
  "SSE2"  = "$base|HWY_AVX2|HWY_SSE4|HWY_SSSE3"
}

# Dependency order, dd cores first: a foundation regression should fail
# under its own name, not a consumer's (same order as tests/CMakeLists.txt
# and the CI ULP report steps).
$gates = @(
  @{ exe = "test_exp_dd";   data = "exp_dd_reference.txt" },
  @{ exe = "test_log_dd";   data = "log_dd_reference.txt" },
  @{ exe = "test_dd_special"; data = "dd_special_reference.txt" },
  @{ exe = "test_erf_ulp";  data = "erf_reference.txt" },
  @{ exe = "test_erfc_ulp"; data = "erfc_reference.txt" },
  @{ exe = "test_lgamma_ulp"; data = "lgamma_reference.txt" },
  @{ exe = "test_erfinv_ulp"; data = "erfinv_reference.txt", "erfcinv_reference.txt" },
  @{ exe = "test_gamma_ulp"; data = "gamma_p_reference.txt", "gamma_q_reference.txt" },
  @{ exe = "test_beta_ulp"; data = "beta_p_reference.txt", "beta_q_reference.txt" },
  @{ exe = "test_digamma_ulp"; data = "digamma_reference.txt" },
  @{ exe = "test_trigamma_ulp"; data = "trigamma_reference.txt" },
  @{ exe = "test_gammainv_ulp"; data = "gammainv_p_reference.txt", "gammainv_q_reference.txt" },
  @{ exe = "test_betainv_ulp"; data = "betainv_p_reference.txt", "betainv_q_reference.txt" },
  @{ exe = "test_bessel_ulp"; data = "i0_reference.txt", "i1_reference.txt", "i0e_reference.txt", "i1e_reference.txt" }
)

$selected = if ($Tier) { $Tier } else { $caps.Keys }

foreach ($t in $selected) {
  Write-Host "`n===== tier $t =====" -ForegroundColor Cyan

  # Built as an array, not a backtick-continued command line: PowerShell does
  # NOT expand $variables on a continuation line in native-argument mode, so
  # the compiler would silently be configured as the literal "$CxxCompiler".
  $cfg = @(
    "-S", $repo, "-B", $bd, "-G", "Ninja",
    "-DCMAKE_BUILD_TYPE=Release",
    "-DCMAKE_C_COMPILER=$CCompiler",
    "-DCMAKE_CXX_COMPILER=$CxxCompiler",
    "-DCORVUS_DISABLED_TARGETS=$($caps[$t])"
  )
  cmake @cfg | Out-Null
  if ($LASTEXITCODE -ne 0) { throw "configure failed for tier $t" }

  cmake --build $bd
  if ($LASTEXITCODE -ne 0) { throw "build failed for tier $t" }

  $env:CORVUS_EXPECT_TARGET = $t
  try {
    foreach ($g in $gates) {
      # @(...) around the WHOLE pipeline, not just the input: PowerShell
      # unwraps a single-item pipeline result to a scalar STRING, and
      # splatting a string with @ below would then iterate its CHARACTERS
      # instead of passing it as one argument.
      $dataArgs = @(@($g.data) | ForEach-Object { Join-Path $repo "tests\data\$_" })
      & (Join-Path $bd "tests\$($g.exe).exe") @dataArgs
      if ($LASTEXITCODE -ne 0) { throw "$($g.exe) failed at tier $t (exit $LASTEXITCODE)" }
    }
  } finally {
    Remove-Item Env:\CORVUS_EXPECT_TARGET -ErrorAction SilentlyContinue
  }
}

Write-Host "`nall tiers passed" -ForegroundColor Green
