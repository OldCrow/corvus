# Configure-time guard for known-bad Windows toolchain/tier combinations (#36).
#
# Two outcomes used to configure and build cleanly while being wrong:
#
# - mingw-w64 GCC miscompiles by-value 256/512-bit vector arguments under the
#   Windows ABI (GCC PR 126741): the build is clean and the tests segfault.
#   GCC stays qualified for the 128-bit tiers, where the ABI's 16-byte stack
#   guarantee covers the temporaries. No fixed GCC release exists yet, so the
#   guard has no version condition; issue #29 is the trigger to add one.
# - Real MSVC (cl.exe) is capped at AVX2. Highway blocklists every AVX-512
#   target under it (HWY_BROKEN_MSVC, no compiler-version floor). corvus states
#   the same cap itself so that an upstream change to that blocklist cannot
#   silently widen an MSVC build nobody validated; CORVUS_MSVC_UNBLOCK_AVX512
#   lifts both.
#
# clang-cl and GNU-driver clang are affected by neither and must pass through
# untouched: the verdict keys on the compiler ID (GNU / MSVC), never on the
# MSVC or MINGW convenience variables, which are also true for clang flavours.
#
# The decision is a pure function of its arguments so that it runs in script
# mode (tests/cmake/toolchain_guard_test.cmake) on every CI leg; only the
# probe below needs a real compiler.

# Every target at and above a tier, in Highway's own idiom: x86 target bits
# descend with capability, so (T | (T - 1)) is T plus everything more advanced
# — including targets a future Highway adds. Same form as HWY_BROKEN_MSVC.
set(CORVUS_TARGETS_AVX2_AND_WIDER "(HWY_AVX2|(HWY_AVX2-1))")
set(CORVUS_TARGETS_AVX3_AND_WIDER "(HWY_AVX3|(HWY_AVX3-1))")

# corvus_toolchain_verdict(<out>
#   SYSTEM_NAME <CMAKE_SYSTEM_NAME>  COMPILER_ID <CMAKE_CXX_COMPILER_ID>
#   WIDE_TARGETS <bool>        256-bit+ x86 targets are compiled in (GNU only)
#   ALLOW_UNSUPPORTED <bool>   CORVUS_ALLOW_UNSUPPORTED_TOOLCHAIN
#   MSVC_UNBLOCK <bool>        CORVUS_MSVC_UNBLOCK_AVX512)
#
# <out> is one of:
#   OK               nothing to do
#   MINGW_REJECT     mingw GCC with 256-bit+ targets: fatal
#   MINGW_OVERRIDDEN the same, explicitly allowed: warn
#   MSVC_CAP         real MSVC: cap at AVX2 and say so
#   MSVC_UNBLOCKED   real MSVC with the unblock option: no cap
function(corvus_toolchain_verdict out)
  cmake_parse_arguments(PARSE_ARGV 1 A ""
    "SYSTEM_NAME;COMPILER_ID;WIDE_TARGETS;ALLOW_UNSUPPORTED;MSVC_UNBLOCK" "")

  set(verdict OK)
  if(A_COMPILER_ID STREQUAL "MSVC")
    if(A_MSVC_UNBLOCK)
      set(verdict MSVC_UNBLOCKED)
    else()
      set(verdict MSVC_CAP)
    endif()
  elseif(A_COMPILER_ID STREQUAL "GNU"
         AND A_SYSTEM_NAME MATCHES "^(Windows|CYGWIN|MSYS)$"
         AND A_WIDE_TARGETS)
    # CYGWIN/MSYS GCC use the same ms_abi argument passing; untested here, and
    # a false reject there is one cache variable away from a build.
    if(A_ALLOW_UNSUPPORTED)
      set(verdict MINGW_OVERRIDDEN)
    else()
      set(verdict MINGW_REJECT)
    endif()
  endif()
  set(${out} ${verdict} PARENT_SCOPE)
endfunction()

# corvus_probe_wide_x86_targets(<out> <disabled-targets-expr> <hwy-include-dirs>)
#
# Asks Highway itself whether any AVX2-or-wider target survives in HWY_TARGETS
# once <disabled-targets-expr> is applied — the compiled set is what PR 126741
# bites on, not the host CPU and not what dispatch would pick. Evaluating the
# real macro avoids parsing CORVUS_DISABLED_TARGETS, where a parsing slip would
# give false confidence in exactly the case the guard exists for.
#
# Two-sided: "capped" needs the no-wide-targets TU to compile, "wide" needs the
# opposite TU to compile. If neither does the probe itself is broken, and
# <out> is set to PROBE_FAILED rather than to either answer.
function(corvus_probe_wide_x86_targets out disabled_expr hwy_includes)
  set(prelude "")
  if(disabled_expr)
    set(prelude "#define HWY_DISABLED_TARGETS (${disabled_expr})\n")
  endif()
  string(APPEND prelude
    "#include \"hwy/detect_targets.h\"\n"
    "#define CORVUS_WIDE (HWY_ARCH_X86 && (HWY_TARGETS & ${CORVUS_TARGETS_AVX2_AND_WIDER}) != 0)\n")

  set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)  # compile only, no link
  try_compile(capped_compiles
    SOURCE_FROM_CONTENT corvus_probe_capped.cpp
      "${prelude}#if CORVUS_WIDE\n#error wide\n#endif\nint corvus_probe_capped;\n"
    CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${hwy_includes}"
    CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON
    NO_CACHE)
  try_compile(wide_compiles
    SOURCE_FROM_CONTENT corvus_probe_wide.cpp
      "${prelude}#if !CORVUS_WIDE\n#error capped\n#endif\nint corvus_probe_wide;\n"
    CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${hwy_includes}"
    CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON
    NO_CACHE)

  if(capped_compiles AND NOT wide_compiles)
    set(${out} FALSE PARENT_SCOPE)
  elseif(wide_compiles AND NOT capped_compiles)
    set(${out} TRUE PARENT_SCOPE)
  else()
    set(${out} PROBE_FAILED PARENT_SCOPE)
  endif()
endfunction()
