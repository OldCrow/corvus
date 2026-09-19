# Script-mode test of the #36 toolchain verdict: cmake -P, no compiler needed,
# so every CI leg exercises the Windows decisions. The probe that feeds
# WIDE_TARGETS needs a real compiler and is covered by the configure-only
# mingw steps in the Windows CI job.
include("${CMAKE_CURRENT_LIST_DIR}/../../cmake/ToolchainGuard.cmake")

set(failures 0)
function(expect want)
  corvus_toolchain_verdict(got ${ARGN})
  if(NOT got STREQUAL want)
    message(SEND_ERROR "verdict ${got}, expected ${want}, for: ${ARGN}")
    math(EXPR failures "${failures} + 1")
    set(failures ${failures} PARENT_SCOPE)
  endif()
endfunction()

# mingw GCC: rejected only while 256-bit+ targets are compiled in.
expect(MINGW_REJECT     SYSTEM_NAME Windows COMPILER_ID GNU WIDE_TARGETS TRUE)
expect(OK               SYSTEM_NAME Windows COMPILER_ID GNU WIDE_TARGETS FALSE)
expect(MINGW_OVERRIDDEN SYSTEM_NAME Windows COMPILER_ID GNU WIDE_TARGETS TRUE
                        ALLOW_UNSUPPORTED ON)
expect(OK               SYSTEM_NAME Windows COMPILER_ID GNU WIDE_TARGETS FALSE
                        ALLOW_UNSUPPORTED ON)
expect(MINGW_REJECT     SYSTEM_NAME CYGWIN  COMPILER_ID GNU WIDE_TARGETS TRUE)
expect(MINGW_REJECT     SYSTEM_NAME MSYS    COMPILER_ID GNU WIDE_TARGETS TRUE)
# A failed probe is not a cap: PROBE_FAILED is truthy and must reject.
expect(MINGW_REJECT     SYSTEM_NAME Windows COMPILER_ID GNU WIDE_TARGETS PROBE_FAILED)

# GCC off Windows is the primary supported toolchain, wide targets and all.
expect(OK SYSTEM_NAME Linux  COMPILER_ID GNU WIDE_TARGETS TRUE)
expect(OK SYSTEM_NAME Darwin COMPILER_ID GNU WIDE_TARGETS TRUE)

# Clang on Windows — clang-cl and the GNU driver both report "Clang" — is the
# supported path and must never be caught, whatever else is set.
expect(OK SYSTEM_NAME Windows COMPILER_ID Clang WIDE_TARGETS TRUE)
expect(OK SYSTEM_NAME Windows COMPILER_ID Clang WIDE_TARGETS TRUE MSVC_UNBLOCK ON)
expect(OK SYSTEM_NAME Darwin  COMPILER_ID AppleClang WIDE_TARGETS TRUE)

# Real MSVC: capped unless explicitly unblocked; the mingw override is unrelated.
expect(MSVC_CAP       SYSTEM_NAME Windows COMPILER_ID MSVC)
expect(MSVC_CAP       SYSTEM_NAME Windows COMPILER_ID MSVC ALLOW_UNSUPPORTED ON)
expect(MSVC_UNBLOCKED SYSTEM_NAME Windows COMPILER_ID MSVC MSVC_UNBLOCK ON)

if(failures)
  message(FATAL_ERROR "toolchain guard: ${failures} verdict(s) wrong")
endif()
message(STATUS "toolchain guard: all verdicts as expected")
