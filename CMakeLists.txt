cmake_minimum_required(VERSION 3.25)  # 3.25: FetchContent SYSTEM keyword
project(corvus VERSION 0.9.0 LANGUAGES CXX
  DESCRIPTION "SIMD-vectorized statistical special functions"
  HOMEPAGE_URL "https://github.com/OldCrow/corvus")

# The public header carries its own kVersion{Major,Minor,Patch} constants;
# they sat at 0.1.0 through the v0.2.0 release unnoticed, so a mismatch is
# now a configure-time error rather than a release-checklist item.
file(READ include/corvus/corvus.h _corvus_hdr)
foreach(_part Major Minor Patch)
  string(TOUPPER "${_part}" _PART)
  if(NOT _corvus_hdr MATCHES "kVersion${_part} = ([0-9]+);")
    message(FATAL_ERROR "corvus.h: kVersion${_part} not found")
  endif()
  if(NOT CMAKE_MATCH_1 STREQUAL PROJECT_VERSION_${_PART})
    message(FATAL_ERROR "corvus.h kVersion${_part} = ${CMAKE_MATCH_1} but "
      "project(VERSION) says ${PROJECT_VERSION_${_PART}} — update the header "
      "constants alongside the CMake version")
  endif()
endforeach()
unset(_corvus_hdr)

# --- Build type ---------------------------------------------------------
# Single-config generators default to Release: this is a numerical library
# and Debug-by-accident produces misleading benchmarks (house rule: perf
# numbers from optimized builds only). Intended roles:
#   Release        -O3, no debug info: benchmarks, accuracy validation, dist
#   RelWithDebInfo -O2 -g: profiling (Instruments/perf need symbols)
#   Debug          -O0 -g: debugger sessions; pair with CORVUS_SANITIZE
# Multi-config generators (Visual Studio, Xcode, Ninja Multi-Config) pick the
# config at build time and never create the CMAKE_BUILD_TYPE cache entry, so
# both the default and the STRINGS property are meaningless there — and
# set_property(CACHE) on a nonexistent entry is a hard configure error.
if(NOT CMAKE_CONFIGURATION_TYPES)
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  endif()
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY
    STRINGS Release RelWithDebInfo Debug)
endif()

# compile_commands.json for clangd/IDEs; harmless for consumers.
if(PROJECT_IS_TOP_LEVEL)
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif()

# --- Options (all CORVUS_-prefixed) -------------------------------------
option(CORVUS_BUILD_TESTS "Build corvus tests" ${PROJECT_IS_TOP_LEVEL})
option(CORVUS_BUILD_EXAMPLES "Build corvus examples" ${PROJECT_IS_TOP_LEVEL})
option(CORVUS_DEV_WARNINGS "Enable -Wall -Wextra -Wpedantic (dev builds)"
       ${PROJECT_IS_TOP_LEVEL})
option(CORVUS_WERROR "Treat warnings as errors (CI)" OFF)

# Compile-time target restriction, e.g. -DCORVUS_DISABLED_TARGETS="HWY_AVX3|HWY_AVX2"
# forces dispatch to stop at the next-lower tier. Used for native per-tier
# validation on a machine whose best tier is higher (Ryzen AVX-512 -> AVX2/SSE2).
set(CORVUS_DISABLED_TARGETS "" CACHE STRING
    "Pipe-separated HWY_* target macros to exclude from compilation/dispatch")

# Highway blocklists every AVX-512 target (HWY_AVX3 and the more advanced
# variants below its bit) under real MSVC: one unconditional entry in
# hwy/detect_targets.h citing a 2016 codegen bug, with no compiler-version
# floor — unlike every sibling entry in that file. The effect is silent, an
# MSVC build simply tops out at AVX2 with no diagnostic.
#
# Measured 2026-07-24 (Ryzen 7 7445HS, MSVC 19.51): unblocking dispatches to
# AVX3_ZEN4, passes every ULP gate with values bit-identical to GCC and
# clang-cl, and runs ~1.6x faster than the AVX2 ceiling. Default is still OFF
# because upstream declares the path broken and does not test it, and corvus
# does not claim what upstream will not stand behind. For AVX-512 coverage on
# Windows prefer clang-cl, which Highway does not blocklist.
option(CORVUS_MSVC_UNBLOCK_AVX512
       "MSVC only: override Highway's AVX-512 blocklist (unsupported upstream)"
       OFF)

# Sanitizers for Debug-type dev builds, e.g. -DCORVUS_SANITIZE=address;undefined
# (UBSan is cheap insurance for numerical kernel code; not for benchmarks).
set(CORVUS_SANITIZE "" CACHE STRING
    "Semicolon-separated sanitizers (address;undefined;thread)")

include(cmake/FindOrFetchHighway.cmake)

add_library(corvus STATIC
  src/bessel.cpp
  src/beta.cpp
  src/betainv.cpp
  src/digamma.cpp
  src/erf.cpp
  src/erf_data.cpp
  src/erfc.cpp
  src/erfinv.cpp
  src/exp.cpp
  src/exp_dd_data.cpp
  src/gamma.cpp
  src/gammainv.cpp
  src/lgamma.cpp
  src/log.cpp
  src/log_dd_data.cpp
  src/target.cpp
  src/trig.cpp
  src/trig_ph_data.cpp
  src/trigamma.cpp
)
add_library(corvus::corvus ALIAS corvus)

# PUBLIC: the header uses std::span, so consumers need C++20 too — the
# feature requirement travels through the installed package config.
target_compile_features(corvus PUBLIC cxx_std_20)
set_target_properties(corvus PROPERTIES
  CXX_EXTENSIONS OFF
  # PIC so the static lib can later be linked into shared objects
  # (pybind11 bindings per the pylibhmm/pylibstats pattern).
  POSITION_INDEPENDENT_CODE ON)

target_include_directories(corvus
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(corvus PRIVATE hwy::hwy)

# Dev warnings: PRIVATE only — never exported to consumers. CORVUS_WARN_FLAGS
# is also applied to test/bench targets in tests/CMakeLists.txt.
set(CORVUS_WARN_FLAGS "")
if(CORVUS_DEV_WARNINGS)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang"
     AND CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
    # clang-cl: -Wall is an ALIAS for MSVC /Wall, which clang-cl maps to
    # -Weverything (thousands of -W*-compat pedantry against our own C++20
    # baseline). /W4 is the spelling that maps to clang's -Wall -Wextra;
    # -Wpedantic passes through as the native clang flag.
    list(APPEND CORVUS_WARN_FLAGS /W4 -Wpedantic)
    if(CORVUS_WERROR)
      list(APPEND CORVUS_WARN_FLAGS /WX)
    endif()
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    list(APPEND CORVUS_WARN_FLAGS -Wall -Wextra -Wpedantic)
    if(CORVUS_WERROR)
      list(APPEND CORVUS_WARN_FLAGS -Werror)
    endif()
  elseif(MSVC)
    list(APPEND CORVUS_WARN_FLAGS /W4)
    if(CORVUS_WERROR)
      list(APPEND CORVUS_WARN_FLAGS /WX)
    endif()
  endif()
endif()
target_compile_options(corvus PRIVATE ${CORVUS_WARN_FLAGS})

# --- Floating-point contraction --------------------------------------------
# OFF, deliberately. The double-double layer (src/dd-inl.h) is a set of exact
# identities — TwoSum, Fast2Sum, Dekker — whose proofs assume each IEEE
# operation is performed and rounded AS WRITTEN. GCC defaults to
# -ffp-contract=fast, which fuses a Mul into an adjacent Add across statement
# boundaries; where that lands inside a TwoSum, s is no longer fl(a+b) and the
# "exact" residual is the error of an operation that never happened.
#
# Measured 2026-07-25 (Ryzen, log_dd reference set): GCC's default gave
# 2^-68.48 while MSVC and clang-cl gave 2^-67.88 for the same source at the
# same tier — contraction was making the result 0.6 bits BETTER by luck, and
# compiler-dependent either way. corvus documents cross-compiler
# reproducibility (docs/ACCURACY.md); leaving this to the optimizer's
# discretion contradicts it. erf and erfc are bit-for-bit unaffected because
# every fusion they want is already an explicit ops::MulAdd — which is the
# point: fusion should be requested in the source, not inferred.
#
# MSVC's default (/fp:precise) already does not contract, so it needs
# nothing — but ONLY from VS 2022 (cl 19.30): before 17.0, /fp:precise DID
# contract a*b+c to FMA whenever /arch:AVX2 or higher was on the command
# line, and MSVC has no contract-off switch to pin the older behavior. A
# VS 2019 build with /arch:AVX2 would silently fuse the dd layer, so the
# floor is fatal rather than advisory (#16; docs/ENVIRONMENT.md records it).
if(MSVC AND NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND MSVC_VERSION LESS 1930)
  message(FATAL_ERROR
    "corvus requires MSVC 19.30+ (Visual Studio 2022): older cl.exe "
    "contracts a*b+c to FMA under /arch:AVX2 with no off switch, which "
    "breaks the double-double layer's exactness proofs (MSVC_VERSION is "
    "${MSVC_VERSION}). Use VS 2022, clang-cl, or GCC/Clang.")
endif()
set(CORVUS_FP_FLAGS "")
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND MSVC)
  list(APPEND CORVUS_FP_FLAGS /clang:-ffp-contract=off)  # clang-cl driver
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU|IntelLLVM")
  list(APPEND CORVUS_FP_FLAGS -ffp-contract=off)
endif()
target_compile_options(corvus PRIVATE ${CORVUS_FP_FLAGS})

# Real MSVC only (clang-cl reports MSVC too but rejects /d2 flags): cl.exe's
# optimizer is superlinear in function size, and gamma.cpp — the heaviest
# per-target TU — stalled its back end past the CI job timeout twice
# (2026-07-29), on an unchanged runner image and compiler. The kernel already
# outlines its region cores and driver (HWY_NOINLINE, see gamma-inl.h); this
# flag is cl.exe's own escape hatch for the pathological cases that remain.
# It changes compile-time behavior only for functions the optimizer deems
# huge, and an unknown /d2 flag on some future MSVC would fail loudly, not
# silently. beta.cpp added 2026-08-07 per the trigger recorded at its
# design ("only if MSVC CI times degrade"): the Windows job went ~7 min →
# ~20–22 min when beta landed and hit its 25-min timeout on a slow runner.
# gammainv.cpp carries it from day one, per its detail design: the inverse
# instantiates all four forward region cores plus erfcinv and the dd
# transcendentals, twice per target (one per export), so it is the heaviest
# TU in the library by construction rather than by measurement.
# betainv.cpp likewise carries it from day one, and takes gammainv's title:
# it instantiates all four of BETA's region cores plus the shared cpsi
# machinery, gammainv's three seed functions, erfcinv, digamma, trigamma and
# the dd transcendentals — twice per target.
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  set_source_files_properties(src/gamma.cpp src/beta.cpp src/gammainv.cpp
    src/betainv.cpp
    PROPERTIES COMPILE_OPTIONS "/d2ReducedOptimizeHugeFunctions")
endif()

if(CORVUS_SANITIZE)
  # MSVC-ABI toolchains (cl AND clang-cl) get a hard configure error, not
  # GNU-style flags: cl emits D9002 and IGNORES -fsanitize, producing a
  # "sanitizer" build that sanitizes nothing (the silent pre-2026-08-12
  # behavior this guard replaces), and the /fsanitize=address runtime
  # plumbing on that ABI is unvalidated for corvus. Sanitizer builds run
  # on the Linux CI job (docs/ENVIRONMENT.md).
  if(MSVC)
    message(FATAL_ERROR
      "CORVUS_SANITIZE is supported on GNU-driver toolchains only "
      "(sanitizer builds are the Linux CI job's); the MSVC ABI would "
      "silently ignore -fsanitize=${CORVUS_SANITIZE}.")
  endif()
  list(JOIN CORVUS_SANITIZE "," _corvus_san)
  target_compile_options(corvus PRIVATE -fsanitize=${_corvus_san})
  target_link_options(corvus PUBLIC -fsanitize=${_corvus_san})
endif()

# Defines that decide which Highway targets a TU compiles. Collected in a
# variable because tests/ has a target that runs foreach_target over an
# internal kernel header (test_exp_dd) and must see the SAME target set — a
# capped library measured by an uncapped test is the failure mode
# CORVUS_EXPECT_TARGET exists to catch.
set(CORVUS_HWY_TARGET_DEFS "")
if(CORVUS_DISABLED_TARGETS)
  list(APPEND CORVUS_HWY_TARGET_DEFS
    "HWY_DISABLED_TARGETS=(${CORVUS_DISABLED_TARGETS})")
endif()

if(CORVUS_MSVC_UNBLOCK_AVX512)
  # CMake's MSVC is true for any cl-compatible driver, clang-cl included; only
  # real cl.exe is blocklisted, matching Highway's HWY_COMPILER_MSVC
  # (_MSC_VER && !__clang__).
  if(MSVC AND NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    # PRIVATE is sufficient: Highway derives the dispatch index from
    # HWY_CHOSEN_TARGET_MASK_TARGETS in the *calling* TU by design (see the
    # comment on ChosenTarget::GetIndex), so corvus's compiled target set may
    # legitimately be wider than that of a prebuilt libhwy. No hwy rebuild.
    list(APPEND CORVUS_HWY_TARGET_DEFS HWY_BROKEN_MSVC=0)
    message(WARNING
      "corvus: AVX-512 unblocked under MSVC (HWY_BROKEN_MSVC=0). Upstream "
      "Highway declares this path broken and does not test it — re-run the ULP "
      "gates before trusting numbers from this build.")
  else()
    message(WARNING
      "corvus: CORVUS_MSVC_UNBLOCK_AVX512 ignored — Highway's blocklist applies "
      "only to real MSVC (cl.exe); this compiler is ${CMAKE_CXX_COMPILER_ID}.")
  endif()
endif()

target_compile_definitions(corvus PRIVATE ${CORVUS_HWY_TARGET_DEFS})

if(CORVUS_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

if(CORVUS_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

# Install/export is only supported against a system (find_package) Highway.
# A FetchContent-built Highway is not in our export set; installing corvus
# would then produce a broken package. See PLAN.md open items.
if(CORVUS_HWY_PROVIDER STREQUAL "system")
  include(GNUInstallDirs)
  include(CMakePackageConfigHelpers)

  install(TARGETS corvus EXPORT corvus-targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
  install(DIRECTORY include/corvus DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  install(EXPORT corvus-targets
    NAMESPACE corvus::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/corvus)

  configure_package_config_file(
    cmake/corvus-config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/corvus-config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/corvus)
  write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/corvus-config-version.cmake
    COMPATIBILITY SameMajorVersion)
  install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/corvus-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/corvus-config-version.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/corvus)

  # pkg-config file. Requires: libhwy matches the module name Highway's
  # own install() emits (lib/pkgconfig/libhwy.pc) — verified by building
  # Highway 1.4.0 locally and inspecting its installed .pc file. It is
  # Requires, not Requires.private, on purpose: corvus is a static archive,
  # so a consumer's link line always needs Highway; Requires.private would
  # make plain `pkg-config --libs corvus` fail to link (only --static
  # would work). CI compiles and links a consumer with the plain flags.
  # Path from <prefix>/<libdir>/pkgconfig back to <prefix> (usually ../..;
  # handles multiarch libdirs too). Consumed by cmake/corvus.pc.in.
  file(RELATIVE_PATH CORVUS_PC_PREFIX_REL
       "${CMAKE_INSTALL_FULL_LIBDIR}/pkgconfig" "${CMAKE_INSTALL_PREFIX}")
  string(REGEX REPLACE "/$" "" CORVUS_PC_PREFIX_REL "${CORVUS_PC_PREFIX_REL}")
  configure_file(
    cmake/corvus.pc.in
    ${CMAKE_CURRENT_BINARY_DIR}/corvus.pc
    @ONLY)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/corvus.pc
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
else()
  message(STATUS "corvus: Highway was fetched, not found — 'install' target disabled")
endif()
