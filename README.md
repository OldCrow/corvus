# corvus

[![CI](https://github.com/OldCrow/corvus/actions/workflows/ci.yml/badge.svg)](https://github.com/OldCrow/corvus/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

SIMD-vectorized statistical special functions for C++20, with runtime
multi-target dispatch. Fills the gap between basic-transcendental SIMD
libraries (SLEEF, Highway's contrib math) and SciPy-level special-function
coverage: erf/erfc, lgamma, regularized incomplete gamma and beta, and their
inverses — the functions that gate vectorized statistical CDFs, quantiles,
and maximum-likelihood fitting.

**Status: early development.** `erf` and `erfc` are production-quality
clean-room kernels validated against an mpmath oracle on every SIMD tier
available across the development fleet — AVX-512 (`AVX3`, `AVX3_DL`,
`AVX3_ZEN4`), AVX2, SSE4, SSSE3, SSE2 and NEON, each on native silicon.
`erf`: max 1 ULP over the full domain. `erfc`: max 1 ULP for |x| <= 6 and
for subnormal results; max 5 ULP in the exp-bound tail (the backend `Exp`'s
contribution). API not yet stable.

## Design

- **Public API is std-only.** `std::span` in, `std::span` out. The SIMD
  backend (Google Highway) is an implementation detail hidden behind a
  ~20-op internal facade, sized so it can later be reimplemented on
  `std::simd` without touching kernel code.
- **Runtime dispatch.** One binary serves SSE2 through AVX-512 and NEON;
  the best available tier is selected at runtime.
- **Audited accuracy.** Every kernel documents its approximation source and
  accuracy bound; claims are made per SIMD tier only after validation on
  native silicon (not emulation).
- **Clean provenance.** Clean-room implementations only; MIT licensed.

## Build

```sh
cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

Uses an installed Highway if found, otherwise fetches a pinned copy at
configure time.

### Platforms and compilers

Built and tested in CI on Linux x86-64 (GCC), macOS arm64 (Apple Clang),
and Windows x86-64 (MSVC). Two Windows-specific points are worth knowing:

- **Without `-G Ninja` you get the Visual Studio generator**, which is
  multi-config: it ignores `CMAKE_BUILD_TYPE`, builds `Debug` by default,
  and needs `--config Release` to build and `-C Release` for `ctest` —
  without the latter, ctest runs no tests at all. Accuracy and performance
  claims only mean anything from an optimized build, so either pass
  `-G Ninja` as above or supply the config explicitly.
- **MSVC cannot reach AVX-512.** Highway places every AVX-512 target on its
  broken list under MSVC, so an MSVC build silently tops out at AVX2. It
  still passes every accuracy gate — the bounds hold on all tiers — but the
  widest vectors go unused. For AVX-512 on Windows, build with `clang-cl`
  (which keeps the MSVC ABI) or mingw-w64 GCC. `corvus::active_target()`
  reports the tier runtime dispatch actually selected, and is the only
  reliable way to know.

Accuracy is independent of optimization level and of compiler FP-contraction
settings: the kernels use explicit, capability-guarded FMA rather than
relying on the compiler to contract, so Debug and Release produce
bit-identical results.

> **Naming note:** Highway calls AVX-512 "AVX3" — so `HWY_AVX3`,
> `AVX3_DL`, `AVX3_ZEN4`, and `AVX3_SPR` in build output, target lists, and
> `CORVUS_DISABLED_TARGETS` values all refer to AVX-512 feature sets
> (baseline, VL/BW/DQ+VNNI, Zen 4, Sapphire Rapids), not some post-AVX2
> Intel extension of that name. `corvus::active_target()` reports these
> Highway names verbatim.

```cpp
#include <corvus/corvus.h>

std::vector<double> x = ..., y(x.size());
corvus::erf(x, y);
corvus::erfc(x, y);
```

Per-function methods, measured ULP bounds, and the validation matrix live
in [docs/ACCURACY.md](docs/ACCURACY.md).
