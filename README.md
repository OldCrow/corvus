# corvus

[![CI](https://github.com/OldCrow/corvus/actions/workflows/ci.yml/badge.svg)](https://github.com/OldCrow/corvus/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

SIMD-vectorized statistical special functions for C++20, with runtime
multi-target dispatch. Fills the gap between basic-transcendental SIMD
libraries (SLEEF, Highway's contrib math) and SciPy-level special-function
coverage: erf/erfc, lgamma, regularized incomplete gamma and beta, and their
inverses — the functions that gate vectorized statistical CDFs, quantiles,
and maximum-likelihood fitting.

**Status: early development.** `erf`, `erfc`, `lgamma`, `erfinv`,
`erfcinv`, `gamma_p` and `gamma_q` are production-quality clean-room
kernels validated against an mpmath oracle on every SIMD tier available
across the development fleet — AVX-512 (`AVX3`, `AVX3_DL`, `AVX3_ZEN4`),
AVX2, SSE4, SSSE3, SSE2 and NEON, each on native silicon (the gamma pair's
NEON and AVX-512 validation is pending its first CI run and next Ryzen
session; see docs/ACCURACY.md). API not yet stable.

- `erf`: max 1 ULP over the full domain.
- `erfc`: max 1 ULP for |x| <= 6 and for subnormal results; max 2 ULP in
  the tail, where what remains is the tail polynomial's fit, not the
  exponential.
- `lgamma`: max 1 ULP across the positive axis — including arbitrarily
  close to the zeros at x = 1 and x = 2, which are exact — and correctly
  rounded throughout the Stirling region. On the negative axis the bound is
  1 ULP where |lgamma| >= 1 and 2^-53 absolute below that, because lgamma
  has infinitely many zeros there with no closed form.
- `erfinv` / `erfcinv`: max 1 ULP everywhere, including subnormal results
  down to the far tail (`erfcinv` reaches x up to ~27.2; `erfinv` never
  leaves x < 6). Useful directly as the normal quantile:
  `probit(p) = -sqrt(2)*erfcinv(2p)`.
- `gamma_p` / `gamma_q` (regularized incomplete gamma): max 2 ULP on the
  directly computed (smaller) side over the whole (a, x) plane, and every
  bound is relative — the routing always computes the smaller of P/Q
  directly, so tiny values keep full relative accuracy down to (and
  through) the subnormals, including Q for arbitrarily small a.

Both transcendental cores the kernels need (`exp_dd`, `log_dd`) are
corvus's own, so no accuracy-critical path depends on the backend's math
library.

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
corvus::lgamma(x, y);
corvus::erfinv(x, y);
corvus::erfcinv(x, y);
```

Per-function methods, measured ULP bounds, and the validation matrix live
in [docs/ACCURACY.md](docs/ACCURACY.md).
