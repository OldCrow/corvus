# corvus

[![CI](https://github.com/OldCrow/corvus/actions/workflows/ci.yml/badge.svg)](https://github.com/OldCrow/corvus/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

SIMD-vectorized statistical special functions for C++20, with runtime
multi-target dispatch. Fills the gap between basic-transcendental SIMD
libraries (SLEEF, Highway's contrib math) and SciPy-level special-function
coverage: erf/erfc, lgamma, regularized incomplete gamma and beta, and their
inverses — the functions that gate vectorized statistical CDFs, quantiles,
and maximum-likelihood fitting.

**Status: early development.** `Erf` and `Erfc` are production-quality
clean-room kernels validated against an mpmath oracle on AVX2, SSE4, SSSE3,
and SSE2 (NEON and AVX-512 validation pending). `Erf`: max 1 ULP over the
full domain. `Erfc`: max 1 ULP for |x| <= 6 and for subnormal results;
max 5 ULP in the exp-bound tail (the backend `Exp`'s contribution). API not
yet stable.

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
cmake -B build && cmake --build build
ctest --test-dir build --output-on-failure
```

Uses an installed Highway if found, otherwise fetches a pinned copy at
configure time.

> **Naming note:** Highway calls AVX-512 "AVX3" — so `HWY_AVX3`,
> `AVX3_DL`, `AVX3_ZEN4`, and `AVX3_SPR` in build output, target lists, and
> `CORVUS_DISABLED_TARGETS` values all refer to AVX-512 feature sets
> (baseline, VL/BW/DQ+VNNI, Zen 4, Sapphire Rapids), not some post-AVX2
> Intel extension of that name. `corvus::ActiveTarget()` reports these
> Highway names verbatim.

```cpp
#include <corvus/corvus.h>

std::vector<double> x = ..., y(x.size());
corvus::erf(x, y);
corvus::erfc(x, y);
```

Per-function methods, measured ULP bounds, and the validation matrix live
in [docs/ACCURACY.md](docs/ACCURACY.md).
