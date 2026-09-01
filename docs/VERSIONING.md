# Versioning and Stability Policy

corvus follows semantic versioning: MAJOR.MINOR.PATCH. This document
states exactly what each component promises, because a numerics
library has two stability surfaces — the API and the accuracy bounds —
and standard SemVer only speaks to the first.

## What 1.0 promises

**API stability.** Within a major version, every public entry point in
`corvus.h` keeps its name, signature, parameter order, `noexcept`
status, and documented special-value behavior. New functions may
appear (minor); existing ones do not change or disappear.

**Accuracy bounds are contractual floors.** The per-tier bounds in
[ACCURACY.md](ACCURACY.md) are part of the interface. Within a major
version a released bound never gets worse: a function documented at
1 ULP on a validated tier stays at or under 1 ULP on that tier in
every later 1.x release. Bounds may tighten at any time.

**Special-value behavior is contractual.** The documented handling of
NaN (propagation, payload), infinities, signed zero, and domain edges
is part of the interface, and the smoke tests pin it.

## What is deliberately NOT promised

**Bit-for-bit reproducibility across versions.** Two releases may
return different correctly-behaving results for the same input, as
long as both sit inside the documented bound. A kernel improvement
that moves results by a fraction of an ulp is a patch, not a break.
Consumers that need run-to-run reproducibility should pin a corvus
version, a compiler, and a SIMD tier together.

**ABI.** corvus ships as source; build it with your toolchain. No
binary compatibility is promised between releases, compilers, or
build configurations.

**Performance.** [PERFORMANCE.md](PERFORMANCE.md) reports
measurements, not commitments.

**Internal headers.** Everything under `src/` is implementation.
Only `include/corvus/corvus.h` is public.

## What each version component means

| Change | Component |
|---|---|
| Remove/rename a function, change a signature or parameter order | MAJOR |
| Loosen a published accuracy bound on a validated tier | MAJOR |
| Change documented special-value behavior | MAJOR |
| Raise the required C++ standard | MAJOR |
| Add a function or overload | MINOR |
| Tighten a published bound; validate a new tier or platform | MINOR |
| Fix a result that violated its documented bound | PATCH |
| Result movement within documented bounds (kernel improvements) | PATCH |
| Documentation, build system, CI, test changes | PATCH |

Boundary case: fixing a function whose actual behavior violated its
documentation is a PATCH even when results change visibly — the
documentation was the contract, and the fix restores it.

## Per-tier claims

Accuracy bounds are stated per SIMD tier and become claims only after
validation on native silicon. A tier the validation matrix in
ACCURACY.md marks open carries no bound for that cell; treating the
bound of a validated tier as if it covered an unvalidated one is a
consumer error, not a corvus regression.

## Version identity

The version lives in two places that a configure-time check forces to
agree: `corvus::kVersion{Major,Minor,Patch}` in the public header and
`project(VERSION)` in CMakeLists.txt. Release tags are `vX.Y.Z` and
must match both; `tools/check_release_version.py` verifies all three
and the release checklist ([RELEASING.md](RELEASING.md)) runs it
before and after tagging.
