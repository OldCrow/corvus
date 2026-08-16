# corvus examples

Working programs that show what the library is *for*: building statistical
CDFs, quantiles, and fits out of the special functions corvus provides.

Each example verifies its own output and exits non-zero on failure, so they can
be run as smoke tests. They are deliberately **not** registered with ctest —
the gates under `tests/` are the accuracy contract, and mixing demonstrations
into that list would blur what a green `ctest` run means.

For the separate question of whether an *installed* corvus can be found and
linked from an outside project, see [`../consumer_example`](../consumer_example),
which exists to test the packaging rather than to teach the API.

## Building

Examples build by default when corvus is the top-level project:

```sh
cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
./build/examples/normal_distribution
```

Turn them off with `-DCORVUS_BUILD_EXAMPLES=OFF`. They are off automatically
when corvus is consumed via `add_subdirectory`.

> **On Windows**, the compiler that happens to be first on `PATH` decides which
> SIMD tier you get and, with one toolchain, whether the binary runs at all.
> See the platform notes in the top-level README: use `clang-cl` for AVX-512,
> and note that mingw-w64 GCC is unsafe above the 128-bit tiers.

## The examples

| example | shows | functions |
|---|---|---|
| `normal_distribution` | CDF, quantile, and why the far tail needs the other identity | `erfc`, `erfcinv` |
| `chi_squared_test` | critical values and p-values, and how subtraction goes wrong quietly before it goes wrong visibly | `gamma_p`, `gamma_q`, `gamma_p_inv` |

More are being added; this table grows with them.

## What they are trying to teach

Beyond "here is the call", each example is built around a decision a user
actually has to make:

- **Batch, don't loop.** Every example calls corvus once for a whole span. The
  kernels dispatch to the widest SIMD target the CPU supports and handle the
  ragged end with a masked tail. Calling them one element at a time is correct
  and pointless.
- **Pick the formulation that computes the small quantity.** `1 - Phi(z)` and
  `0.5*erfc(z/sqrt2)` are the same function on paper and completely different
  in floating point. corvus's bounds are relative, so choosing the right form
  is what lets you keep them.
- **The bounds are the product.** Where an example asserts something exact —
  `Phi(0) == 0.5`, a tail probability surviving where the naive route
  underflows to zero — that assertion is the point of the example, not
  decoration around it.
