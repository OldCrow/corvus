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
| `gamma_mle_fit` | a batched Newton fit, and an accuracy floor that is arithmetic rather than a defect | `digamma`, `trigamma` |
| `students_t` | p-values and critical values, and significance far past anything a table lists | `beta_p`, `beta_p_inv` |
| `von_mises_density` | working past the point I₀ overflows, and a ratio whose scalings cancel exactly | `i0`, `i1`, `i0e`, `i1e` |
| `log_space_counting` | binomial coefficients and Beta normalisers, and a dedicated function beating a correct assembly | `lgamma`, `lbeta` |

More are being added; this table grows with them.

## Read this before concluding corvus is inaccurate

Every example above prints a place where the answer is worse than corvus's
documented bound. That is deliberate, and it is not corvus missing its bound.

**The bounds are on what corvus returns, not on what you do with it.** They do
not survive an operation that throws the accuracy away, and that operation is
nearly always a subtraction of two nearly equal numbers. Computing `1 - F` when
`F` is near 1 costs you about `log2(1 / (1 - F))` bits — forty of them by the
time `1 - F` is 1e-12, all of them by 1e-16.

That is why the library ships in pairs: `erf`/`erfc`, `gamma_p`/`gamma_q`,
`beta_p`/`beta_q`, `i0`/`i0e` exist so you can obtain whichever side you need
*directly*, rather than subtracting the other from one.

> A 1-ULP function does not give you a 1-ULP pipeline. Choose the formulation
> that computes your quantity directly; where you cannot, work out the floor.

**[docs/USER-GUIDE.md](../docs/USER-GUIDE.md) is the full treatment** — the
arithmetic behind that rule, rules of thumb for picking between each pair based
on your parameters, and why some bounds are stated as absolute rather than
relative. Read it before writing against the library; these examples are the
same lessons in runnable form.

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
  is what lets you keep them. This is the big one — see the section above.
- **The bounds are the product.** Where an example asserts something exact —
  `Phi(0) == 0.5`, a tail probability surviving where the naive route
  underflows to zero — that assertion is the point of the example, not
  decoration around it.
