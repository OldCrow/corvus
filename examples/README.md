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
not survive an operation that throws the accuracy away, and the operation that
does so is nearly always a subtraction of two nearly equal numbers.

The rule is quantitative, so you can predict it rather than discover it. If a
function returns `F` correct to about half an ulp, and you then compute `1 - F`
where `F` is close to 1, the absolute error stays around `ulp(1)/2 ≈ 1.1e-16`
while the answer itself shrinks. The relative error of the result is therefore
roughly

```
1.1e-16 / (1 - F)
```

which is to say: you lose about `log2(1 / (1 - F))` bits. At `1 - F = 1e-12`
that is forty bits — around four surviving digits — and at `1 - F = 1e-16` it
is all of them. `chi_squared_test` shows exactly this: at x = 50 the two routes
differ in the fifth significant digit, matching the prediction, while both
still look perfectly healthy on the page.

**This is why the library ships in pairs.** `erf`/`erfc`, `gamma_p`/`gamma_q`,
`beta_p`/`beta_q`, `i0`/`i0e` are not conveniences. Each pair exists so that
whichever side of the distribution you need, you can obtain it *directly*
instead of by subtracting the other from one. corvus's routing always evaluates
whichever of P/Q is the smaller and reports it at a relative bound, which is
what keeps a p-value of 1e-107 as trustworthy as one of 0.05 — but you have to
call the member of the pair that computes the quantity you actually want.

**Sometimes the remedy is a dedicated function rather than a pair.** `lbeta`
is not `lgamma(a) + lgamma(b) - lgamma(a+b)` evaluated carefully — the `a+b`
cancellation is removed analytically before any floating point happens, which
is why it is correctly rounded on every measured row while the assembly of
three correctly-rounded lgammas is not. `log_space_counting` puts the two side
by side: at a = 1e15 the assembly is subtracting numbers near 3.35e16, where
one ulp is 4, to produce an answer near 11.5.

Sometimes the cancellation is unavoidable, and then the honest move is to know
where the floor is rather than to iterate below it. `gamma_mle_fit` is that
case: the MLE condition *is* a difference of two nearly equal quantities, the
floor works out to a few parts in 1e13 at shape 300, and the example derives
that number and then stops there rather than reporting digits it does not have.

## Relative bounds, and where they stop applying

A related thing to know before reading `docs/ACCURACY.md`: some bounds are
stated as **absolute** rather than relative, and that is not hedging.

Where a function has a zero, no relative bound is possible — every nonzero
double is infinitely far from zero in relative terms. So corvus documents an
absolute bound in those bands and the relative bound elsewhere: `lgamma` on the
negative axis, `digamma` near its negative-axis zeros, `lbeta` around the zero
curve of ln B through (1,1). `log_space_counting` shows this directly —
`lbeta(1,1)` returns −8.1e-23 rather than exactly 0, which is comfortably
inside the documented `0.5·2⁻⁵³` and would be a meaningless number to hold to a
relative standard.

Reading a bound of that shape as a weakness gets it backwards. A library
claiming 1 ULP across a zero would be claiming something that cannot be
measured, never mind met.

So the short version, and the thing these examples are really teaching:

> A 1-ULP function does not give you a 1-ULP pipeline. Choose the formulation
> that computes your quantity directly; where you cannot, work out the floor.

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
