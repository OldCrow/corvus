# corvus — User Guide

This guide is for people writing code against corvus. It covers what the
library gives you, what it deliberately does not, and the one thing that
catches nearly everybody: the difference between how accurate corvus's answers
are and how accurate *your* answer ends up being.

`docs/ACCURACY.md` is the reference for the measured bounds. This is the guide
for using them.

---

## 1. What corvus is

corvus computes about twenty **special functions** — the awkward mathematical
functions that statistical work is built out of. It computes them on whole
arrays at a time, using whatever vector instructions your CPU turns out to
have, and every function comes with a measured accuracy bound.

That is the whole library. It is not a statistics package. It gives you the
pieces that are hard to get right, and you assemble the distribution you
actually want out of them.

---

## 2. A few terms

Skip this if it is familiar. Everything below uses these words precisely.

**Double.** A 64-bit floating point number, about 15–17 significant decimal
digits. Every number corvus takes or returns is a double.

**ULP** — *unit in the last place.* Doubles are not evenly spaced: the gap
between neighbouring representable values grows with magnitude. One ULP is that
gap, wherever you happen to be. Near 1.0 it is about `2.2e-16`; near 1000 it is
about `1.1e-13`; near `1e16` it is `2` — by that magnitude, doubles cannot
represent odd numbers at all.

So "**max 1 ULP**" means the returned value is at most one neighbouring-double
step away from the true answer. That is very close to the best possible — you
cannot do better than landing on one of the two doubles nearest the truth.

**Correctly rounded.** Stronger than 1 ULP: the value returned is *the* nearest
double to the true answer, with no better choice available. Sometimes written
as 0.5 ULP. corvus achieves this for `lbeta` and in several regions elsewhere.

**Relative vs absolute error.** Relative error is the error as a fraction of
the answer — "correct to 15 digits" is a relative statement, and it holds
whether the answer is 1e-300 or 1e300. Absolute error is the plain difference,
regardless of size. Relative bounds are much stronger for small answers, which
is why corvus states them wherever it can. §8 covers where it cannot.

**Subnormal.** Doubles smaller than about `2.2e-308`, where the format starts
trading precision for range. corvus's bounds are stated to hold down through
this range for several functions, which is unusual and deliberate.

**Cancellation.** What happens when you subtract two nearly equal numbers: the
leading digits agree and destroy each other, and only the low-order digits —
the least reliable ones — survive into the answer. This is §6, and it is the
main thing this guide is about.

**Regularized incomplete gamma / beta.** "Incomplete" means the integral is
taken only part of the way, up to some `x`, instead of over the whole domain.
"Regularized" means it is then divided by the complete integral, so the result
runs from 0 to 1 and is directly a probability. `gamma_p(a, x)` is the fraction
of a Gamma(a) distribution lying below `x`; `gamma_q` is the fraction above.
The two always sum to 1.

**Lower and upper tail.** The lower tail at `x` is the probability of landing
below `x`; the upper tail is the probability of landing above. A p-value is
almost always an upper tail. The `_p` functions give lower tails, the `_q`
functions upper.

**SIMD, and tier.** Your CPU can apply one instruction to several numbers at
once. How many depends on which instruction set it supports — SSE2 does 2
doubles, AVX2 does 4, AVX-512 does 8. Those are the *tiers*. corvus picks the
best one available when your program starts; you do not configure anything.

---

## 3. What you get

Grouped by what you would reach for them to do.

**Normal distribution work**

| function | computes | typical bound |
|---|---|---|
| `erf`, `erfc` | the error function and its complement | 1 ULP (2 in erfc's far tail) |
| `erfinv`, `erfcinv` | their inverses | 1 ULP |

**Gamma, chi-squared, Poisson work**

| function | computes | typical bound |
|---|---|---|
| `gamma_p`, `gamma_q` | regularized incomplete gamma, lower and upper | 2 ULP |
| `gamma_p_inv`, `gamma_q_inv` | their inverses | 1 ULP |
| `lgamma` | log of the gamma function | 1 ULP |
| `digamma`, `trigamma` | first and second derivatives of `lgamma` | 1 ULP |

**Beta, Student's t, F, binomial work**

| function | computes | typical bound |
|---|---|---|
| `beta_p`, `beta_q` | regularized incomplete beta, lower and upper | 3 ULP |
| `beta_p_inv`, `beta_q_inv` | their inverses | 1 ULP |
| `lbeta` | log of the beta function | correctly rounded |

**Circular / von Mises work**

| function | computes | typical bound |
|---|---|---|
| `i0`, `i1` | modified Bessel functions, orders 0 and 1 | 1 ULP |
| `i0e`, `i1e` | the same, scaled by `exp(-|x|)` | 1 ULP |

Plus `active_target()`, which tells you which vector instruction set runtime
dispatch actually chose. That is more useful than it sounds — see §9.

The bounds above are the headline figures. Several functions have bands where
the guarantee changes shape; §8 explains why, and `docs/ACCURACY.md` gives the
exact statement for each.

---

## 4. What you don't get

Worth being blunt about, so you can plan around it.

- **No distributions.** There is no `normal_cdf`. You build it from `erfc`, in
  one line. The `examples/` directory shows how for the common ones.
- **No basic transcendentals.** No `exp`, `log`, `sin`, `pow`. Use `<cmath>`,
  or your SIMD library's own math functions. corvus covers the special
  functions specifically.
- **No random number generation, sampling, or fitting.** `digamma` and
  `trigamma` are what a maximum-likelihood fit needs, but you write the fit.
- **Real arguments only.** No complex arguments, no arbitrary precision, no
  `float` overloads — `double` throughout.
- **Not a SciPy port.** Names and conventions are similar in places, but
  nothing is guaranteed to match call for call.

---

## 5. Calling it

Every function takes input spans and an output span, all the same length:

```cpp
#include <corvus/corvus.h>

std::vector<double> x = { /* ... */ };
std::vector<double> y(x.size());

corvus::erfc(x, y);
```

Functions with more than one parameter take one span per parameter, and they
are evaluated elementwise — there is no broadcasting. To evaluate the
incomplete gamma at many `x` for a single shape `a`, you materialise `a` as an
array of equal values:

```cpp
std::vector<double> a(x.size(), 3.5);
corvus::gamma_p(a, x, out);
```

Things worth knowing:

- **Batch, don't loop.** One call for the whole array. The kernel handles the
  ragged end itself. Calling it per element is correct and pointless.
- **Writing into the input is fine**, as long as the spans are exactly the same
  (`in.data() == out.data()`). Partial overlap is undefined behaviour.
- **No allocation, no exceptions, thread-safe.** Safe to call from anywhere,
  including a hot loop or a worker pool. Every public function is declared
  `noexcept` (since v0.6.0), so the promise is compile-checked, not just
  documented.
- **One binary serves every CPU.** Dispatch happens on the first call.

---

## 6. The thing that will bite you

Here is the most important idea in this guide, stated plainly:

> **corvus's accuracy bounds describe the number it hands back to you. They say
> nothing about what happens to that number afterwards.**

A 1-ULP function does not give you a 1-ULP pipeline. If you take an excellent
answer and subtract it from another excellent answer, you can be left with
almost nothing, and no amount of care inside the library can prevent that.

### Why subtraction is the villain

Suppose you want the probability that a standard normal exceeds 10. You could
compute the CDF and subtract:

```cpp
double p = 1.0 - normal_cdf(10.0);     // don't
```

`normal_cdf(10.0)` is 0.9999999999999999999999924. Correct to the last bit a
double can hold — but doubles just below 1 are spaced about `1.1e-16` apart.
Everything that made this number interesting sits far below that spacing, so it
rounds to exactly `1.0`, and your answer is exactly `0`. The true value is
7.6e-24.

Computed the other way it is fine:

```cpp
double p = 0.5 * erfc(10.0 / std::sqrt(2.0));    // 7.6198530241605945e-24
```

Same function, same library, same input. The difference is entirely in which
expression you wrote.

### The rule, with numbers

When you compute `1 - F` and `F` is close to 1, the error stays around
`1.1e-16` while the answer shrinks. So the relative error of your result is
about

```
1.1e-16 / (1 - F)
```

Equivalently: you lose about **log₂(1 / (1 − F)) bits**.

| if `1 - F` is about | bits lost | digits you keep |
|---|---|---|
| 1e-3 | 10 | ~12 |
| 1e-8 | 27 | ~7 |
| 1e-12 | 40 | ~4 |
| 1e-16 | 53 | none |

More generally, whenever you subtract two numbers and the answer is much
smaller than the things you subtracted, you lose about `log₂(size of the
inputs / size of the answer)` bits.

### The dangerous part is the middle of that table

Total failure is easy to notice — a zero where you expected a small number
tends to get investigated. The rows around 1e-12 are the problem: you get a
plausible-looking answer with four good digits and twelve bad ones, and nothing
anywhere signals it.

The `chi_squared_test` example shows exactly this. At one point the two routes
differ in the fifth significant figure while both still look perfectly healthy
on the page. A p-value that collapses to zero announces itself; one that is
quietly wrong in the fourth digit gets published.

---

## 7. Choosing the right function

The library gives you pairs — `erf`/`erfc`, `gamma_p`/`gamma_q`,
`beta_p`/`beta_q`, `i0`/`i0e` — and this is why. Each pair exists so that
whichever quantity you need, you can get it **directly**, instead of computing
the other one and subtracting.

The general rule fits on one line:

> **Compute the small thing directly. Never subtract your way to it.**

Applying it means knowing roughly where your answer sits before you compute it.
That is not corvus being demanding; it is what working in floating point costs.
The rules of thumb below are meant to make it quick.

### Normal distribution

`Φ(z) = 0.5 * erfc(-z / √2)` handles **both** tails correctly, because the sign
of the argument flips with `z`. Just use it everywhere and don't think about
it.

For the upper tail, don't compute `1 - Φ(z)` — use `0.5 * erfc(z / √2)`.

For quantiles, `Φ⁻¹(p) = -√2 * erfcinv(2p)`.

### Incomplete gamma — `gamma_p` or `gamma_q`?

`Gamma(a)` has mean `a`, so:

| your situation | use |
|---|---|
| `x < a` — you are below the bulk | `gamma_p` |
| `x > a` — you are above the bulk | `gamma_q` |

For chi-squared with `k` degrees of freedom, `a = k/2` and `x = X²/2`, so the
comparison is just **X² against k**. A significance test is asking about the
upper tail, so it wants `gamma_q`.

### Incomplete beta — `beta_p` or `beta_q`?

`Beta(a, b)` has mean `a / (a + b)`:

| your situation | use |
|---|---|
| `x < a/(a+b)` | `beta_p` |
| `x > a/(a+b)` | `beta_q` |

For Student's t the substitution `x = ν/(ν + t²)` sends large `|t|` to small
`x`, which is convenient: the interesting tail is already the small side, and
`beta_p` computes it directly.

### The inverses — choose by the number you have

This one is different, and it is the case people get wrong most often.

**Pick the inverse that accepts the probability you already hold.** Not the one
you can convert your way into.

If you have a small upper-tail probability `q` and you want the value it
corresponds to, call `gamma_q_inv` / `beta_q_inv` / `erfcinv` with `q`
directly. Do not compute `1 - q` and pass that to the lower-tail inverse.

The reason is stark. At `q = 1e-30`, `1 - q` is **exactly 1.0** as a double.
Not approximately — exactly. The question you were asking no longer exists in
the number you passed. The `students_t` example asks for the t-statistic at a
one-sided p of 1e-100 and gets it right; there is no way to pose that question
through the other inverse at all.

### Bessel — `i0` or `i0e`?

| your situation | use |
|---|---|
| `x` might exceed ~700 | `i0e` / `i1e` — `i0` overflows at 713.99 |
| you want a ratio like `I₁/I₀` | `i1e/i0e` — the scalings cancel exactly |
| you want `log I₀(x)` | `log(i0e(x)) + x` |
| small `x`, and you want `I₀` itself | `i0` is fine |

A concentration of 700 is not exotic — for a von Mises distribution it is an
angular spread of about two degrees. Reach for the scaled forms by default in
that work.

### Log-gamma and log-beta

If you want `log B(a,b)`, call `lbeta`. Do not assemble it from three
`lgamma` calls. `lbeta` removes the cancellation as part of its derivation, not
afterwards, which is why it is correctly rounded where the assembly is not.

For binomial coefficients the same identity helps:

```
log C(n,k) = -log(n+1) - lbeta(k+1, n-k+1)
```

At `n = 1e15` the three-`lgamma` route is subtracting numbers near 3.35e16 to
produce an answer near 135, and it is wrong in the first decimal place. The
identity above never forms anything larger than 170.

---

## 8. Why some bounds are absolute instead of relative

Reading `docs/ACCURACY.md` you will find bounds stated two ways — some as a
number of ULP, some as an absolute figure like `2⁻⁵³`. That is not hedging.

A relative bound is meaningless where a function has a zero. Every nonzero
double is infinitely far from zero in relative terms, so "1 ULP" cannot be
stated, let alone met. In those bands corvus documents an absolute bound
instead:

- `lgamma` on the negative axis, which has infinitely many zeros
- `digamma` near its negative-axis zeros
- `lbeta` around the curve where `ln B` passes through zero, including the
  point (1, 1)

So `lbeta(1,1)` returns about `-8e-23` rather than exactly zero. That is inside
the documented `0.5 × 2⁻⁵³` and is the correct behaviour. A library claiming
1 ULP across a zero would be claiming something nobody can measure.

Practical consequence: if your work lives near one of those zeros, think in
absolute terms rather than relative ones, and expect the same from any library
you compare against.

---

## 9. Checking your own work

Three habits worth having.

**Estimate the size of your answer before you compute it.** If it is far
smaller than the intermediates you subtracted, count the bits you lost with the
rule in §6. If that number is near 53, you have nothing left, and you need a
different formulation rather than a more careful implementation.

**Where cancellation is genuinely unavoidable, find the floor and stop there.**
Some quantities really are small differences of large numbers — a maximum
likelihood condition, a log-probability assembled from log-factorials. You
cannot fix that by iterating harder. Work out the floor, report to it, and do
not present digits below it. The `gamma_mle_fit` example does this explicitly.

**Print `active_target()` when a result surprises you.** It tells you which
instruction set is actually running. Bounds are validated per tier, so knowing
the tier is the first step in reproducing anything — and on Windows the
compiler you built with silently decides the answer (see the platform notes in
the README).

---

## 10. Where to go next

- **`examples/`** — six worked programs, each verifying its own output. Start
  with `normal_distribution`, then `chi_squared_test`.
- **`docs/ACCURACY.md`** — the measured bound for every function, per SIMD tier,
  with the domains they apply on.
- **`docs/ARCHITECTURE.md`** — how the library is put together, if you are
  curious or intend to contribute.
- **`README.md`** — build instructions and platform notes. The Windows section
  is worth reading before trusting any number from a Windows build.
