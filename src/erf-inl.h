// erf kernel assembly (#8: hoisted from the dispatch TU so it is
// reachable by a corvus_kernel_test_target() test and reusable from
// another TU, like every other family; the method itself lives in
// src/erf_core-inl.h). Per-target include guard (Highway -inl.h idiom).
#if defined(CORVUS_ERF_INL_H_) == defined(HWY_TARGET_TOGGLE)
#ifdef CORVUS_ERF_INL_H_
#undef CORVUS_ERF_INL_H_
#else
#define CORVUS_ERF_INL_H_
#endif

#include "src/erf_core-inl.h"
#include "src/erf_data.h"
#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {
namespace op = ops;

// See ErfTableCore (src/erf_core-inl.h) for the method. Max 1 ULP on
// validated tiers. |x| >= kErfAMax saturates to +/-1 (the bound is the
// smallest double whose erf rounds to 1, so this is correctly rounded);
// NaN is masked explicitly rather than relying on per-target Min/convert
// NaN semantics.
template <class D>
static HWY_INLINE op::V<D> ErfVec(D d, op::V<D> x) {
  const auto one = op::Set(d, 1.0);
  const auto amax = op::Set(d, detail::kErfAMax);

  const auto ax = op::Abs(x);
  const auto nan = op::IsNaN(x);
  auto ac = op::Min(ax, amax);
  ac = op::IfThenElse(nan, op::Zero(d), ac);  // safe table index for NaN lanes

  const auto core = ErfTableCore(d, ac);
  auto res = op::Add(core.e_hi, core.small);

  res = op::IfThenElse(op::Ge(ax, amax), one, res);
  res = op::IfThenElse(nan, x, res);      // propagate NaN (payload preserved)
  return op::CopySign(res, x);            // erf is odd; also -0 -> -0
}

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#endif  // include guard
