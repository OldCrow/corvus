#include "corvus/corvus.h"
#include "src/erf_data.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "src/erf.cpp"
#include "hwy/foreach_target.h"

#include "src/ops-inl.h"

HWY_BEFORE_NAMESPACE();
namespace corvus {
namespace HWY_NAMESPACE {

namespace op = ops;

// Table + local Taylor correction, ported from libstats vector_erf_neon
// (independently derived; see tools/gen_erf_table.py header and libstats
// docs/NEON_ERF_DERIVATION.md). Max 1 ULP on validated tiers.
//
// For a = |x| clamped to the grid and r = round(a*256)/256 (exact: 256 is a
// power of two, round-to-nearest-even), with residual d = a - r (exact by
// Sterbenz), the table supplies E = erf(r) compensated as E_hi + E_lo and
// S = erf'(r); the kernel evaluates
//   erf(a) = E + S*(d + c2 d^2 + c3 d^3 + c4 d^4 + c5 d^5)
// with c_k derived at run time from w = r^2:
//   c2 = -r, c3 = (2w-1)/3, c4 = r(3-2w)/6, c5 = (4w^2-12w+3)/30.
// |x| >= kErfAMax saturates to +/-1 (bound chosen so this is correctly
// rounded); NaN is masked explicitly rather than relying on per-target
// Min/convert NaN semantics (unlike the NEON original, this kernel must be
// portable across targets where those semantics differ).
template <class D>
HWY_INLINE op::V<D> ErfVec(D d, op::V<D> x) {
  const op::SignedTag<D> di;
  const auto one = op::Set(d, 1.0);
  const auto amax = op::Set(d, detail::kErfAMax);

  const auto ax = op::Abs(x);
  const auto nan = op::IsNaN(x);
  auto ac = op::Min(ax, amax);
  ac = op::IfThenElse(nan, op::Zero(d), ac);  // safe table index for NaN lanes

  const auto rf = op::Round(op::Mul(ac, op::Set(d, 256.0)));  // integral, <= 1536
  const auto ji = op::ConvertToInt(di, rf);                   // exact truncation
  const auto j3 = op::Add(op::ShiftLeft<1>(ji), ji);          // 3*j: field stride

  const auto E  = op::GatherIndex(d, detail::kErfTable + 0, j3);
  const auto S  = op::GatherIndex(d, detail::kErfTable + 1, j3);
  const auto El = op::GatherIndex(d, detail::kErfTable + 2, j3);

  const auto r = op::Mul(rf, op::Set(d, 1.0 / 256.0));  // exact
  const auto dd = op::Sub(ac, r);                       // exact residual
  const auto w = op::Mul(r, r);

  const auto c2 = op::Neg(r);
  const auto c3 = op::MulAdd(w, op::Set(d, 2.0 / 3.0), op::Set(d, -1.0 / 3.0));
  const auto c4 = op::Mul(r, op::MulAdd(w, op::Set(d, -1.0 / 3.0), op::Set(d, 0.5)));
  const auto c5 = op::MulAdd(
      w, op::MulAdd(w, op::Set(d, 2.0 / 15.0), op::Set(d, -2.0 / 5.0)),
      op::Set(d, 0.1));

  // t = c2 + d(c3 + d(c4 + d*c5)); p = d + d^2*t; erf = E_hi + (E_lo + S*p)
  auto t = op::MulAdd(dd, c5, c4);
  t = op::MulAdd(dd, t, c3);
  t = op::MulAdd(dd, t, c2);
  const auto p = op::MulAdd(op::Mul(dd, dd), t, dd);
  auto res = op::Add(E, op::MulAdd(S, p, El));

  res = op::IfThenElse(op::Ge(ax, amax), one, res);
  res = op::IfThenElse(nan, x, res);      // propagate NaN (payload preserved)
  return op::CopySign(res, x);            // erf is odd; also -0 -> -0
}

void ErfImpl(const double* in, double* out, size_t n) {
  const op::ScalableTag<double> d;
  const size_t N = op::Lanes(d);
  size_t i = 0;
  for (; i + N <= n; i += N) {
    op::Store(ErfVec(d, op::Load(d, in + i)), d, out + i);
  }
  if (i < n) {
    // Masked tail: same code path as the full lanes, no scalar fallback.
    op::StoreN(ErfVec(d, op::LoadN(d, in + i, n - i)), d, out + i, n - i);
  }
}

const char* TargetNameImpl() { return hwy::TargetName(HWY_TARGET); }

}  // namespace HWY_NAMESPACE
}  // namespace corvus
HWY_AFTER_NAMESPACE();

#if HWY_ONCE
namespace corvus {

HWY_EXPORT(ErfImpl);
HWY_EXPORT(TargetNameImpl);

void Erf(std::span<const double> in, std::span<double> out) {
  HWY_DYNAMIC_DISPATCH(ErfImpl)(in.data(), out.data(), in.size());
}

const char* ActiveTarget() {
  return HWY_DYNAMIC_DISPATCH(TargetNameImpl)();
}

}  // namespace corvus
#endif  // HWY_ONCE
