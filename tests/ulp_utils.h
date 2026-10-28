// Shared ULP-measurement helpers and the strict reference-file loader for
// every test in this directory (#15). Before this header existed the three
// helpers below were hand-copied into 21 test files with four divergent
// NaN/Inf treatments, and every loader silently zero-filled malformed rows.
//
// THE ONE NaN/Inf POLICY (deliberate, applies to every gate):
//   * both NaN            -> 0 (equal). Payload identity is a SEPARATE,
//     explicit SameBits assertion where a contract requires it; a ULP
//     metric has no meaningful distance between two NaNs.
//   * exactly one NaN     -> UINT64_MAX (never inside any gate).
//   * any infinity        -> SameBits ? 0 : UINT64_MAX. An infinity in a
//     reference row is an exact saturation boundary, not a value with a
//     neighbourhood; matching sign and finiteness is all-or-nothing.
//   * finite vs finite    -> |OrderedBits(a) - OrderedBits(b)|. Note this
//     maps +0 and -0 to the SAME point (distance 0): a signed-zero
//     contract must be asserted with SameBits, never through UlpDiff --
//     no ULP gate can see a signed-zero regression (#14 N6).
// Tests that pre-filter non-finite rows at the call site keep working
// unchanged: the guards here only fire on rows such filters already
// removed.
#ifndef CORVUS_TESTS_ULP_UTILS_H_
#define CORVUS_TESTS_ULP_UTILS_H_

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include "corvus/corvus.h"

namespace corvus_test {

// Monotonic sign-magnitude mapping: adjacent doubles differ by 1 everywhere,
// including across +/-0 (which share one point -- see the policy above).
inline int64_t OrderedBits(double x) {
  int64_t b;
  std::memcpy(&b, &x, sizeof(b));
  return b < 0 ? (INT64_MIN - b) : b;
}

inline bool SameBits(double a, double b) {
  uint64_t ba, bb;
  std::memcpy(&ba, &a, sizeof(ba));
  std::memcpy(&bb, &b, sizeof(bb));
  return ba == bb;
}

// Every lane position (#79): each probe alone, then the whole probe list at
// offsets 0..8 behind filler drawn from the far end of the list. On every
// tier up to AVX-512 each probe then lands in every lane, next to a
// different region's point each time. A result that depends on its
// neighbours -- a vector-wide branch or early exit, a lane-crossing op --
// differs from the probe alone. `fn(cols, out)` makes one call over the K
// argument columns. Returns the number of drifted results; each is reported.
template <size_t K, class Fn>
int LaneMix(const char* name, const std::vector<std::array<double, K>>& probes,
            Fn fn) {
  const size_t n = probes.size();
  auto call = [&](size_t off, size_t first, size_t count) {
    std::array<std::vector<double>, K> cols;
    for (size_t k = 0; k < off; ++k) {
      for (size_t j = 0; j < K; ++j) cols[j].push_back(probes[(n - 1 - k) % n][j]);
    }
    for (size_t i = first; i < first + count; ++i) {
      for (size_t j = 0; j < K; ++j) cols[j].push_back(probes[i][j]);
    }
    std::vector<double> out(cols[0].size());
    fn(cols, out);
    return out;
  };
  std::vector<double> alone(n);
  for (size_t i = 0; i < n; ++i) alone[i] = call(0, i, 1)[0];
  int drift = 0;
  for (size_t off = 0; off < 9; ++off) {
    const std::vector<double> got = call(off, 0, n);
    for (size_t i = 0; i < n; ++i) {
      if (SameBits(got[off + i], alone[i])) continue;
      std::fprintf(stderr, "FAIL: %s lane-mix drift at (", name);
      for (size_t j = 0; j < K; ++j) {
        std::fprintf(stderr, "%s%a", j ? ", " : "", probes[i][j]);
      }
      std::fprintf(stderr, ") offset %zu: batched %a vs alone %a\n", off,
                   got[off + i], alone[i]);
      ++drift;
    }
  }
  return drift;
}

inline uint64_t UlpDiff(double a, double b) {
  const bool na = std::isnan(a), nb = std::isnan(b);
  if (na && nb) return 0;
  if (na || nb) return UINT64_MAX;
  if (std::isinf(a) || std::isinf(b)) return SameBits(a, b) ? 0 : UINT64_MAX;
  // Subtract in uint64 (#53): for opposite signs the int64 difference
  // overflows once |a|, |b| >= 2. Flipping the sign bit maps OrderedBits
  // onto uint64 in the same order, and the larger minus the smaller cannot
  // wrap.
  const uint64_t ua = static_cast<uint64_t>(OrderedBits(a)) ^ (1ULL << 63);
  const uint64_t ub = static_cast<uint64_t>(OrderedBits(b)) ^ (1ULL << 63);
  return ua > ub ? ua - ub : ub - ua;
}

// Strict double parse: the WHOLE token must be consumed. strtod without the
// endptr check turns a malformed token into 0.0 silently, which is exactly
// the failure a reference-driven gate cannot be allowed to absorb.
inline double ParseDouble(const std::string& tok, const char* path,
                          size_t line) {
  const char* s = tok.c_str();
  char* end = nullptr;
  const double v = std::strtod(s, &end);
  if (end == s || *end != '\0') {
    std::fprintf(stderr, "%s:%zu: malformed numeric token '%s'\n", path, line,
                 tok.c_str());
    std::exit(2);
  }
  return v;
}

// One reference row: its whitespace-separated tokens plus the 1-based source
// line, kept so a later ParseDouble failure can name where it came from.
// `res` holds the row's residuals (#82 T1) when the file carries them: the
// optional last token `r=<r0>[,<r1>...]`, one per output column, where
// r = (exact - stored CR value) / ulp(stored value's binade), |r| <= 1/2.
// tools/hunt/annotate.py writes it. Empty when the row has none.
// `stratum` names the `# stratum <name>` section the row sits in (#82 T3);
// rows before the first such line are stratum "base".
struct RefRow {
  std::vector<std::string> tok;
  size_t line;
  std::vector<double> res;
  std::string stratum;
};

// True when environment variable `name` is set to something other than
// empty or "0".
inline bool EnvFlag(const char* name) {
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4996)  // one read, see expect_target.h
#endif
  const char* v = std::getenv(name);
#if defined(_MSC_VER)
#pragma warning(pop)
#endif
  return v != nullptr && v[0] != '\0' && v[0] != '0';
}

// A file beside the reference file `ref_path` (tests/data/<name>).
inline std::string SiblingPath(const char* ref_path, const char* name) {
  std::string p(ref_path);
  const size_t slash = p.find_last_of("/\\");
  return (slash == std::string::npos ? std::string() : p.substr(0, slash + 1)) +
         name;
}

// Exact row counts (#82 T3): tests/data/ref_counts.txt holds
// `file|stratum|rows` for every reference file a gate loads. A file whose
// strata or counts differ -- a truncated file, a dropped or renamed stratum,
// a regenerated population of another size -- is exit 2, never a quietly
// different gate. tools/hunt/strata_append.py writes the entries.
inline void CheckCounts(const char* ref_path,
                        const std::vector<std::pair<std::string, size_t>>& got) {
  const std::string p = SiblingPath(ref_path, "ref_counts.txt");
  std::string base(ref_path);
  const size_t slash = base.find_last_of("/\\");
  if (slash != std::string::npos) base = base.substr(slash + 1);
  std::ifstream f(p);
  if (!f) {
    std::fprintf(stderr, "cannot open row-count file: %s\n", p.c_str());
    std::exit(2);
  }
  std::vector<std::pair<std::string, size_t>> want;
  std::string line;
  while (std::getline(f, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty() || line[0] == '#') continue;
    const size_t b1 = line.find('|');
    const size_t b2 = line.rfind('|');
    if (b1 == std::string::npos || b1 == b2) {
      std::fprintf(stderr, "%s: malformed line '%s'\n", p.c_str(), line.c_str());
      std::exit(2);
    }
    if (line.substr(0, b1) != base) continue;
    want.emplace_back(line.substr(b1 + 1, b2 - b1 - 1),
                      std::strtoull(line.c_str() + b2 + 1, nullptr, 10));
  }
  if (want != got) {
    std::fprintf(stderr, "%s: row counts differ from %s\n", base.c_str(), p.c_str());
    for (const auto& w : want) {
      std::fprintf(stderr, "  want %-40s %zu\n", w.first.c_str(), w.second);
    }
    for (const auto& g : got) {
      std::fprintf(stderr, "  got  %-40s %zu\n", g.first.c_str(), g.second);
    }
    std::exit(2);
  }
}

// Residual k of a row, NaN when the file is not annotated.
inline double Residual(const RefRow& row, size_t k) {
  return k < row.res.size() ? row.res[k] : std::nan("");
}

// Line-based reference loader. Whole-line tokenization with a per-row field
// count check replaces the old `while (f >> a >> b)` pattern, under which a
// truncated row silently shifted every following column by one. Blank lines
// and '#' comment lines are skipped (the dd reference files carry a comment
// header). Any malformed row -- wrong field count -- or a suspiciously small
// file is exit 2, never a silently shorter gate. A trailing `r=` token is
// split off into RefRow::res and does not count toward ncols.
// A `# stratum <name>` line starts a stratum section (#82 T3); the per-stratum
// counts must match ref_counts.txt exactly. Sections named `deferred:#NN`
// hold rows that break a published bound until issue #NN's fix lands: they
// are counted but not returned unless CORVUS_REF_DEFERRED=1, which is how
// the fix row shows its gate failing first. `deferred:?` (no issue named)
// is exit 2.
inline std::vector<RefRow> LoadRef(const char* path, size_t ncols,
                                   size_t min_rows) {
  std::ifstream f(path);
  if (!f) {
    std::fprintf(stderr, "cannot open reference file: %s\n", path);
    std::exit(2);
  }
  std::vector<RefRow> rows;
  std::string linebuf;
  size_t lineno = 0;
  std::string stratum = "base";
  std::vector<std::pair<std::string, size_t>> counts{{stratum, 0}};
  const bool with_deferred = EnvFlag("CORVUS_REF_DEFERRED");
  size_t deferred = 0;
  while (std::getline(f, linebuf)) {
    ++lineno;
    if (linebuf.rfind("# stratum ", 0) == 0) {
      stratum = linebuf.substr(10);
      while (!stratum.empty() && (stratum.back() == '\r' || stratum.back() == ' ')) {
        stratum.pop_back();
      }
      if (stratum == "deferred:?") {
        std::fprintf(stderr, "%s:%zu: deferred rows name no issue\n", path, lineno);
        std::exit(2);
      }
      counts.emplace_back(stratum, 0);
      continue;
    }
    RefRow row;
    row.line = lineno;
    row.stratum = stratum;
    std::string tok;
    for (size_t i = 0; i <= linebuf.size(); ++i) {
      const char c = i < linebuf.size() ? linebuf[i] : ' ';
      if (c == ' ' || c == '\t' || c == '\r') {
        if (!tok.empty()) {
          row.tok.push_back(tok);
          tok.clear();
        }
      } else {
        tok.push_back(c);
      }
    }
    if (row.tok.empty()) continue;           // blank line
    if (row.tok[0][0] == '#') continue;      // comment header
    if (row.tok.back().rfind("r=", 0) == 0) {
      std::string list = row.tok.back().substr(2) + ",";
      row.tok.pop_back();
      for (size_t at = 0, comma; (comma = list.find(',', at)) != std::string::npos;
           at = comma + 1) {
        row.res.push_back(ParseDouble(list.substr(at, comma - at), path, lineno));
      }
    }
    if (row.tok.size() != ncols) {
      std::fprintf(stderr, "%s:%zu: expected %zu fields, found %zu\n", path,
                   lineno, ncols, row.tok.size());
      std::exit(2);
    }
    ++counts.back().second;
    if (stratum.rfind("deferred:", 0) == 0) {
      ++deferred;
      if (!with_deferred) continue;
    }
    rows.push_back(std::move(row));
  }
  if (counts.front().second == 0) counts.erase(counts.begin());
  CheckCounts(path, counts);
  if (counts.size() > 1 || deferred) {
    std::printf("%s: %zu rows gated in %zu strata, %zu deferred%s\n", path,
                rows.size(), counts.size(), deferred,
                with_deferred ? " (included: CORVUS_REF_DEFERRED=1)" : "");
  }
  if (rows.size() < min_rows) {
    std::fprintf(stderr, "%s: suspiciously small: %zu rows (expected >= %zu)\n",
                 path, rows.size(), min_rows);
    std::exit(2);
  }
  return rows;
}

// ---- real error (#82 T1) ----------------------------------------------------
// An integer bound B fails only once the real error |out - exact| passes
// B + 1/2, which a finite reference set reaches rarely; the real error's
// maximum shows the erosion long before. With the residual r of the stored
// correctly rounded value cr, exact = cr + r * u, u = ulp(cr's binade), so
//   out - exact = (out - cr) / u - r     ulps of cr's binade,
// rescaled to ulps of exact's binade: exact sits one binade lower when cr
// is a power of two and r points toward zero (above the subnormal range).
// Returns NaN when the row carries no residual, +inf for a non-finite
// mismatch, 0 for a matching non-finite pair.
inline double RealErr(double out, double cr, double r) {
  if (std::isnan(r)) return r;
  if (!std::isfinite(cr) || !std::isfinite(out)) {
    return (SameBits(out, cr) || (std::isnan(out) && std::isnan(cr)))
               ? 0.0
               : HUGE_VAL;
  }
  int e;
  const double m = std::frexp(cr, &e);  // |m| in [1/2, 1), or 0
  const int ue = (cr == 0.0 ? -1074 : (e - 1 < -1022 ? -1022 : e - 1) - 52);
  // Power-of-two scaling is exact; the subtraction is exact for nearby out.
  const double d = std::ldexp(out, -ue) - std::ldexp(cr, -ue) - r;
  const bool lower = std::fabs(m) == 0.5 && e - 1 > -1022 && r * cr < 0.0;
  return std::fabs(d) * (lower ? 2.0 : 1.0);
}

// Running maximum of RealErr over a region; rows without a residual count
// as missing, which fails the region (a gate that opted in must see every
// row).
struct RealMax {
  double max = 0.0;
  size_t n = 0, missing = 0;
  size_t worst = 0;  // caller's row index of the maximum
  void Add(double e, size_t row) {
    if (std::isnan(e)) {
      ++missing;
      return;
    }
    ++n;
    if (e > max) {
      max = e;
      worst = row;
    }
  }
};

// Tier class for the pins: the 128-bit x86 tiers and the emulated ones have
// no FMA, so their dd products take the Dekker path with its own error.
inline const char* TierClass() {
  const char* t = corvus::active_target();
  static const char* const kNoFma[] = {"SSE2", "SSSE3", "SSE4", "EMU128",
                                       "SCALAR"};
  for (const char* n : kNoFma) {
    if (std::strcmp(t, n) == 0) return "nofma";
  }
  return "fma";
}

// The two real-error ceilings of a region (#82 T1), from tests/data/
// real_pins.txt, rows `gate|region|class|pin[|issue]`:
//   * contract ceiling: B + 1/2 - 1/4 for an integer bound B >= 1 (a quarter
//     ulp of safety below where the integer gate itself would fail); 1/2 for
//     a correctly rounded cell, where the integer gate is already exact.
//   * regression pin: the measured maximum over the tier class plus a small
//     margin (tools/hunt/pins.py). It catches drift long before the contract.
// A pin above the contract ceiling must name its issue: that cell is a known
// defect, and the fix row lowers the pin (the fail-first hook for 12d-12k).
// Every call prints `PIN gate|region|class|measured|contract`, so a local run
// or a CI log feeds pins.py. CORVUS_REAL_PIN_REPORT=1 (pins.py sets it) also
// stops a missing or exceeded pin from failing the gate.
struct Pin {
  double pin;
  std::string issue;
};

inline const std::vector<std::pair<std::string, Pin>>& Pins(
    const char* ref_path) {
  static std::vector<std::pair<std::string, Pin>> pins;
  static bool loaded = false;
  if (loaded) return pins;
  loaded = true;
  const std::string p = SiblingPath(ref_path, "real_pins.txt");
  std::ifstream f(p);
  if (!f) {
    std::fprintf(stderr, "cannot open pin file: %s\n", p.c_str());
    std::exit(2);
  }
  std::string line;
  size_t lineno = 0;
  while (std::getline(f, line)) {
    ++lineno;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty() || line[0] == '#') continue;
    std::vector<std::string> fld;
    for (size_t at = 0;;) {
      const size_t bar = line.find('|', at);
      fld.push_back(line.substr(at, bar - at));
      if (bar == std::string::npos) break;
      at = bar + 1;
    }
    if (fld.size() < 4 || fld.size() > 5) {
      std::fprintf(stderr, "%s:%zu: expected gate|region|class|pin[|issue]\n",
                   p.c_str(), lineno);
      std::exit(2);
    }
    pins.push_back({fld[0] + "|" + fld[1] + "|" + fld[2],
                    {ParseDouble(fld[3], p.c_str(), lineno),
                     fld.size() == 5 ? fld[4] : std::string()}});
  }
  return pins;
}

inline bool PinReport() { return EnvFlag("CORVUS_REAL_PIN_REPORT"); }

// Prints the region's real-error line and enforces both ceilings. Returns
// 0 on pass, 1 on failure. `ref_path` locates real_pins.txt beside it.
// ---- bound table (#82 T6) ---------------------------------------------------
// tests/data/bounds.txt, `gate|region|unit|gate|published|claim[|issue]`, is
// the one statement of every gated bound: the gates must pass exactly its
// gate field (CheckBound), and tools/claims_check.py holds the docs to its
// published field. A gate change that does not edit the table fails here.
struct BoundRow {
  std::string unit;
  double gate;
};

inline const std::vector<std::pair<std::string, BoundRow>>& Bounds(
    const char* ref_path) {
  static std::vector<std::pair<std::string, BoundRow>> rows;
  static bool loaded = false;
  if (loaded) return rows;
  loaded = true;
  const std::string p = SiblingPath(ref_path, "bounds.txt");
  std::ifstream f(p);
  if (!f) {
    std::fprintf(stderr, "cannot open bound table: %s\n", p.c_str());
    std::exit(2);
  }
  std::string line;
  size_t lineno = 0;
  while (std::getline(f, line)) {
    ++lineno;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty() || line[0] == '#') continue;
    std::vector<std::string> fld;
    for (size_t at = 0;;) {
      const size_t bar = line.find('|', at);
      fld.push_back(line.substr(at, bar - at));
      if (bar == std::string::npos) break;
      at = bar + 1;
    }
    if (fld.size() < 6 || fld.size() > 7) {
      std::fprintf(stderr, "%s:%zu: expected gate|region|unit|gate|published|claim[|issue]\n",
                   p.c_str(), lineno);
      std::exit(2);
    }
    rows.emplace_back(fld[0] + "|" + fld[1],
                      BoundRow{fld[2], ParseDouble(fld[3], p.c_str(), lineno)});
  }
  return rows;
}

// 0 when (gate, region) is in the bound table with this unit and gate value.
inline int CheckBound(const char* ref_path, const char* gate, const char* region,
                      const char* unit, double bound) {
  const std::string key = std::string(gate) + "|" + region;
  std::printf("BOUND %s|%s|%g\n", key.c_str(), unit, bound);
  for (const auto& kv : Bounds(ref_path)) {
    if (kv.first != key) continue;
    if (kv.second.unit == unit && kv.second.gate == bound) return 0;
    std::fprintf(stderr,
                 "FAIL: %s: the gate enforces %g %s, tests/data/bounds.txt says "
                 "%g %s\n",
                 key.c_str(), bound, unit, kv.second.gate, kv.second.unit.c_str());
    return 1;
  }
  std::fprintf(stderr, "FAIL: %s: no row in tests/data/bounds.txt\n", key.c_str());
  return 1;
}

inline int CheckReal(const char* ref_path, const char* gate, const char* region,
                     uint64_t bound, const RealMax& m) {
  const int bound_rc =
      CheckBound(ref_path, gate, region, "ulp", static_cast<double>(bound));
  const char* cls = TierClass();
  if (std::strchr(gate, '|') != nullptr || std::strchr(region, '|') != nullptr) {
    std::fprintf(stderr, "%s / %s: '|' separates pin-file fields\n", gate,
                 region);
    std::exit(2);
  }
  const std::string key = std::string(gate) + "|" + region + "|" + cls;
  const double contract =
      bound == 0 ? 0.5 : static_cast<double>(bound) + 0.25;
  const bool report = PinReport();
  std::printf("PIN %s|%.6f|%.2f\n", key.c_str(), m.max, contract);
  const Pin* pin = nullptr;
  for (const auto& kv : Pins(ref_path)) {
    if (kv.first == key) pin = &kv.second;
  }
  char pin_s[32] = "none";
  if (pin) std::snprintf(pin_s, sizeof(pin_s), "%.2f", pin->pin);
  std::printf("  real %-32s max=%.4f  contract %.2f  pin %s%s%s\n", region,
              m.max, contract, pin_s,
              pin && !pin->issue.empty() ? "  known: " : "",
              pin ? pin->issue.c_str() : "");
  int rc = bound_rc;
  if (m.missing) {
    std::fprintf(stderr, "FAIL: %s %s: %zu rows carry no residual (r=)\n",
                 gate, region, m.missing);
    rc = 1;
  }
  if (report) return rc;
  if (pin == nullptr) {
    std::fprintf(stderr, "FAIL: %s: no real-error pin in real_pins.txt\n",
                 key.c_str());
    return 1;
  }
  // pins.py writes "?" where a person still has to name the issue.
  const bool named = !pin->issue.empty() && pin->issue != "?";
  if (pin->pin > contract && !named) {
    std::fprintf(stderr,
                 "FAIL: %s: pin %.4f is above the contract ceiling %.2f and "
                 "names no issue\n",
                 key.c_str(), pin->pin, contract);
    rc = 1;
  }
  if (m.max > pin->pin) {
    std::fprintf(stderr, "FAIL: %s: real error %.4f exceeds pin %.4f\n",
                 key.c_str(), m.max, pin->pin);
    rc = 1;
  }
  if (m.max > contract && !named) {
    std::fprintf(stderr,
                 "FAIL: %s: real error %.4f exceeds contract ceiling %.2f\n",
                 key.c_str(), m.max, contract);
    rc = 1;
  }
  return rc;
}

}  // namespace corvus_test

#endif  // CORVUS_TESTS_ULP_UTILS_H_
