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

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

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

inline uint64_t UlpDiff(double a, double b) {
  const bool na = std::isnan(a), nb = std::isnan(b);
  if (na && nb) return 0;
  if (na || nb) return UINT64_MAX;
  if (std::isinf(a) || std::isinf(b)) return SameBits(a, b) ? 0 : UINT64_MAX;
  return static_cast<uint64_t>(std::llabs(OrderedBits(a) - OrderedBits(b)));
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
struct RefRow {
  std::vector<std::string> tok;
  size_t line;
};

// Line-based reference loader. Whole-line tokenization with a per-row field
// count check replaces the old `while (f >> a >> b)` pattern, under which a
// truncated row silently shifted every following column by one. Blank lines
// and '#' comment lines are skipped (the dd reference files carry a comment
// header). Any malformed row -- wrong field count -- or a suspiciously small
// file is exit 2, never a silently shorter gate.
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
  while (std::getline(f, linebuf)) {
    ++lineno;
    RefRow row;
    row.line = lineno;
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
    if (row.tok.size() != ncols) {
      std::fprintf(stderr, "%s:%zu: expected %zu fields, found %zu\n", path,
                   lineno, ncols, row.tok.size());
      std::exit(2);
    }
    rows.push_back(std::move(row));
  }
  if (rows.size() < min_rows) {
    std::fprintf(stderr, "%s: suspiciously small: %zu rows (expected >= %zu)\n",
                 path, rows.size(), min_rows);
    std::exit(2);
  }
  return rows;
}

}  // namespace corvus_test

#endif  // CORVUS_TESTS_ULP_UTILS_H_
