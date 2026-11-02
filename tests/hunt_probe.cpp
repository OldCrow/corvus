// Evaluator for tools/hunt (#82, T5). Not a test: reads lines
//   <fn> <arg1> [<arg2> [<arg3>]]
// on stdin (hex or decimal doubles) and writes one line per input,
//   <result bits as 16 hex digits>
// in the same order. Unknown functions and malformed lines write "?".
// Each line is one public span call of length 1, so the masked-tail path
// evaluates every row; lane independence is guarded by the smoke tests.
// A census build (CORVUS_ROUTE_CENSUS, #82 T7) appends the routes lane 0
// took, ` <route>[,<route>...]` or ` -`, for tools/hunt/route_census.py.
#include <corvus/corvus.h>

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <span>
#include <string_view>

#if defined(CORVUS_ROUTE_CENSUS)
#include "src/route_census.h"
#endif

namespace {

using Unary = void (*)(std::span<const double>, std::span<double>) noexcept;
using Binary = void (*)(std::span<const double>, std::span<const double>,
                        std::span<double>) noexcept;
using Ternary = void (*)(std::span<const double>, std::span<const double>,
                         std::span<const double>, std::span<double>) noexcept;

struct Entry {
  const char* name;
  Unary u;
  Binary b;
  Ternary t;
};

const Entry kTable[] = {
    {"erf", corvus::erf, nullptr, nullptr},
    {"erfc", corvus::erfc, nullptr, nullptr},
    {"erfinv", corvus::erfinv, nullptr, nullptr},
    {"erfcinv", corvus::erfcinv, nullptr, nullptr},
    {"lgamma", corvus::lgamma, nullptr, nullptr},
    {"digamma", corvus::digamma, nullptr, nullptr},
    {"trigamma", corvus::trigamma, nullptr, nullptr},
    {"i0", corvus::i0, nullptr, nullptr},
    {"i1", corvus::i1, nullptr, nullptr},
    {"i0e", corvus::i0e, nullptr, nullptr},
    {"i1e", corvus::i1e, nullptr, nullptr},
    {"exp", corvus::exp, nullptr, nullptr},
    {"log", corvus::log, nullptr, nullptr},
    {"log1p", corvus::log1p, nullptr, nullptr},
    {"expm1", corvus::expm1, nullptr, nullptr},
    {"log1pmx", corvus::log1pmx, nullptr, nullptr},
    {"cos", corvus::cos, nullptr, nullptr},
    {"sin", corvus::sin, nullptr, nullptr},
    {"gamma_p", nullptr, corvus::gamma_p, nullptr},
    {"gamma_q", nullptr, corvus::gamma_q, nullptr},
    {"gamma_p_inv", nullptr, corvus::gamma_p_inv, nullptr},
    {"gamma_q_inv", nullptr, corvus::gamma_q_inv, nullptr},
    {"lbeta", nullptr, corvus::lbeta, nullptr},
    {"beta_p", nullptr, nullptr, corvus::beta_p},
    {"beta_q", nullptr, nullptr, corvus::beta_q},
    {"beta_p_inv", nullptr, nullptr, corvus::beta_p_inv},
    {"beta_q_inv", nullptr, nullptr, corvus::beta_q_inv},
};

const Entry* Find(std::string_view name) {
  for (const auto& e : kTable) {
    if (name == e.name) return &e;
  }
  return nullptr;
}

}  // namespace

int main() {
  char line[512];
  while (std::fgets(line, sizeof line, stdin)) {
    char* p = line;
    while (*p == ' ') ++p;
    char* name_end = p;
    while (*name_end && *name_end != ' ' && *name_end != '\n') ++name_end;
    const Entry* e = Find(std::string_view(p, static_cast<size_t>(name_end - p)));
    double arg[3] = {0, 0, 0};
    int n = 0;
    char* q = name_end;
    while (n < 3) {
      char* end = nullptr;
      const double v = std::strtod(q, &end);
      if (end == q) break;
      arg[n++] = v;
      q = end;
    }
    const int want = e == nullptr ? -1 : (e->u ? 1 : e->b ? 2 : 3);
    if (n != want) {
      std::fputs("?\n", stdout);
      continue;
    }
    double r = 0;
    const std::span<double> out(&r, 1);
    const std::span<const double> a0(&arg[0], 1), a1(&arg[1], 1), a2(&arg[2], 1);
    if (e->u) e->u(a0, out);
    else if (e->b) e->b(a0, a1, out);
    else e->t(a0, a1, a2, out);
    uint64_t bits;
    std::memcpy(&bits, &r, sizeof bits);
#if defined(CORVUS_ROUTE_CENSUS)
    std::printf("%016llx ", static_cast<unsigned long long>(bits));
    const auto routes = corvus::census::Take();
    for (size_t i = 0; i < routes.size(); ++i) {
      std::printf("%s%s", i ? "," : "", routes[i]);
    }
    std::fputs(routes.empty() ? "-\n" : "\n", stdout);
#else
    std::printf("%016llx\n", static_cast<unsigned long long>(bits));
#endif
  }
  return 0;
}
