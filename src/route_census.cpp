// Route census record (#82 T7); compiled only with CORVUS_ROUTE_CENSUS.
#include "src/route_census.h"

#include <utility>
#include <vector>

namespace corvus::census {
namespace {
thread_local std::vector<const char*> g_hits;
thread_local int g_outside = 0;  // open scopes lane 0 is not in
}  // namespace

void Hit(const char* name, bool lane0) noexcept {
  if (!lane0 || g_outside != 0) return;
  if (!g_hits.empty() && g_hits.back() == name) return;
  try {
    g_hits.push_back(name);
  } catch (...) {  // a census build only; a lost record is not a result
  }
}

std::vector<const char*> Take() {
  std::vector<const char*> out;
  std::swap(out, g_hits);
  return out;
}

void Enter(bool lane0) noexcept {
  if (!lane0) ++g_outside;
}

void Leave(bool lane0) noexcept {
  if (!lane0) --g_outside;
}

}  // namespace corvus::census
