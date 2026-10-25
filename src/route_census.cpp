// Route census record (#82 T7); compiled only with CORVUS_ROUTE_CENSUS.
#include "src/route_census.h"

#include <set>
#include <string>
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

void Floor(const char* file, unsigned line, bool lane0) noexcept {
  if (!lane0 || g_outside != 0) return;
  try {
    // Interned, so the record holds stable pointers like the route names.
    static std::set<std::string> names;
    std::string f(file);
    const size_t slash = f.find_last_of("/\\");
    if (slash != std::string::npos) f = f.substr(slash + 1);
    const auto it = names.insert("floor:" + f + ":" + std::to_string(line)).first;
    Hit(it->c_str(), true);
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
