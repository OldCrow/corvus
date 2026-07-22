# Consumer Example: find_package

Demonstrates consuming corvus from an external project after installation.

Only the `find_package` variant exists — the FetchContent path makes no
sense while `install` is gated on a system (find_package) Highway (see the
CMakeLists.txt comment above the install block, and PLAN.md).

## Prerequisites

A system Highway install (`find_package(hwy)` must succeed) and corvus
built and installed against it:

```bash
cd /path/to/corvus
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/path/to/hwy-prefix
cmake --build build --parallel
cmake --install build --prefix /tmp/corvus-install
```

## Build this example

```bash
cd consumer_example
cmake -S . -B build -DCMAKE_PREFIX_PATH="/tmp/corvus-install;/path/to/hwy-prefix"
cmake --build build
./build/consumer_demo
```

`CMAKE_PREFIX_PATH` needs both prefixes: corvus's own install, and the
Highway install `find_dependency(hwy)` resolves at `find_package(corvus)`
time (the exported `corvus-targets.cmake` references `hwy::hwy` via
`$<LINK_ONLY:...>`).

## What it tests

- `find_package(corvus REQUIRED)` locates the installed package
- `find_dependency(hwy CONFIG)` inside `corvus-config.cmake` resolves Highway
- `corvus::corvus` target provides headers and the static library, with the
  transitive `hwy::hwy` link resolving correctly
- `#include "corvus/corvus.h"` resolves correctly from the install tree
- `corvus::erf` computes correct values on a small span
