#!/usr/bin/env python3
"""Verify the three places a corvus version lives agree.

Sources compared:
  1. kVersion{Major,Minor,Patch} in include/corvus/corvus.h
  2. project(VERSION) in CMakeLists.txt
  3. the release tag: --tag vX.Y.Z before tagging, or the exact
     vX.Y.Z tag on HEAD after (the default when --tag is omitted).

Exit 0 when all agree, 1 on any mismatch or missing source. The
release checklist (docs/RELEASING.md) runs this before and after
tagging; the post-tag run exists because two historical releases were
tagged while the header/CMake pair still read a stale number.
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def header_version() -> str:
    text = (ROOT / "include" / "corvus" / "corvus.h").read_text(
        encoding="utf-8"
    )
    parts = {}
    for name in ("Major", "Minor", "Patch"):
        m = re.search(
            rf"inline constexpr int kVersion{name} = (\d+);", text
        )
        if not m:
            sys.exit(f"FAIL: kVersion{name} not found in corvus.h")
        parts[name] = m.group(1)
    return f"{parts['Major']}.{parts['Minor']}.{parts['Patch']}"


def cmake_version() -> str:
    text = (ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
    m = re.search(r"project\(corvus VERSION (\d+\.\d+\.\d+)", text)
    if not m:
        sys.exit("FAIL: project(corvus VERSION ...) not found")
    return m.group(1)


def head_tag() -> str:
    r = subprocess.run(
        ["git", "describe", "--tags", "--exact-match", "HEAD"],
        cwd=ROOT, capture_output=True, text=True,
    )
    if r.returncode != 0:
        sys.exit(
            "FAIL: HEAD carries no tag. Pass --tag vX.Y.Z for the "
            "pre-tag check, or tag first for the post-tag check."
        )
    return r.stdout.strip()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--tag",
        help="intended release tag (vX.Y.Z); omit to check the "
        "exact tag on HEAD instead",
    )
    args = ap.parse_args()

    tag = args.tag if args.tag else head_tag()
    m = re.fullmatch(r"v(\d+\.\d+\.\d+)", tag)
    if not m:
        sys.exit(f"FAIL: tag '{tag}' is not of the form vX.Y.Z")
    tag_version = m.group(1)

    versions = {
        "corvus.h": header_version(),
        "CMakeLists.txt": cmake_version(),
        f"tag {tag}": tag_version,
    }
    for source, version in versions.items():
        print(f"  {source}: {version}")

    if len(set(versions.values())) != 1:
        print("FAIL: version sources disagree")
        return 1
    print(f"OK: all sources agree on {tag_version}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
