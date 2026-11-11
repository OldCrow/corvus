"""Check an extracted release tarball for the two things a build cannot see.

A release tarball is `git archive` output: the tree minus everything
`.gitattributes` marks export-ignore (docs/RELEASING.md step 14). Given the
directory such an archive was extracted into, this script fails if

1. a repo-only path is present, which means export-ignore was not applied
   and whatever is built from that directory is not the release tree; or
2. a shipped Markdown file links to a relative path the tarball lacks.
   Shipped docs may *name* repo-only files as provenance, but a link that
   resolves to nothing is a defect.

Usage:
    git archive HEAD | tar -x -C <dir>
    python3 tools/check_release_tarball.py <dir>

CI runs this, then configures, builds and tests <dir> (ci.yml, Linux job).
Standard library only.
"""

import re
import sys
from pathlib import Path

# One representative per export-ignore group in .gitattributes. Their absence
# shows the filter ran; the full list lives in .gitattributes, not here.
REPO_ONLY = (
    "AGENTS.md",
    "PLAN.md",
    ".gitattributes",
    ".github",
    ".clang-format",
    ".pre-commit-config.yaml",
    "tools",
    "docs/RELEASING.md",
    "docs/ENVIRONMENT.md",
    "docs/bench-evidence",
)

# [text](target) with an optional #fragment; inline code and URLs are skipped.
LINK = re.compile(r"\]\(([^)#\s]+)(?:#[^)]*)?\)")
EXTERNAL = re.compile(r"^(?:[a-z][a-z0-9+.-]*:|//)", re.IGNORECASE)


def main(argv: list[str]) -> int:
    if len(argv) != 1:
        print(__doc__)
        return 2
    root = Path(argv[0]).resolve()
    if not root.is_dir():
        print(f"error: {root} is not a directory")
        return 2

    problems = 0

    for name in REPO_ONLY:
        if (root / name).exists():
            problems += 1
            print(f"repo-only path is in the tarball: {name}")

    for doc in sorted(root.rglob("*.md")):
        text = doc.read_text(encoding="utf-8", errors="replace")
        for line_no, line in enumerate(text.split("\n"), 1):
            for match in LINK.finditer(line):
                target = match.group(1)
                if EXTERNAL.match(target):
                    continue
                if not (doc.parent / target).exists():
                    problems += 1
                    rel = doc.relative_to(root).as_posix()
                    print(f"{rel}:{line_no}: link to a file not in the tarball: {target}")

    if problems:
        print(f"{problems} problem(s) in the release tarball")
        return 1
    print("release tarball: no repo-only paths, no dangling Markdown links")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
