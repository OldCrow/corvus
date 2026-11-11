# Release Checklist

Maintainer ritual for tagging a corvus release. Work top to bottom; a
failed step stops the release until fixed.

## Preconditions

1. The milestone is clean: every issue closed or explicitly moved.
2. CI is green on `main` — check per-job conclusions, not just the
   run banner.
3. The validation matrix in [ACCURACY.md](ACCURACY.md) is current:
   every cell either validated or explicitly marked open, and no open
   cell contradicts a bound the release notes will claim.
4. [PERFORMANCE.md](PERFORMANCE.md) status headers reflect what has
   and has not been reproduced.
5. [PLAN.md](../PLAN.md) holds live state only: finished rows collapsed to
   one line, superseded status blocks and closed issues' working notes cut
   (issues, commit messages and git history keep them). Keep the
   per-release errata list.
6. A docs pass is done: every document checked against the code and
   `tests/data/bounds.txt` for stale, inconsistent or wrong statements,
   and trimmed; `python tools/claims_check.py` passes.
7. The release branch (`dev/vX.Y.Z`) has `main`'s own commits (Dependabot
   bumps, say) merged in and is then merged into `main`; the working tree
   is clean and local `main` matches `origin/main`. Precondition 2's CI
   run is on that merge.

## Version ladder

8. Bump `kVersion{Major,Minor,Patch}` in `include/corvus/corvus.h`
   and `project(VERSION)` in `CMakeLists.txt` in the same commit. The
   configure-time gate fails the build if they disagree, but only a
   configure run exercises it — so run one.
9. Run the consistency check against the intended tag:

   ```sh
   python tools/check_release_version.py --tag vX.Y.Z
   ```

   It verifies header constants == CMake version == intended tag.
10. Fresh release configure + build + full ctest on this machine.

## Tag and publish

11. Signed annotated tag: `git tag -s vX.Y.Z -m "corvus vX.Y.Z"`,
   then push the tag. Verify the signature with `git tag -v vX.Y.Z`.
12. Re-run `python tools/check_release_version.py` (no `--tag`): it
    now also compares against `git describe` and catches a tag placed
    on the wrong commit. This step exists because two historical
    releases were tagged while the version pair still read a stale
    number.
13. Confirm CI green on the tagged commit.
14. GitHub release from the tag. The source tarball GitHub generates
    is the release artifact — corvus distributes source only. It ships
    what a user or packager needs to build, use and verify corvus;
    `export-ignore` in `.gitattributes` keeps contributor and
    maintainer files (agent and plan files, CI, `tools/`, maintainer
    docs, bench evidence) in the repo only. A new top-level file or
    doc gets classified there in the change that adds it. CI's Linux
    job extracts an archive, runs `tools/check_release_tarball.py` on it
    (no repo-only path present, no Markdown link to an unshipped file)
    and builds and tests it, so a green run on the tagged commit covers
    this. Read the listing with `git archive HEAD | tar -t` before
    tagging all the same. Release
    notes cover: surface changes, bound changes (tightened or newly
    validated cells), notable fixes, and **errata**. List every published
    bound corrected under VERSIONING.md's erratum clause, with the cell,
    the old and new bound, and the witness input. ACCURACY.md's
    "Erratum" paragraphs and PLAN.md's per-release list are the
    sources.

## After

15. Close the milestone.
16. Run any downstream handshakes (consumer libraries waiting on the
    release).

## Binary-channel prerequisites (recorded for a future decision)

corvus currently ships no binaries. Any future binary channel
(Homebrew bottles, vcpkg, conan, GitHub release binaries) must first
settle, in order:

- **NOTICE propagation**: binary artifacts include the NOTICE file
  and Highway's Apache-2.0 license text (see NOTICE at the repo
  root — the source tree carries no Highway code, binaries do).
- **macOS**: codesign + notarization for distributed dylibs/archives.
- **Windows**: Authenticode signing for distributed DLLs.
- **ABI/SONAME policy** for any shared-library variant — today no ABI
  is promised at all ([VERSIONING.md](VERSIONING.md)), which is
  incompatible with a shared-library channel as-is.
- **Per-channel metadata**: license fields must state MIT +
  Apache-2.0 (Highway) for binary artifacts, MIT alone for source.

The CMake plumbing a channel needs (stable tag tarballs, system-
Highway resolution, installed package config, pkg-config) already
exists and is CI-tested.
