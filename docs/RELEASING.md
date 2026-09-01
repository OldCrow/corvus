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
5. The working tree is clean and local `main` matches `origin/main`.

## Version ladder

6. Bump `kVersion{Major,Minor,Patch}` in `include/corvus/corvus.h`
   and `project(VERSION)` in `CMakeLists.txt` in the same commit. The
   configure-time gate fails the build if they disagree, but only a
   configure run exercises it — so run one.
7. Run the consistency check against the intended tag:

   ```sh
   python tools/check_release_version.py --tag vX.Y.Z
   ```

   It verifies header constants == CMake version == intended tag.
8. Fresh release configure + build + full ctest on this machine.

## Tag and publish

9. Signed annotated tag: `git tag -s vX.Y.Z -m "corvus vX.Y.Z"`,
   then push the tag. Verify the signature with `git tag -v vX.Y.Z`.
10. Re-run `python tools/check_release_version.py` (no `--tag`): it
    now also compares against `git describe` and catches a tag placed
    on the wrong commit. This step exists because two historical
    releases were tagged while the version pair still read a stale
    number.
11. Confirm CI green on the tagged commit.
12. GitHub release from the tag. The source tarball GitHub generates
    is the release artifact — corvus distributes source only. Release
    notes cover: surface changes, bound changes (tightened or newly
    validated cells), and notable fixes.

## After

13. Close the milestone.
14. Run any downstream handshakes (consumer libraries waiting on the
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
