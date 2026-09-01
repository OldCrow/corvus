# Security Policy

## Reporting a vulnerability

Report vulnerabilities privately through GitHub's security advisory
form: [Report a vulnerability](https://github.com/OldCrow/corvus/security/advisories/new).
Please do not open a public issue for a suspected vulnerability.

Reports receive an acknowledgment within a week. Fixes land as a patch
release with credit to the reporter, unless the reporter prefers
anonymity.

## Scope

corvus is a computational library with a small attack surface: it
allocates nothing, throws nothing, and touches no I/O, network, or
files. The security-relevant surface is memory safety at the API
boundary — span handling, table indexing, and the SIMD tail path.
Reports most likely to matter:

- Out-of-bounds reads or writes reachable through the documented API
  contract (matched span lengths, exact aliasing).
- Table-gather indexing escapes on any input, including NaN and
  infinity.
- Miscompilation reports for a supported compiler/target pair that
  produce memory-unsafe code.

Numerical accuracy defects — a result outside its documented bound —
are quality bugs, not vulnerabilities. Report those as ordinary GitHub
issues.

## Supported versions

The latest release receives security fixes. Older releases do not;
upgrading is the fix.
