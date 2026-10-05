# 2026-10-02 — M1 quiet-bench gate abort (macOS 27 Golden Gate)

Side task C attempt. `tools/quiet_bench.sh -b build-m1-gg -t NEON -m 10 -s 10 -w 30`
at corvus `4534d30` (code `a78eddd`), 20:55–21:25 local, screen locked,
Claude desktop minimized, Backblaze paused (20:25–22:25), Cloudflare WARP
killed, the in-flight Time Machine backup skipped before arming.

Result: ABORT after 30 min. 173 ambient samples, min 11.88%, mean 28.3%,
none below 10%. Top consumer in 134 of 173 failed-gate snapshots:
`mediaanalysisd` (Photos video analysis, 85–180% of one core), which was
idle while the user was active and launched once the machine went idle —
the same consumer that defeated every 5% window in August. Also present:
`corespotlightd` / `spotlightknowledged.updater` (post-upgrade indexing
not yet settled), `backupd` (Time Machine's hourly schedule restarted
mid-window), `deleted_helper`, and `bztransmit` three times despite the
pause. `quiet_bench.log` is the full gate history for the run.
