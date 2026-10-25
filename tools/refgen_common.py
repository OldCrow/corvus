"""Shared reference-writer helpers (#13 N14).

float(mpf) double-rounds subnormal results: mpmath's to_float rounds the
mantissa to 53 bits, then ldexp rounds AGAIN onto the 2^-1074 grid — two
roundings, up to 1 ulp off correctly rounded in the subnormal band.
round_to_double rounds once.
"""

import math
import random
import struct
from typing import Optional

import mpmath as mp

_MIN_NORMAL = 2.0 ** -1022


def round_to_double(v):
    """Nearest double to mpf/int/float v, single rounding, ties-to-even.

    Normal range: float(mpf) IS a single correct rounding (53-bit mantissa,
    exponent representable), so it is used directly. Subnormal range: scale
    by 2^1074 (exact — pure exponent shift in binary fp), round to an
    integer once (ties-to-even, matching IEEE), ldexp back (exact: the
    integer is <= 2^52).

    Precision-independent (#82 T1): v keeps every bit it carries, whatever
    the caller's ambient precision. mp.mpf(v) on an mpf re-rounds it to the
    working precision, so a caller back at dps 15 got a 53-bit v here and
    the subnormal path rounded twice -- 13 trigamma rows near x = 1e308
    were stored 1 ulp off.
    """
    if not isinstance(v, mp.mpf):
        with mp.workprec(2200):  # ints and floats convert exactly
            v = mp.mpf(v)
    if not mp.isfinite(v):
        return float(v)
    if v == 0:
        # An mpf zero carries no usable sign; callers owning a signed-zero
        # contract (erfinv(-0) = -0) must emit that sign explicitly.
        return 0.0
    if abs(v) >= _MIN_NORMAL:
        return float(v)
    # t = |v| * 2^1074 < 2^52. frac = t - n is exact at this workprec:
    # t's mantissa spans <= prec bits below an MSB <= 2^51, so the
    # fractional part needs < prec bits — no bits are dropped, and the
    # tie comparison against exactly 0.5 is therefore decisive.
    with mp.workprec(max(mp.mp.prec, int(v.man).bit_length()) + 64):
        t = abs(v) * mp.mpf(2) ** 1074
        n = int(mp.floor(t))
        frac = t - n
        if frac > 0.5 or (frac == 0.5 and (n & 1)):
            n += 1
    return math.copysign(math.ldexp(float(n), -1074), 1.0 if v > 0 else -1.0)


# ---- point selection (#82 T3) -------------------------------------------------
# Shared by the reference generators, the hunt strata (tools/hunt) and the
# strata appendix writer (tools/hunt/strata_append.py), so a stratum means
# the same population wherever it is drawn.
def bits(d: float) -> int:
    return struct.unpack("<Q", struct.pack("<d", d))[0]


def from_bits(u: int) -> float:
    return struct.unpack("<d", struct.pack("<Q", u))[0]


def ordered(d: float) -> int:
    u = struct.unpack("<q", struct.pack("<d", d))[0]
    return u if u >= 0 else -(u & 0x7FFFFFFFFFFFFFFF)


def neighbours(x: float, k: int) -> list[float]:
    """x and its k bit-neighbours on each side (finite, same sign)."""
    out = []
    b = bits(x)
    for j in range(-k, k + 1):
        nb = b + j
        if 0 <= nb < 0x7FF0000000000000 or 0x8000000000000000 <= nb < 0xFFF0000000000000:
            out.append(from_bits(nb))
    return out


def full_mantissa(rng: random.Random, lo: float, hi: float) -> float:
    """Uniform on [lo, hi] then re-randomised in its low 30 bits, so rows
    carry full mantissas rather than the sampler's short ones."""
    x = rng.uniform(lo, hi)
    if x == 0 or not math.isfinite(x):
        return x
    return from_bits((bits(x) & ~0x3FFFFFFF) | rng.getrandbits(30))


def log_uniform(rng: random.Random, lo: float, hi: float) -> float:
    return math.exp(rng.uniform(math.log(lo), math.log(hi)))


def subnormal(rng: random.Random, top_binades: Optional[int] = None) -> float:
    """Positive subnormal; top_binades=k restricts to the top k binades."""
    if top_binades is None:
        return from_bits(rng.getrandbits(52) | 1)
    e = rng.randint(52 - top_binades, 51)
    return from_bits((1 << e) | rng.getrandbits(e))


# ---- oracle worker memory cap -----------------------------------------------
# One mpmath evaluation at an extreme argument (huge shapes, deep series) can
# allocate gigabytes in seconds -- faster than a parent polling free memory
# can react. Z (16 GB) froze on 2026-10-22 and had to be hard-rebooted on
# 2026-10-25. So every oracle worker caps ITSELF: an allocation past the cap
# fails with MemoryError inside that one evaluation, which the caller records
# as an oracle failure, and the host keeps its headroom.
DEFAULT_WORKER_MB = 1024


def worker_mb() -> int:
    """Per-worker cap in MB (CORVUS_ORACLE_WORKER_MB overrides the default)."""
    import os
    try:
        return int(os.environ.get("CORVUS_ORACLE_WORKER_MB", DEFAULT_WORKER_MB))
    except ValueError:
        return DEFAULT_WORKER_MB


def cap_process_memory(mb: Optional[int] = None) -> bool:
    """Cap this process's committed memory at mb MB. Windows: a Job Object
    with a per-process memory limit; elsewhere: RLIMIT_AS. Returns whether
    the cap was applied. Call it first in a pool worker's initializer."""
    mb = worker_mb() if mb is None else mb
    limit = int(mb) << 20
    import sys
    if sys.platform != "win32":
        try:
            import resource
            resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
            return True
        except (ImportError, ValueError, OSError):
            return False
    import ctypes
    from ctypes import wintypes

    class BasicLimit(ctypes.Structure):
        _fields_ = [("PerProcessUserTimeLimit", ctypes.c_int64),
                    ("PerJobUserTimeLimit", ctypes.c_int64),
                    ("LimitFlags", wintypes.DWORD),
                    ("MinimumWorkingSetSize", ctypes.c_size_t),
                    ("MaximumWorkingSetSize", ctypes.c_size_t),
                    ("ActiveProcessLimit", wintypes.DWORD),
                    ("Affinity", ctypes.c_size_t),
                    ("PriorityClass", wintypes.DWORD),
                    ("SchedulingClass", wintypes.DWORD)]

    class IoCounters(ctypes.Structure):
        _fields_ = [(n, ctypes.c_uint64) for n in
                    ("ReadOps", "WriteOps", "OtherOps",
                     "ReadBytes", "WriteBytes", "OtherBytes")]

    class ExtendedLimit(ctypes.Structure):
        _fields_ = [("Basic", BasicLimit), ("Io", IoCounters),
                    ("ProcessMemoryLimit", ctypes.c_size_t),
                    ("JobMemoryLimit", ctypes.c_size_t),
                    ("PeakProcessMemoryUsed", ctypes.c_size_t),
                    ("PeakJobMemoryUsed", ctypes.c_size_t)]

    # Never raise: this runs in pool initializers, and a raising initializer
    # makes multiprocessing respawn dying workers forever.
    try:
        k32 = ctypes.WinDLL("kernel32", use_last_error=True)
        k32.CreateJobObjectW.restype = wintypes.HANDLE
        k32.CreateJobObjectW.argtypes = [ctypes.c_void_p, wintypes.LPCWSTR]
        k32.GetCurrentProcess.restype = wintypes.HANDLE
        k32.GetCurrentProcess.argtypes = []
        k32.SetInformationJobObject.restype = wintypes.BOOL
        k32.SetInformationJobObject.argtypes = [wintypes.HANDLE, ctypes.c_int,
                                                ctypes.c_void_p, wintypes.DWORD]
        k32.AssignProcessToJobObject.restype = wintypes.BOOL
        k32.AssignProcessToJobObject.argtypes = [wintypes.HANDLE, wintypes.HANDLE]
        job = k32.CreateJobObjectW(None, None)
        if not job:
            return False
        info = ExtendedLimit()
        info.Basic.LimitFlags = 0x100  # JOB_OBJECT_LIMIT_PROCESS_MEMORY
        info.ProcessMemoryLimit = limit
        ok = k32.SetInformationJobObject(job, 9,  # JobObjectExtendedLimitInformation
                                         ctypes.byref(info), ctypes.sizeof(info))
        ok = ok and k32.AssignProcessToJobObject(job, k32.GetCurrentProcess())
        return bool(ok)  # the job handle stays open for the process's life
    except Exception:
        return False


def _self_test():
    # Tie at the +0/min-subnormal boundary: 2^-1075 is exactly half the
    # smallest subnormal; ties-to-even rounds to 0 (even).
    assert round_to_double(mp.mpf(2) ** -1075) == 0.0
    # Just above the tie rounds up to the min subnormal.
    assert round_to_double(mp.mpf(2) ** -1075 * (1 + mp.mpf(2) ** -40)) == 2.0 ** -1074
    # Just below rounds down to 0.
    assert round_to_double(mp.mpf(2) ** -1075 * (1 - mp.mpf(2) ** -40)) == 0.0
    # Odd-multiple tie rounds to the EVEN neighbor (up): 3 * 2^-1075 -> 2 * 2^-1074.
    assert round_to_double(mp.mpf(3) * mp.mpf(2) ** -1075) == 2 * 2.0 ** -1074
    # Even-multiple neighborhood: 5 * 2^-1075 -> tie between 2 and 3 ulps -> 2 (even...
    # 5/2 = 2.5 -> even = 2).
    assert round_to_double(mp.mpf(5) * mp.mpf(2) ** -1075) == 2 * 2.0 ** -1074
    # Subnormal/normal boundary: rounding up across it is exact.
    below = mp.mpf(2) ** -1022 * (1 - mp.mpf(2) ** -70)
    assert round_to_double(below) == 2.0 ** -1022
    # Normal range unchanged versus float().
    assert round_to_double(mp.mpf(1) / 3) == float(mp.mpf(1) / 3)
    # Sign preserved through the subnormal path.
    assert round_to_double(-(mp.mpf(3) * mp.mpf(2) ** -1075)) == -(2 * 2.0 ** -1074)
    # The double-rounding defect this module exists to fix, with a
    # construction that actually DISCRIMINATES (#35 M1: the previous one,
    # 2^-1070*(1 + 2^-53 + 2^-90), scales to 16 + 2^-49 + 2^-86 -- nowhere
    # near a half-integer of the 2^-1074 grid, so float() passes it too,
    # and at import-time default precision it collapsed to exactly 2^-1070
    # besides). u = 33*2^-1075 + 2^-1130 scales to 16.5 + 2^-55: single
    # rounding gives 17 (just above the tie); mpmath's float() first
    # rounds the mantissa to 53 bits -- landing EXACTLY on the 16.5 tie --
    # then ties-to-even DOWN to 16. Verified against pinned mpmath 1.4.1:
    # float(u) returns 16*2^-1074, so this assert fails on any regression
    # to float() in the subnormal branch.
    u = mp.mpf(33) * mp.mpf(2) ** -1075 + mp.mpf(2) ** -1130
    assert round_to_double(u) == math.ldexp(17.0, -1074)
    assert float(u) == math.ldexp(16.0, -1074), \
        "float(mpf) no longer double-rounds; round_to_double may be droppable"
    # #82 T1: the same value handed over at the default 53-bit precision
    # (a generator that restored mp.dps before rounding) must not re-round.
    with mp.workprec(53):
        assert round_to_double(u) == math.ldexp(17.0, -1074)


# #35 M1: run under explicit precision -- module import happens BEFORE the
# generators set their working dps, and the old ambient-precision run let
# every discriminating construction collapse to an exactly-representable
# value at prec 53.
with mp.workdps(60):
    _self_test()
