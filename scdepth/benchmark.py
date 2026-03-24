#SPDX-License-Identifier: MIT
#Copyright (c) Gavin W. Wilson

import os, time, typing, resource, psutil, sys

def _ru_maxrss_bytes(ru_maxrss: int) -> int:
    #Convert ru_maxrss to bytes (linux: KB, darwin: bytes)
    if sys.platform == "darwin":
        return int(ru_maxrss)
    else:
        return int(ru_maxrss) * 1024 

def run(fn: typing.Callable, *args, **kwargs):
    proc = psutil.Process(os.getpid())

    rss_pre = proc.memory_info().rss
    r_pre = resource.getrusage(resource.RUSAGE_SELF)
    peak_pre = _ru_maxrss_bytes(r_pre.ru_maxrss)
    wall_pre = time.perf_counter()

    result = fn(*args, **kwargs)

    wall_post = time.perf_counter()
    r_post = resource.getrusage(resource.RUSAGE_SELF)
    rss_post = proc.memory_info().rss
    peak_post = _ru_maxrss_bytes(r_post.ru_maxrss)

    user = r_post.ru_utime - r_pre.ru_utime
    sysc = r_post.ru_stime - r_pre.ru_stime

    stats = {
        "wall_m": (wall_post - wall_pre)/60,
        "cpu_m": (user + sysc) / 60,
        "rss_after_bytes": rss_post,
        "rss_delta_bytes": max(rss_post - rss_pre, 0),
        "peak_rss_after_bytes": peak_post,
        "peak_rss_delta_bytes": max(peak_post - peak_pre, 0),
    }
    return result, stats
