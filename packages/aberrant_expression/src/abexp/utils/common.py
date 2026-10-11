import os
import sys
import warnings


def url_encode(raw: str):
    """
    Returns a url-encoded version of a raw string
    """
    from urllib.parse import quote

    return quote(raw.strip())


def url_decode(raw: str):
    """
    Returns a url-decoded version of a raw string
    """
    from urllib.parse import quote, unquote

    return unquote(raw.strip())


def set_polars_memory_env(limit: int):
    """Tell polars the memory limit (in bytes) of this process.

    polars reads the memory limit only at the root of /sys/fs/cgroup, which is the own cgroup of a container.
    A SLURM job sits in a nested cgroup, so polars would plan with the RAM of the node.
    polars reads the limit once, so call this function before the first import of polars.
    It warns if polars is imported already.

    polars 1.x reads the limit in bytes from POLARS_OVERRIDE_TOTAL_MEMORY.
    polars 2 reads it in units of 10**6 bytes from POLARS_OVERRIDE_TOTAL_MEMORY_MB.
    This function sets both with os.environ.setdefault, so a variable that is set already wins.
    """
    if "polars" in sys.modules:
        warnings.warn("set_polars_memory_env() runs after the import of polars, which may ignore it")
    os.environ.setdefault("POLARS_OVERRIDE_TOTAL_MEMORY", str(limit))
    os.environ.setdefault("POLARS_OVERRIDE_TOTAL_MEMORY_MB", str(limit // 10**6))
