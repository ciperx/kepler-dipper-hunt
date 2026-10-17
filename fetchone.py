#!/usr/bin/env python3
"""Download one or more Kepler light curves by KIC id into HUNT_HOME/fits_cache.

  fetchone.py 8364969 5255552 8953257
"""
import os
import sys

BASE = os.environ.get("HUNT_HOME", "/opt/agent/repo")
CACHE = os.path.join(BASE, "fits_cache")
MANIFEST = "/opt/agent/astro/manifest.tsv"


def main(kepids):
    os.makedirs(CACHE, exist_ok=True)
    want = {int(k): None for k in kepids}
    with open(MANIFEST) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 2:
                continue
            try:
                k = int(parts[0])
            except ValueError:
                continue
            if k in want and want[k] is None:
                want[k] = parts[1]
                if all(v for v in want.values()):
                    break

    import subprocess

    for k, uri in want.items():
        if not uri:
            print(f"{k}: not in manifest")
            continue
        out = os.path.join(CACHE, f"{k}.fits")
        if os.path.exists(out) and os.path.getsize(out) > 50_000:
            print(f"{k}: cached")
            continue
        url = f"https://mast.stsci.edu/api/v0.1/Download/file?uri={uri}"
        r = subprocess.run(["curl", "-sL", "-m", "90", "-o", out + ".part", url])
        if r.returncode == 0 and os.path.getsize(out + ".part") > 50_000:
            os.replace(out + ".part", out)
            print(f"{k}: downloaded {os.path.getsize(out)} B")
        else:
            print(f"{k}: FAILED")


if __name__ == "__main__":
    main(sys.argv[1:])
