"""Null-test the detector: how often does it fire on random quiet-ish stars?"""
import os
import random
import sys

import numpy as np

sys.path.insert(0, "/opt/agent/astro")
import hunt4  # noqa: E402

BASE = "/opt/agent/astro"
REPO = "/opt/agent/repo"

rows = []
with open(f"{BASE}/manifest.tsv") as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) == 2:
            try:
                rows.append((int(p[0]), p[1]))
            except ValueError:
                pass

random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 11)
random.shuffle(rows)

n = 0
stats = []
for kepid, uri in rows[:25]:
    path = os.path.join(REPO, "fits_cache", f"{kepid}.fits")
    if not os.path.exists(path) or os.path.getsize(path) < 50000:
        if not hunt4.fetch(uri, kepid):
            continue
    try:
        from astropy.io import fits
        with fits.open(path, memmap=False) as h:
            d = h[1].data
            t = np.asarray(d["TIME"], float)
            f = np.asarray(d["PDCSAP_FLUX"], float)
            q = np.asarray(d["SAP_QUALITY"], int)
    except Exception:
        continue
    ok = np.isfinite(t) & np.isfinite(f) & (q == 0)
    t, f = t[ok], f[ok]
    if len(f) < 500:
        continue
    med = np.median(f)
    if not np.isfinite(med) or med <= 0:
        continue
    f = f / med

    dips, sc = hunt4.detect(t, f / med, min_depth=0.02)
    strong = [x for x in dips if x["depth"] >= 0.05 and x["snr"] >= 50]
    deep = max((x["depth"] for x in dips), default=0.0)
    stats.append((kepid, len(f), sc, len(dips), len(strong), deep))
    n += 1
    print(f"{kepid} n={len(f)} scatter={sc:.4f} dips>2%={len(dips)} "
          f"strong={len(strong)} deepest={deep:.3f}", flush=True)

print(f"\n{ n } stars; with >=1 strong dip: "
      f"{sum(1 for s in stats if s[4] > 0)}")
