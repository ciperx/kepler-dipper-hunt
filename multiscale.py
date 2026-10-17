"""Multi-scale re-scan of cached Kepler light curves.

The single-window detector uses a 61-cadence rolling median (~1.25 d) as its
baseline, so any dip longer than that window is absorbed into the trend and
disappears. This script reruns the same detection at several windows and
unions the results, which restores sensitivity to 1-3 day dips.

Only cached FITS are used, so a full pass over the cache costs no downloads.

  multiscale.py --windows 61,121,241,481 --min-depth 0.05
"""
import argparse
import glob
import json
import os
import sys
import warnings

import numpy as np
from astropy.io import fits

warnings.filterwarnings("ignore")
sys.path.insert(0, "/opt/agent/astro")
import hunt4  # noqa: E402

CACHE = hunt4.CACHE
OUT = "/opt/agent/astro/multiscale.jsonl"


def load(path):
    try:
        with fits.open(path, memmap=False) as h:
            d = h[1].data
            t = np.asarray(d["TIME"], float)
            f = np.asarray(d["PDCSAP_FLUX"], float)
            q = np.asarray(d["SAP_QUALITY"], int)
            names = set(d.columns.names)
            m1 = np.asarray(d["MOM_CENTR1"], float) if "MOM_CENTR1" in names else None
            m2 = np.asarray(d["MOM_CENTR2"], float) if "MOM_CENTR2" in names else None
    except Exception:
        return None
    ok = np.isfinite(t) & np.isfinite(f) & (q == 0)
    t, f = t[ok], f[ok]
    if m1 is not None:
        m1, m2 = m1[ok], m2[ok]
    if len(f) < 500:
        return None
    med = np.median(f)
    if not np.isfinite(med) or med <= 0:
        return None
    return t, f / med, m1, m2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--windows", default="61,121,241,481")
    ap.add_argument("--min-depth", type=float, default=0.05)
    ap.add_argument("--min-snr", type=float, default=50.0)
    ap.add_argument("--out", default=OUT)
    args = ap.parse_args()

    wins = [int(w) for w in args.windows.split(",")]
    files = sorted(glob.glob(os.path.join(CACHE, "*.fits")))
    print(f"{len(files)} cached curves, windows={wins}, "
          f"min_depth={args.min_depth} min_snr={args.min_snr}", flush=True)

    n_hit = 0
    n_done = 0
    fh = open(args.out, "w")
    for i, path in enumerate(files, 1):
        kepid = int(os.path.basename(path).replace(".fits", ""))
        lc = load(path)
        if lc is None:
            continue
        t, f, m1, m2 = lc
        n_done += 1

        best = None
        for w in wins:
            dips, sc = hunt4.detect(t, f, win=w, min_depth=args.min_depth,
                                    min_snr=args.min_snr)
            for d in dips:
                if d["depth"] < args.min_depth or d["snr"] < args.min_snr:
                    continue
                if not (0.0 < d["depth"] < 1.0):
                    continue
                if best is None or d["depth"] > best["depth"]:
                    best = dict(d)
                    best["window"] = w
                    best["scatter"] = sc

        if best is None:
            if i % 2000 == 0:
                print(f"{i}/{len(files)} scanned, hits={n_hit}", flush=True)
            continue

        cen = hunt4.centroid_shift(m1, m2, best)
        rec = dict(kepid=kepid, depth=best["depth"], snr=best["snr"],
                   ratio=best["ratio"], ratio_half=best["ratio_half"],
                   ratio_tenth=best["ratio_tenth"], dur_d=best["dur_d"],
                   window=best["window"], scatter=best["scatter"],
                   t_in=best["t_in"], centroid_millipix=cen)
        lopsided = best["ratio"] <= 0.67 or best["ratio"] >= 1.5
        rec["lopsided"] = bool(lopsided)
        rec["centroid_ok"] = bool(cen is not None and cen <= 1.0)
        fh.write(json.dumps(rec) + "\n")
        fh.flush()
        n_hit += 1
        print(f"  HIT {kepid} d={best['depth']:.3f} dur={best['dur_d']:.2f}d "
              f"w={best['window']} r={best['ratio']:.2f} cen={cen if cen is None else round(cen, 2)} "
              f"lop={lopsided}", flush=True)

    fh.close()
    print(f"done: {n_done} curves, {n_hit} hits -> {args.out}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
