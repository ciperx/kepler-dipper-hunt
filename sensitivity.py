"""Detection completeness + false-positive rate for the dipper hunter.

Two measurements on the same cached Kepler light curves:

1. FALSE POSITIVE RATE - run the detector on the untouched curve, count how
   often it reports a dip that passes the pipeline's own thresholds.
2. RECOVERY (completeness) - inject a synthetic trapezoidal dip of known
   depth / duration / asymmetry at a random phase, rerun, count how often the
   detector recovers it. This is the sensitivity curve: what a null result
   does and does not rule out.

Usage: sensitivity.py --stars 60 --out /opt/agent/astro/sensitivity.json
"""
import argparse
import glob
import json
import os
import random
import sys

import numpy as np
from astropy.io import fits

sys.path.insert(0, "/opt/agent/astro")
import hunt4  # noqa: E402

CACHE = hunt4.CACHE          # /opt/agent/astro/f4
MIN_DEPTH = 0.05             # the pipeline's search threshold


def load(path):
    try:
        with fits.open(path, memmap=False) as h:
            d = h[1].data
            t = np.asarray(d["TIME"], float)
            f = np.asarray(d["PDCSAP_FLUX"], float)
            q = np.asarray(d["SAP_QUALITY"], int)
    except Exception:
        return None
    ok = np.isfinite(t) & np.isfinite(f) & (q == 0)
    t, f = t[ok], f[ok]
    if len(f) < 500:
        return None
    med = np.median(f)
    if not np.isfinite(med) or med <= 0:
        return None
    return t, f / med


def inject(t, f, depth, dur_d, ingress_frac, phase):
    """Multiply in a trapezoidal dip; ingress_frac<0.5 makes it lopsided."""
    g = f.copy()
    cad = float(np.median(np.diff(t)))
    t_in = t[0] + (t[-1] - t[0]) * phase
    t_out = t_in + dur_d
    rise_end = t_in + dur_d * ingress_frac
    fall_start = t_out - dur_d * ingress_frac
    m = (t >= t_in) & (t <= t_out)
    if not m.any():
        return None, t_in, dur_d
    tt = t[m]
    shape = np.ones_like(tt)
    r = tt < rise_end
    shape[r] = (tt[r] - t_in) / max(dur_d * ingress_frac, cad)
    fl = tt > fall_start
    shape[fl] = (t_out - tt[fl]) / max(dur_d * ingress_frac, cad)
    g[m] = g[m] * (1.0 - depth * np.clip(shape, 0.0, 1.0))
    return g, t_in + dur_d / 2.0, dur_d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stars", type=int, default=60)
    ap.add_argument("--seed", type=int, default=3)
    ap.add_argument("--out", default="/opt/agent/astro/sensitivity.json")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    files = sorted(glob.glob(os.path.join(CACHE, "*.fits")))
    rng.shuffle(files)

    curves = []
    for p in files:
        if len(curves) >= args.stars:
            break
        lc = load(p)
        if lc is None:
            continue
        curves.append((os.path.basename(p).replace(".fits", ""), *lc))
        if len(curves) % 20 == 0:
            print(f"loaded {len(curves)}/{args.stars}", flush=True)

    if not curves:
        print("no usable curves")
        return 1
    print(f"loaded {len(curves)} curves", flush=True)

    # ---- 1. false positives on the untouched curves
    fp = 0
    depths_fp = []
    for kepid, t, f in curves:
        dips, _ = hunt4.detect(t, f, min_depth=MIN_DEPTH)
        strong = [d for d in dips if d["depth"] >= MIN_DEPTH and d["snr"] >= 50]
        if strong:
            fp += 1
            depths_fp.append(max(d["depth"] for d in strong))
    print(f"\nFALSE POSITIVE: {fp}/{len(curves)} untouched curves fire "
          f"({100.0 * fp / len(curves):.1f}%)", flush=True)

    # ---- 2. recovery grid
    grid = []
    for depth in (0.02, 0.03, 0.05, 0.08, 0.12, 0.20):
        for dur in (0.3, 0.7, 1.5, 3.0):
            for a in (0.30, 0.50):
                ndet = nlop = 0
                ntry = 0
                for kepid, t, f in curves:
                    phase = rng.uniform(0.08, 0.85)
                    g, centre, dd = inject(t, f, depth, dur, a, phase)
                    if g is None:
                        continue
                    ntry += 1
                    dips, _ = hunt4.detect(t, g, min_depth=MIN_DEPTH)
                    for d in dips:
                        c = (d["t_in"] + d["t_out"]) / 2.0
                        if abs(c - centre) <= max(0.6, dd) and \
                           abs(d["depth"] - depth) / depth <= 0.6:
                            ndet += 1
                            if d["ratio"] <= 0.67 or d["ratio"] >= 1.5:
                                nlop += 1
                            break
                grid.append(dict(depth=depth, dur_d=dur, ingress_frac=a,
                                 n=ntry, detected=ndet, lopsided=nlop))
                print(f"  d={depth:<5} dur={dur:<4} in={a}  "
                      f"det={ndet:>3}/{ntry:<3} lop={nlop:>3}/{ntry:<3}", flush=True)

    json.dump(dict(n_curves=len(curves), min_depth=MIN_DEPTH,
                   false_positive=fp, fp_depths=depths_fp, grid=grid),
              open(args.out, "w"), indent=1)
    print("\nwrote", args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
