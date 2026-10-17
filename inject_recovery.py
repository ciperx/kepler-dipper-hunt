"""Injection-recovery completeness test for the dipper hunter.

Takes quiet Kepler light curves (no dip at the search threshold), injects
synthetic Tabby-style dips of known depth / duration / asymmetry, and reruns
the same detector. The recovery fraction is what the pipeline *can* see; the
rest is what an honest null result does NOT rule out.

Output: JSON with a recovery matrix + a text table.
"""
import argparse
import json
import os
import random
import sys

import numpy as np
from astropy.io import fits

sys.path.insert(0, "/opt/agent/astro")
import hunt4  # noqa: E402

BASE = "/opt/agent/astro"
REPO = "/opt/agent/repo"
CACHE = os.path.join(REPO, "fits_cache")


def load_manifest():
    rows = []
    with open(f"{BASE}/manifest.tsv") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) == 2:
                try:
                    rows.append((int(p[0]), p[1]))
                except ValueError:
                    pass
    return rows


def get_lc(kepid, uri):
    path = os.path.join(CACHE, f"{kepid}.fits")
    if not os.path.exists(path) or os.path.getsize(path) < 50000:
        if not hunt4.fetch(uri, kepid):
            return None
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


def inject(t, f, depth, dur_d, ingress_frac, phase=0.37):
    """Multiply in a trapezoidal dip; ingress_frac<0.5 makes it lopsided."""
    g = f.copy()
    cad = float(np.median(np.diff(t)))
    span = t[-1] - t[0]
    centre = t[0] + span * phase
    half = dur_d / 2.0
    t_in = centre - half
    t_out = centre + half
    flat = max(dur_d * (1.0 - 2.0 * ingress_frac), 0.0)
    t_flat_in = t_in + dur_d * ingress_frac
    t_flat_out = t_out - dur_d * ingress_frac

    m = (t >= t_in) & (t <= t_out)
    tt = t[m]
    shape = np.ones_like(tt)
    rise = tt < t_flat_in
    shape[rise] = (tt[rise] - t_in) / max(dur_d * ingress_frac, cad)
    fall = tt > t_flat_out
    shape[fall] = (t_out - tt[fall]) / max(dur_d * ingress_frac, cad)
    shape = np.clip(shape, 0.0, 1.0)
    g[m] = g[m] * (1.0 - depth * shape)
    return g, centre, dur_d


def recovered(dips, centre, dur_d, depth):
    for d in dips:
        c = (d["t_in"] + d["t_out"]) / 2.0
        if abs(c - centre) <= max(0.6, dur_d) and abs(d["depth"] - depth) / depth <= 0.5:
            return d
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stars", type=int, default=25)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--out", default="/opt/agent/astro/injection.json")
    args = ap.parse_args()

    random.seed(args.seed)
    rows = load_manifest()
    random.shuffle(rows)

    depths = [0.02, 0.03, 0.05, 0.08, 0.12, 0.20]
    durs = [0.3, 0.7, 1.5]
    asyms = [0.30, 0.50]   # 0.30 -> ratio ~0.43 (lopsided), 0.50 -> 1.00 (symmetric)

    quiet = []
    tried = 0
    for kepid, uri in rows:
        if len(quiet) >= args.stars or tried >= args.stars * 12:
            break
        tried += 1
        lc = get_lc(kepid, uri)
        if lc is None:
            continue
        t, f = lc
        # "quiet" = no dip that would survive the pipeline's own thresholds
        dips, _ = hunt4.detect(t, f, min_depth=0.05)
        strong = [d for d in dips if d["depth"] >= 0.10 and d["snr"] >= 100]
        if strong:
            continue
        quiet.append((kepid, t, f))
        print(f"quiet {len(quiet):>2}/{args.stars} kepid={kepid} n={len(f)}", flush=True)

    if not quiet:
        print("no quiet stars found")
        return 1

    res = []
    for depth in depths:
        for dur in durs:
            for a in asyms:
                n_det = n_lop = 0
                for kepid, t, f in quiet:
                    g, centre, dd = inject(t, f, depth, dur, a)
                    dips, _ = hunt4.detect(t, g, min_depth=0.03)
                    hit = recovered(dips, centre, dd, depth)
                    if hit:
                        n_det += 1
                        if hit["ratio"] <= 0.67 or hit["ratio"] >= 1.5:
                            n_lop += 1
                res.append(dict(depth=depth, dur_d=dur, ingress_frac=a,
                                n=len(quiet), detected=n_det, lopsided=n_lop))
                print(f"d={depth:<5} dur={dur:<4} in={a}  detected={n_det}/{len(quiet)} "
                      f"lopsided={n_lop}/{len(quiet)}", flush=True)

    json.dump(dict(stars=[k for k, _, _ in quiet], results=res),
              open(args.out, "w"), indent=1)
    print("wrote", args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
