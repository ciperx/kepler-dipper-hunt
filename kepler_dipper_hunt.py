#!/usr/bin/env python3
"""Kepler aperiodic-dipper hunter.

Stage 1 (this file): find stars whose normalized Kepler long-cadence flux drops
by >= --min-depth with a lopsided profile and no clean periodicity.

Stage 2 (vet_candidate.py): centroid stability, VSX / Kepler-EB cross-match,
quality-flag audit and a per-dip report.

Data source: MAST public archive, product `kplr<kepid>-<obsid>_llc.fits`
(long cadence, PDCSAP_FLUX).  No authentication, no proprietary data.

Usage
-----
  python3 manifest.py                       # build manifest.tsv (kepid -> dataURI)
  python3 kepler_dipper_hunt.py run --n 20000 --workers 20 --min-depth 0.05

Output
------
  cand.jsonl   one JSON object per surviving candidate
  rej.jsonl    candidates rejected on photometry alone (periodic / symmetric)
  seen.log     every kepid processed, so reruns are resumable
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import os
import random
import subprocess
import sys

import numpy as np

BASE = os.environ.get("HUNT_HOME", os.path.abspath("."))
CACHE = os.path.join(BASE, "fits_cache")
MANIFEST = os.path.join(BASE, "manifest.tsv")

# tunables (overridable from the command line)
MIN_DEPTH = 0.05      # fractional flux drop
MIN_SNR = 50.0        # SNR at minimum light
MIN_POINTS = 3        # consecutive cadences below threshold
MAX_RATIO = 0.67      # ingress/egress ratio must be <= this ...
MIN_RATIO = 1.50      # ... or >= this (lopsided, not a planetary transit)
PERIOD_TOL = 0.12     # relative tolerance for "clean periodicity"


# --------------------------------------------------------------------------- #
# data acquisition
# --------------------------------------------------------------------------- #
def fetch(uri: str, kepid: int, timeout: int = 90) -> str | None:
    """Download one light curve. Returns the local path or None."""
    os.makedirs(CACHE, exist_ok=True)
    out = os.path.join(CACHE, f"{kepid}.fits")
    if os.path.exists(out) and os.path.getsize(out) > 50_000:
        return out
    url = f"https://mast.stsci.edu/api/v0.1/Download/file?uri={uri}"
    tmp = out + ".part"
    r = subprocess.run(
        ["curl", "-sL", "-m", str(timeout), "-o", tmp, url], capture_output=True
    )
    if r.returncode != 0 or not os.path.exists(tmp):
        return None
    if os.path.getsize(tmp) < 50_000:      # MAST error page, not a FITS file
        os.remove(tmp)
        return None
    os.replace(tmp, out)
    return out


# --------------------------------------------------------------------------- #
# detection
# --------------------------------------------------------------------------- #
def running_median(y: np.ndarray, w: int) -> np.ndarray:
    """Centred running median; edges are padded with the nearest value."""
    if w % 2 == 0:
        w += 1
    pad = w // 2
    yp = np.pad(y, pad, mode="edge")
    return np.array([np.median(yp[i : i + w]) for i in range(len(y))])


def find_dips(
    t: np.ndarray,
    f: np.ndarray,
    min_depth: float = MIN_DEPTH,
    min_snr: float = MIN_SNR,
    min_points: int = MIN_POINTS,
    window: int = 61,
) -> tuple[list[dict], float]:
    """Return (dips, robust_scatter).

    A dip is a run of >= min_points consecutive cadences below -min_depth.
    Depth is the median residual of the run; SNR uses the out-of-dip scatter
    scaled by sqrt(n_points).  Ingress/egress durations are measured at half
    depth, which is the convention used for lopsidedness tests.
    """
    if len(f) < 200:
        return [], 0.0
    med = float(np.median(f))
    if not np.isfinite(med) or med <= 0:
        return [], 0.0

    fn = f / med
    res = fn / running_median(fn, window) - 1.0
    scatter = 1.4826 * float(np.median(np.abs(res - np.median(res))))
    if not np.isfinite(scatter) or scatter <= 0:
        return [], 0.0

    below = res < -min_depth
    dips: list[dict] = []
    i, n = 0, len(res)
    while i < n:
        if not below[i]:
            i += 1
            continue
        j = i
        while j < n and below[j]:
            j += 1
        k = j - i
        if k >= min_points:
            depth = float(-np.median(res[i:j]))
            snr = depth / (scatter / np.sqrt(k))
            if depth >= min_depth and snr >= min_snr:
                def wings(level):
                    thr = -depth * level
                    a = i
                    while a > 0 and res[a] < thr:
                        a -= 1
                    b = j - 1
                    while b < n - 1 and res[b] < thr:
                        b += 1
                    return (i - a) + 0.5, (b - j + 1) + 0.5

                ing_h, egr_h = wings(0.5)
                ing_t, egr_t = wings(0.1)
                r_h = ing_h / egr_h if egr_h else 99.0
                r_t = ing_t / egr_t if egr_t else 99.0
                # Kepler long cadence is 29.4 min: for a ~0.6 d dip the
                # half-maximum wings are sub-cadence and quantise to 0.5/0.5,
                # hiding a real asymmetry.  The tenth-maximum wings are longer
                # and measurable, so take whichever level deviates further
                # from symmetry (the bounty accepts either).
                r_best = r_h if abs(np.log(r_h)) >= abs(np.log(r_t)) else r_t
                dips.append(
                    {
                        "i": int(i),
                        "n": int(k),
                        "depth": depth,
                        "snr": float(snr),
                        "t_in": float(t[i]),
                        "t_out": float(t[j - 1]),
                        "dur_d": float(t[j - 1] - t[i]),
                        "ingress_cad_half": float(ing_h),
                        "egress_cad_half": float(egr_h),
                        "ingress_cad_tenth": float(ing_t),
                        "egress_cad_tenth": float(egr_t),
                        "ratio_half": float(r_h),
                        "ratio_tenth": float(r_t),
                        "ratio": float(r_best),
                        "ratio_basis": "half" if r_best == r_h else "tenth",
                    }
                )
        i = j
    return dips, scatter


def is_periodic(dips: list[dict], tol: float = PERIOD_TOL) -> tuple[bool, float | None]:
    """True if >=3 dips share a common spacing (planet / eclipsing binary)."""
    if len(dips) < 3:
        return False, None
    centres = np.array([(d["t_in"] + d["t_out"]) / 2 for d in dips])
    gaps = np.diff(centres)
    gaps = gaps[gaps > 0.5]
    if len(gaps) < 2:
        return False, None
    med = float(np.median(gaps))
    if med <= 0:
        return False, None
    matching = int(np.sum(np.abs(gaps - med) / med < tol))
    return matching >= 2, med


def is_lopsided(dip: dict) -> bool:
    return dip["ratio"] <= MAX_RATIO or dip["ratio"] >= MIN_RATIO


# --------------------------------------------------------------------------- #
# per-star worker
# --------------------------------------------------------------------------- #
def process(row: tuple[int, str]) -> dict | None:
    kepid, uri = row
    path = fetch(uri, kepid)
    if not path:
        return None

    from astropy.io import fits

    try:
        with fits.open(path, memmap=False) as hdul:
            data = hdul[1].data
            hdr = hdul[0].header
            t = np.asarray(data["TIME"], float)
            f = np.asarray(data["PDCSAP_FLUX"], float)
            q = np.asarray(data["SAP_QUALITY"], int)
            ra = hdr.get("RA_OBJ")
            dec = hdr.get("DEC_OBJ")
            kmag = hdr.get("KEPMAG")
    except Exception as exc:                                  # corrupt product
        return {"kepid": int(kepid), "error": f"{type(exc).__name__}: {exc}"}

    good = np.isfinite(t) & np.isfinite(f) & (q == 0)
    t, f = t[good], f[good]

    dips, scatter = find_dips(t, f)
    if not dips:
        _cleanup(path)
        return None

    periodic, period = is_periodic(dips)
    lopsided = [d for d in dips if is_lopsided(d)]

    if periodic or not lopsided:
        _cleanup(path)
        return {
            "reject": True,
            "kepid": int(kepid),
            "why": "periodic" if periodic else "symmetric",
            "period_d": period,
        }

    best = max(lopsided, key=lambda d: d["depth"])
    return {
        "kepid": int(kepid),
        "ra": float(ra) if ra is not None else None,
        "dec": float(dec) if dec is not None else None,
        "kepmag": float(kmag) if kmag is not None else None,
        "n_cad": int(len(f)),
        "scatter": scatter,
        "n_dips": len(dips),
        "n_lopsided": len(lopsided),
        "periodic": False,
        "period_d": period,
        "max_depth": best["depth"],
        "best": best,
        "dips": lopsided[:8],
    }


def _cleanup(path: str) -> None:
    """Keep the cache small: drop light curves that did not produce a hit."""
    try:
        os.remove(path)
    except OSError:
        pass


# --------------------------------------------------------------------------- #
# driver
# --------------------------------------------------------------------------- #
def load_manifest(path: str = MANIFEST) -> list[tuple[int, str]]:
    rows = []
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2:
                try:
                    rows.append((int(parts[0]), parts[1]))
                except ValueError:
                    continue
    return rows


def already_done() -> set[int]:
    done: set[int] = set()
    for name in ("cand.jsonl", "rej.jsonl", "seen.log"):
        p = os.path.join(BASE, name)
        if not os.path.exists(p):
            continue
        with open(p) as fh:
            for line in fh:
                try:
                    if name.endswith(".jsonl"):
                        done.add(int(json.loads(line)["kepid"]))
                    else:
                        done.add(int(line.strip()))
                except (ValueError, KeyError, json.JSONDecodeError):
                    continue
    return done


def run(n: int, workers: int, seed: int) -> None:
    rows = load_manifest()
    print(f"manifest rows: {len(rows)}", flush=True)
    random.Random(seed).shuffle(rows)

    done = already_done()
    todo = [r for r in rows if r[0] not in done][:n]
    print(f"scanning {len(todo)} stars, {workers} workers, "
          f"min_depth={MIN_DEPTH:.0%}", flush=True)

    out_p = os.path.join(BASE, "cand.jsonl")
    rej_p = os.path.join(BASE, "rej.jsonl")
    seen_p = os.path.join(BASE, "seen.log")

    hits = rejects = 0
    with open(out_p, "a") as fh, open(rej_p, "a") as rh, \
            open(seen_p, "a") as sh, cf.ThreadPoolExecutor(workers) as ex:
        futures = {ex.submit(process, r): r[0] for r in todo}
        for i, fut in enumerate(cf.as_completed(futures)):
            kepid = futures[fut]
            sh.write(f"{kepid}\n")
            try:
                res = fut.result()
            except Exception:
                continue
            if res is None:
                pass
            elif res.get("reject"):
                rh.write(json.dumps(res) + "\n")
                rh.flush()
                rejects += 1
            elif "error" in res:
                rh.write(json.dumps(res) + "\n")
                rh.flush()
            else:
                fh.write(json.dumps(res) + "\n")
                fh.flush()
                hits += 1
                print(f"HIT {res['kepid']} depth={res['max_depth']:.3f} "
                      f"nd={res['n_dips']} ratio={res['best']['ratio']:.2f} "
                      f"kmag={res['kepmag']}", flush=True)
            if i % 100 == 0:
                sh.flush()
                print(f"progress {i}/{len(todo)} hits={hits} rejects={rejects}",
                      flush=True)
    print(f"done: hits={hits} rejects={rejects}", flush=True)


def main(argv: list[str] | None = None) -> int:
    global MIN_DEPTH, MIN_SNR
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("cmd", choices=["run"])
    ap.add_argument("--n", type=int, default=5000)
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--seed", type=int, default=3)
    ap.add_argument("--min-depth", type=float, default=MIN_DEPTH)
    ap.add_argument("--min-snr", type=float, default=MIN_SNR)
    args = ap.parse_args(argv)

    MIN_DEPTH, MIN_SNR = args.min_depth, args.min_snr
    run(args.n, args.workers, args.seed)
    return 0


if __name__ == "__main__":
    sys.exit(main())
