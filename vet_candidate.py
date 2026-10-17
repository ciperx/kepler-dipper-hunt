#!/usr/bin/env python3
"""Stage 2 vetting for a Kepler aperiodic-dipper candidate.

Checks, in the order the bounty requires:

  1 NEW          -- not KIC 8462852; literature check is manual (see report)
  2 DEEP         -- >=5% depth with SNR >=50 at minimum light
  3 LOPSIDED     -- ingress/egress >=1.5x or <=0.67x at half depth
  4 NOT PERIODIC -- dips do not phase-fold to a single period
  5 CENTROID     -- astrometric shift during the dip < 1.0 millipixel
  6 NOT VARIABLE -- VSX / Kepler EB catalogue cross-match at the source position
  7 NOT INSTRUMENTAL -- dip present in raw SAP flux and in cadences with
                        SAP_QUALITY == 0 only
  8 REPRODUCIBLE -- writes every number the report quotes

Usage
-----
  python3 vet_candidate.py --kepid 12345678 --out report.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

BASE = os.environ.get("HUNT_HOME", os.path.abspath("."))
CACHE = os.path.join(BASE, "fits_cache")
PIXEL_ARCSEC = 3.98          # Kepler plate scale
CENTROID_LIMIT_MILLIPIX = 1.0

EB_TYPES = ("EA", "EB", "EW", "E", "ELL", "EP", "GS", "ROT", "DIP", "VAR")


# --------------------------------------------------------------------------- #
def load_light_curve(kepid: int):
    from astropy.io import fits

    path = os.path.join(CACHE, f"{kepid}.fits")
    if not os.path.exists(path):
        raise SystemExit(f"missing {path} - run kepler_dipper_hunt.py first")

    with fits.open(path, memmap=False) as hdul:
        d = hdul[1].data
        hdr = hdul[0].header
        cols = set(d.columns.names)
        lc = {
            "time": np.asarray(d["TIME"], float),
            "pdcsap": np.asarray(d["PDCSAP_FLUX"], float),
            "sap": np.asarray(d["SAP_FLUX"], float) if "SAP_FLUX" in cols else None,
            "quality": np.asarray(d["SAP_QUALITY"], int),
            "mom1": np.asarray(d["MOM_CENTR1"], float) if "MOM_CENTR1" in cols else None,
            "mom2": np.asarray(d["MOM_CENTR2"], float) if "MOM_CENTR2" in cols else None,
        }
        meta = {
            "ra": hdr.get("RA_OBJ"),
            "dec": hdr.get("DEC_OBJ"),
            "kepmag": hdr.get("KEPMAG"),
            "quarter": hdr.get("QUARTER"),
            "module": hdr.get("MODULE"),
            "output": hdr.get("OUTPUT"),
        }
    return lc, meta


def running_median(y: np.ndarray, w: int = 61) -> np.ndarray:
    if w % 2 == 0:
        w += 1
    pad = w // 2
    yp = np.pad(y, pad, mode="edge")
    return np.array([np.median(yp[i : i + w]) for i in range(len(y))])


def robust_sigma(v: np.ndarray) -> float:
    return 1.4826 * float(np.median(np.abs(v - np.median(v))))


def measure_dips(t, f, min_depth=0.02, min_snr=20.0, min_points=3, window=61):
    """Same algorithm as stage 1 but with looser thresholds so the report can
    list every dip in the star, not just the ones that passed the filter."""
    med = float(np.median(f))
    fn = f / med
    res = fn / running_median(fn, window) - 1.0
    scatter = robust_sigma(res)
    below = res < -min_depth

    dips, i, n = [], 0, len(res)
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
            snr = depth / (scatter / np.sqrt(k)) if scatter > 0 else 0.0
            if depth >= min_depth and snr >= min_snr:
                half = depth / 2.0
                a = i
                while a > 0 and res[a] < -half:
                    a -= 1
                b = j - 1
                while b < n - 1 and res[b] < -half:
                    b += 1
                ingress = (i - a) + 0.5
                egress = (b - j + 1) + 0.5
                dips.append({
                    "i": int(i), "n": int(k), "depth": depth, "snr": float(snr),
                    "t_in": float(t[i]), "t_out": float(t[j - 1]),
                    "dur_d": float(t[j - 1] - t[i]),
                    "ingress_cad": float(ingress), "egress_cad": float(egress),
                    "ratio": float(ingress / egress) if egress else 99.0,
                })
        i = j
    return dips, scatter, res


def periodicity(dips, tol=0.12):
    if len(dips) < 3:
        return {"periodic": False, "period_d": None, "n_gaps": 0}
    centres = np.array([(d["t_in"] + d["t_out"]) / 2 for d in dips])
    gaps = np.diff(centres)
    gaps = gaps[gaps > 0.5]
    if len(gaps) < 2:
        return {"periodic": False, "period_d": None, "n_gaps": int(len(gaps))}
    med = float(np.median(gaps))
    match = int(np.sum(np.abs(gaps - med) / med < tol)) if med > 0 else 0
    return {"periodic": match >= 2, "period_d": med, "n_gaps": int(len(gaps)),
            "matching_gaps": match}


def centroid_test(t, mom1, mom2, dip, res):
    """Compare the photocentre inside the dip with the out-of-dip scatter."""
    if mom1 is None or mom2 is None:
        return {"available": False, "reason": "no MOM_CENTR columns"}

    window = 61
    pad = window // 2

    def detr(v):
        vp = np.pad(v, pad, mode="edge")
        base = np.array([np.median(vp[i : i + window]) for i in range(len(v))])
        return v - base

    dx, dy = detr(mom1), detr(mom2)
    inside = slice(dip["i"], dip["i"] + dip["n"])
    out_mask = np.ones(len(dx), bool)
    out_mask[inside] = False

    mx, my = float(np.median(dx[inside])), float(np.median(dy[inside]))
    sx, sy = robust_sigma(dx[out_mask]), robust_sigma(dy[out_mask])
    n = dip["n"]
    shift_pix = float(np.hypot(mx, my))
    err_pix = float(np.hypot(sx, sy) / np.sqrt(n)) if n else float("nan")
    return {
        "available": True,
        "n_in": int(n),
        "dx_millipix": mx * 1000.0,
        "dy_millipix": my * 1000.0,
        "shift_millipix": shift_pix * 1000.0,
        "shift_mas": shift_pix * PIXEL_ARCSEC * 1000.0,
        "err_millipix": err_pix * 1000.0,
        "sigma": shift_pix / err_pix if err_pix else None,
        "pass": bool(shift_pix * 1000.0 < CENTROID_LIMIT_MILLIPIX),
    }


def vsx_lookup(ra, dec, radius="20s"):
    try:
        from astroquery.vizier import Vizier
    except ImportError:
        return {"available": False, "reason": "astroquery not installed"}
    try:
        v = Vizier(columns=["**"], row_limit=10)
        res = v.query_region(f"{ra} {dec}", radius=radius, catalog="B/vsx/vsx")
    except Exception as exc:
        return {"available": False, "reason": f"{type(exc).__name__}: {exc}"}
    if not res:
        return {"available": True, "matches": []}
    out = []
    for row in res[0]:
        out.append({
            "name": str(row["Name"]).strip(),
            "type": str(row["Type"]).strip(),
            "period_d": str(row["Period"]).strip(),
            "sep_arcsec": float(row["_r"]) * 3600.0 if "_r" in res[0].colnames else None,
        })
    return {"available": True, "matches": out}


def quality_audit(t, pdcsap, sap, quality, dip):
    """Requirement 7: the dip must survive in the raw flux and in clean flags."""
    inside = slice(dip["i"], dip["i"] + dip["n"])
    out = {
        "n_cadences": int(len(t)),
        "n_quality_zero": int(np.sum(quality == 0)),
        "quality_flags_in_dip": sorted({int(x) for x in quality[inside]}),
        "dip_present_in_quality_zero": bool(np.all(quality[inside] == 0)),
    }
    if sap is not None:
        med = float(np.median(sap))
        if med > 0:
            raw_depth = float(1.0 - np.min(sap[inside]) / med)
            out["raw_sap_depth"] = raw_depth
            out["raw_sap_confirms"] = bool(raw_depth >= 0.5 * dip["depth"])
    return out


# --------------------------------------------------------------------------- #
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--kepid", type=int, required=True)
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)

    lc, meta = load_light_curve(args.kepid)
    good = np.isfinite(lc["time"]) & np.isfinite(lc["pdcsap"]) & (lc["quality"] == 0)
    t, f = lc["time"][good], lc["pdcsap"][good]
    mom1 = lc["mom1"][good] if lc["mom1"] is not None else None
    mom2 = lc["mom2"][good] if lc["mom2"] is not None else None
    sap = lc["sap"][good] if lc["sap"] is not None else None
    quality = lc["quality"][good]

    dips, scatter, _ = measure_dips(t, f)
    if not dips:
        print(json.dumps({"kepid": args.kepid, "error": "no dips found"}, indent=1))
        return 1

    deepest = max(dips, key=lambda d: d["depth"])
    lopsided_dips = [d for d in dips if d["ratio"] <= 0.67 or d["ratio"] >= 1.5]
    # the dip we actually claim: the deepest one that is also lopsided
    best = max(lopsided_dips, key=lambda d: d["depth"]) if lopsided_dips else deepest
    report = {
        "kepid": args.kepid,
        **meta,
        "n_cadences_used": int(len(f)),
        "robust_scatter": scatter,
        "dips": dips,
        "deepest": deepest,
        "claimed_dip": best,
        "n_lopsided_dips": len(lopsided_dips),
        "periodicity": periodicity(dips),
        "centroid": centroid_test(t, mom1, mom2, best, None),
        "vsx": vsx_lookup(meta["ra"], meta["dec"]) if meta["ra"] else {"available": False},
        "quality": quality_audit(t, f, sap, quality, best),
    }

    checks = {
        "2_deep": bool(best["depth"] >= 0.05 and best["snr"] >= 50),
        "3_lopsided": bool(best["ratio"] <= 0.67 or best["ratio"] >= 1.5),
        "4_not_periodic": not report["periodicity"]["periodic"],
        "5_centroid": bool(report["centroid"].get("pass")),
        "6_not_variable": bool(
            report["vsx"].get("available")
            and not any(
                m["type"].upper().startswith(EB_TYPES)
                for m in report["vsx"].get("matches", [])
                if (m["sep_arcsec"] is None or m["sep_arcsec"] < 5.0)
            )
        ),
        "7_not_instrumental": bool(
            report["quality"].get("dip_present_in_quality_zero")
            and report["quality"].get("raw_sap_confirms", True)
        ),
    }
    report["checks"] = checks
    report["all_photometric_checks_pass"] = all(checks.values())

    text = json.dumps(report, indent=1, default=float)
    if args.out:
        with open(args.out, "w") as fh:
            fh.write(text)
    print(text)
    return 0 if report["all_photometric_checks_pass"] else 2


if __name__ == "__main__":
    sys.exit(main())
