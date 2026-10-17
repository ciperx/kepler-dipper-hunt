#!/usr/bin/env python3
"""Build manifest.tsv: kepid -> MAST dataURI for Kepler long-cadence light curves.

Queries MAST once for the whole Kepler timeseries collection, then pages
through the product list in chunks so a 200k-row query does not time out.
The result is resumable: rerunning continues from manifest.state.

Usage
-----
  python3 manifest.py            # build / resume
  python3 manifest.py --reset    # start over
"""
from __future__ import annotations

import argparse
import os
import warnings

warnings.filterwarnings("ignore")

BASE = os.environ.get("HUNT_HOME", os.path.abspath("."))
OUT = os.path.join(BASE, "manifest.tsv")
STATE = os.path.join(BASE, "manifest.state")
CHUNK = 2000


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reset", action="store_true")
    args = ap.parse_args()

    if args.reset:
        for p in (OUT, STATE):
            if os.path.exists(p):
                os.remove(p)

    from astroquery.mast import Observations

    obs = Observations.query_criteria(
        obs_collection="Kepler", dataproduct_type="timeseries"
    )
    total = len(obs)
    start = int(open(STATE).read()) if os.path.exists(STATE) else 0
    print(f"{total} observations, resuming at {start}", flush=True)

    seen: set[int] = set()
    if os.path.exists(OUT):
        with open(OUT) as fh:
            for line in fh:
                try:
                    seen.add(int(line.split("\t")[0]))
                except ValueError:
                    continue
    print(f"{len(seen)} kepids already mapped", flush=True)

    with open(OUT, "a") as fh:
        for i in range(start, total, CHUNK):
            chunk = obs[i : i + CHUNK]
            try:
                products = Observations.get_product_list(chunk)
            except Exception as exc:
                print(f"chunk {i}: {type(exc).__name__}: {exc}", flush=True)
                continue
            if products is None or len(products) == 0:
                continue
            if "dataURI" not in set(products.colnames):
                print("unexpected product columns:", products.colnames, flush=True)
                return

            added = 0
            for row in products:
                uri = row["dataURI"] or ""
                if "_llc.fits" not in uri:      # long cadence only
                    continue
                obs_id = str(row["obs_id"])
                if not obs_id.startswith("kplr"):
                    continue
                try:
                    kepid = int(obs_id[4:13])
                except ValueError:
                    continue
                if kepid in seen:
                    continue
                seen.add(kepid)
                fh.write(f"{kepid}\t{uri}\n")
                added += 1
            fh.flush()
            with open(STATE, "w") as sh:
                sh.write(str(i + CHUNK))
            print(f"chunk {i}: +{added} (total {len(seen)})", flush=True)


if __name__ == "__main__":
    main()
