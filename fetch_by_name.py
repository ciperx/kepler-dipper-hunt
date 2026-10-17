"""Resolve the MAST dataURI for one KIC id and download it into fits_cache."""
import os
import sys

BASE = os.environ.get("HUNT_HOME", "/opt/agent/repo")
CACHE = os.path.join(BASE, "fits_cache")

kepid = int(sys.argv[1])
os.makedirs(CACHE, exist_ok=True)
out = os.path.join(CACHE, f"{kepid}.fits")

from astroquery.mast import Observations

obs = Observations.query_criteria(
    obs_collection="Kepler", target_name=f"KIC {kepid}", dataproduct_type="timeseries"
)
print("obs rows:", 0 if obs is None else len(obs))
if obs is None or len(obs) == 0:
    obs = Observations.query_criteria(
        obs_collection="Kepler", objectname=f"KIC {kepid}", dataproduct_type="timeseries"
    )
    print("retry rows:", 0 if obs is None else len(obs))
if obs is None or len(obs) == 0:
    sys.exit("no observations")

prods = Observations.get_product_list(obs)
uris = [u for u in prods["dataURI"] if u and "_llc.fits" in u]
print("llc products:", len(uris))
if not uris:
    sys.exit("no llc product")

url = f"https://mast.stsci.edu/api/v0.1/Download/file?uri={uris[0]}"
import subprocess

r = subprocess.run(["curl", "-sL", "-m", "90", "-o", out + ".part", url])
if r.returncode == 0 and os.path.getsize(out + ".part") > 50000:
    os.replace(out + ".part", out)
    print("saved", out, os.path.getsize(out), "B")
else:
    print("download failed")
