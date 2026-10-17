"""Download a specific Kepler quarter for a KIC id (for calibration)."""
import os
import subprocess
import sys

BASE = os.environ.get("HUNT_HOME", "/opt/agent/repo")
CACHE = os.path.join(BASE, "fits_cache")

kepid = int(sys.argv[1])
quarter = int(sys.argv[2])

from astroquery.mast import Observations

obs = Observations.query_criteria(
    obs_collection="Kepler", objectname=f"KIC {kepid}", dataproduct_type="timeseries"
)
prods = Observations.get_product_list(obs)
rows = []
for i, uri in enumerate(prods["dataURI"]):
    if not uri or "_llc.fits" not in uri:
        continue
    oid = str(prods["obs_id"][i])
    if not oid.startswith(f"kplr{kepid:09d}"):
        continue
    rows.append((uri, oid))
print("llc products for target:", len(rows))

# obs_id format kplr008462852-2009166043257 ; quarter is in the product's
# t_min/t_max, so just try each until the quarter header matches
os.makedirs(CACHE, exist_ok=True)
target = os.path.join(CACHE, f"{kepid}_q{quarter}.fits")
if os.path.exists(target) and os.path.getsize(target) > 50000:
    print("cached", target)
    sys.exit()

from astropy.io import fits

for uri, obs_id in rows:
    tmp = target + ".part"
    url = f"https://mast.stsci.edu/api/v0.1/Download/file?uri={uri}"
    subprocess.run(["curl", "-sL", "-m", "60", "-o", tmp, url], capture_output=True)
    if not os.path.exists(tmp) or os.path.getsize(tmp) < 50000:
        continue
    try:
        with fits.open(tmp, memmap=False) as h:
            qh = h[0].header.get("QUARTER")
    except Exception:
        os.remove(tmp)
        continue
    if qh == quarter:
        os.replace(tmp, target)
        print("saved", target, os.path.getsize(target), "B from", obs_id)
        sys.exit()
    os.remove(tmp)

print("quarter not found")
