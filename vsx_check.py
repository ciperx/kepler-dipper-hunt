"""Cross-check a candidate against VSX with explicit coordinates and separations."""
import sys

from astroquery.vizier import Vizier

ra, dec = float(sys.argv[1]), float(sys.argv[2])
print(f"target RA={ra} DEC={dec}")

v = Vizier(columns=["**"], row_limit=20)
res = v.query_region(f"{ra} {dec}", radius="60s", catalog="B/vsx/vsx")
if not res:
    print("no VSX object within 60 arcsec")
    sys.exit()

t = res[0]
print("columns:", t.colnames)
for row in t:
    name = str(row["Name"]).strip()
    typ = str(row["Type"]).strip()
    per = str(row["Period"]).strip()
    vra = row["RAJ2000"] if "RAJ2000" in t.colnames else None
    vdec = row["DEJ2000"] if "DEJ2000" in t.colnames else None
    sep = None
    if vra is not None and vdec is not None:
        import numpy as np
        sep = float(np.hypot((float(vra) - ra) * np.cos(np.radians(dec)),
                             float(vdec) - dec) * 3600.0)
    print(f"  {name:<20} type={typ:<16} P={per:<10} RA={vra} DEC={vdec} sep={sep} arcsec")
