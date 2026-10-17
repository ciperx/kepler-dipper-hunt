import numpy as np
from astropy.io import fits

p = "/opt/agent/repo/fits_cache/8462852.fits"
with fits.open(p, memmap=False) as h:
    print("HDUs:", len(h))
    hdr = h[0].header
    print("QUARTER:", hdr.get("QUARTER"), "MODULE:", hdr.get("MODULE"), "OUTPUT:", hdr.get("OUTPUT"))
    print("KEPMAG:", hdr.get("KEPMAG"), "RA:", hdr.get("RA_OBJ"), "DEC:", hdr.get("DEC_OBJ"))
    d = h[1].data
    print("cols:", d.columns.names)
    t = np.asarray(d["TIME"], float)
    f = np.asarray(d["PDCSAP_FLUX"], float)
    q = np.asarray(d["SAP_QUALITY"], int)
    print("n:", len(t), "t range:", np.nanmin(t), np.nanmax(t))
    print("finite flux:", int(np.isfinite(f).sum()), "q==0:", int((q == 0).sum()))
    med = np.nanmedian(f)
    print("median flux:", med)
    g = np.isfinite(f) & (q == 0)
    print("min rel flux (q0):", np.nanmin(f[g]) / med - 1)
    print("min rel flux (all):", np.nanmin(f[np.isfinite(f)]) / med - 1)
