"""Sanity-check the scanner's centroid filter on one known-positive and one blend."""
import sys

sys.path.insert(0, "/opt/agent/astro")
import hunt4

hunt4.MIN_DEPTH = 0.05
hunt4.USE_VSX = False          # isolate the centroid logic

rows = {}
with open("/opt/agent/astro/manifest.tsv") as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) == 2:
            rows[int(p[0])] = p[1]

for kepid in (8462852, 9775475, 9896435, 9702891):
    uri = rows.get(kepid)
    if not uri:
        print(kepid, "not in manifest")
        continue
    r = hunt4.process((kepid, uri))
    if r is None:
        print(kepid, "no dips at min_depth 5%")
    elif r.get("reject"):
        print(kepid, "REJECT", r.get("why"))
    else:
        print(kepid, f"PASS depth={r['max_depth']:.3f} cen={r['centroid_millipix']:.3f} mp")
