"""Vet every candidate in cand4.jsonl and emit a table for the report."""
import json
import os
import subprocess
import sys

BASE = "/opt/agent/astro"
REPO = "/opt/agent/repo"
PY = f"{BASE}/venv/bin/python"

rows = {}
with open(f"{BASE}/manifest.tsv") as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) == 2:
            try:
                rows[int(p[0])] = p[1]
            except ValueError:
                pass

cands = []
seen = set()
for line in open(f"{BASE}/cand4.jsonl"):
    try:
        d = json.loads(line)
    except json.JSONDecodeError:
        continue
    if d["kepid"] in seen:
        continue
    seen.add(d["kepid"])
    cands.append(d)

print(f"{len(cands)} unique candidates", flush=True)
out = []
for i, c in enumerate(cands, 1):
    k = c["kepid"]
    fits = f"{REPO}/fits_cache/{k}.fits"
    if not os.path.exists(fits):
        uri = rows.get(k)
        if uri:
            subprocess.run(
                ["curl", "-sL", "-m", "90", "-o", fits,
                 f"https://mast.stsci.edu/api/v0.1/Download/file?uri={uri}"],
                capture_output=True,
            )
    if not os.path.exists(fits) or os.path.getsize(fits) < 50000:
        print(f"{i:>2} {k} FETCH_FAIL", flush=True)
        continue

    jf = f"/tmp/vv_{k}.json"
    subprocess.run(
        [PY, f"{REPO}/vet_candidate.py", "--kepid", str(k), "--out", jf],
        capture_output=True, timeout=300,
    )
    try:
        d = json.load(open(jf))
    except Exception:
        print(f"{i:>2} {k} VET_FAIL", flush=True)
        continue
    if "error" in d:
        print(f"{i:>2} {k} {d['error']}", flush=True)
        continue

    b = d["claimed_dip"]
    failed = [kk for kk, v in d["checks"].items() if not v]
    vsx = [(m["name"], m["type"]) for m in d["vsx"]["matches"]][:2]
    rec = dict(kepid=k, depth=b["depth"], ratio=b["ratio"], snr=b["snr"],
               ndips=len(d["dips"]), periodic=d["periodicity"]["periodic"],
               cen=d["centroid"]["shift_millipix"], vsx=vsx,
               failed=failed, passed=d["all_photometric_checks_pass"])
    out.append(rec)
    print(f"{i:>2} {k} d={b['depth']:.3f} r={b['ratio']:.2f} snr={b['snr']:.0f} "
          f"nd={len(d['dips'])} cen={d['centroid']['shift_millipix']:.2f} "
          f"vsx={vsx} fail={failed}", flush=True)

json.dump(out, open(f"{BASE}/vet_all.json", "w"), indent=1)
print("wrote", f"{BASE}/vet_all.json", flush=True)
