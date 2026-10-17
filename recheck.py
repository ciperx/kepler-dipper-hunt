"""Re-run a list of KIC ids through hunt4's detector with the current logic.

Used after a filter change: stars rejected by the old (half-maximum only)
asymmetry test are re-evaluated with the half/tenth-maximum test.

  recheck.py /tmp/sym.txt
"""
import concurrent.futures as cf
import json
import sys

sys.path.insert(0, "/opt/agent/astro")
import hunt4

hunt4.MIN_DEPTH = 0.05
hunt4.USE_VSX = True

BASE = "/opt/agent/astro"

rows = []
with open(f"{BASE}/manifest.tsv") as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) == 2:
            try:
                rows.append((int(p[0]), p[1]))
            except ValueError:
                pass
index = {k: u for k, u in rows}

want = [int(x) for x in open(sys.argv[1]).read().split() if x.strip()]
todo = [(k, index[k]) for k in want if k in index]
print(f"rechecking {len(todo)} of {len(want)} requested", flush=True)

hits = 0
with cf.ThreadPoolExecutor(12) as ex:
    for fut in cf.as_completed([ex.submit(hunt4.process, r) for r in todo]):
        try:
            res = fut.result()
        except Exception as exc:
            print("  exc", type(exc).__name__, flush=True)
            continue
        if not res:
            continue
        if res.get("reject"):
            print(f"  rej {res['kepid']} {res.get('why')}", flush=True)
            continue
        if "error" in res:
            continue
        hits += 1
        b = res["best"]
        print(f"  HIT {res['kepid']} depth={b['depth']:.3f} ratio={b['ratio']:.2f} "
              f"({b.get('ratio_basis','?')}) nd={res['n_dips']}", flush=True)
        with open(f"{BASE}/cand4.jsonl", "a") as fh:
            fh.write(json.dumps(res) + "\n")
print(f"done, {hits} new candidates", flush=True)
