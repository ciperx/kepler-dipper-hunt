#!/usr/bin/env python3
"""Push a local directory to a GitHub repo using the Contents API (no gh CLI)."""
import base64
import json
import os
import sys
import urllib.error
import urllib.request

TOKEN = open(
    "/opt/data/profiles/idea/cache/scratch/autonomous-income/secrets/gh_token.txt"
).read().strip()
OWNER = "ciperx"
REPO = "kepler-dipper-hunt"
BRANCH = "main"
API = "https://api.github.com"


def call(method, path, payload=None):
    url = f"{API}/repos/{OWNER}/{REPO}/{path}" if path else f"{API}/repos/{OWNER}/{REPO}"
    data = json.dumps(payload).encode() if payload is not None else None
    req = urllib.request.Request(url, data=data, method=method)
    req.add_header("Authorization", f"token {TOKEN}")
    req.add_header("Accept", "application/vnd.github+json")
    if data:
        req.add_header("Content-Type", "application/json")
    try:
        with urllib.request.urlopen(req, timeout=40) as r:
            return r.status, json.load(r)
    except urllib.error.HTTPError as e:
        return e.code, json.loads(e.read() or b"{}")


def sha_of(path):
    code, body = call("GET", f"contents/{path}?ref={BRANCH}")
    if code == 200:
        return body["sha"]
    return None


def put(path, content, message):
    payload = {
        "message": message,
        "content": base64.b64encode(content.encode()).decode(),
        "branch": BRANCH,
    }
    existing = sha_of(path)
    if existing:
        payload["sha"] = existing
    code, body = call("PUT", f"contents/{path}", payload)
    status = "updated" if existing else "created"
    if code in (200, 201):
        print(f"  {status}: {path}")
    else:
        print(f"  FAILED {path}: {code} {body.get('message')}")
    return code in (200, 201)


def main(local_dir):
    ok = True
    for name in sorted(os.listdir(local_dir)):
        full = os.path.join(local_dir, name)
        if not os.path.isfile(full):
            continue
        if name.startswith(".") or name.endswith(".pyc"):
            continue
        with open(full, encoding="utf-8", errors="replace") as fh:
            content = fh.read()
        ok &= put(name, content, f"Add {name}")
    print("done, all ok:", ok)


if __name__ == "__main__":
    main(sys.argv[1])
