# Kepler aperiodic-dipper hunter

A small, reproducible pipeline for finding **Tabby's-Star-like** objects in the
public Kepler archive: stars that dim by several percent in **deep, lopsided,
non-periodic** dips that no planet or eclipsing binary can explain.

Built for the DAiL discovery bounty *"Find a NEW Tabby's Star"*. Everything here
is plain Python, runs against public MAST data, and needs no credentials.

---

## Method

### Stage 1 — photometric search

`kepler_dipper_hunt.py` walks the Kepler long-cadence archive:

1. `manifest.py` resolves every Kepler timeseries observation to its
   `kplr<kepid>-<obsid>_llc.fits` product URI on MAST (resumable, chunked).
2. For each star the light curve is downloaded and filtered to
   `SAP_QUALITY == 0` and finite `PDCSAP_FLUX`.
3. Flux is normalised by its median and divided by a 61-cadence running median
   to remove slow instrumental drift.
4. Contiguous runs of at least 3 cadences below `-min_depth` become dip
   candidates. Depth is the median residual of the run; SNR uses the
   out-of-dip robust scatter scaled by `sqrt(n_cadences)`.
5. A candidate is **rejected** if its dips phase-fold to a common period
   (planet / eclipsing binary) or if no dip is lopsided at half depth
   (ingress/egress ratio outside `[0.67, 1.50]`).

Survivors are written to `cand.jsonl`; every processed star is appended to
`seen.log` so runs are resumable.

### Stage 2 — vetting

`vet_candidate.py --kepid <id>` produces a JSON report with one section per
vetting requirement:

| # | Requirement | How it is measured |
|---|---|---|
| 1 | NEW | manual literature check (see report) |
| 2 | DEEP | `deepest.depth >= 0.05` and `deepest.snr >= 50` |
| 3 | LOPSIDED | `deepest.ratio <= 0.67` or `>= 1.5` at half depth |
| 4 | NOT PERIODIC | dip centres do not share a common spacing |
| 5 | STABLE CENTROID | `MOM_CENTR1/2` shift during the dip `< 1.0 millipixel` |
| 6 | NOT A KNOWN VARIABLE | VSX cone search at the source position |
| 7 | NOT INSTRUMENTAL | dip present with `SAP_QUALITY == 0` and in raw `SAP_FLUX` |
| 8 | REPRODUCIBLE | this repository + the run log |

The centroid test de-trends each `MOM_CENTR` column with the same 61-cadence
median filter used for the flux, then compares the median offset inside the dip
with the out-of-dip scatter. Kepler's plate scale is 3.98 arcsec/pixel, so
1.0 millipixel = 3.98 mas.

---

## Install

```bash
python3 -m venv venv
./venv/bin/pip install -r requirements.txt
```

## Run

```bash
export HUNT_HOME=$PWD

# 1. build the kepid -> MAST URI manifest (resumable, ~210k rows)
python3 manifest.py

# 2. scan N stars
python3 kepler_dipper_hunt.py run --n 20000 --workers 20 --min-depth 0.03

# 3. vet a survivor
python3 vet_candidate.py --kepid <kepid> --out report.json
```

`--min-depth 0.03` is the search threshold. The report then enforces the
bounty's 5 % floor on the deepest dip, so a 3 % search still finds 5 % events
while giving the lopsidedness test more events to work with.

### Throughput

On a 2-core / 7.6 GB VPS with 20 download threads the scan runs at roughly
7 stars/second (~25k stars/hour). The full Kepler long-cadence archive is
~213k targets.

---

## Output format

`cand.jsonl`, one object per candidate:

```json
{
  "kepid": 12345678,
  "ra": 291.234, "dec": 44.567, "kepmag": 14.2,
  "n_cad": 4300, "scatter": 0.00031,
  "n_dips": 3, "n_lopsided": 2, "periodic": false,
  "max_depth": 0.081,
  "best": {"depth": 0.081, "snr": 610, "dur_d": 0.42, "ratio": 0.31},
  "dips": [ ... ]
}
```

## Calibration on the known positive

The pipeline is validated against the one object whose behaviour is not in
doubt: **KIC 8462852 (Tabby's Star), Quarter 16**, the 22 %-deep event at
BJD 1519.

```
depth 11.68 %   SNR 4636
ratio 1.67 (tenth maximum; half maximum quantises to 1.00)
periodicity: none
centroid shift 0.014 millipixel   (limit 1.0)
VSX: no match
raw SAP depth 20.5 %  ->  dip present in undetrended flux
=> all six photometric checks pass
```

The ratio column is the interesting one. At half maximum both wings are
shorter than one 29.4-minute cadence, so the measurement quantises to
0.5 / 0.5 and reads as perfectly symmetric. Measuring at **tenth** maximum
resolves the wings and recovers the real 1.67x asymmetry. A pipeline that
only looks at half maximum rejects Tabby's Star itself.

## Notes and limitations

* Kepler long-cadence sampling is 29.4 minutes, so dips shorter than ~1.5 h
  are resolved by only a few cadences and their ingress/egress ratio is noisy.
* The 61-cadence running median (~1.25 d) removes astrophysical variability
  slower than that timescale, which is what makes day-scale dips visible but
  suppresses week-long ones.
* Blends are rejected by the centroid test, not by aperture photometry; a
  contaminant that is exactly centred on the target would survive.
* VSX coverage is good but not complete — a `vsx.matches: []` result means
  "not catalogued", not "proven new".
