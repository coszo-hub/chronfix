
End-to-end pipeline that detects the HYS14 (OOI Hydrate Ridge OBS) clock
error from ambient-noise cross-correlation and applies a per-sample
timestamp correction to the raw MiniSEED data, with closed-loop
validation. Companion document to HYS14 — Timing Diagnostic [[05-01-26 Notes]]
(diagnosis-only) and `Chronos — Implementation Plan.md` (project
context).

---

## 1. Pipeline overview

```
   raw MiniSEED                     chronos                                chronfix
   (EarthScope)        ──>       (clock measurement)         ──>       (clock correction)        ──> corrected MiniSEED
                                  diagnostics/                          chronfix/

  scripts/download_hys.py   →   diagnostics/hys_ccf.py
                                diagnostics/peak_lag_hourly.py
                                diagnostics/combine_clock_hourly.py
                                                                 chronfix/scripts/filter_and_triggers.py
                                                                 chronfix/clock_model.py
                                                                 chronfix/correct.py
                                                                 chronfix/scripts/correct_hys14.py
```

Two packages:

- **chronos** — measures Δt(t), the clock error of HYS14 vs UTC, from
  cross-correlations against a reference station with known good timing
  (HYS12).
- **chronfix** — consumes a segment-modeled hourly Δt time series and a list of
  resync events ("triggers"), and applies the timestamp shift to raw
  HYS14 MiniSEED, splitting the output at each trigger.

---

## 2. Data acquisition

`scripts/download_hys.py` pulls daily MHZ (8 Hz vertical) MiniSEED for
`OO.HYS{12,14,B1}` from EarthScope FDSN over 2022-01-01 → present, plus
the corresponding StationXML response files. Files land in the layout
the diagnostic pipeline expects:

```
/data/wsd02/maleen_data/OOI-Data/{sta}/{yr}/{doy:03d}/{sta}.OO.{yr}.{doy:03d}.MHZ
/data/wsd02/maleen_data/OOI-Data/StationXML/OO.{sta}..MHZ.xml
```

Downloads are idempotent — a re-run skips files already on disk.

---

## 3. Cross-correlation (chronos)

`diagnostics/hys_ccf.py` computes daily inter-station ZZ
cross-correlations for the HYS12-HYS14 pair, using the recipe adapted
from Earthnote's single-station SC pipeline.

### Per-day preprocessing

1. Read daily mseed; cast to float64; merge with `merge(method=1, fill_value=0.0)`.
2. Cosine-taper around any zero-fill gaps (100 s on each side).
3. `detrend("demean")`, `detrend("linear")`.
4. High-pass at 0.4 Hz; remove instrument response to velocity with a
   Nyquist-clamped pre-filter.
5. Resample / interpolate onto a uniform 8 Hz grid.
6. Trim/pad to span exactly the canonical UTC day so sample 0 of each
   trace always corresponds to the same UTC second across stations.

### Cross-correlation per day

7. 30-min windows with 7.5-min step (75 % overlap), per-window cosine
   edge taper, raised-cosine phase whitening over 0.5–3.8 Hz, one-bit
   amplitude normalisation.
8. Frequency-domain CC `irfft(conj(FFT(a)) * FFT(b))`, truncated to ±60 s.
9. Per-day median across all 30-min windows ⇒ robust daily stack.

Outputs under `data/ccf/HYS12-HYS14/`:

```
cc_30min.npy        (N_seg, n_lags)     all 30-min CCs concatenated
cc_30min_times.npy  fractional day index of each segment midpoint
cc_daily.npy        (N_days, n_lags)    per-lag median per day
cc_dates.npy        datetime64[D] day axis
cc_ref.npy          long-term mean reference stack
lags.npy            lag axis in seconds
```

---

## 4. Hourly peak-lag tracking (chronos)

`diagnostics/peak_lag_hourly.py` re-stacks the per-30-min CCs into
hourly bins (median across the ~8 windows whose midpoints fall in each
hour) and computes the lag of the global maximum of the envelope
of `CC²` per hour. Hilbert-envelope-of-squared-CC sharpens the
ballistic peak; the result is one peak-lag estimate per UTC hour.

For HYS12-HYS14 the post-fix anchor lag is ≈ 0 s (the stations are
near-co-located), so the peak-lag track is a direct reading of the
HYS14 clock error trajectory: drift episodes followed by sharp
resyncs back to ≈ 0.

---

## 5. Combining into a canonical hourly Δt (chronos)

`diagnostics/combine_clock_hourly.py` converts the peak-lag track to a
clock-error time series:

- Anchor lag = late-window median of valid hourly picks (≈ 0 s for
  HYS12-HYS14, after the post-fix stable period).
- Per-hour shift = `anchor − peak_lag`.
- Sign convention: HYS14 is the B-side of HYS12-HYS14, so
  `Δt_HYS14 = −shift`. Δt > 0 means HYS14 timestamps are *late*
  relative to true UTC.
- Cross-validation: the same procedure on HYS14-HYSB1 (lowband, 0.1–0.3
  Hz) is reported as a residual but does not contribute to the canonical
  estimate. HYS12-HYS14 has higher SNR and better coverage.

Outputs under `data/clock_estimate/HYS14/`:

- `delta_t_hourly.npy` — primary hourly Δt (s), NaN where unknown
- `hour_times.npy` — datetime64[h] master axis
- `dt_hourly_from_HYS12_HYS14.npy` — same series, stored separately
- `dt_hourly_from_HYS14_HYSB1.npy` — cross-check series
- `residual_hourly.npy` — primary − cross-check where both valid

---

## 6. Outlier filtering and trigger detection (chronfix)

`chronfix/scripts/filter_and_triggers.py` is the single bridge that
turns the chronos hourly Δt into a chronfix-ready clock model. It
combines outlier rejection and resync detection in one pass.

### Stage 1 — strict multi-pass Hampel filter

Three rolling Hampel passes are applied in sequence. Each pass marks
outliers and replaces them with NaN before the next runs, so a wider-
window pass isn't biased by single-hour artifacts:

| Pass | Window | σ-threshold | Min |residual| |
|---|---|---|---|
| 1 | 168 h (7 d) | 4.5 | 0.75 s |
| 2 | 72 h (3 d)  | 4.0 | 0.60 s |
| 3 | 25 h (~1 d) | 3.5 | 0.50 s |

A point is flagged when `|x − rolling_median(x)| > max(σ_threshold ×
1.4826 × rolling_MAD, min_abs_deviation)`.

### Stage 1b — 3-point continuity check

For each hour, compare the value with the median of its two immediate
neighbours. Flag if the residual exceeds 3 s. Catches isolated single-
hour spikes the wider Hampel passes can miss.

### Stage 2 — final "bit less-strict" 10-day pass

One more Hampel pass at a 240-hour window with `σ=6.5` and
`min_abs=1.30` removes any remaining moderate outliers without
re-flagging the legitimate drift ramp interiors.

For the HYS12-HYS14 record this filter masks **673 outlier hours** out
of 31,169 valid (2.2 %), leaving 30,496 retained samples. The retained
series is preserved as `delta_t_hourly_filtered_raw.npy` (used only
for QC); the file chronfix actually consumes is the segment-modeled
series produced in Stage 3 below.

### Trigger detection

After filtering, compute the difference between consecutive **retained**
hourly samples (NaNs are skipped, so a retained-sample-to-retained-
sample diff naturally ignores intervening masked hours):

```
Δdt_k = retained[k] − retained[k − 1]
```

If `|Δdt_k| > 1.0 s`, the interval `[t(k−1), t(k)]` is flagged as a
**trigger** — a region where Δt changed too fast to be the result of
ordinary slow drift. Touching/overlapping intervals are merged. For the
HYS12-HYS14 record, **32 trigger periods** are detected.

### Stage 3 — per-segment robust smoothing of Δt

The Hampel-filtered series above is a picker-quantum staircase: the
0.125 s picker quantum (one sample at fs = 8 Hz) is small relative to
the drift signal but not relative to one CCF window's coherence budget,
and chronfix's `_resample` linearly interpolates whatever it sees,
turning the staircase into sub-second false time-warp. To eliminate
that, each inter-trigger segment of Δt is replaced with a robust
smooth model before being written as the file chronfix consumes:

- **Slope test (Theil-Sen).** Robust slope of the segment's finite
  Δt samples. If `|slope| < 0.05 s/day`, the segment is treated as
  stable.
- **Stable segments → robust median.** Replace the entire segment
  with `nanmedian(seg)`. Constant; no sub-quantum artifacts; no
  end-effect wobble.
- **Drift segments → rolling robust median + light moving average.**
  Centered rolling robust median with a 24-hour window, then a
  centered 6-hour moving average on top of that. No assumption that
  the drift is linear, so curved drifts (e.g. the 2022-08→2023-01
  ramp, the 2023-06→2023-09 ramp) are tracked rather than
  approximated by a straight line.

Trigger intervals stay NaN; chronfix splits its corrected output at
those boundaries. Segments with fewer than 5 finite samples are
skipped (left NaN).

Validated on the HYS12-HYS14 record: per-segment residual = (raw
cleaned Δt) − (modeled Δt) has median ≈ 0.000 s and MAD ≤ 0.115 s on
every segment, including the longest ones (126 days / 56.7 s of total
drift). No detectable timing signal remains in the residuals.

Outputs:

- `data/clock_estimate/HYS14/delta_t_hourly_clean.npy` — **modeled** Δt
  (this is the file chronfix consumes)
- `data/clock_estimate/HYS14/delta_t_hourly_filtered_raw.npy` — raw
  Hampel-filtered Δt (QC only, not consumed by chronfix)
- `data/clock_estimate/HYS14/delta_t_hourly_outlier_mask.npy` — boolean mask
- `data/clock_estimate/HYS14/trigger_periods.csv` — merged trigger intervals
- `data/clock_estimate/HYS14/filter_and_triggers.png` — 2-panel diagnostic

---

## 7. Clock model (chronfix)

`chronfix.clock_model.ClockModel` wraps the filter outputs into a query
interface used by the correction step:

- **`interp_delta_t(t)`** — linear interpolation of Δt at any UTC time
  using the segment-modeled hourly samples (see Stage 3 of section 6).
  Returns NaN if `t` falls inside any trigger interval (Δt is undefined
  there).
- **`stable_intervals(t0, t1)`** — UTC ranges between consecutive
  triggers within `[t0, t1]`. The correction processes one stable
  interval at a time and never spans a trigger.

Δt within a stable interval is treated as smoothly drifting between
hourly samples. The hourly samples themselves are no longer the raw
picker output — they are the segment-modeled rolling-median + MA
series, so the linear interpolation chronfix performs is interpolating
between points on a smooth curve rather than between staircase
quanta. This matters because a 0.125 s picker quantum interpolated
linearly across one hour produces sub-second false time-warp inside
each CCF window, which is what the 2026-05-04 fix eliminates.

---

## 8. Per-sample correction (chronfix)

`chronfix.correct.correct_trace` applies the clock model to one trace.

For each input trace:

1. Identify all stable intervals overlapping the trace's apparent time
   range.
2. For each overlap, slice the input trace to the overlap and either:
   - **`method="resample"`** (default): build a regular UTC sample grid
     covering the corrected interval; for each output sample at UTC
     time `t_utc`, compute `t_apparent ≈ t_utc + Δt(t_utc)` and read
     the input trace value via linear interpolation. Output trace has
     starttime aligned to true UTC and the same nominal sampling rate.
   - **`method="shift_only"`**: subtract `Δt(starttime)` from the
     trace's starttime; data samples are bit-identical. Within-segment
     drift is approximated as a constant offset; useful as a sanity
     check.
3. Drop slices where Δt is NaN (no chronos measurement available).
4. Return one corrected sub-trace per stable-segment overlap.

`chronfix.correct.correct_stream` runs the per-trace logic over a
Stream and returns a flat list of corrected sub-traces.

### Sign convention

Δt > 0 ⇒ HYS14 timestamps are *late* relative to true UTC.
Therefore a sample with apparent timestamp `T` was actually recorded at
true UTC time `T − Δt(T)`. The correction relabels the sample to that
true UTC time. Order of samples is preserved.

### Trigger handling

A trigger interval is treated as a hard cut. The trace is split there
and one separate output mseed is written per stable-segment overlap.
Across-trigger discontinuities are an unavoidable physical consequence
of the clock having been wrong before the resync:

- Δt decreasing at the jump → corrected output has a UTC **gap** equal
  to the jump magnitude.
- Δt increasing at the jump → corrected output has a UTC **overlap**.

chronfix surfaces these gaps/overlaps explicitly between segment
output files. Downstream tools handle them as ordinary gaps or
duplicates.

---

## 9. Daily driver (chronfix)

`chronfix/scripts/correct_hys14.py` walks the HYS14 input tree day by
day. For each day:

1. Read the input mseed.
2. Apply `correct_stream(..., method="resample")`.
3. Write all corrected sub-traces to a single output mseed file under
   `/data/wsd02/maleen_data/OOI-Data-corrected/HYS14/{yr}/{doy:03d}/HYS14.OO.{yr}.{doy:03d}.MHZ`
   (one record per stable segment).
4. Append a manifest row recording input path, output path, segment
   index, true UTC bounds, sample count, and method.

For the HYS14 record (1582 days input), 951 days produced corrected
output (some days were skipped because the input file was missing or
because Δt was unavailable for the entire day after outlier
filtering). 1059 corrected mseed segments were written across those
days.

---

## 10. Closed-loop validation

The diagnostic pipeline is re-run on the corrected HYS14 against the
unchanged HYS12 to confirm the correction succeeded:

```bash
python diagnostics/hys_ccf.py --pairs HYS12-HYS14 --tag corrected --workers 8 \
    --input-root-override HYS14=/data/wsd02/maleen_data/OOI-Data-corrected
python diagnostics/peak_lag_hourly.py --pair HYS12-HYS14_corrected
```

Reading the outputs at `data/ccf/HYS12-HYS14_corrected/` and
`data/peak_lag_hourly/HYS12-HYS14_corrected/`:

| | Before correction | After correction |
|---|---|---|
| Reference CC RMS | 12.8 | **40.6** (≈ 3.2× stronger) |
| Reference shape | wave packet smeared across ±5 s | sharp peak at lag 0 |
| Daily 2D stack | drift bands wandering across lag | tight vertical stripe at lag 0 |
| Hourly peak-lag time series | clear drift ramps + 23 visible resync resets | flat at 0 across all 4 years |

The closed loop confirms that the per-sample correction has eliminated
the clock error visible in the original cross-correlations to the
limit of measurement precision. The remaining sparse off-trend hours
in the corrected peak-lag plot are noise-driven picks at low-SNR hours
(uniformly scattered, no time structure), not residual clock error.

Comparison plot: `data/peak_lag_hourly/HYS12-HYS14_corrected/before_after_ccf.png`.

### Note on the corrected hourly peak-lag scatter

An earlier draft of this doc explained the corrected hourly peak-lag
scatter as the per-hour "noise floor" being unmasked once the dominant
clock-drift signal was removed. **That explanation was wrong.** A
diagnostic pass in 2026-05-04 (see `issue1.md`) showed the corrected
hourly track had genuinely degraded CCF coherence — the ballistic peak
amplitude in single-hour CCFs was halved, and the global picker was
locking onto far noise lobes (±60 s) on a substantial fraction of
hours, including periods that had been clean monotonic ramps in the
uncorrected plot.

**Root cause:** the cleaned hourly Δt fed to chronfix was a picker-
quantum staircase (picker resolution = 0.125 s = 1 sample at 8 Hz, so
hour-to-hour Δt zigzagged by ±0.5 s on top of the slow drift).
chronfix's `_resample` linearly interpolated this staircase, injecting
sub-second false time-warp into every CCF window. Within a 30-min
window that's enough relative timing jitter to scramble phase
coherence at 1–3 Hz — the band the ballistic peak lives in.

**Fix (2026-05-04):**

1. Per-segment robust smoothing of Δt at the chronos clock-model stage
   so chronfix consumes a shape-preserving smooth function instead of
   the staircase. Implemented in
   `chronos/scripts/filter_and_triggers.py::model_segments`: rolling
   robust median (24 h window) + light moving average (6 h) for drift
   segments; robust segment median for flat segments. Driven by the
   existing trigger-detection output — no station-specific knobs.
2. Leading-NaN boundary snap in `chronfix/correct.py::_resample` so
   float-precision overshoots near the apparent-grid boundaries no
   longer drop the entire output. Recovered ~7,800 corrected hours
   that were silently being discarded.

**Validation after the fix** (full-timeline hourly outlier rates):

| | n | \|x\|>5s | \|x\|>20s | \|x\|>40s | med\|x\| |
|---|---|---|---|---|---|
| uncorrected | 31,563 | 34.76 % | 12.76 % | 3.55 % | 1.875 s |
| pre-fix corrected (jittery Δt + boundary bug) | 22,244 | 0.81 % | 0.57 % | 0.30 % | 0.125 s |
| post-fix, linear per-segment fit | 30,084 | 6.33 % | 0.19 % | 0.08 % | 0.375 s |
| **post-fix, rolling-median per-segment (deployed)** | **30,084** | **0.30 %** | **0.20 %** | **0.10 %** | **0.125 s** |

A residuals diagnostic on the deployed model (raw cleaned Δt minus
modeled Δt, per inter-trigger segment) shows median ≈ 0.000 s and
MAD ≤ 0.115 s on every segment, including the longest and most curved
drift segments (e.g. 126 days / 56.7 s of drift). There is no
detectable systematic timing signal left in the residuals — what
remains is consistent with the picker's quantization noise.

The stacked-product evidence (long-term reference RMS up from 12.8 to
40.6, daily-stack heatmap collapsed to a tight vertical stripe at
lag 0) is also stronger than in the pre-fix run, where it had risen
to 33.1 / 2.6×.

Archived versioned figures and arrays (uncorrected, pre-fix,
linear-per-segment, rolling-per-segment) live under
`data/results_archive/` with a per-version README that pins each
iteration's outlier rates and figures.

### Validation bug uncovered and fixed

The first attempt at this validation step produced a corrected hourly
peak-lag track that *still* showed the original drift ramps — i.e. the
correction appeared to have done nothing. Investigation revealed that
`hys_ccf.py`'s `load_day_z` helper returned the trace data as a plain
NumPy array, dropping the trace's UTC starttime in the process. The
existing pipeline used per-station daily input files that always
started at exactly 00:00:00 UTC, so the absent timestamp was a hidden
assumption (sample 0 of HYS12 and sample 0 of HYS14 implicitly aligned
to the same UTC second).

The chronfix-corrected HYS14 daily files break that assumption: their
starttime is offset from midnight by `-Δt(midnight)` (e.g.
`23:59:33.375` on the previous day when Δt ≈ +27 s). Reading them with
the original loader implicitly re-introduced an apparent-vs-true
misalignment of size Δt, which is exactly the same offset we had just
removed via correction — net zero correction in the cross-correlation
output, hence the unchanged drift ramps.

Fix: at the end of `load_day_z`, trim and pad each trace to span
exactly the canonical UTC day before extracting the data array:

```python
day_start = UTCDateTime(d.year, d.month, d.day)
tr.trim(starttime=day_start, endtime=day_start + 86400.0,
        pad=True, fill_value=0.0, nearest_sample=True)
```

This forces sample 0 of every loaded trace to correspond to the same
UTC second across stations, regardless of where the input file
actually begins or ends. Original (midnight-aligned) input files are
unaffected; corrected input files are aligned to the canonical day
(losing the first/last few tens of seconds in exchange). After the
fix, the corrected hourly peak-lag track is flat at zero across the
full record — the validation passes and the table above reflects the
corrected behaviour.

---

## 11. Reproducing the full pipeline

From `/home/seismic/chronos/`, with conda env `noisepy2` (obspy ≥ 1.5):

```bash
# 1. Acquire data
python scripts/download_hys.py --workers 4

# 2. Daily cross-correlation
python diagnostics/hys_ccf.py --pairs HYS12-HYS14 --workers 8

# 3. Hourly peak-lag
python diagnostics/peak_lag_hourly.py --pair HYS12-HYS14

# 4. Combine into canonical hourly Δt
python diagnostics/combine_clock_hourly.py

# 5. chronfix outlier filter + trigger detection
python -m chronfix.scripts.filter_and_triggers

# 6. Apply correction to all HYS14 daily mseed files
python -m chronfix.scripts.correct_hys14 --workers 8

# 7. Validation: re-run the diagnostic on the corrected output
python diagnostics/hys_ccf.py --pairs HYS12-HYS14 --tag corrected --workers 8 \
    --input-root-override HYS14=/data/wsd02/maleen_data/OOI-Data-corrected
python diagnostics/peak_lag_hourly.py --pair HYS12-HYS14_corrected
```

End-to-end runtime on cascadia is roughly:
- Download: 15 min (one-off)
- CCF: ~40 min per run (×2: original and corrected)
- Hourly peak-lag: seconds
- chronfix filter: seconds
- chronfix correction: ~5 min (8 workers)
- Total clock-on time: ~1.5 hours

---

## 12. Caveats and limitations

- **Δt is sampled at one value per UTC hour and the segment-smoothed
  model has a 24-hour rolling window.** Within-hour drift is modelled
  as the linear interpolation between consecutive segment-smoothed
  hourly samples. This is accurate when drift rates are smooth on
  multi-hour timescales (the regime the data live in) but cannot
  resolve sub-hour events. The smoother also has small (~0.1 s)
  end-effects within ~½ window of each trigger boundary; below the
  noise floor of CCF-coherence sensitivity, but worth noting.
- **Trigger localisation is to ~1 hour** because each 30-min CC window
  contributes data spanning ~30 min around the hour bucket center.
  Real resyncs that happen at, say, 14:23 UTC are reported with hour-
  level precision (14:00 or 15:00), not minute-level.
- **Days with NaN Δt yield no output.** The user can re-run with a
  relaxed outlier filter or accept the gap. Skipping is safer than
  fabricating a Δt across long unmeasured stretches.
- **Single channel (MHZ).** The same Δt(t) applies to all HYS14
  channels — clock errors are per-instrument — so extending to BHZ /
  HHZ / etc. is mechanical: just point `correct_hys14.py` at a
  different channel via a flag (not yet implemented).
- **Peak-lag picker is single-channel ZZ.** Three-component or
  whitening alternatives could in principle improve sub-second
  precision but are not needed for the 1–60 s drift magnitudes we
  observe here.
- **High-rate channel reuse.** chronfix currently resamples the trace
  data using linear interpolation, which is fine at 8 Hz with signals
  ≤ 3 Hz. Re-using the same Δt(t) on HHZ (200 Hz) to correct higher-
  frequency body-wave content would need a sinc / lanczos
  interpolator and a sub-second Δt model.

---

## 13. Output products

| Path | Contents |
|---|---|
| `data/ccf/HYS12-HYS14/` | Original (uncorrected) daily/30-min CCs and reference |
| `data/ccf/HYS12-HYS14_corrected/` | Corrected daily/30-min CCs and reference (validation) |
| `data/peak_lag_hourly/HYS12-HYS14/` | Hourly peak-lag track (uncorrected) |
| `data/peak_lag_hourly/HYS12-HYS14_corrected/` | Hourly peak-lag track + before/after plot |
| `data/clock_estimate/HYS14/delta_t_hourly_clean.npy` | Cleaned hourly Δt fed to chronfix |
| `data/clock_estimate/HYS14/trigger_periods.csv` | 32 trigger intervals |
| `data/clock_estimate/HYS14/filter_and_triggers.png` | Filter + trigger diagnostic |
| `data/clock_estimate/HYS14/correction_function.png` | Δt(t) actually applied to timestamps |
| `/data/wsd02/maleen_data/OOI-Data-corrected/HYS14/` | Corrected MiniSEED, mirroring input layout |
| `/data/wsd02/maleen_data/OOI-Data-corrected/HYS14/manifest.csv` | Per-output-segment manifest |
