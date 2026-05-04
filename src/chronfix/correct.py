"""Apply a ClockModel to ObsPy traces.

Splits each input trace at trigger boundaries (real clock discontinuities)
and writes one corrected sub-trace per stable segment overlap. Two methods:

- ``method="resample"`` (default, recommended for hourly Δt input):
  resamples the trace data onto a regular UTC grid, accounting for
  within-segment drift via linear interp of Δt(t).

- ``method="shift_only"``: subtracts Δt at the chunk start from the
  starttime; data is bit-identical. Within-chunk drift is approximated
  as constant. Useful for spot-check / validation.
"""
from __future__ import annotations

from typing import Literal

import numpy as np
from obspy import Trace, Stream, UTCDateTime

from chronfix.clock_model import ClockModel


def _to_utc(t: np.datetime64) -> UTCDateTime:
    return UTCDateTime(t.astype("datetime64[us]").astype(str))


def correct_trace(
    tr: Trace,
    model: ClockModel,
    method: Literal["resample", "shift_only"] = "resample",
    interpolator: Literal["linear"] = "linear",
) -> list[Trace]:
    """Return one corrected sub-trace per stable-segment overlap.

    The output traces have:
        - starttime aligned to true UTC (resample) or apparent-start - Δt
          (shift_only).
        - Same nominal sampling_rate as input.
        - Data dtype preserved.

    Days where Δt is unavailable (NaN throughout the segment) yield no
    output.
    """
    if tr.stats.npts == 0:
        return []

    fs = float(tr.stats.sampling_rate)
    t_apparent_start = tr.stats.starttime
    t_apparent_end = tr.stats.endtime

    # Find stable intervals overlapping the trace's apparent time range.
    stables = model.stable_intervals(t_apparent_start.datetime, t_apparent_end.datetime)
    out: list[Trace] = []

    for stable_t0, stable_t1 in stables:
        s0 = _to_utc(stable_t0)
        s1 = _to_utc(stable_t1)
        sub_t0 = max(s0, t_apparent_start)
        sub_t1 = min(s1, t_apparent_end)
        if sub_t1 <= sub_t0:
            continue

        try:
            sub_tr = tr.slice(sub_t0, sub_t1, nearest_sample=False)
        except Exception:
            continue
        if sub_tr.stats.npts < 2:
            continue

        if method == "shift_only":
            corrected = _shift_only(sub_tr, model)
        elif method == "resample":
            corrected = _resample(sub_tr, model)
        else:
            raise ValueError(f"unknown method: {method!r}")

        if corrected is not None:
            out.append(corrected)

    return out


def correct_stream(
    st: Stream,
    model: ClockModel,
    method: Literal["resample", "shift_only"] = "resample",
) -> Stream:
    """Apply correct_trace to every trace and return a flat Stream of chunks."""
    chunks: list[Trace] = []
    for tr in st:
        chunks.extend(correct_trace(tr, model, method=method))
    return Stream(chunks)


# --------------------------- internals ---------------------------

def _shift_only(tr: Trace, model: ClockModel) -> Trace | None:
    """Subtract Δt at the trace start from starttime; leave data alone."""
    dt0 = model.interp_delta_t(np.datetime64(tr.stats.starttime.datetime, "s"))
    if not np.isfinite(dt0):
        return None
    new = tr.copy()
    new.stats.starttime = tr.stats.starttime - float(dt0)
    return new


def _resample(tr: Trace, model: ClockModel) -> Trace | None:
    """Resample the trace onto a regular UTC grid using linear interp.

    For each output sample at true UTC time `t_utc`, look up the input
    trace at apparent time `t_apparent ≈ t_utc + Δt(t_utc)`. Write the
    interpolated value at `t_utc` on a regular sample grid.
    """
    fs = float(tr.stats.sampling_rate)
    n_in = tr.stats.npts
    if n_in < 2:
        return None

    # Apparent time axis of the input trace (UTCDateTime relative to its starttime).
    apparent_start = tr.stats.starttime
    apparent_end = tr.stats.endtime
    # True UTC bounds of this chunk:
    dt_start = model.interp_delta_t(np.datetime64(apparent_start.datetime, "s"))
    dt_end = model.interp_delta_t(np.datetime64(apparent_end.datetime, "s"))
    if not (np.isfinite(dt_start) and np.isfinite(dt_end)):
        return None
    utc_start = apparent_start - float(dt_start)
    utc_end = apparent_end - float(dt_end)
    if utc_end <= utc_start:
        return None

    # Output UTC sample grid at the same nominal fs.
    n_out = int(np.floor((utc_end - utc_start) * fs)) + 1
    if n_out < 2:
        return None
    t_utc_offsets = np.arange(n_out, dtype=np.float64) / fs  # seconds since utc_start

    # Compute Δt at each output sample's UTC time. Vectorized: build datetime64
    # array of UTC times in seconds.
    utc_start_dt64 = np.datetime64(utc_start.datetime, "us")
    utc_query_us = utc_start_dt64 + (t_utc_offsets * 1e6).astype("timedelta64[us]")
    dt_at_utc = model.interp_delta_t(utc_query_us)
    if np.isnan(dt_at_utc).any():
        # Mask: only keep contiguous valid run starting from sample 0.
        valid_mask = np.isfinite(dt_at_utc)
        if not valid_mask[0]:
            return None
        # Trim to first NaN
        first_nan = int(np.argmax(~valid_mask)) if (~valid_mask).any() else len(valid_mask)
        if first_nan < 2:
            return None
        n_out = first_nan
        t_utc_offsets = t_utc_offsets[:n_out]
        dt_at_utc = dt_at_utc[:n_out]

    # Apparent-time targets at each output sample.
    apparent_targets_offsets = t_utc_offsets + dt_at_utc + (utc_start - apparent_start)
    # ^ offsets relative to apparent_start

    # Input data on apparent-time grid relative to apparent_start.
    apparent_grid = np.arange(n_in, dtype=np.float64) / fs

    data_in = np.asarray(tr.data, dtype=np.float64)
    data_out = np.interp(apparent_targets_offsets, apparent_grid, data_in,
                         left=np.nan, right=np.nan)

    # Drop trailing NaNs (target fell off the input range).
    valid_out = np.isfinite(data_out)
    if not valid_out.all():
        last = int(np.argmax(~valid_out)) if (~valid_out).any() else len(valid_out)
        if last < 2:
            return None
        data_out = data_out[:last]

    new = Trace(data=data_out.astype(tr.data.dtype, copy=False))
    new.stats = tr.stats.copy()
    new.stats.starttime = utc_start
    new.stats.npts = len(data_out)
    return new
