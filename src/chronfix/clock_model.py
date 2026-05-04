"""ClockModel: the Δt(t) function chronfix applies to mseed traces.

Loads the cleaned hourly Δt and trigger intervals produced by
`filter_and_triggers.py` and exposes:

- `interp_delta_t(t)` — linear interpolation of Δt at any UTC time.
- `is_in_trigger(t)` — True if `t` falls within a trigger interval.
- `stable_intervals(t0, t1)` — UTC ranges between triggers within [t0, t1].

Δt convention: Δt > 0 means the station's clock is late vs UTC, so a
sample timestamped `t_apparent` was actually recorded at true UTC time
`t_apparent - Δt(t_apparent)`.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_CHRONOS_ROOT = Path("/home/seismic/chronos/data/clock_estimate/HYS14")


def _to_dt64(t) -> np.datetime64:
    """Coerce ObsPy UTCDateTime / datetime / np.datetime64 / str to datetime64[s]."""
    if isinstance(t, np.datetime64):
        return t.astype("datetime64[s]")
    if hasattr(t, "datetime"):  # ObsPy UTCDateTime
        return np.datetime64(t.datetime, "s")
    return np.datetime64(t, "s")


@dataclass
class ClockModel:
    """Cleaned hourly Δt + trigger intervals.

    Attributes
    ----------
    hour_times : np.ndarray, datetime64[h]
        UTC hour timestamps. Strictly increasing.
    delta_t : np.ndarray, float64
        Cleaned hourly Δt in seconds. NaN at outliers / missing.
    trigger_starts, trigger_ends : np.ndarray, datetime64[s]
        UTC bounds of each trigger interval (inclusive on both ends).
        Inside a trigger, Δt is undefined; chronfix must split here.
    station : str
    """
    hour_times: np.ndarray
    delta_t: np.ndarray
    trigger_starts: np.ndarray
    trigger_ends: np.ndarray
    station: str = "HYS14"

    def __post_init__(self) -> None:
        if self.hour_times.shape != self.delta_t.shape:
            raise ValueError("hour_times and delta_t shape mismatch")
        if len(self.trigger_starts) != len(self.trigger_ends):
            raise ValueError("trigger arrays length mismatch")

    @classmethod
    def from_chronos(
        cls,
        root: Path | str = DEFAULT_CHRONOS_ROOT,
        delta_t_file: str = "delta_t_hourly_clean.npy",
        triggers_file: str = "trigger_periods.csv",
        hour_times_file: str = "hour_times.npy",
        station: str = "HYS14",
    ) -> "ClockModel":
        root = Path(root)
        hour_times = np.load(root / hour_times_file)
        delta_t = np.load(root / delta_t_file).astype(np.float64, copy=False)

        df = pd.read_csv(root / triggers_file)
        # Map start_index / end_index from filter_and_triggers into UTC bounds.
        if len(df) == 0:
            t_starts = np.array([], dtype="datetime64[s]")
            t_ends = np.array([], dtype="datetime64[s]")
        else:
            si = df["start_index"].to_numpy(dtype=int)
            ei = df["end_index"].to_numpy(dtype=int)
            t_starts = hour_times[si].astype("datetime64[s]")
            t_ends = hour_times[ei].astype("datetime64[s]")

        return cls(
            hour_times=hour_times,
            delta_t=delta_t,
            trigger_starts=t_starts,
            trigger_ends=t_ends,
            station=station,
        )

    # -------------------------- queries --------------------------

    def is_in_trigger(self, t) -> bool:
        """True if `t` is inside any trigger interval (inclusive endpoints)."""
        ts = _to_dt64(t)
        if len(self.trigger_starts) == 0:
            return False
        idx = np.searchsorted(self.trigger_ends, ts, side="left")
        if idx >= len(self.trigger_ends):
            return False
        return bool(self.trigger_starts[idx] <= ts <= self.trigger_ends[idx])

    def interp_delta_t(self, t) -> float | np.ndarray:
        """Linear interp of Δt at UTC time(s). NaN inside triggers / out-of-range."""
        scalar = not isinstance(t, np.ndarray)
        if scalar:
            ts = np.array([_to_dt64(t)], dtype="datetime64[s]")
        else:
            ts = t.astype("datetime64[s]", copy=False)

        # Convert the hour-time grid and the query times to seconds-since-epoch
        # for linear interpolation.
        h_s = self.hour_times.astype("datetime64[s]").astype(np.int64)
        q_s = ts.astype(np.int64)

        valid = ~np.isnan(self.delta_t)
        if valid.sum() < 2:
            out = np.full(len(ts), np.nan, dtype=np.float64)
            return float(out[0]) if scalar else out

        h_valid = h_s[valid]
        d_valid = self.delta_t[valid]
        out = np.interp(q_s, h_valid, d_valid, left=np.nan, right=np.nan)

        # NaN where the query falls inside a trigger interval.
        if len(self.trigger_starts):
            ts_starts = self.trigger_starts.astype(np.int64)
            ts_ends = self.trigger_ends.astype(np.int64)
            for s, e in zip(ts_starts, ts_ends):
                in_range = (q_s >= s) & (q_s <= e)
                out[in_range] = np.nan

        return float(out[0]) if scalar else out

    def stable_intervals(self, t0=None, t1=None) -> list[tuple[np.datetime64, np.datetime64]]:
        """UTC intervals between triggers within [t0, t1].

        If t0 / t1 are None, defaults to the full hour_times span.
        """
        full_t0 = self.hour_times[0].astype("datetime64[s]") if len(self.hour_times) else None
        full_t1 = self.hour_times[-1].astype("datetime64[s]") if len(self.hour_times) else None
        a = _to_dt64(t0) if t0 is not None else full_t0
        b = _to_dt64(t1) if t1 is not None else full_t1
        if a is None or b is None or b <= a:
            return []

        # Build segment cuts at trigger ends (a stable interval begins right
        # after a trigger ends) and trigger starts (a stable interval ends just
        # before a trigger starts).
        cuts = [a]
        for s, e in zip(self.trigger_starts.astype("datetime64[s]"),
                        self.trigger_ends.astype("datetime64[s]")):
            if e < a or s > b:
                continue
            cuts.append(min(b, s))   # end of preceding stable interval
            cuts.append(max(a, e))   # start of next stable interval
        cuts.append(b)

        out: list[tuple[np.datetime64, np.datetime64]] = []
        for i in range(0, len(cuts), 2):
            t_lo = cuts[i]
            t_hi = cuts[i + 1] if i + 1 < len(cuts) else b
            if t_hi > t_lo:
                out.append((t_lo, t_hi))
        return out

    def __repr__(self) -> str:
        n = len(self.hour_times)
        nv = int((~np.isnan(self.delta_t)).sum())
        ntr = len(self.trigger_starts)
        return (f"ClockModel(station={self.station!r}, "
                f"n_hours={n}, n_valid={nv}, n_triggers={ntr})")
