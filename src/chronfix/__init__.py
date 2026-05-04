"""chronfix — apply chronos-measured timing corrections to MiniSEED data.

Inputs: a chronos correction directory (delta_t_hourly_clean.npy +
hour_times.npy + trigger_periods.csv) and raw MiniSEED files.
Outputs: corrected MiniSEED, split at trigger boundaries.

Public API:
    chronfix.ClockModel       — load + query the correction file
    chronfix.correct_trace    — apply correction to one trace
    chronfix.correct_stream   — apply correction to a Stream
"""
from chronfix.clock_model import ClockModel
from chronfix.correct import correct_trace, correct_stream

__all__ = ["ClockModel", "correct_trace", "correct_stream"]
__version__ = "0.1.0"
