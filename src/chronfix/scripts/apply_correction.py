#!/usr/bin/env python
"""Apply a chronos correction file to one station's MiniSEED daily files.

For each input day file, splits the trace at trigger boundaries (real
clock discontinuities) and writes one corrected output per stable
segment overlap. Output mirrors the input layout under a parallel root.

Run from anywhere:

    python -m chronfix.scripts.apply_correction \
        --correction-dir /home/seismic/chronos/data/clock_estimate/HYS14 \
        --network OO --station HYS14 --channel MHZ \
        --input-root  /data/wsd02/maleen_data/OOI-Data \
        --output-root /data/wsd02/maleen_data/OOI-Data-corrected \
        --start 2022-01-01 --workers 8

The correction directory is the output of chronos.scripts.filter_and_triggers
and must contain:
    delta_t_hourly_clean.npy
    hour_times.npy
    trigger_periods.csv
"""
from __future__ import annotations

import argparse
import csv
import logging
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import date, timedelta
from pathlib import Path

from obspy import read

from chronfix.clock_model import ClockModel
from chronfix.correct import correct_stream

LOG = logging.getLogger("apply_correction")


def daily_input(root: Path, network: str, sta: str, ch: str, d: date) -> Path:
    doy = d.timetuple().tm_yday
    return (root / sta / f"{d.year}" / f"{doy:03d}"
            / f"{sta}.{network}.{d.year}.{doy:03d}.{ch}")


def daily_output_dir(root: Path, sta: str, d: date) -> Path:
    doy = d.timetuple().tm_yday
    return root / sta / f"{d.year}" / f"{doy:03d}"


def correct_day(
    d: date, correction_dir: str, method: str,
    network: str, station: str, channel: str,
    input_root: str, output_root: str,
) -> dict:
    """Worker: correct one day's mseed. Returns a manifest result."""
    in_path = daily_input(Path(input_root), network, station, channel, d)
    if not in_path.exists():
        return {"date": str(d), "input": str(in_path), "status": "missing", "n_chunks": 0}

    try:
        st = read(str(in_path))
    except Exception as ex:
        return {"date": str(d), "input": str(in_path), "status": f"read_err:{ex}",
                "n_chunks": 0}

    model = ClockModel.from_chronos(correction_dir, station=station)
    corrected = correct_stream(st, model, method=method)
    if len(corrected) == 0:
        return {"date": str(d), "input": str(in_path),
                "status": "no_overlap_or_dt_unavailable", "n_chunks": 0}

    out_dir = daily_output_dir(Path(output_root), station, d)
    out_dir.mkdir(parents=True, exist_ok=True)
    doy = d.timetuple().tm_yday
    out_path = out_dir / f"{station}.{network}.{d.year}.{doy:03d}.{channel}"
    corrected.write(str(out_path), format="MSEED")
    rows = [{
        "date": str(d),
        "input": str(in_path),
        "output": str(out_path),
        "segment": i,
        "utc_start": str(tr.stats.starttime),
        "utc_end": str(tr.stats.endtime),
        "n_samples": int(tr.stats.npts),
        "method": method,
    } for i, tr in enumerate(corrected, start=1)]
    return {"date": str(d), "n_chunks": len(rows), "rows": rows}


def daterange(start: date, end: date):
    d = start
    while d <= end:
        yield d
        d += timedelta(days=1)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--correction-dir", required=True,
                   help="Directory containing delta_t_hourly_clean.npy, "
                        "hour_times.npy, trigger_periods.csv (output of "
                        "chronos.scripts.filter_and_triggers).")
    p.add_argument("--input-root", required=True,
                   help="Raw MiniSEED tree root, e.g. /data/.../OOI-Data.")
    p.add_argument("--output-root", required=True,
                   help="Where to write corrected files (parallel layout).")
    p.add_argument("--network", required=True, help="Network code, e.g. OO.")
    p.add_argument("--station", required=True, help="Station code being corrected.")
    p.add_argument("--channel", required=True, help="Channel code, e.g. MHZ.")
    p.add_argument("--start", type=date.fromisoformat, default=date(2022, 1, 1))
    p.add_argument("--end", type=date.fromisoformat, default=date.today())
    p.add_argument("--method", choices=["resample", "shift_only"], default="resample")
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--manifest", default=None,
                   help="Manifest CSV path. Default: <output-root>/<station>/manifest.csv")
    p.add_argument("--verbose", "-v", action="store_true")
    args = p.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    days = list(daterange(args.start, args.end))
    LOG.info("processing %d days  station=%s.%s.%s  method=%s  workers=%d",
             len(days), args.network, args.station, args.channel,
             args.method, args.workers)

    manifest_rows: list[dict] = []
    counts = {"ok": 0, "missing": 0, "no_overlap": 0, "err": 0, "chunks": 0}

    fixed = dict(
        correction_dir=args.correction_dir, method=args.method,
        network=args.network, station=args.station, channel=args.channel,
        input_root=args.input_root, output_root=args.output_root,
    )

    if args.workers <= 1:
        for d in days:
            res = correct_day(d, **fixed)
            _tally(res, counts, manifest_rows)
            LOG.info("[%s] %s n_chunks=%d", d, res.get("status", "ok"), res["n_chunks"])
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futs = {ex.submit(correct_day, d, **fixed): d for d in days}
            for fut in as_completed(futs):
                d = futs[fut]
                res = fut.result()
                _tally(res, counts, manifest_rows)
                LOG.info("[%s] %s n_chunks=%d", d, res.get("status", "ok"), res["n_chunks"])

    if manifest_rows:
        manifest_path = Path(args.manifest) if args.manifest else (
            Path(args.output_root) / args.station / "manifest.csv"
        )
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        write_header = not manifest_path.exists()
        with open(manifest_path, "a", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(manifest_rows[0].keys()))
            if write_header:
                w.writeheader()
            w.writerows(manifest_rows)
        LOG.info("manifest -> %s (%d rows)", manifest_path, len(manifest_rows))

    LOG.info("done: ok=%d missing=%d no_overlap/no_dt=%d err=%d chunks=%d",
             counts["ok"], counts["missing"], counts["no_overlap"],
             counts["err"], counts["chunks"])
    return 0


def _tally(res: dict, counts: dict, manifest_rows: list[dict]) -> None:
    status = res.get("status")
    if status == "missing":
        counts["missing"] += 1
    elif status == "no_overlap_or_dt_unavailable":
        counts["no_overlap"] += 1
    elif status and status.startswith("read_err"):
        counts["err"] += 1
    else:
        counts["ok"] += 1
        counts["chunks"] += res["n_chunks"]
        manifest_rows.extend(res.get("rows", []))


if __name__ == "__main__":
    sys.exit(main())
