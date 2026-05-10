"""Round-trip: generate synthetic case → write to disk → load → analyze."""
from __future__ import annotations

import shutil
import zipfile
from pathlib import Path

import numpy as np
import pytest

from kinetics_lems.config import load_config
from kinetics_lems.io import load_case, read_wave
from kinetics_lems.runner import run_analysis
from kinetics_lems.synthetic import generate_case, write_case


def test_wave_reader_handles_comma_decimals(tmp_path: Path):
    p = tmp_path / "wave.txt"
    p.write_text("300,0\t0,1\n400,5\t0,5\n500,0\t0,9\n", encoding="utf-8")
    w = read_wave(p)
    assert w.x.shape == (3,)
    assert pytest.approx(w.x[1]) == 400.5
    assert pytest.approx(w.y[2]) == 0.9


def test_round_trip_folder(tmp_path: Path):
    case_in = generate_case([2.5, 5.0, 10.0], Ea_J_per_mol=110_000.0, n_points=1500, seed=0)
    case_dir = tmp_path / "case"
    write_case(case_in, case_dir)
    case_out = load_case(case_dir)
    assert sorted(case_out.waves.keys()) == [2.5, 5.0, 10.0]


def test_round_trip_zip(tmp_path: Path):
    case_in = generate_case([5.0, 10.0, 20.0], Ea_J_per_mol=110_000.0, n_points=1500, seed=0)
    case_dir = tmp_path / "case"
    write_case(case_in, case_dir)

    zip_path = tmp_path / "case.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        for f in case_dir.iterdir():
            zf.write(f, arcname=f.name)
    shutil.rmtree(case_dir)

    case_out = load_case(zip_path)
    cfg = load_config()
    res = run_analysis(case_out, cfg)
    assert abs(np.nanmean(res.isoconversional["friedman"].Ea_kJ_per_mol) - 110.0) < 2.0
