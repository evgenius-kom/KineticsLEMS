"""3-column (t, T, y) wave reader + conversion path."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from kinetics_lems.conversion import build_runs
from kinetics_lems.io.wave_reader import read_wave
from kinetics_lems.models import (
    CaseData,
    CaseParams,
    ExperimentType,
    Method,
    Wave,
)


def _write_lines(path: Path, lines: list[str]) -> None:
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_two_column_path_still_works(tmp_path: Path) -> None:
    path = tmp_path / "two_col.txt"
    _write_lines(path, ["300\t0.1", "350\t0.5", "400\t0.9"])
    wave = read_wave(path)
    assert wave.has_recorded_time is False
    np.testing.assert_array_equal(wave.x, [300.0, 350.0, 400.0])
    np.testing.assert_array_equal(wave.y, [0.1, 0.5, 0.9])


def test_three_column_detected_and_time_preserved(tmp_path: Path) -> None:
    path = tmp_path / "three_col.txt"
    _write_lines(path, ["0\t300\t0.1", "60\t350\t0.5", "120\t400\t0.9"])
    wave = read_wave(path)
    assert wave.has_recorded_time is True
    np.testing.assert_array_equal(wave.t_seconds, [0.0, 60.0, 120.0])
    np.testing.assert_array_equal(wave.x, [300.0, 350.0, 400.0])
    np.testing.assert_array_equal(wave.y, [0.1, 0.5, 0.9])


def test_comma_decimal_handled_in_three_column(tmp_path: Path) -> None:
    path = tmp_path / "comma.txt"
    _write_lines(path, ["0\t300,0\t0,10", "60\t350,0\t0,50", "120\t400,0\t0,90"])
    wave = read_wave(path)
    assert wave.has_recorded_time is True
    np.testing.assert_array_almost_equal(wave.y, [0.10, 0.50, 0.90])


def test_three_column_conversion_matches_linear_baseline_case() -> None:
    """A linearly-heated 3-column wave must produce the same α(T) as a
    matching 2-column wave with the same nominal β."""
    beta_K_per_min = 10.0
    beta_K_per_sec = beta_K_per_min / 60.0
    T = np.linspace(400.0, 600.0, 201)
    t = (T - T[0]) / beta_K_per_sec  # seconds since start of heating

    # Single Gaussian peak in dα/dT.
    y = np.exp(-((T - 500.0) ** 2) / (2 * 20.0**2))

    wave_2col = Wave(x=T, y=y)
    wave_3col = Wave(x=T, y=y, t_seconds=t)

    case_2 = CaseData(
        params=CaseParams(
            material="m",
            experiment_type=ExperimentType.HEATING,
            method=Method.DSC,
            file_to_condition={"a.txt": beta_K_per_min},
        ),
        waves={beta_K_per_min: wave_2col},
    )
    case_3 = CaseData(
        params=case_2.params,
        waves={beta_K_per_min: wave_3col},
    )
    run2 = build_runs(case_2)[0]
    run3 = build_runs(case_3)[0]
    # Same α(T) curve up to numerical precision.
    np.testing.assert_allclose(run2.alpha, run3.alpha, atol=1e-12)
    # dα/dt may differ in absolute scale because the 3-column path uses
    # ∫y dt while the 2-col path uses ∫y dT, but under linear heating both
    # should produce a similar peak shape.
    assert run3.dalpha_dt.max() > 0


def test_strictly_monotone_time_required() -> None:
    T = np.array([300.0, 350.0, 400.0])
    y = np.array([0.1, 0.5, 0.9])
    bad_t = np.array([0.0, 60.0, 60.0])  # not strictly increasing
    import pytest

    with pytest.raises(ValueError):
        Wave(x=T, y=y, t_seconds=bad_t)
