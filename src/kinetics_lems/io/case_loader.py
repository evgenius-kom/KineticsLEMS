"""Load a CaseData (settings + waves) from a folder or .zip archive.

Case layout (folder or root of a zip):
    settings.json  — see schema below
    *.txt          — two-column experimental waves referenced from settings

settings.json (legacy schema, kept for compatibility):
{
  "Settings": [
    {"name": "Material",        "value": "..."},
    {"name": "ExperimentType",  "value": "heating"|"cooling"|"isothermal"},
    {"name": "Method",          "value": "DSC"|"TGA"|"FSC"|"POM"},
    {"name": "Conditions",      "value": {"file1.txt": "5", "file2.txt": "10", ...}}
  ]
}
"""
from __future__ import annotations

import json
import tempfile
import zipfile
from pathlib import Path

from ..models import CaseData, CaseParams, ExperimentType, Method, Wave
from .wave_reader import read_wave

CASE_SETTINGS_FILE = "settings.json"


def load_case(path: Path | str) -> CaseData:
    """Load a case from either a directory or a .zip archive."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    if path.is_dir():
        return _load_from_dir(path)

    if path.suffix.lower() == ".zip":
        return _load_from_zip(path)

    raise ValueError(f"Expected a directory or .zip, got {path}")


def _load_from_zip(zip_path: Path) -> CaseData:
    with tempfile.TemporaryDirectory(prefix="kinetics_lems_") as tmp:
        tmp_dir = Path(tmp)
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(tmp_dir)
        # Some zips wrap content in a single subfolder; unwrap if so.
        entries = list(tmp_dir.iterdir())
        if len(entries) == 1 and entries[0].is_dir():
            tmp_dir = entries[0]
        return _load_from_dir(tmp_dir)


def _load_from_dir(folder: Path) -> CaseData:
    settings_path = folder / CASE_SETTINGS_FILE
    if not settings_path.is_file():
        raise FileNotFoundError(f"Missing {CASE_SETTINGS_FILE} in {folder}")
    with settings_path.open("r", encoding="utf-8") as f:
        raw = json.load(f)

    params = _parse_params(raw)
    waves: dict[float, Wave] = {}
    for filename, condition in params.file_to_condition.items():
        waves[float(condition)] = read_wave(folder / filename)
    return CaseData(params=params, waves=waves)


def _parse_params(raw: dict) -> CaseParams:
    if "Settings" not in raw or not isinstance(raw["Settings"], list):
        raise ValueError("settings.json: 'Settings' must be an array")

    items: dict[str, object] = {}
    for item in raw["Settings"]:
        items[item["name"]] = item["value"]

    try:
        material = str(items["Material"])
        exp_type = ExperimentType(str(items["ExperimentType"]).lower())
        method = Method(str(items["Method"]).upper())
        conditions = items["Conditions"]
    except KeyError as e:
        raise ValueError(f"settings.json: missing field {e}") from None

    if not isinstance(conditions, dict) or not conditions:
        raise ValueError("settings.json: 'Conditions' must be a non-empty mapping")

    file_to_condition = {str(k): float(v) for k, v in conditions.items()}
    return CaseParams(
        material=material,
        experiment_type=exp_type,
        method=method,
        file_to_condition=file_to_condition,
    )


__all__ = ["load_case", "CASE_SETTINGS_FILE"]
