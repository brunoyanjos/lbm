from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import xml.etree.ElementTree as ET

import numpy as np


@dataclass(frozen=True)
class VtiVelocity:
    path: Path
    nx: int
    ny: int
    step: int | None
    ux: np.ndarray
    uy: np.ndarray

    @property
    def speed(self) -> np.ndarray:
        return np.hypot(self.ux, self.uy)


def step_from_name(path: Path) -> int | None:
    match = re.search(r"output_(\d+)\.vti$", path.name)
    return int(match.group(1)) if match else None


def list_vti_files(vtk_dir: Path) -> list[Path]:
    if not vtk_dir.exists():
        raise FileNotFoundError(f"VTK directory not found: {vtk_dir}")

    files = sorted(vtk_dir.glob("output_*.vti"), key=lambda item: step_from_name(item) or -1)
    if not files:
        raise FileNotFoundError(f"no output_*.vti files found in: {vtk_dir}")

    return files


def select_vti_files(vtk_dir: Path, which: str) -> list[Path]:
    files = list_vti_files(vtk_dir)

    if which == "all":
        return files
    if which == "first":
        return [files[0]]
    if which == "latest":
        return [files[-1]]

    raise ValueError(f"unsupported VTI selection: {which}")


def _extent_size(piece: ET.Element) -> tuple[int, int]:
    extent = piece.attrib.get("Extent")
    if not extent:
        raise ValueError("VTI Piece is missing Extent")

    x0, x1, y0, y1, _z0, _z1 = (int(value) for value in extent.split())
    return x1 - x0 + 1, y1 - y0 + 1


def _data_array(point_data: ET.Element, name: str) -> ET.Element:
    for child in point_data.findall("DataArray"):
        if child.attrib.get("Name") == name:
            return child
    raise ValueError(f"VTI file is missing PointData array: {name}")


def _read_scalar_array(point_data: ET.Element, name: str, nx: int, ny: int) -> np.ndarray:
    data_array = _data_array(point_data, name)
    values = np.fromstring(data_array.text or "", sep=" ", dtype=np.float64)

    expected = nx * ny
    if values.size != expected:
        raise ValueError(f"array {name!r} has {values.size} values, expected {expected}")

    return values.reshape((ny, nx))


def read_vti_velocity(path: Path) -> VtiVelocity:
    root = ET.parse(path).getroot()
    piece = root.find("./ImageData/Piece")
    if piece is None:
        raise ValueError(f"invalid VTI ImageData/Piece structure: {path}")

    nx, ny = _extent_size(piece)
    point_data = piece.find("PointData")
    if point_data is None:
        raise ValueError(f"VTI file is missing PointData: {path}")

    ux = _read_scalar_array(point_data, "ux", nx, ny)
    uy = _read_scalar_array(point_data, "uy", nx, ny)

    if np.any(~np.isfinite(ux)) or np.any(~np.isfinite(uy)):
        raise ValueError(f"non-finite velocity values found in: {path}")

    return VtiVelocity(path=path, nx=nx, ny=ny, step=step_from_name(path), ux=ux, uy=uy)
