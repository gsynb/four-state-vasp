#!/usr/bin/env python3
"""
Prepare and postprocess VASP four-state magnetic-interaction calculations.

This script is intentionally dependency-light so it can run on login nodes where
only Python and the usual VASP/Slurm tools are available.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
import textwrap
from dataclasses import dataclass, replace
from itertools import combinations, permutations, product
from pathlib import Path
from typing import Iterable


JANI_COMPONENTS = ["Jxx", "Jxy", "Jyy", "Jyz", "Jzz", "Jzx", "Jxz", "Jyx", "Jzy"]
SIA_COMPONENTS = ["Axy", "Axz", "Ayz", "Ayy_minus_Axx", "Azz_minus_Axx"]
JISO_SPINS = ["upup", "updn", "dnup", "dndn"]
STAGE_SINGLE = [{"name": "single", "base": "", "suffix": "single"}]
STAGE_PBE_HSE = [
    {"name": "pbe_pre", "base": "pbe", "suffix": "pbe_pre"},
    {"name": "hse_no_u", "base": "", "suffix": "hse_no_u"},
]
REQUIRED_INPUTS = ["INCAR", "KPOINTS", "POSCAR", "POTCAR"]
AXIS = {
    "x": (1.0, 0.0, 0.0),
    "y": (0.0, 1.0, 0.0),
    "z": (0.0, 0.0, 1.0),
}
KITAEV_AXIS_LABELS = ("x", "y", "z")
IDEAL_KITAEV_BASIS = (
    (-(6.0**0.5) / 3.0, 0.0, (3.0**0.5) / 3.0),
    ((6.0**0.5) / 6.0, (2.0**0.5) / 2.0, (3.0**0.5) / 3.0),
    ((6.0**0.5) / 6.0, -(2.0**0.5) / 2.0, (3.0**0.5) / 3.0),
)
DEFAULT_U = {
    "Ag": 1.5,
    "Co": 3.4,
    "Cr": 3.5,
    "Cu": 4.0,
    "Fe": 4.0,
    "Mn": 3.9,
    "Mo": 3.5,
    "Nb": 1.5,
    "Ni": 6.0,
    "Rh": 3.0,
    "V": 3.1,
}


@dataclass(frozen=True)
class PoscarInfo:
    path: Path
    elements: list[str]
    counts: list[int]
    atom_symbols: list[str]
    lattice: list[tuple[float, float, float]]
    frac_coords: list[tuple[float, float, float]]

    @property
    def natoms(self) -> int:
        return len(self.atom_symbols)


@dataclass(frozen=True)
class PairInfo:
    label_i: int
    label_j: int
    global_i: int
    global_j: int
    magnetic_i: int
    magnetic_j: int

    @property
    def label(self) -> str:
        return f"pair_{self.label_i}_{self.label_j}"


@dataclass(frozen=True)
class PairBondContext:
    image_shift: tuple[int, int, int]
    distance: float
    nearest_image_shift: tuple[int, int, int]
    nearest_distance: float
    multiplicity: int
    equivalent_shifts: list[tuple[int, int, int]]


@dataclass(frozen=True)
class NeighborContact:
    shell: int
    distance: float
    i_global: int
    j_global: int
    i_mag: int
    j_mag: int
    image_shift: tuple[int, int, int]
    disp_frac: tuple[float, float, float]
    disp_cart: tuple[float, float, float]

    @property
    def crosses_boundary(self) -> bool:
        return any(self.image_shift)

    @property
    def pair_label(self) -> str:
        return f"pair_{self.i_global + 1}_{self.j_global + 1}"


@dataclass(frozen=True)
class KitaevFrame:
    pair: PairInfo
    gamma_axis: tuple[float, float, float]
    local_x: tuple[float, float, float]
    local_y: tuple[float, float, float]
    gamma_axis_j: tuple[float, float, float]
    local_x_j: tuple[float, float, float]
    local_y_j: tuple[float, float, float]
    bond_axis: tuple[float, float, float]
    pair_shift: tuple[int, int, int]
    shared_ligands: list[tuple[int, tuple[int, int, int], float, float]]
    gamma_label: str
    alpha_label: str
    beta_label: str
    kitaev_axes: list[tuple[str, tuple[float, float, float]]]
    kitaev_axes_j: list[tuple[str, tuple[float, float, float]]]
    reference_basis: list[tuple[float, float, float]]
    reference_basis_j: list[tuple[float, float, float]]
    axis_overlaps: tuple[float, float, float]
    axis_overlaps_j: tuple[float, float, float]
    axis_bond_dots: tuple[float, float, float]
    axis_consistency_degrees: tuple[float, float, float]
    ligand_method: str
    axis_match: str
    axis_match_j: str
    octahedral_ligands: list[tuple[int, tuple[int, int, int], float]]
    octahedral_ligands_j: list[tuple[int, tuple[int, int, int], float]]
    octahedral_angle_error: float | None
    octahedral_angle_error_j: float | None


def is_int_token(token: str) -> bool:
    try:
        int(token)
        return True
    except ValueError:
        return False


def parse_vec3(line: str, context: str) -> tuple[float, float, float]:
    parts = line.split()
    if len(parts) < 3:
        raise ValueError(f"Expected three numbers for {context}: {line!r}")
    try:
        return (float(parts[0]), float(parts[1]), float(parts[2]))
    except ValueError as exc:
        raise ValueError(f"Bad numeric vector for {context}: {line!r}") from exc


def vec_add(a: Iterable[float], b: Iterable[float]) -> tuple[float, float, float]:
    ax, ay, az = a
    bx, by, bz = b
    return (float(ax) + float(bx), float(ay) + float(by), float(az) + float(bz))


def vec_sub(a: Iterable[float], b: Iterable[float]) -> tuple[float, float, float]:
    ax, ay, az = a
    bx, by, bz = b
    return (float(ax) - float(bx), float(ay) - float(by), float(az) - float(bz))


def vec_scale(a: Iterable[float], scale: float) -> tuple[float, float, float]:
    ax, ay, az = a
    return (float(ax) * scale, float(ay) * scale, float(az) * scale)


def vec_dot(a: Iterable[float], b: Iterable[float]) -> float:
    ax, ay, az = a
    bx, by, bz = b
    return float(ax) * float(bx) + float(ay) * float(by) + float(az) * float(bz)


def vec_cross(a: Iterable[float], b: Iterable[float]) -> tuple[float, float, float]:
    ax, ay, az = a
    bx, by, bz = b
    return (
        float(ay) * float(bz) - float(az) * float(by),
        float(az) * float(bx) - float(ax) * float(bz),
        float(ax) * float(by) - float(ay) * float(bx),
    )


def vec_norm(a: Iterable[float]) -> float:
    return math.sqrt(vec_dot(a, a))


def vec_normalize(a: Iterable[float], label: str = "vector") -> tuple[float, float, float]:
    length = vec_norm(a)
    if length < 1e-12:
        raise ValueError(f"Cannot normalize near-zero {label}")
    return tuple(float(x) / length for x in a)  # type: ignore[return-value]


def vec_angle_degrees(a: Iterable[float], b: Iterable[float]) -> float:
    ua = vec_normalize(a, "angle vector a")
    ub = vec_normalize(b, "angle vector b")
    cos_theta = max(-1.0, min(1.0, vec_dot(ua, ub)))
    return math.degrees(math.acos(cos_theta))


def basis_linear_combination(
    coeffs: Iterable[float],
    basis: list[tuple[float, float, float]],
    label: str,
) -> tuple[float, float, float]:
    cx, cy, cz = coeffs
    vec = (0.0, 0.0, 0.0)
    for coeff, axis in zip((cx, cy, cz), basis):
        vec = vec_add(vec, vec_scale(axis, float(coeff)))
    return vec_normalize(vec, label)


def det3(rows: list[tuple[float, float, float]]) -> float:
    a, b, c = rows
    return vec_dot(a, vec_cross(b, c))


def inverse_lattice_rows(lattice: list[tuple[float, float, float]]) -> list[tuple[float, float, float]]:
    a, b, c = lattice
    volume = det3(lattice)
    if abs(volume) < 1e-12:
        raise ValueError("POSCAR lattice is singular")
    # Columns of the inverse row-lattice transform are reciprocal dual vectors.
    dual_a = vec_scale(vec_cross(b, c), 1.0 / volume)
    dual_b = vec_scale(vec_cross(c, a), 1.0 / volume)
    dual_c = vec_scale(vec_cross(a, b), 1.0 / volume)
    return [dual_a, dual_b, dual_c]


def frac_to_cart(frac: Iterable[float], lattice: list[tuple[float, float, float]]) -> tuple[float, float, float]:
    f0, f1, f2 = frac
    a, b, c = lattice
    return (
        float(f0) * a[0] + float(f1) * b[0] + float(f2) * c[0],
        float(f0) * a[1] + float(f1) * b[1] + float(f2) * c[1],
        float(f0) * a[2] + float(f1) * b[2] + float(f2) * c[2],
    )


def cart_to_frac(cart: Iterable[float], lattice: list[tuple[float, float, float]]) -> tuple[float, float, float]:
    dual = inverse_lattice_rows(lattice)
    return (vec_dot(cart, dual[0]), vec_dot(cart, dual[1]), vec_dot(cart, dual[2]))


def nearest_int(value: float) -> int:
    return int(math.floor(value + 0.5)) if value >= 0 else int(math.ceil(value - 0.5))


def nearest_image_delta(
    frac_i: tuple[float, float, float],
    frac_j: tuple[float, float, float],
    lattice: list[tuple[float, float, float]],
    search_radius: int = 2,
) -> tuple[tuple[float, float, float], tuple[int, int, int]]:
    raw = vec_sub(frac_j, frac_i)
    center = tuple(-nearest_int(x) for x in raw)
    best_shift, best_frac, best_distance, best_tie = find_nearest_shift_in_box(
        raw,
        lattice,
        [
            (center[axis] - max(0, search_radius), center[axis] + max(0, search_radius))
            for axis in range(3)
        ],
    )
    if best_shift is None or best_frac is None:
        raise ValueError("Could not determine nearest periodic image")

    heights = lattice_heights(lattice)
    if any(height <= 1e-12 for height in heights):
        raise ValueError("Cannot determine nearest periodic image for a singular lattice")
    eps = 1e-10
    bounds: list[tuple[int, int]] = []
    for axis, height in enumerate(heights):
        limit = best_distance / height + eps
        low = math.ceil(-raw[axis] - limit - eps)
        high = math.floor(-raw[axis] + limit + eps)
        bounds.append((low, high))

    strict_shift, strict_frac, strict_distance, strict_tie = find_nearest_shift_in_box(raw, lattice, bounds)
    if strict_shift is not None and strict_frac is not None:
        best_shift, best_frac, best_distance, best_tie = strict_shift, strict_frac, strict_distance, strict_tie
    if best_shift is None or best_frac is None:
        raise ValueError("Could not determine nearest periodic image")
    return best_frac, best_shift


def find_nearest_shift_in_box(
    raw_frac_delta: tuple[float, float, float],
    lattice: list[tuple[float, float, float]],
    bounds: list[tuple[int, int]],
) -> tuple[tuple[int, int, int] | None, tuple[float, float, float] | None, float, tuple[int, int, int, int]]:
    best_shift: tuple[int, int, int] | None = None
    best_frac: tuple[float, float, float] | None = None
    best_distance = float("inf")
    best_tie = (10**9, 10**9, 10**9, 10**9)
    for tx in range(bounds[0][0], bounds[0][1] + 1):
        for ty in range(bounds[1][0], bounds[1][1] + 1):
            for tz in range(bounds[2][0], bounds[2][1] + 1):
                shift = (tx, ty, tz)
                disp_frac = (
                    raw_frac_delta[0] + tx,
                    raw_frac_delta[1] + ty,
                    raw_frac_delta[2] + tz,
                )
                distance = vec_norm(frac_to_cart(disp_frac, lattice))
                tie = (abs(tx) + abs(ty) + abs(tz), abs(tx), abs(ty), abs(tz))
                if distance < best_distance - 1e-10 or (abs(distance - best_distance) <= 1e-10 and tie < best_tie):
                    best_shift = shift
                    best_frac = disp_frac
                    best_distance = distance
                    best_tie = tie
    return best_shift, best_frac, best_distance, best_tie


def lattice_heights(lattice: list[tuple[float, float, float]]) -> list[float]:
    a, b, c = lattice
    volume = abs(det3(lattice))
    areas = [vec_norm(vec_cross(b, c)), vec_norm(vec_cross(c, a)), vec_norm(vec_cross(a, b))]
    return [volume / area if area > 1e-12 else 0.0 for area in areas]


def parse_lattice_scale(
    scale_line: str,
    raw_lattice: list[tuple[float, float, float]],
) -> tuple[list[tuple[float, float, float]], tuple[float, float, float]]:
    tokens = [float(tok) for tok in scale_line.split()]
    if not tokens:
        raise ValueError("POSCAR scale line is empty")
    if len(tokens) == 1:
        scale = tokens[0]
        if scale < 0:
            raw_volume = abs(det3(raw_lattice))
            if raw_volume < 1e-12:
                raise ValueError("Cannot use negative POSCAR volume scale with singular lattice")
            scale = ((-scale) / raw_volume) ** (1.0 / 3.0)
        lattice = [vec_scale(vec, scale) for vec in raw_lattice]
        return lattice, (scale, scale, scale)
    if len(tokens) == 3:
        sx, sy, sz = tokens
        lattice = [(x * sx, y * sy, z * sz) for x, y, z in raw_lattice]
        return lattice, (sx, sy, sz)
    raise ValueError(f"Unsupported POSCAR scale line: {scale_line!r}")


def read_poscar(path: Path) -> PoscarInfo:
    lines = path.read_text().splitlines()
    if len(lines) < 8:
        raise ValueError(f"POSCAR is too short: {path}")
    raw_lattice = [parse_vec3(lines[idx], f"lattice vector {idx - 1}") for idx in range(2, 5)]
    lattice, cart_scale = parse_lattice_scale(lines[1], raw_lattice)
    line5 = lines[5].split()
    line6 = lines[6].split()
    if not line5:
        raise ValueError("POSCAR element/count line is empty")
    if all(is_int_token(tok) for tok in line5):
        counts = [int(tok) for tok in line5]
        counts_idx = 5
        comment_symbols = lines[0].split()
        if len(comment_symbols) == len(counts) and all(re.fullmatch(r"[A-Z][a-z]?", tok) for tok in comment_symbols):
            elements = comment_symbols
        else:
            elements = [f"X{i + 1}" for i in range(len(counts))]
    else:
        elements = line5
        if not line6 or not all(is_int_token(tok) for tok in line6):
            raise ValueError("POSCAR does not look like VASP 5 format with element names")
        counts = [int(tok) for tok in line6]
        counts_idx = 6
    if len(elements) != len(counts):
        raise ValueError(f"Element/count mismatch in POSCAR: {elements} vs {counts}")
    atom_symbols: list[str] = []
    for elem, count in zip(elements, counts):
        atom_symbols.extend([elem] * count)
    coord_idx = counts_idx + 1
    if coord_idx >= len(lines):
        raise ValueError("POSCAR is missing coordinate mode line")
    if lines[coord_idx].strip().lower().startswith("s"):
        coord_idx += 1
    if coord_idx >= len(lines):
        raise ValueError("POSCAR is missing coordinate mode line after Selective dynamics")
    coord_mode = lines[coord_idx].strip().lower()
    coord_start = coord_idx + 1
    if len(lines) < coord_start + len(atom_symbols):
        raise ValueError("POSCAR does not contain enough coordinate rows")
    frac_coords: list[tuple[float, float, float]] = []
    for atom_idx in range(len(atom_symbols)):
        raw = parse_vec3(lines[coord_start + atom_idx], f"atom {atom_idx + 1} coordinate")
        if coord_mode.startswith("d"):
            frac = raw
        elif coord_mode.startswith("c") or coord_mode.startswith("k"):
            cart = (raw[0] * cart_scale[0], raw[1] * cart_scale[1], raw[2] * cart_scale[2])
            frac = cart_to_frac(cart, lattice)
        else:
            raise ValueError(f"Unsupported POSCAR coordinate mode: {lines[coord_idx]!r}")
        frac_coords.append(tuple(x % 1.0 for x in frac))  # type: ignore[arg-type]
    return PoscarInfo(
        path=path,
        elements=elements,
        counts=counts,
        atom_symbols=atom_symbols,
        lattice=lattice,
        frac_coords=frac_coords,
    )


def magnetic_indices(info: PoscarInfo, elems: list[str] | None) -> list[int]:
    if not elems:
        elems = [info.elements[0]]
        print(f"[WARN] --magnetic-elements not set; using first POSCAR element: {elems[0]}", file=sys.stderr)
    missing = [elem for elem in elems if elem not in info.elements]
    if missing:
        raise ValueError(f"Magnetic element(s) not in POSCAR: {', '.join(missing)}")
    return [idx for idx, elem in enumerate(info.atom_symbols) if elem in elems]


def one_based_to_global(idx: int, mode: str, mag: list[int], natoms: int) -> tuple[int, int]:
    if idx <= 0:
        raise ValueError("Atom indices are 1-based and must be positive")
    if mode == "global":
        global_idx = idx - 1
        if global_idx < 0 or global_idx >= natoms:
            raise ValueError(f"Global atom index {idx} outside 1..{natoms}")
        if global_idx not in mag:
            raise ValueError(f"Global atom index {idx} is not in the selected magnetic atom list")
        mag_ord = mag.index(global_idx) + 1
        return global_idx, mag_ord
    if mode == "magnetic":
        if idx > len(mag):
            raise ValueError(f"Magnetic atom index {idx} outside 1..{len(mag)}")
        return mag[idx - 1], idx
    raise ValueError(f"Unknown index mode: {mode}")


def make_pair(i: int, j: int, mode: str, mag: list[int], natoms: int) -> PairInfo:
    gi, mi = one_based_to_global(i, mode, mag, natoms)
    gj, mj = one_based_to_global(j, mode, mag, natoms)
    if gi == gj:
        raise ValueError("Pair atoms must be different")
    return PairInfo(label_i=i, label_j=j, global_i=gi, global_j=gj, magnetic_i=mi, magnetic_j=mj)


def parse_pairs(raw_pairs: list[str], mode: str, mag: list[int], natoms: int) -> list[PairInfo]:
    pairs: list[PairInfo] = []
    for raw in raw_pairs:
        for item in raw.split(","):
            item = item.strip()
            if not item:
                continue
            match = re.fullmatch(r"(\d+)\s*[-:]\s*(\d+)", item)
            if not match:
                raise ValueError(f"Bad pair format '{item}'. Use I-J or repeat --pair I-J.")
            pairs.append(make_pair(int(match.group(1)), int(match.group(2)), mode, mag, natoms))
    return pairs


def atom_label(info: PoscarInfo, global_idx: int) -> str:
    return f"{info.atom_symbols[global_idx]}{global_idx + 1}"


def magnetic_ord(mag: list[int], global_idx: int) -> int:
    return mag.index(global_idx) + 1 if global_idx in mag else 0


def translation_search_extent(info: PoscarInfo, cutoff: float) -> int:
    lengths = [vec_norm(vec) for vec in info.lattice]
    heights = [height for height in lattice_heights(info.lattice) if height > 1e-12]
    scale = min(lengths + heights)
    if scale <= 1e-12:
        raise ValueError("Cannot determine POSCAR cell size")
    return max(1, int(math.ceil(cutoff / scale)) + 1)


def pair_distance_for_shift(info: PoscarInfo, pair: PairInfo, shift: tuple[int, int, int]) -> float:
    disp_frac = (
        info.frac_coords[pair.global_j][0] + shift[0] - info.frac_coords[pair.global_i][0],
        info.frac_coords[pair.global_j][1] + shift[1] - info.frac_coords[pair.global_i][1],
        info.frac_coords[pair.global_j][2] + shift[2] - info.frac_coords[pair.global_i][2],
    )
    return vec_norm(frac_to_cart(disp_frac, info.lattice))


def pair_bond_context(info: PoscarInfo, pair: PairInfo, args: argparse.Namespace) -> PairBondContext:
    requested_shift = tuple(getattr(args, "pair_image_shift", None) or ())
    nearest_shift, _, _, nearest_distance = nearest_pair_image(info, pair)
    if requested_shift:
        if len(requested_shift) != 3:
            raise ValueError("--pair-image-shift requires exactly three integers")
        image_shift = requested_shift  # type: ignore[assignment]
        distance = pair_distance_for_shift(info, pair, image_shift)
    else:
        image_shift = nearest_shift
        distance = nearest_distance
    tol = getattr(args, "bond_distance_tol", 0.02)
    if requested_shift and distance > nearest_distance + tol and not getattr(args, "allow_periodic_bond_sum", False):
        raise ValueError(
            f"--pair-image-shift {format_shift(image_shift)} is not the nearest image for {pair.label}; "
            f"nearest is {format_shift(nearest_shift)} at {nearest_distance:.6f} A. Four-state MAGMOM/M_CONSTR "
            "selects POSCAR atom indices, not a specific farther periodic image. Use the nearest explicit pair "
            "or enable --allow-periodic-bond-sum only for exploratory summed/geometry analysis."
        )
    extent = translation_search_extent(info, distance + tol + 0.5)
    equivalent_shifts: list[tuple[int, int, int]] = []
    for tx in range(-extent, extent + 1):
        for ty in range(-extent, extent + 1):
            for tz in range(-extent, extent + 1):
                shift = (tx, ty, tz)
                cand_distance = pair_distance_for_shift(info, pair, shift)
                if abs(cand_distance - distance) <= tol:
                    equivalent_shifts.append(shift)
    equivalent_shifts.sort(key=lambda item: (abs(item[0]) + abs(item[1]) + abs(item[2]), item))
    if image_shift not in equivalent_shifts:
        equivalent_shifts.insert(0, image_shift)
    return PairBondContext(
        image_shift=image_shift,
        distance=distance,
        nearest_image_shift=nearest_shift,
        nearest_distance=nearest_distance,
        multiplicity=max(1, len(equivalent_shifts)),
        equivalent_shifts=equivalent_shifts,
    )


def validate_pair_bond_context(kind: str, pair: PairInfo, bond: PairBondContext, args: argparse.Namespace) -> None:
    if bond.multiplicity <= 1:
        return
    if getattr(args, "allow_periodic_bond_sum", False):
        print(
            f"[WARN] {kind} {pair.label}: POSCAR atom pair has {bond.multiplicity} equal-distance periodic images. "
            "The extracted value is the sum over these periodic translations; no multiplicity division is applied.",
            file=sys.stderr,
        )
        return
    shifts = ", ".join(format_shift(shift) for shift in bond.equivalent_shifts)
    raise ValueError(
        f"{kind} {pair.label} is not a unique explicit pair: {bond.multiplicity} periodic images are within "
        f"--bond-distance-tol at {bond.distance:.6f} A ({shifts}). Four-state MAGMOM/M_CONSTR selects POSCAR "
        "atom indices, not a specific image shift, so the energy maps to a sum over periodic translations. "
        "Build a larger supercell with a unique explicit target pair, or rerun with --allow-periodic-bond-sum "
        "only for exploratory summed couplings."
    )


def bond_formula_metadata(bond: PairBondContext, args: argparse.Namespace) -> dict[str, object]:
    if bond.multiplicity > 1 and getattr(args, "allow_periodic_bond_sum", False):
        note = (
            "Selected POSCAR atom pair has multiple equal-distance periodic images. "
            "Four-state extraction reports their translation sum; denominator is not divided by multiplicity."
        )
    else:
        note = "Selected POSCAR atom pair is unique within --bond-distance-tol."
    return {
        "bond_multiplicity": bond.multiplicity,
        "pair_image_shift": bond.image_shift,
        "nearest_pair_image_shift": bond.nearest_image_shift,
        "nearest_pair_distance_A": bond.nearest_distance,
        "equivalent_pair_shifts": bond.equivalent_shifts,
        "pair_distance_A": bond.distance,
        "pair_extraction_note": note,
    }


def spin_scale(args: argparse.Namespace, power: int) -> float:
    if getattr(args, "spin_convention", "unit_vector") == "unit_vector":
        return 1.0
    spin_length = float(getattr(args, "spin_length_S", 1.0))
    if spin_length <= 0:
        raise ValueError("--spin-length-S must be positive")
    return spin_length**power


def hamiltonian_prefactor(args: argparse.Namespace) -> float:
    return -1.0 if getattr(args, "hamiltonian_sign", "plus") == "minus" else 1.0


def energy_denominator(
    args: argparse.Namespace,
    base_denominator: float,
    spin_power: int,
    bond_multiplicity: int = 1,
) -> float:
    if bond_multiplicity != 1:
        raise ValueError("periodic bond multiplicity is not a valid denominator factor")
    return base_denominator * spin_scale(args, spin_power)


def sia_energy_denominator(component: str) -> float:
    return 2.0 if component in {"Ayy_minus_Axx", "Azz_minus_Axx"} else 4.0


def format_shifts(shifts: list[tuple[int, int, int]]) -> str:
    return ";".join(format_shift(shift) for shift in shifts)


def cluster_contacts(contacts: list[NeighborContact], shell_tol: float) -> list[NeighborContact]:
    sorted_contacts = sorted(contacts, key=lambda c: (c.distance, c.i_global, c.j_global, c.image_shift))
    shell_refs: list[float] = []
    assigned: list[NeighborContact] = []
    for contact in sorted_contacts:
        shell = None
        for idx, ref in enumerate(shell_refs, start=1):
            if abs(contact.distance - ref) <= shell_tol:
                shell = idx
                break
        if shell is None:
            shell_refs.append(contact.distance)
            shell = len(shell_refs)
        assigned.append(replace(contact, shell=shell))
    return assigned


def find_neighbor_contacts(
    info: PoscarInfo,
    mag: list[int],
    cutoff: float,
    shell_tol: float,
) -> list[NeighborContact]:
    extent = translation_search_extent(info, cutoff)
    contacts: list[NeighborContact] = []
    for pos_i, i_global in enumerate(mag):
        for j_global in mag[pos_i + 1 :]:
            for tx in range(-extent, extent + 1):
                for ty in range(-extent, extent + 1):
                    for tz in range(-extent, extent + 1):
                        shift = (tx, ty, tz)
                        disp_frac = (
                            info.frac_coords[j_global][0] + tx - info.frac_coords[i_global][0],
                            info.frac_coords[j_global][1] + ty - info.frac_coords[i_global][1],
                            info.frac_coords[j_global][2] + tz - info.frac_coords[i_global][2],
                        )
                        disp_cart = frac_to_cart(disp_frac, info.lattice)
                        distance = vec_norm(disp_cart)
                        if distance < 1e-8 or distance > cutoff:
                            continue
                        contacts.append(
                            NeighborContact(
                                shell=0,
                                distance=distance,
                                i_global=i_global,
                                j_global=j_global,
                                i_mag=magnetic_ord(mag, i_global),
                                j_mag=magnetic_ord(mag, j_global),
                                image_shift=shift,
                                disp_frac=disp_frac,
                                disp_cart=disp_cart,
                            )
                        )
    return cluster_contacts(contacts, shell_tol)


def find_self_image_contacts(info: PoscarInfo, mag: list[int], cutoff: float) -> list[tuple[int, tuple[int, int, int], float]]:
    extent = translation_search_extent(info, cutoff)
    contacts: list[tuple[int, tuple[int, int, int], float]] = []
    for global_idx in mag:
        for tx in range(-extent, extent + 1):
            for ty in range(-extent, extent + 1):
                for tz in range(-extent, extent + 1):
                    shift = (tx, ty, tz)
                    if shift == (0, 0, 0):
                        continue
                    disp_cart = frac_to_cart(shift, info.lattice)
                    distance = vec_norm(disp_cart)
                    if 1e-8 < distance <= cutoff:
                        contacts.append((global_idx, shift, distance))
    return sorted(contacts, key=lambda item: item[2])


def default_center_atom(info: PoscarInfo, mag: list[int]) -> int:
    center = (0.5, 0.5, 0.5)
    return min(mag, key=lambda idx: vec_dot(vec_sub(info.frac_coords[idx], center), vec_sub(info.frac_coords[idx], center)))


def representative_contacts(contacts: list[NeighborContact], center_global: int | None) -> list[NeighborContact]:
    reps: list[NeighborContact] = []
    for shell in sorted({contact.shell for contact in contacts}):
        group = [contact for contact in contacts if contact.shell == shell]
        candidates = group
        if center_global is not None:
            centered = [c for c in group if c.i_global == center_global or c.j_global == center_global]
            if centered:
                candidates = centered
        nonboundary = [c for c in candidates if not c.crosses_boundary]
        reps.append(sorted(nonboundary or candidates, key=lambda c: (c.distance, c.crosses_boundary, c.i_global, c.j_global))[0])
    return reps


def suggested_expansion_for_radius(info: PoscarInfo, radius: float, margin: float) -> list[int]:
    needed = 2.0 * radius + margin
    mults: list[int] = []
    for height in lattice_heights(info.lattice):
        if height <= 1e-12:
            mults.append(1)
        else:
            mults.append(max(1, int(math.ceil(needed / height))))
    return mults


def neighbor_warnings(
    info: PoscarInfo,
    contacts: list[NeighborContact],
    self_contacts: list[tuple[int, tuple[int, int, int], float]],
    center_global: int | None,
    margin: float,
) -> list[str]:
    warnings: list[str] = []
    if self_contacts:
        atom, shift, distance = self_contacts[0]
        warnings.append(
            f"Periodic self-image within cutoff: {atom_label(info, atom)} shift={shift} distance={distance:.4f} A. "
            "This is a finite-size warning for long-range exchange."
        )
    for shell in sorted({contact.shell for contact in contacts}):
        group = [contact for contact in contacts if contact.shell == shell]
        radius = sum(c.distance for c in group) / len(group)
        mults = suggested_expansion_for_radius(info, radius, margin)
        if any(mult > 1 for mult in mults):
            warnings.append(
                f"Shell {shell} at {radius:.4f} A may need a larger supercell for a centered four-state pair; "
                f"heuristic minimum expansion is {mults[0]}x{mults[1]}x{mults[2]}."
            )
        if center_global is not None:
            centered = [c for c in group if c.i_global == center_global or c.j_global == center_global]
            centered_boundary = [c for c in centered if c.crosses_boundary]
            if centered_boundary:
                warnings.append(
                    f"Shell {shell}: {len(centered_boundary)}/{len(centered)} contacts around center "
                    f"{atom_label(info, center_global)} cross the periodic boundary. "
                    "Do not skip this shell by choosing only in-cell pairs; enlarge the cell or choose a more central equivalent."
                )
        if group and all(c.crosses_boundary for c in group):
            warnings.append(
                f"Shell {shell}: every contact uses a periodic image. The current cell is too small for an in-cell representative pair."
            )
    return warnings


def format_shift(shift: tuple[int, int, int]) -> str:
    return f"{shift[0]},{shift[1]},{shift[2]}"


def write_neighbor_outputs(
    out_dir: Path,
    info: PoscarInfo,
    contacts: list[NeighborContact],
    reps: list[NeighborContact],
    warnings: list[str],
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    contact_lines = [
        "shell\tdistance_A\tpair\tglobal_i\tglobal_j\tmagnetic_i\tmagnetic_j\timage_shift_j\tdx_A\tdy_A\tdz_A\tcrosses_boundary\n"
    ]
    for c in contacts:
        contact_lines.append(
            f"{c.shell}\t{c.distance:.8f}\t{c.pair_label}\t{c.i_global + 1}\t{c.j_global + 1}\t"
            f"{c.i_mag}\t{c.j_mag}\t{format_shift(c.image_shift)}\t"
            f"{c.disp_cart[0]:.8f}\t{c.disp_cart[1]:.8f}\t{c.disp_cart[2]:.8f}\t{int(c.crosses_boundary)}\n"
        )
    (out_dir / "neighbor_contacts.tsv").write_text("".join(contact_lines))
    rep_lines = ["shell\tdistance_A\trepresentative_pair\tglobal_i\tglobal_j\timage_shift_j\tcrosses_boundary\n"]
    for c in reps:
        rep_lines.append(
            f"{c.shell}\t{c.distance:.8f}\t{c.pair_label}\t{c.i_global + 1}\t{c.j_global + 1}\t"
            f"{format_shift(c.image_shift)}\t{int(c.crosses_boundary)}\n"
        )
    (out_dir / "neighbor_representatives.tsv").write_text("".join(rep_lines))
    (out_dir / "supercell_warnings.txt").write_text("\n".join(warnings) + ("\n" if warnings else ""))


def neighbors(args: argparse.Namespace) -> None:
    info = read_poscar(Path(args.poscar))
    mag = magnetic_indices(info, args.magnetic_elements)
    if args.center_atom:
        center_global, _ = one_based_to_global(args.center_atom, args.index_mode, mag, info.natoms)
    else:
        center_global = default_center_atom(info, mag)
    contacts = find_neighbor_contacts(info, mag, args.cutoff, args.shell_tol)
    self_contacts = find_self_image_contacts(info, mag, args.cutoff)
    reps = representative_contacts(contacts, center_global)
    warnings = neighbor_warnings(info, contacts, self_contacts, center_global, args.boundary_margin)

    print(f"# POSCAR: {Path(args.poscar).resolve()}")
    print(f"# Magnetic elements: {', '.join(args.magnetic_elements or [info.elements[0]])}")
    print(f"# Magnetic atoms global: {' '.join(str(idx + 1) for idx in mag)}")
    print(f"# Cell heights A: {' '.join(f'{height:.4f}' for height in lattice_heights(info.lattice))}")
    print(f"# Center atom for representatives: {atom_label(info, center_global)} (global {center_global + 1}, magnetic {magnetic_ord(mag, center_global)})")
    print(f"\nRepresentative shells within {args.cutoff:.2f} A:")
    print("shell distance_A count boundary_count representative_pair image_shift_j suggested_expansion")
    for rep in reps:
        group = [c for c in contacts if c.shell == rep.shell]
        boundary_count = sum(1 for c in group if c.crosses_boundary)
        mults = suggested_expansion_for_radius(info, rep.distance, args.boundary_margin)
        print(
            f"{rep.shell} {rep.distance:.4f} {len(group)} {boundary_count} {rep.pair_label} "
            f"{format_shift(rep.image_shift)} {mults[0]}x{mults[1]}x{mults[2]}"
        )
    print(f"\nAll magnetic-pair contacts within {args.cutoff:.2f} A:")
    print("shell distance_A pair global_i global_j magnetic_i magnetic_j image_shift_j boundary")
    for row_no, c in enumerate(contacts, start=1):
        if row_no > args.max_rows:
            print(f"... truncated at {args.max_rows} rows; use --out to write full TSV.")
            break
        print(
            f"{c.shell} {c.distance:.4f} {c.pair_label} {c.i_global + 1} {c.j_global + 1} "
            f"{c.i_mag} {c.j_mag} {format_shift(c.image_shift)} {int(c.crosses_boundary)}"
        )
    if warnings:
        print("\nSupercell warnings:")
        for warning in warnings:
            print(f"- {warning}")
    else:
        print("\nSupercell warnings: none from the current heuristic.")
    if args.out:
        write_neighbor_outputs(Path(args.out), info, contacts, reps, warnings)
        print(f"\n[OK] Wrote neighbor analysis to {Path(args.out).resolve()}")


def ligand_indices(info: PoscarInfo, mag: list[int], ligand_elements: list[str] | None) -> list[int]:
    if ligand_elements:
        missing = [elem for elem in ligand_elements if elem not in info.elements]
        if missing:
            raise ValueError(f"Ligand element(s) not in POSCAR: {', '.join(missing)}")
        return [idx for idx, elem in enumerate(info.atom_symbols) if elem in ligand_elements]
    mag_set = set(mag)
    ligands = [idx for idx in range(info.natoms) if idx not in mag_set]
    if not ligands:
        raise ValueError("No ligand atoms inferred; pass --ligand-elements explicitly")
    return ligands


def canonical_direction(vec: Iterable[float]) -> tuple[float, float, float]:
    unit = vec_normalize(vec)
    max_idx = max(range(3), key=lambda idx: abs(unit[idx]))
    if unit[max_idx] < 0:
        unit = neg_vec(unit)
    return unit


def cartesian_axis_label(vec: Iterable[float]) -> str:
    unit = vec_normalize(vec)
    labels = ["x", "y", "z"]
    return labels[max(range(3), key=lambda idx: abs(unit[idx]))]


def nearest_pair_image(
    info: PoscarInfo,
    pair: PairInfo,
) -> tuple[tuple[int, int, int], tuple[float, float, float], tuple[float, float, float], float]:
    disp_frac, shift = nearest_image_delta(info.frac_coords[pair.global_i], info.frac_coords[pair.global_j], info.lattice)
    disp_cart = frac_to_cart(disp_frac, info.lattice)
    return shift, disp_frac, disp_cart, vec_norm(disp_cart)


def find_shared_ligands(
    info: PoscarInfo,
    pair: PairInfo,
    ligand_global_indices: list[int],
    cutoff: float,
    pair_shift: tuple[int, int, int],
) -> list[tuple[int, tuple[int, int, int], float, float, tuple[float, float, float]]]:
    extent = translation_search_extent(info, cutoff)
    i_frac = info.frac_coords[pair.global_i]
    j_frac = vec_add(info.frac_coords[pair.global_j], pair_shift)
    shared: list[tuple[int, tuple[int, int, int], float, float, tuple[float, float, float]]] = []
    for lig_idx in ligand_global_indices:
        for tx in range(-extent, extent + 1):
            for ty in range(-extent, extent + 1):
                for tz in range(-extent, extent + 1):
                    lig_shift = (tx, ty, tz)
                    lig_frac = vec_add(info.frac_coords[lig_idx], lig_shift)
                    vi_frac = vec_sub(lig_frac, i_frac)
                    vj_frac = vec_sub(lig_frac, j_frac)
                    vi_cart = frac_to_cart(vi_frac, info.lattice)
                    vj_cart = frac_to_cart(vj_frac, info.lattice)
                    di = vec_norm(vi_cart)
                    dj = vec_norm(vj_cart)
                    if di <= cutoff and dj <= cutoff:
                        shared.append((lig_idx, lig_shift, di, dj, lig_frac))
    shared.sort(key=lambda item: (item[2] + item[3], max(item[2], item[3]), item[0], item[1]))
    return shared


def nearest_ligand_images(
    info: PoscarInfo,
    center_frac: tuple[float, float, float],
    ligand_global_indices: list[int],
    count: int = 6,
) -> list[tuple[int, tuple[int, int, int], float, tuple[float, float, float], tuple[float, float, float]]]:
    ranked: list[tuple[int, tuple[int, int, int], float, tuple[float, float, float], tuple[float, float, float]]] = []
    for lig_idx in ligand_global_indices:
        disp_frac, shift = nearest_image_delta(center_frac, info.frac_coords[lig_idx], info.lattice)
        lig_frac = vec_add(info.frac_coords[lig_idx], shift)
        disp_cart = frac_to_cart(disp_frac, info.lattice)
        ranked.append((lig_idx, shift, vec_norm(disp_cart), lig_frac, disp_cart))
    ranked.sort(key=lambda item: (item[2], item[0], item[1]))
    return ranked[:count]


def rank_ligand_references(
    info: PoscarInfo,
    i_frac: tuple[float, float, float],
    j_frac: tuple[float, float, float],
    ligand_global_indices: list[int],
    search_cutoff: float,
) -> list[tuple[int, tuple[int, int, int], float, float, tuple[float, float, float]]]:
    extent = translation_search_extent(info, search_cutoff)
    ranked: list[tuple[int, tuple[int, int, int], float, float, tuple[float, float, float]]] = []
    for lig_idx in ligand_global_indices:
        for tx in range(-extent, extent + 1):
            for ty in range(-extent, extent + 1):
                for tz in range(-extent, extent + 1):
                    lig_shift = (tx, ty, tz)
                    lig_frac = vec_add(info.frac_coords[lig_idx], lig_shift)
                    vi_cart = frac_to_cart(vec_sub(lig_frac, i_frac), info.lattice)
                    vj_cart = frac_to_cart(vec_sub(lig_frac, j_frac), info.lattice)
                    di = vec_norm(vi_cart)
                    dj = vec_norm(vj_cart)
                    ranked.append((lig_idx, lig_shift, di, dj, lig_frac))
    ranked.sort(key=lambda item: (item[2] + item[3], max(item[2], item[3]), item[0], item[1]))
    return ranked


def choose_kitaev_ligands(
    info: PoscarInfo,
    pair: PairInfo,
    ligand_global_indices: list[int],
    cutoff: float,
    pair_shift: tuple[int, int, int],
    pair_distance: float,
    allow_fallback: bool = False,
) -> tuple[list[tuple[int, tuple[int, int, int], float, float, tuple[float, float, float]]], str]:
    shared = find_shared_ligands(info, pair, ligand_global_indices, cutoff, pair_shift)
    if len(shared) >= 2:
        return shared[:2], "shared_within_cutoff"

    if not allow_fallback:
        raise ValueError(
            f"{pair.label} cannot confirm two shared ligand images within {cutoff:.3f} A. "
            "Pass --ligand-elements/--metal-ligand-cutoff explicitly, or use "
            "--allow-kitaev-ligand-fallback only for exploratory geometry checks."
        )

    i_frac = info.frac_coords[pair.global_i]
    j_frac = vec_add(info.frac_coords[pair.global_j], pair_shift)
    nearest_i = nearest_ligand_images(info, i_frac, ligand_global_indices, 6)
    nearest_j = nearest_ligand_images(info, j_frac, ligand_global_indices, 6)
    shared_nearest = sorted({item[0] for item in nearest_i} & {item[0] for item in nearest_j})
    fallback_cutoff = max(cutoff, pair_distance + cutoff, 6.0)
    if len(shared_nearest) >= 2:
        ranked = rank_ligand_references(info, i_frac, j_frac, shared_nearest, fallback_cutoff)
        if len(ranked) >= 2:
            return ranked[:2], "nearest_six_shared"

    ranked = rank_ligand_references(info, i_frac, j_frac, ligand_global_indices, fallback_cutoff)
    if len(ranked) >= 2:
        return ranked[:2], "distance_sum_fallback"
    raise ValueError(f"Metal pair {pair.label} cannot find two ligand references")


def angle_error_from_orthogonal(vec_a: tuple[float, float, float], vec_b: tuple[float, float, float]) -> float:
    cos_theta = max(-1.0, min(1.0, vec_dot(vec_normalize(vec_a), vec_normalize(vec_b))))
    return abs(math.degrees(math.acos(cos_theta)) - 90.0)


def gram_schmidt_rows(rows: list[tuple[float, float, float]]) -> list[tuple[float, float, float]]:
    basis: list[tuple[float, float, float]] = []
    for row_no, row in enumerate(rows, start=1):
        vec = row
        for prev in basis:
            vec = vec_sub(vec, vec_scale(prev, vec_dot(vec, prev)))
        basis.append(vec_normalize(vec, f"octahedral basis row {row_no}"))
    return basis


def find_octahedral_basis(
    info: PoscarInfo,
    center_global: int,
    ligand_global_indices: list[int],
    cutoff: float,
    center_frac_override: tuple[float, float, float] | None = None,
) -> tuple[list[tuple[float, float, float]], list[tuple[int, tuple[int, int, int], float]], float] | None:
    center_frac = center_frac_override or info.frac_coords[center_global]
    shell = nearest_ligand_images(info, center_frac, ligand_global_indices, 12)
    local_shell = [item for item in shell if item[2] <= cutoff]
    if len(local_shell) < 6:
        local_shell = shell[:6]
    if len(local_shell) < 3:
        return None

    best_combo = None
    best_score = float("inf")
    best_max_error = float("inf")
    for combo in combinations(local_shell, 3):
        vecs = [item[4] for item in combo]
        errors = [angle_error_from_orthogonal(a, b) for a, b in combinations(vecs, 2)]
        score = sum(errors)
        max_error = max(errors)
        if (score, max_error) < (best_score, best_max_error):
            best_combo = combo
            best_score = score
            best_max_error = max_error

    if best_combo is None:
        return None
    raw_basis = [item[4] for item in best_combo]
    basis = gram_schmidt_rows(raw_basis)
    ligand_summary = [(item[0], item[1], item[2]) for item in best_combo]
    return basis, ligand_summary, best_max_error


def construct_shared_ligand_reference_basis(
    info: PoscarInfo,
    center_frac: tuple[float, float, float],
    chosen_ligands: list[tuple[int, tuple[int, int, int], float, float, tuple[float, float, float]]],
    octahedral_basis: list[tuple[float, float, float]] | None,
) -> list[tuple[float, float, float]]:
    ligand_vecs = [frac_to_cart(vec_sub(item[4], center_frac), info.lattice) for item in chosen_ligands[:2]]
    x_axis = vec_normalize(ligand_vecs[0], "first shared-ligand vector")
    second_axis = vec_normalize(ligand_vecs[1], "second shared-ligand vector")
    y_axis = vec_sub(second_axis, vec_scale(x_axis, vec_dot(second_axis, x_axis)))
    y_axis = vec_normalize(y_axis, "orthogonalized second shared-ligand vector")
    z_axis = vec_normalize(vec_cross(x_axis, y_axis), "shared-ligand plane normal")

    if octahedral_basis:
        oct_axis = max(octahedral_basis, key=lambda axis: abs(vec_dot(axis, z_axis)))
        if vec_dot(oct_axis, z_axis) < 0:
            y_axis = neg_vec(y_axis)
            z_axis = neg_vec(z_axis)
    return [x_axis, y_axis, z_axis]


def align_ideal_kitaev_basis(
    reference_basis: list[tuple[float, float, float]],
    allow_component_permutation: bool = True,
) -> tuple[list[tuple[float, float, float]], tuple[float, float, float], str]:
    reference = [vec_normalize(row, f"reference basis row {idx + 1}") for idx, row in enumerate(reference_basis)]
    ideal = [vec_normalize(row, f"ideal Kitaev basis row {idx + 1}") for idx, row in enumerate(IDEAL_KITAEV_BASIS)]
    if allow_component_permutation:
        basis = [
            basis_linear_combination(row, reference, f"Kitaev {label} axis")
            for label, row in zip(KITAEV_AXIS_LABELS, ideal)
        ]
        overlaps = tuple(abs(vec_dot(basis[idx], reference[idx])) for idx in range(3))
        return basis, overlaps, "continuous_reference_projection"

    component_perms = list(permutations(range(3))) if allow_component_permutation else [(0, 1, 2)]
    component_sign_sets = list(product((-1.0, 1.0), repeat=3)) if allow_component_permutation else [(1.0, 1.0, 1.0)]

    best_score = -float("inf")
    best_basis: list[tuple[float, float, float]] | None = None
    best_overlaps: tuple[float, float, float] | None = None
    best_row_perm: tuple[int, int, int] | None = None
    best_row_signs: tuple[float, float, float] | None = None
    best_component_perm: tuple[int, int, int] | None = None
    best_component_signs: tuple[float, float, float] | None = None

    for row_perm in permutations(range(3)):
        permuted = [ideal[idx] for idx in row_perm]
        for row_signs in product((-1.0, 1.0), repeat=3):
            row_signed = [vec_scale(permuted[idx], row_signs[idx]) for idx in range(3)]
            for component_perm in component_perms:
                for component_signs in component_sign_sets:
                    candidate = [
                        tuple(row[component_perm[col]] * component_signs[col] for col in range(3))
                        for row in row_signed
                    ]
                    overlaps = tuple(vec_dot(candidate[idx], reference[idx]) for idx in range(3))
                    score = sum(overlaps)
                    if score > best_score:
                        best_score = score
                        best_basis = candidate
                        best_overlaps = overlaps  # type: ignore[assignment]
                        best_row_perm = row_perm  # type: ignore[assignment]
                        best_row_signs = row_signs  # type: ignore[assignment]
                        best_component_perm = component_perm  # type: ignore[assignment]
                        best_component_signs = component_signs  # type: ignore[assignment]

    if best_basis is None or best_overlaps is None:
        raise ValueError("Could not align ideal Kitaev basis to reference basis")
    match = (
        f"ideal_rows_1based={','.join(str(idx + 1) for idx in best_row_perm or ())};"
        f"row_signs={','.join(str(int(sign)) for sign in best_row_signs or ())};"
        f"component_order_1based={','.join(str(idx + 1) for idx in best_component_perm or ())};"
        f"component_signs={','.join(str(int(sign)) for sign in best_component_signs or ())}"
    )
    return best_basis, best_overlaps, match


def choose_gamma_from_bond(
    kitaev_axes: list[tuple[str, tuple[float, float, float]]],
    bond_axis: tuple[float, float, float],
) -> tuple[int, int, int, tuple[float, float, float]]:
    axis_bond_dots = tuple(abs(vec_dot(axis, bond_axis)) for _, axis in kitaev_axes)
    gamma_idx = min(range(3), key=lambda idx: axis_bond_dots[idx])
    alpha_idx, beta_idx = [idx for idx in range(3) if idx != gamma_idx]
    alpha_axis = kitaev_axes[alpha_idx][1]
    beta_axis = kitaev_axes[beta_idx][1]
    gamma_axis = kitaev_axes[gamma_idx][1]
    if vec_dot(vec_cross(alpha_axis, beta_axis), gamma_axis) < 0:
        alpha_idx, beta_idx = beta_idx, alpha_idx
    return alpha_idx, beta_idx, gamma_idx, axis_bond_dots  # type: ignore[return-value]


def fallback_perpendicular_axis(axis: tuple[float, float, float]) -> tuple[float, float, float]:
    trial = min(AXIS.values(), key=lambda basis: abs(vec_dot(axis, basis)))
    projected = vec_sub(trial, vec_scale(axis, vec_dot(trial, axis)))
    return vec_normalize(projected, "fallback perpendicular axis")


def detect_kitaev_frame(info: PoscarInfo, mag: list[int], pair: PairInfo, args: argparse.Namespace) -> KitaevFrame:
    ligands = ligand_indices(info, mag, getattr(args, "ligand_elements", None))
    bond = pair_bond_context(info, pair, args)
    pair_shift = bond.image_shift
    pair_disp_frac = vec_sub(vec_add(info.frac_coords[pair.global_j], pair_shift), info.frac_coords[pair.global_i])
    pair_disp_cart = frac_to_cart(pair_disp_frac, info.lattice)
    pair_distance = bond.distance
    if pair_distance > getattr(args, "kitaev_pair_cutoff", 10.0):
        raise ValueError(
            f"{pair.label} nearest image distance is {pair_distance:.4f} A, larger than --kitaev-pair-cutoff"
        )
    chosen, ligand_method = choose_kitaev_ligands(
        info,
        pair,
        ligands,
        args.metal_ligand_cutoff,
        pair_shift,
        pair_distance,
        getattr(args, "allow_kitaev_ligand_fallback", False),
    )
    i_frac = info.frac_coords[pair.global_i]
    j_frac = vec_add(info.frac_coords[pair.global_j], pair_shift)
    octahedral_i = find_octahedral_basis(info, pair.global_i, ligands, args.metal_ligand_cutoff, i_frac)
    octahedral_j = find_octahedral_basis(info, pair.global_j, ligands, args.metal_ligand_cutoff, j_frac)
    if octahedral_i:
        octahedral_basis, octahedral_ligands, octahedral_angle_error = octahedral_i
    else:
        octahedral_basis = []
        octahedral_ligands = []
        octahedral_angle_error = None
    if octahedral_j:
        octahedral_basis_j, octahedral_ligands_j, octahedral_angle_error_j = octahedral_j
    else:
        octahedral_basis_j = []
        octahedral_ligands_j = []
        octahedral_angle_error_j = None
    reference_basis_i = construct_shared_ligand_reference_basis(
        info,
        i_frac,
        chosen,
        octahedral_basis if octahedral_basis else None,
    )
    reference_basis_j = construct_shared_ligand_reference_basis(
        info,
        j_frac,
        chosen,
        octahedral_basis_j if octahedral_basis_j else None,
    )
    kitaev_basis_i, axis_overlaps_i, axis_match_i = align_ideal_kitaev_basis(
        reference_basis_i,
        allow_component_permutation=not getattr(args, "kitaev_no_component_permutation", False),
    )
    kitaev_basis_j, axis_overlaps_j, axis_match_j = align_ideal_kitaev_basis(
        reference_basis_j,
        allow_component_permutation=not getattr(args, "kitaev_no_component_permutation", False),
    )
    kitaev_axes_i = list(zip(KITAEV_AXIS_LABELS, kitaev_basis_i))
    kitaev_axes_j = list(zip(KITAEV_AXIS_LABELS, kitaev_basis_j))
    bond_axis = vec_normalize(pair_disp_cart, "metal-metal bond")
    alpha_idx, beta_idx, gamma_idx, axis_bond_dots = choose_gamma_from_bond(kitaev_axes_i, bond_axis)
    alpha_label, local_x = kitaev_axes_i[alpha_idx]
    beta_label, local_y = kitaev_axes_i[beta_idx]
    gamma_label, gamma = kitaev_axes_i[gamma_idx]
    axis_map_j = {label: axis for label, axis in kitaev_axes_j}
    local_x_j = axis_map_j[alpha_label]
    local_y_j = axis_map_j[beta_label]
    gamma_j = axis_map_j[gamma_label]
    axis_consistency = []
    for label, axis_i in kitaev_axes_i:
        angle = vec_angle_degrees(axis_i, axis_map_j[label])
        axis_consistency.append(min(angle, 180.0 - angle))
    max_axis_consistency = getattr(args, "max_kitaev_axis_consistency_deg", 25.0)
    if max_axis_consistency is not None and max(axis_consistency) > float(max_axis_consistency):
        raise ValueError(
            f"Kitaev local axes differ too much between the two sites: max axis consistency angle "
            f"{max(axis_consistency):.3f} deg > --max-kitaev-axis-consistency-deg {float(max_axis_consistency):.3f}. "
            "Use full Jani tensor reporting or inspect the octahedral geometry before compressing to a standard model."
        )
    shared_summary = [(idx, shift, di, dj) for idx, shift, di, dj, _ in chosen]
    return KitaevFrame(
        pair=pair,
        gamma_axis=gamma,
        local_x=local_x,
        local_y=local_y,
        gamma_axis_j=gamma_j,
        local_x_j=local_x_j,
        local_y_j=local_y_j,
        bond_axis=bond_axis,
        pair_shift=pair_shift,
        shared_ligands=shared_summary,
        gamma_label=gamma_label,
        alpha_label=alpha_label,
        beta_label=beta_label,
        kitaev_axes=kitaev_axes_i,
        kitaev_axes_j=kitaev_axes_j,
        reference_basis=reference_basis_i,
        reference_basis_j=reference_basis_j,
        axis_overlaps=axis_overlaps_i,
        axis_overlaps_j=axis_overlaps_j,
        axis_bond_dots=axis_bond_dots,
        axis_consistency_degrees=tuple(axis_consistency),  # type: ignore[arg-type]
        ligand_method=ligand_method,
        axis_match=axis_match_i,
        axis_match_j=axis_match_j,
        octahedral_ligands=octahedral_ligands,
        octahedral_ligands_j=octahedral_ligands_j,
        octahedral_angle_error=octahedral_angle_error,
        octahedral_angle_error_j=octahedral_angle_error_j,
    )


def frame_rows(frame: KitaevFrame, info: PoscarInfo) -> list[str]:
    ligands = ";".join(
        f"{atom_label(info, idx)}@{format_shift(shift)}(d_i={di:.4f},d_j={dj:.4f})"
        for idx, shift, di, dj in frame.shared_ligands
    )
    octahedral_ligands = ";".join(
        f"{atom_label(info, idx)}@{format_shift(shift)}(d={dist:.4f})"
        for idx, shift, dist in frame.octahedral_ligands
    )
    octahedral_ligands_j = ";".join(
        f"{atom_label(info, idx)}@{format_shift(shift)}(d={dist:.4f})"
        for idx, shift, dist in frame.octahedral_ligands_j
    )
    axis_rows = [
        f"kitaev_i_{label}_axis_cart\t{axis[0]:.10f}\t{axis[1]:.10f}\t{axis[2]:.10f}"
        for label, axis in frame.kitaev_axes
    ]
    axis_rows_j = [
        f"kitaev_j_{label}_axis_cart\t{axis[0]:.10f}\t{axis[1]:.10f}\t{axis[2]:.10f}"
        for label, axis in frame.kitaev_axes_j
    ]
    reference_rows = [
        f"reference_i_basis_{idx}_cart\t{axis[0]:.10f}\t{axis[1]:.10f}\t{axis[2]:.10f}"
        for idx, axis in enumerate(frame.reference_basis, start=1)
    ]
    reference_rows_j = [
        f"reference_j_basis_{idx}_cart\t{axis[0]:.10f}\t{axis[1]:.10f}\t{axis[2]:.10f}"
        for idx, axis in enumerate(frame.reference_basis_j, start=1)
    ]
    return [
        f"pair\t{frame.pair.label}",
        f"pair_shift_j\t{format_shift(frame.pair_shift)}",
        f"gamma_label\t{frame.gamma_label}",
        f"alpha_label\t{frame.alpha_label}",
        f"beta_label\t{frame.beta_label}",
        f"gamma_axis_cart_i\t{frame.gamma_axis[0]:.10f}\t{frame.gamma_axis[1]:.10f}\t{frame.gamma_axis[2]:.10f}",
        f"gamma_axis_cart_j\t{frame.gamma_axis_j[0]:.10f}\t{frame.gamma_axis_j[1]:.10f}\t{frame.gamma_axis_j[2]:.10f}",
        f"local_x_cart_i\t{frame.local_x[0]:.10f}\t{frame.local_x[1]:.10f}\t{frame.local_x[2]:.10f}",
        f"local_x_cart_j\t{frame.local_x_j[0]:.10f}\t{frame.local_x_j[1]:.10f}\t{frame.local_x_j[2]:.10f}",
        f"local_y_cart_i\t{frame.local_y[0]:.10f}\t{frame.local_y[1]:.10f}\t{frame.local_y[2]:.10f}",
        f"local_y_cart_j\t{frame.local_y_j[0]:.10f}\t{frame.local_y_j[1]:.10f}\t{frame.local_y_j[2]:.10f}",
        f"bond_axis_cart\t{frame.bond_axis[0]:.10f}\t{frame.bond_axis[1]:.10f}\t{frame.bond_axis[2]:.10f}",
        f"axis_overlaps_i\t{frame.axis_overlaps[0]:.8f}\t{frame.axis_overlaps[1]:.8f}\t{frame.axis_overlaps[2]:.8f}",
        f"axis_overlaps_j\t{frame.axis_overlaps_j[0]:.8f}\t{frame.axis_overlaps_j[1]:.8f}\t{frame.axis_overlaps_j[2]:.8f}",
        f"axis_consistency_deg_xyz\t{frame.axis_consistency_degrees[0]:.6f}\t{frame.axis_consistency_degrees[1]:.6f}\t{frame.axis_consistency_degrees[2]:.6f}",
        f"axis_bond_abs_dots_xyz\t{frame.axis_bond_dots[0]:.8f}\t{frame.axis_bond_dots[1]:.8f}\t{frame.axis_bond_dots[2]:.8f}",
        f"ligand_selection_method\t{frame.ligand_method}",
        f"axis_match_i\t{frame.axis_match}",
        f"axis_match_j\t{frame.axis_match_j}",
        f"octahedral_i_max_orthogonality_error_deg\t{frame.octahedral_angle_error:.6f}" if frame.octahedral_angle_error is not None else "octahedral_i_max_orthogonality_error_deg\tNA",
        f"octahedral_j_max_orthogonality_error_deg\t{frame.octahedral_angle_error_j:.6f}" if frame.octahedral_angle_error_j is not None else "octahedral_j_max_orthogonality_error_deg\tNA",
        f"octahedral_ligands_i\t{octahedral_ligands}",
        f"octahedral_ligands_j\t{octahedral_ligands_j}",
        f"shared_ligands\t{ligands}",
        *axis_rows,
        *axis_rows_j,
        *reference_rows,
        *reference_rows_j,
    ]


def write_kitaev_frames(out_root: Path, frames: list[KitaevFrame], info: PoscarInfo) -> None:
    lines = [
        "pair\tgamma_label\talpha_label\tbeta_label\tpair_shift_j\t"
        "gamma_i_x\tgamma_i_y\tgamma_i_z\tgamma_j_x\tgamma_j_y\tgamma_j_z\t"
        "alpha_i_x\talpha_i_y\talpha_i_z\talpha_j_x\talpha_j_y\talpha_j_z\t"
        "beta_i_x\tbeta_i_y\tbeta_i_z\tbeta_j_x\tbeta_j_y\tbeta_j_z\t"
        "kitaev_i_x_x\tkitaev_i_x_y\tkitaev_i_x_z\tkitaev_i_y_x\tkitaev_i_y_y\tkitaev_i_y_z\tkitaev_i_z_x\tkitaev_i_z_y\tkitaev_i_z_z\t"
        "kitaev_j_x_x\tkitaev_j_x_y\tkitaev_j_x_z\tkitaev_j_y_x\tkitaev_j_y_y\tkitaev_j_y_z\tkitaev_j_z_x\tkitaev_j_z_y\tkitaev_j_z_z\t"
        "overlap_i_x\toverlap_i_y\toverlap_i_z\toverlap_j_x\toverlap_j_y\toverlap_j_z\t"
        "axis_consistency_deg_x\taxis_consistency_deg_y\taxis_consistency_deg_z\tbond_dot_x\tbond_dot_y\tbond_dot_z\t"
        "ligand_method\toctahedral_i_max_orthogonality_error_deg\toctahedral_j_max_orthogonality_error_deg\t"
        "shared_ligands\toctahedral_ligands_i\toctahedral_ligands_j\taxis_match_i\taxis_match_j\n"
    ]
    for frame in frames:
        ligands = ";".join(f"{atom_label(info, idx)}@{format_shift(shift)}" for idx, shift, _, _ in frame.shared_ligands)
        octahedral_ligands = ";".join(
            f"{atom_label(info, idx)}@{format_shift(shift)}" for idx, shift, _ in frame.octahedral_ligands
        )
        octahedral_ligands_j = ";".join(
            f"{atom_label(info, idx)}@{format_shift(shift)}" for idx, shift, _ in frame.octahedral_ligands_j
        )
        axis_map = {label: axis for label, axis in frame.kitaev_axes}
        axis_map_j = {label: axis for label, axis in frame.kitaev_axes_j}
        oct_err_i = f"{frame.octahedral_angle_error:.6f}" if frame.octahedral_angle_error is not None else "NA"
        oct_err_j = f"{frame.octahedral_angle_error_j:.6f}" if frame.octahedral_angle_error_j is not None else "NA"
        lines.append(
            f"{frame.pair.label}\t{frame.gamma_label}\t{frame.alpha_label}\t{frame.beta_label}\t{format_shift(frame.pair_shift)}\t"
            f"{frame.gamma_axis[0]:.10f}\t{frame.gamma_axis[1]:.10f}\t{frame.gamma_axis[2]:.10f}\t"
            f"{frame.gamma_axis_j[0]:.10f}\t{frame.gamma_axis_j[1]:.10f}\t{frame.gamma_axis_j[2]:.10f}\t"
            f"{frame.local_x[0]:.10f}\t{frame.local_x[1]:.10f}\t{frame.local_x[2]:.10f}\t"
            f"{frame.local_x_j[0]:.10f}\t{frame.local_x_j[1]:.10f}\t{frame.local_x_j[2]:.10f}\t"
            f"{frame.local_y[0]:.10f}\t{frame.local_y[1]:.10f}\t{frame.local_y[2]:.10f}\t"
            f"{frame.local_y_j[0]:.10f}\t{frame.local_y_j[1]:.10f}\t{frame.local_y_j[2]:.10f}\t"
            f"{axis_map['x'][0]:.10f}\t{axis_map['x'][1]:.10f}\t{axis_map['x'][2]:.10f}\t"
            f"{axis_map['y'][0]:.10f}\t{axis_map['y'][1]:.10f}\t{axis_map['y'][2]:.10f}\t"
            f"{axis_map['z'][0]:.10f}\t{axis_map['z'][1]:.10f}\t{axis_map['z'][2]:.10f}\t"
            f"{axis_map_j['x'][0]:.10f}\t{axis_map_j['x'][1]:.10f}\t{axis_map_j['x'][2]:.10f}\t"
            f"{axis_map_j['y'][0]:.10f}\t{axis_map_j['y'][1]:.10f}\t{axis_map_j['y'][2]:.10f}\t"
            f"{axis_map_j['z'][0]:.10f}\t{axis_map_j['z'][1]:.10f}\t{axis_map_j['z'][2]:.10f}\t"
            f"{frame.axis_overlaps[0]:.8f}\t{frame.axis_overlaps[1]:.8f}\t{frame.axis_overlaps[2]:.8f}\t"
            f"{frame.axis_overlaps_j[0]:.8f}\t{frame.axis_overlaps_j[1]:.8f}\t{frame.axis_overlaps_j[2]:.8f}\t"
            f"{frame.axis_consistency_degrees[0]:.6f}\t{frame.axis_consistency_degrees[1]:.6f}\t{frame.axis_consistency_degrees[2]:.6f}\t"
            f"{frame.axis_bond_dots[0]:.8f}\t{frame.axis_bond_dots[1]:.8f}\t{frame.axis_bond_dots[2]:.8f}\t"
            f"{frame.ligand_method}\t{oct_err_i}\t{oct_err_j}\t{ligands}\t{octahedral_ligands}\t{octahedral_ligands_j}\t"
            f"{frame.axis_match}\t{frame.axis_match_j}\n"
        )
    (out_root / "kitaev_frames.tsv").write_text("".join(lines))


def mat_vec(matrix: list[list[float]], vec: tuple[float, float, float]) -> tuple[float, float, float]:
    return tuple(sum(matrix[row][col] * vec[col] for col in range(3)) for row in range(3))  # type: ignore[return-value]


def rotate_tensor_to_basis(
    matrix: list[list[float]],
    basis: list[tuple[float, float, float]],
) -> list[list[float]]:
    return [[vec_dot(basis[row], mat_vec(matrix, basis[col])) for col in range(3)] for row in range(3)]


def rotate_tensor_between_bases(
    matrix: list[list[float]],
    left_basis: list[tuple[float, float, float]],
    right_basis: list[tuple[float, float, float]],
) -> list[list[float]]:
    return [[vec_dot(left_basis[row], mat_vec(matrix, right_basis[col])) for col in range(3)] for row in range(3)]


def rotate_tensor_to_frame(
    matrix: list[list[float]],
    frame: KitaevFrame,
) -> list[list[float]]:
    return rotate_tensor_between_bases(
        matrix,
        [frame.local_x, frame.local_y, frame.gamma_axis],
        [frame.local_x_j, frame.local_y_j, frame.gamma_axis_j],
    )


def common_kitaev_basis(frame: KitaevFrame) -> list[tuple[float, float, float]]:
    return [frame.local_x, frame.local_y, frame.gamma_axis]


def kitaev_axis_basis(frame: KitaevFrame) -> list[tuple[float, float, float]]:
    axis_map = {label: axis for label, axis in frame.kitaev_axes}
    return [axis_map[label] for label in KITAEV_AXIS_LABELS]


def kitaev_axis_basis_j(frame: KitaevFrame) -> list[tuple[float, float, float]]:
    axis_map = {label: axis for label, axis in frame.kitaev_axes_j}
    return [axis_map[label] for label in KITAEV_AXIS_LABELS]


def symmetric_matrix(matrix: list[list[float]]) -> list[list[float]]:
    return [[0.5 * (matrix[i][j] + matrix[j][i]) for j in range(3)] for i in range(3)]


def antisymmetric_dmi(matrix: list[list[float]]) -> tuple[float, float, float]:
    return (
        0.5 * (matrix[1][2] - matrix[2][1]),
        0.5 * (matrix[2][0] - matrix[0][2]),
        0.5 * (matrix[0][1] - matrix[1][0]),
    )


def decompose_common_exchange(matrix: list[list[float]]) -> dict[str, float]:
    symmetric = symmetric_matrix(matrix)
    dmi = antisymmetric_dmi(matrix)
    trace_iso = sum(symmetric[i][i] for i in range(3)) / 3.0
    j_ab_avg = 0.5 * (symmetric[0][0] + symmetric[1][1])
    return {
        "J_trace_iso_meV": trace_iso,
        "J_alpha_beta_avg_meV": j_ab_avg,
        "J_gamma_gamma_meV": symmetric[2][2],
        "traceless_gamma_anisotropy_meV": symmetric[2][2] - trace_iso,
        "K_gamma_minus_alpha_beta_avg_meV": symmetric[2][2] - j_ab_avg,
        "Gamma_alpha_beta_meV": symmetric[0][1],
        "Gamma_prime_avg_meV": 0.5 * (symmetric[0][2] + symmetric[1][2]),
        "Gamma_prime_split_meV": 0.5 * (symmetric[0][2] - symmetric[1][2]),
        "alpha_beta_diag_split_meV": 0.5 * (symmetric[0][0] - symmetric[1][1]),
        "DMI_alpha_meV": dmi[0],
        "DMI_beta_meV": dmi[1],
        "DMI_gamma_meV": dmi[2],
    }


def kitaev_exchange_matrices(
    global_matrix: list[list[float]],
    frame: KitaevFrame,
) -> tuple[list[list[float]], list[list[float]], list[list[float]]]:
    common_matrix = rotate_tensor_to_basis(global_matrix, common_kitaev_basis(frame))
    local_gauge_matrix = rotate_tensor_to_frame(global_matrix, frame)
    kitaev_gauge_matrix = rotate_tensor_between_bases(
        global_matrix,
        kitaev_axis_basis(frame),
        kitaev_axis_basis_j(frame),
    )
    return common_matrix, local_gauge_matrix, kitaev_gauge_matrix


def unique_summary_stages(summary: Path) -> list[str]:
    stages: list[str] = []
    seen: set[str] = set()
    for raw in summary.read_text().splitlines():
        line = raw.strip()
        if line.startswith("[") and line.endswith("]"):
            stage = line.strip("[]")
            if stage not in seen:
                stages.append(stage)
                seen.add(stage)
    return stages


def parse_jani_summary(root: Path, pair_label: str, stage_name: str | None) -> tuple[list[list[float]], str | None]:
    summary = root / "final_summary.txt"
    if not summary.exists():
        raise FileNotFoundError(f"Jani final_summary.txt not found: {summary}")
    stages = unique_summary_stages(summary)
    selected_stage = stage_name
    if selected_stage is None and len(stages) > 1:
        raise ValueError(
            f"{summary} contains multiple stages ({', '.join(stages)}); pass --stage explicitly, "
            "for example --stage hse_no_u"
        )
    if selected_stage is None and len(stages) == 1:
        selected_stage = stages[0]
    if selected_stage and stages and selected_stage not in stages:
        raise ValueError(f"Stage {selected_stage!r} not found in {summary}; available stages: {', '.join(stages)}")
    current_stage: str | None = None
    values: dict[str, float] = {}
    for raw in summary.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            current_stage = line.strip("[]")
            continue
        if line.startswith("#"):
            continue
        if selected_stage and current_stage != selected_stage:
            continue
        parts = line.split()
        if len(parts) < 7:
            continue
        label, quantity, value_text = parts[0], parts[1], parts[-1]
        for comp in JANI_COMPONENTS:
            if quantity == f"{comp}_meV" and (label == comp or label == f"{pair_label}_{comp}"):
                values[comp] = float(value_text)
    missing = [comp for comp in JANI_COMPONENTS if comp not in values]
    if missing:
        raise ValueError(f"Missing Jani components in {summary}: {', '.join(missing)}")
    return [
        [values["Jxx"], values["Jxy"], values["Jxz"]],
        [values["Jyx"], values["Jyy"], values["Jyz"]],
        [values["Jzx"], values["Jzy"], values["Jzz"]],
    ], selected_stage


def matrix_lines(title: str, matrix: list[list[float]]) -> list[str]:
    lines = [title]
    for row in matrix:
        lines.append("  " + " ".join(f"{value: .8f}" for value in row))
    return lines


def kitaev_report(args: argparse.Namespace) -> None:
    info = read_poscar(Path(args.poscar))
    mag = magnetic_indices(info, args.magnetic_elements)
    if not args.pair:
        raise ValueError("--pair is required for kitaev-report")
    pairs = parse_pairs(args.pair, args.index_mode, mag, info.natoms)
    if len(pairs) != 1:
        raise ValueError("kitaev-report expects exactly one pair")
    report_bond = pair_bond_context(info, pairs[0], args)
    validate_pair_bond_context("kitaev-report", pairs[0], report_bond, args)
    frame = detect_kitaev_frame(info, mag, pairs[0], args)
    lines = ["Kitaev local-frame report", *frame_rows(frame, info)]
    if report_bond.multiplicity > 1:
        lines.append("periodic_bond_sum_warning\tJani tensor is a POSCAR-pair translation sum; no multiplicity division was applied")
    if args.jani_root:
        global_matrix, selected_stage = parse_jani_summary(Path(args.jani_root), frame.pair.label, args.stage)
        common_matrix, local_gauge_matrix, kitaev_gauge_matrix = kitaev_exchange_matrices(global_matrix, frame)
        symmetric_common = symmetric_matrix(common_matrix)
        physical = decompose_common_exchange(common_matrix)
        lines.extend(matrix_lines("global_J_meV", global_matrix))
        lines.extend(matrix_lines("common_J_meV rows=(alpha,beta,gamma) cols=(alpha,beta,gamma)", common_matrix))
        lines.extend(matrix_lines("common_symmetric_J_meV", symmetric_common))
        lines.extend(
            matrix_lines(
                "local_gauge_J_meV rows_i=(alpha,beta,gamma) cols_j=(alpha,beta,gamma)",
                local_gauge_matrix,
            )
        )
        lines.extend(matrix_lines("kitaev_gauge_J_meV rows_i=(x,y,z) cols_j=(x,y,z)", kitaev_gauge_matrix))
        if selected_stage:
            lines.append(f"jani_stage\t{selected_stage}")
        lines.append("physical_decomposition_frame\tcommon alpha,beta,gamma bond frame from site-i octahedral axes")
        lines.append("local_gauge_note\tlocal_gauge/kitaev_gauge matrices use different left/right site frames and are diagnostics, not physical DMI/K/Gamma decompositions")
        for key, value in physical.items():
            lines.append(f"physical_{key}\t{value:.8f}")
    text = "\n".join(lines) + "\n"
    print(text, end="")
    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text)
        print(f"[OK] Wrote {out.resolve()}")


def scale_vec(vec: Iterable[float], moment: float) -> tuple[float, float, float]:
    return tuple(float(x) * moment for x in vec)  # type: ignore[return-value]


def neg_vec(vec: Iterable[float]) -> tuple[float, float, float]:
    return tuple(-float(x) for x in vec)  # type: ignore[return-value]


def axis_vec(axis: str, moment: float) -> tuple[float, float, float]:
    return scale_vec(AXIS[axis], moment)


def base_vectors(natoms: int, mag: list[int], moment: float, background_axis: str) -> list[tuple[float, float, float]]:
    zero = (0.0, 0.0, 0.0)
    background = axis_vec(background_axis, moment)
    vecs = [zero for _ in range(natoms)]
    for idx in mag:
        vecs[idx] = background
    return vecs


def flatten_magmom(vecs: list[tuple[float, float, float]]) -> str:
    return " ".join(f"{value:.8f}" for vec in vecs for value in vec)


def jani_state_units(component: str) -> list[tuple[tuple[float, float, float], tuple[float, float, float]]]:
    if not re.fullmatch(r"J[xyz][xyz]", component):
        raise ValueError(f"Bad Jani component: {component}")
    a = AXIS[component[1]]
    b = AXIS[component[2]]
    return [(a, b), (a, neg_vec(b)), (neg_vec(a), b), (neg_vec(a), neg_vec(b))]


def sia_state_units(component: str) -> tuple[str, list[tuple[float, float, float]]]:
    inv = 1.0 / math.sqrt(2.0)
    if component == "Axy":
        return "z", [(inv, inv, 0.0), (inv, -inv, 0.0), (-inv, inv, 0.0), (-inv, -inv, 0.0)]
    if component == "Axz":
        return "y", [(inv, 0.0, inv), (inv, 0.0, -inv), (-inv, 0.0, inv), (-inv, 0.0, -inv)]
    if component == "Ayz":
        return "x", [(0.0, inv, inv), (0.0, inv, -inv), (0.0, -inv, inv), (0.0, -inv, -inv)]
    if component == "Ayy_minus_Axx":
        return "z", [AXIS["y"], AXIS["x"], neg_vec(AXIS["x"]), neg_vec(AXIS["y"])]
    if component == "Azz_minus_Axx":
        return "y", [AXIS["z"], AXIS["x"], neg_vec(AXIS["x"]), neg_vec(AXIS["z"])]
    raise ValueError(f"Bad SIA component: {component}")


def biqua_state_units() -> list[tuple[tuple[float, float, float], tuple[float, float, float]]]:
    inv = 1.0 / math.sqrt(2.0)
    return [
        (AXIS["x"], AXIS["x"]),
        (AXIS["x"], neg_vec(AXIS["x"])),
        (AXIS["x"], (inv, inv, 0.0)),
        (AXIS["x"], (-inv, -inv, 0.0)),
    ]


def read_incar(path: Path) -> list[str]:
    if path.exists():
        return path.read_text().splitlines(keepends=True)
    return []


def strip_incar_comment(line: str) -> str:
    return re.split(r"[#!]", line, maxsplit=1)[0].strip()


def read_incar_tags(path: Path) -> dict[str, str]:
    tags: dict[str, str] = {}
    if not path.exists():
        return tags
    for raw in path.read_text(errors="ignore").splitlines():
        line = strip_incar_comment(raw)
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        tags[key.strip().upper()] = value.strip()
    return tags


def parse_numeric_values(raw: str) -> list[float]:
    values: list[float] = []
    for token in re.split(r"[\s,;]+", raw.strip()):
        if not token:
            continue
        try:
            values.append(float(token.replace("D", "E").replace("d", "e")))
        except ValueError:
            continue
    return values


def saxis_is_default(raw_values: list[str] | list[float]) -> bool:
    try:
        values = [float(str(value).replace("D", "E").replace("d", "e")) for value in raw_values]
    except ValueError:
        return False
    return (
        len(values) == 3
        and abs(values[0]) <= 1e-12
        and abs(values[1]) <= 1e-12
        and abs(values[2] - 1.0) <= 1e-12
    )


def validate_cli_saxis(args: argparse.Namespace) -> None:
    if not getattr(args, "saxis", None):
        return
    if not saxis_is_default(args.saxis):
        raise ValueError(
            "Only the default SAXIS = 0 0 1 is currently supported. Non-default SAXIS changes the spinor basis, "
            "while generated M_CONSTR vectors are Cartesian; remove --saxis or set --saxis 0 0 1."
        )


def version_tuple(raw: str | None) -> tuple[int, int, int] | None:
    if not raw:
        return None
    nums = [int(match) for match in re.findall(r"\d+", raw)[:3]]
    if not nums:
        return None
    while len(nums) < 3:
        nums.append(0)
    return tuple(nums[:3])


def validate_constraint_mode(args: argparse.Namespace) -> None:
    mode = int(getattr(args, "constraint_mode", 4))
    if mode not in {1, 2, 4}:
        raise ValueError("--constraint-mode must be one of 1, 2, or 4")
    parsed = version_tuple(getattr(args, "vasp_version", None))
    if mode == 4 and parsed is not None and parsed < (6, 4, 0):
        raise ValueError(
            "--constraint-mode 4 requires VASP >= 6.4.0 because it constrains direction and sign. "
            "Use a newer VASP for production four-state runs, or explicitly choose --constraint-mode 1/2 "
            "and rely on strict signed MW_int diagnostics."
        )
    if mode == 1:
        print(
            "[WARN] I_CONSTRAINED_M=1 constrains only the axis, not the +/- sign; "
            "postprocess will fail if MW_int flips sign.",
            file=sys.stderr,
        )
    if mode == 2:
        print(
            "[WARN] I_CONSTRAINED_M=2 constrains direction and moment size; verify moment-length convergence "
            "before using tiny energy differences.",
            file=sys.stderr,
        )


def incar_bool(raw: str) -> bool | None:
    token = raw.strip().split()[0].strip(".").upper() if raw.strip() else ""
    if token in {"TRUE", "T", "1", "YES"}:
        return True
    if token in {"FALSE", "F", "0", "NO"}:
        return False
    return None


def validate_template_incar(incar: Path, info: PoscarInfo, kind: str, label: str) -> None:
    tags = read_incar_tags(incar)
    errors: list[str] = []
    lambda_values = parse_numeric_values(tags.get("LAMBDA", ""))
    if not lambda_values or max(lambda_values) <= 0:
        errors.append("missing positive LAMBDA; VASP defaults LAMBDA to 0, so constrained moments would not be penalized")
    rwigs_values = parse_numeric_values(tags.get("RWIGS", ""))
    if len(rwigs_values) < len(info.elements) or any(value <= 0 for value in rwigs_values[: len(info.elements)]):
        errors.append(f"missing positive RWIGS values for all {len(info.elements)} POSCAR element(s)")
    if tags.get("ISPIN", "").strip().split()[:1] == ["2"]:
        errors.append("ISPIN=2 is incompatible with the generated noncollinear setup; remove ISPIN from the template")
    if "SAXIS" in tags and not saxis_is_default(parse_numeric_values(tags["SAXIS"])):
        errors.append(
            "non-default SAXIS is not supported; generated MAGMOM/M_CONSTR vectors assume the default SAXIS = 0 0 1"
        )
    if kind in {"jani", "sia", "kitaev"} and incar_bool(tags.get("LSORBIT", "")) is not True:
        errors.append(f"{kind} requires LSORBIT=.TRUE. for SOC-driven anisotropic terms")
    if kind == "biqua" and incar_bool(tags.get("LSORBIT", "")) is True:
        errors.append("biquadratic extraction should be run without LSORBIT to avoid SOC anisotropy contamination")
    if errors:
        formatted = "\n  - ".join(errors)
        raise ValueError(f"Unsafe INCAR template for {label}: {incar}\n  - {formatted}")


def validate_input_templates(info: PoscarInfo, args: argparse.Namespace) -> None:
    templates: list[tuple[str, Path]] = []
    if args.workflow == "pbe-hse":
        templates.append(("pbe_pre", Path(args.pbe_input_dir or args.input_dir).resolve() / "INCAR"))
        templates.append(("hse_no_u", Path(args.hse_input_dir or args.input_dir).resolve() / "INCAR"))
    else:
        templates.append(("single", Path(args.input_dir).resolve() / "INCAR"))
    for label, incar in templates:
        validate_template_incar(incar, info, args.kind, label)


def set_incar_tags(path: Path, tags: dict[str, str | None]) -> None:
    lines = read_incar(path)
    seen: set[str] = set()
    out: list[str] = []
    patterns = {tag.upper(): re.compile(rf"^\s*{re.escape(tag)}\s*=", re.I) for tag in tags}
    for line in lines:
        matched = None
        for tag, pattern in patterns.items():
            if pattern.match(line):
                matched = tag
                break
        if matched is None:
            out.append(line)
            continue
        if matched in seen:
            continue
        if tags[matched] is not None:
            out.append(f"{matched} = {tags[matched]}\n")
        seen.add(matched)
    for tag, value in tags.items():
        tag = tag.upper()
        if tag not in seen and value is not None:
            out.append(f"{tag} = {value}\n")
    path.write_text("".join(out))


def copy_required_inputs(src: Path, dst: Path, poscar: Path | None = None) -> None:
    dst.mkdir(parents=True, exist_ok=True)
    for name in REQUIRED_INPUTS:
        source = poscar if name == "POSCAR" and poscar is not None else src / name
        if not source.exists():
            raise FileNotFoundError(f"Required input missing: {source}")
        shutil.copy2(source, dst / name)


def common_incar_tags(magmom: str, args: argparse.Namespace) -> dict[str, str | None]:
    tags: dict[str, str | None] = {
        "MAGMOM": magmom,
        "M_CONSTR": magmom,
        "LNONCOLLINEAR": ".TRUE.",
        "I_CONSTRAINED_M": str(int(getattr(args, "constraint_mode", 4))),
        "ISPIN": None,
        "NSW": "0",
        "IBRION": "-1",
        "LASPH": ".TRUE.",
        "GGA_COMPAT": ".FALSE.",
    }
    if args.saxis:
        tags["SAXIS"] = " ".join(args.saxis)
    return tags


def stage_tags(stage: str, args: argparse.Namespace) -> dict[str, str | None]:
    if args.no_stage_tag_edits:
        return {}
    if stage == "pbe":
        return {
            "ISTART": "0",
            "ICHARG": "2",
            "LWAVE": ".TRUE.",
            "LCHARG": ".TRUE.",
        }
    if stage == "hse":
        return {
            "ISTART": "1",
            "ICHARG": "0",
            "LWAVE": ".TRUE.",
            "LCHARG": ".TRUE.",
            "LHFCALC": ".TRUE.",
            "GGA": "PE",
            "HFSCREEN": "0.2",
            "AEXX": "0.25",
            "ALGO": "All",
            "LDAU": ".FALSE.",
        }
    return {}


def write_state(
    out_root: Path,
    relpath: str,
    magmom: str,
    system_label: str,
    args: argparse.Namespace,
) -> None:
    poscar_override = Path(args.poscar).resolve()
    if args.workflow == "pbe-hse":
        pbe_src = Path(args.pbe_input_dir or args.input_dir).resolve()
        hse_src = Path(args.hse_input_dir or args.input_dir).resolve()
        pbe_dst = out_root / "pbe" / relpath
        hse_dst = out_root / relpath
        copy_required_inputs(pbe_src, pbe_dst, poscar_override)
        copy_required_inputs(hse_src, hse_dst, poscar_override)
        pbe_tags = common_incar_tags(magmom, args) | stage_tags("pbe", args)
        hse_tags = common_incar_tags(magmom, args) | stage_tags("hse", args)
        pbe_tags["SYSTEM"] = f"{system_label} PBE+U pre"
        hse_tags["SYSTEM"] = f"{system_label} HSE06 no-U"
        set_incar_tags(pbe_dst / "INCAR", pbe_tags)
        set_incar_tags(hse_dst / "INCAR", hse_tags)
        return
    src = Path(args.input_dir).resolve()
    dst = out_root / relpath
    copy_required_inputs(src, dst, poscar_override)
    tags = common_incar_tags(magmom, args)
    tags["SYSTEM"] = system_label
    set_incar_tags(dst / "INCAR", tags)


def add_job(jobs: list[dict[str, str]], relpath: str, name: str, description: str) -> None:
    jobs.append({"relpath": relpath, "job_name": name, "description": description})


def generate_jani(
    info: PoscarInfo,
    mag: list[int],
    pairs: list[PairInfo],
    out_root: Path,
    args: argparse.Namespace,
) -> tuple[list[dict[str, str]], list[dict[str, object]]]:
    jobs: list[dict[str, str]] = []
    formulas: list[dict[str, object]] = []
    nest_pairs = len(pairs) > 1
    background_axis = args.background_axis or "z"
    for pair in pairs:
        bond = pair_bond_context(info, pair, args)
        validate_pair_bond_context("jani", pair, bond, args)
        denominator = energy_denominator(args, 4.0, 2)
        prefix = f"{pair.label}/" if nest_pairs else ""
        for comp in JANI_COMPONENTS:
            states = []
            for state_no, (ui, uj) in enumerate(jani_state_units(comp), start=1):
                vecs = base_vectors(info.natoms, mag, args.moment, background_axis)
                vecs[pair.global_i] = scale_vec(ui, args.moment)
                vecs[pair.global_j] = scale_vec(uj, args.moment)
                rel = f"{prefix}{comp}/{state_no}"
                write_state(out_root, rel, flatten_magmom(vecs), f"{pair.label}_{comp}_{state_no}", args)
                add_job(jobs, rel, f"{pair.label}_{comp}_{state_no}", f"{pair.label} {comp} state {state_no}")
                states.append(rel)
            formulas.append(
                {
                    "label": f"{pair.label}_{comp}" if nest_pairs else comp,
                    "kind": "four_state",
                    "quantity": f"{comp}_meV",
                    "states": states,
                    "state_labels": ["E1", "E2", "E3", "E4"],
                    "formula": f"prefactor * (E1 - E2 - E3 + E4) * 1000 / {denominator:g}",
                    "energy_denominator": denominator,
                    "hamiltonian_prefactor": hamiltonian_prefactor(args),
                    "target_global_indices_1based": [pair.global_i + 1, pair.global_j + 1],
                    **bond_formula_metadata(bond, args),
                }
            )
    return jobs, formulas


def generate_jiso(
    info: PoscarInfo,
    mag: list[int],
    pairs: list[PairInfo],
    out_root: Path,
    args: argparse.Namespace,
) -> tuple[list[dict[str, str]], list[dict[str, object]]]:
    jobs: list[dict[str, str]] = []
    formulas: list[dict[str, object]] = []
    background_axis = args.background_axis or "x"
    spin_units = {
        "upup": (AXIS[args.pair_axis], AXIS[args.pair_axis]),
        "updn": (AXIS[args.pair_axis], neg_vec(AXIS[args.pair_axis])),
        "dnup": (neg_vec(AXIS[args.pair_axis]), AXIS[args.pair_axis]),
        "dndn": (neg_vec(AXIS[args.pair_axis]), neg_vec(AXIS[args.pair_axis])),
    }
    for pair in pairs:
        bond = pair_bond_context(info, pair, args)
        validate_pair_bond_context("jiso", pair, bond, args)
        denominator = energy_denominator(args, 4.0, 2)
        states = []
        for spin in JISO_SPINS:
            ui, uj = spin_units[spin]
            vecs = base_vectors(info.natoms, mag, args.moment, background_axis)
            vecs[pair.global_i] = scale_vec(ui, args.moment)
            vecs[pair.global_j] = scale_vec(uj, args.moment)
            rel = f"{pair.label}/{spin}"
            write_state(out_root, rel, flatten_magmom(vecs), f"{pair.label}_{spin}", args)
            add_job(jobs, rel, f"{pair.label}_{spin}", f"{pair.label} {spin}")
            states.append(rel)
        formulas.append(
            {
                "label": pair.label,
                "kind": "four_state",
                "quantity": f"J{args.pair_axis}{args.pair_axis}_meV",
                "states": states,
                "state_labels": JISO_SPINS,
                    "formula": f"prefactor * (E_upup - E_updn - E_dnup + E_dndn) * 1000 / {denominator:g}",
                    "energy_denominator": denominator,
                    "hamiltonian_prefactor": hamiltonian_prefactor(args),
                    "target_global_indices_1based": [pair.global_i + 1, pair.global_j + 1],
                    **bond_formula_metadata(bond, args),
                    "note": "Axis-resolved four-state exchange along --pair-axis; use Jani trace for a rotational isotropic average.",
                }
            )
    return jobs, formulas


def generate_kitaev(
    info: PoscarInfo,
    mag: list[int],
    pairs: list[PairInfo],
    out_root: Path,
    args: argparse.Namespace,
) -> tuple[list[dict[str, str]], list[dict[str, object]]]:
    jobs: list[dict[str, str]] = []
    formulas: list[dict[str, object]] = []
    frames: list[KitaevFrame] = []
    background_axis = args.background_axis or "z"
    spin_names = ["pp", "pm", "mp", "mm"]
    nest_pairs = len(pairs) > 1
    for pair in pairs:
        bond = pair_bond_context(info, pair, args)
        validate_pair_bond_context("kitaev", pair, bond, args)
        denominator = energy_denominator(args, 4.0, 2)
        frame = detect_kitaev_frame(info, mag, pair, args)
        frames.append(frame)
        spin_units = {
            "pp": (frame.gamma_axis, frame.gamma_axis_j),
            "pm": (frame.gamma_axis, neg_vec(frame.gamma_axis_j)),
            "mp": (neg_vec(frame.gamma_axis), frame.gamma_axis_j),
            "mm": (neg_vec(frame.gamma_axis), neg_vec(frame.gamma_axis_j)),
        }
        prefix = f"{pair.label}/" if nest_pairs else ""
        states: list[str] = []
        for spin in spin_names:
            ui, uj = spin_units[spin]
            vecs = base_vectors(info.natoms, mag, args.moment, background_axis)
            vecs[pair.global_i] = scale_vec(ui, args.moment)
            vecs[pair.global_j] = scale_vec(uj, args.moment)
            rel = f"{prefix}Kitaev_{frame.gamma_label}/{spin}"
            write_state(out_root, rel, flatten_magmom(vecs), f"{pair.label}_Kitaev_{frame.gamma_label}_{spin}", args)
            add_job(jobs, rel, f"{pair.label}_K_{frame.gamma_label}_{spin}", f"{pair.label} Kitaev {frame.gamma_label} {spin}")
            states.append(rel)
        formulas.append(
            {
                "label": f"{pair.label}_Kitaev_{frame.gamma_label}" if nest_pairs else f"Kitaev_{frame.gamma_label}",
                "kind": "four_state",
                "quantity": "local_gauge_J_gamma_i_gamma_j_meV",
                "states": states,
                "state_labels": spin_names,
                "formula": f"prefactor * (E_pp - E_pm - E_mp + E_mm) * 1000 / {denominator:g}",
                "energy_denominator": denominator,
                "hamiltonian_prefactor": hamiltonian_prefactor(args),
                "target_global_indices_1based": [pair.global_i + 1, pair.global_j + 1],
                **bond_formula_metadata(bond, args),
                "note": "Direct Kitaev prepare constrains site-i gamma_i and site-j gamma_j, so this is a local-gauge projection gamma_i^T J gamma_j. Use full Jani plus kitaev-report common-frame decomposition for physical K/Gamma/DMI.",
            }
        )
    write_kitaev_frames(out_root, frames, info)
    return jobs, formulas


def generate_sia(
    info: PoscarInfo,
    mag: list[int],
    out_root: Path,
    args: argparse.Namespace,
) -> tuple[list[dict[str, str]], list[dict[str, object]]]:
    if args.atom is None:
        raise ValueError("--atom is required for --kind sia")
    target_global, target_mag = one_based_to_global(args.atom, args.index_mode, mag, info.natoms)
    label = f"{info.atom_symbols[target_global]}{args.atom}"
    jobs: list[dict[str, str]] = []
    formulas: list[dict[str, object]] = []
    for comp in SIA_COMPONENTS:
        bg_axis, units = sia_state_units(comp)
        states = []
        for state_no, unit in enumerate(units, start=1):
            state_label = f"E{state_no}"
            vecs = base_vectors(info.natoms, mag, args.moment, bg_axis)
            vecs[target_global] = scale_vec(unit, args.moment)
            rel = f"{comp}/{state_label}"
            write_state(out_root, rel, flatten_magmom(vecs), f"SIA_{label}_{comp}_{state_label}", args)
            add_job(jobs, rel, f"SIA_{label}_{comp}_{state_label}", f"SIA {label} {comp} {state_label}")
            states.append(rel)
        formulas.append(
            {
                "label": comp,
                "kind": "four_state",
                "quantity": f"{comp}_meV",
                "states": states,
                "state_labels": ["E1", "E2", "E3", "E4"],
                "formula": (
                    f"prefactor * (E1 - E2 - E3 + E4) * 1000 / "
                    f"{energy_denominator(args, sia_energy_denominator(comp), 2):g}"
                ),
                "energy_denominator": energy_denominator(args, sia_energy_denominator(comp), 2),
                "hamiltonian_prefactor": hamiltonian_prefactor(args),
                "target_global_indices_1based": [target_global + 1],
            }
        )
    (out_root / "sia_target.tsv").write_text(
        "input_index\tindex_mode\tglobal_index\tmagnetic_index\telement\n"
        f"{args.atom}\t{args.index_mode}\t{target_global + 1}\t{target_mag}\t{info.atom_symbols[target_global]}\n"
    )
    return jobs, formulas


def generate_biqua(
    info: PoscarInfo,
    mag: list[int],
    pairs: list[PairInfo],
    out_root: Path,
    args: argparse.Namespace,
) -> tuple[list[dict[str, str]], list[dict[str, object]]]:
    if len(pairs) != 1:
        raise ValueError("biquadratic generation expects exactly one pair")
    pair = pairs[0]
    bond = pair_bond_context(info, pair, args)
    validate_pair_bond_context("biqua", pair, bond, args)
    j_denominator = energy_denominator(args, 2.0, 2)
    b_denominator = energy_denominator(args, 1.0, 4)
    jobs: list[dict[str, str]] = []
    states = []
    background_axis = args.background_axis or "z"
    for state_no, (ui, uj) in enumerate(biqua_state_units(), start=1):
        vecs = base_vectors(info.natoms, mag, args.moment, background_axis)
        vecs[pair.global_i] = scale_vec(ui, args.moment)
        vecs[pair.global_j] = scale_vec(uj, args.moment)
        rel = f"Biquadratic/{state_no}"
        write_state(out_root, rel, flatten_magmom(vecs), f"{pair.label}_Biquadratic_{state_no}", args)
        add_job(jobs, rel, f"{pair.label}_Bq_{state_no}", f"{pair.label} biquadratic state {state_no}")
        states.append(rel)
    (out_root / "state_map.tsv").write_text(
        "state\tS1_i\tS2_j\told_script_label\n"
        "1\t(1,0,0)\t(1,0,0)\tE1\n"
        "2\t(1,0,0)\t(-1,0,0)\tE2\n"
        "3\t(1,0,0)\t(1/sqrt2,1/sqrt2,0)\tE3\n"
        "4\t(1,0,0)\t(-1/sqrt2,-1/sqrt2,0)\tE4\n"
    )
    return jobs, [
        {
            "label": pair.label,
            "kind": "biquadratic",
            "quantity": "Biquadratic_B_meV",
            "states": states,
            "state_labels": ["E1", "E2", "E3", "E4"],
            "formula": (
                f"B = prefactor * (E1 + E2 - E3 - E4) * 1000 / {b_denominator:g}; "
                f"J = prefactor * (E1 - E2) * 1000 / {j_denominator:g}"
            ),
            "j_energy_denominator": j_denominator,
            "b_energy_denominator": b_denominator,
            "hamiltonian_prefactor": hamiltonian_prefactor(args),
            "target_global_indices_1based": [pair.global_i + 1, pair.global_j + 1],
            **bond_formula_metadata(bond, args),
        }
    ]


def write_jobs_tsv(out_root: Path, jobs: list[dict[str, str]]) -> None:
    lines = ["relpath\tjob_name\tdescription\n"]
    for job in jobs:
        lines.append(f"{job['relpath']}\t{job['job_name']}\t{job['description']}\n")
    (out_root / "state_jobs.tsv").write_text("".join(lines))


def write_pair_indexing(out_root: Path, pairs: list[PairInfo]) -> None:
    if not pairs:
        return
    lines = ["pair\tinput_i\tinput_j\tglobal_i\tglobal_j\tmagnetic_i\tmagnetic_j\n"]
    for pair in pairs:
        lines.append(
            f"{pair.label}\t{pair.label_i}\t{pair.label_j}\t"
            f"{pair.global_i + 1}\t{pair.global_j + 1}\t{pair.magnetic_i}\t{pair.magnetic_j}\n"
        )
    (out_root / "pair_indexing.tsv").write_text("".join(lines))


def write_run_scripts(out_root: Path, workflow: str) -> None:
    run_state = r"""#!/usr/bin/env bash
#SBATCH --partition=queue1-1,queue3-1,queue3-2,queue3-3
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --time=10-00:00:00
#SBATCH --output=%x_%j.log
#SBATCH --error=%x_%j.err

set -euo pipefail
REL=${1:?relative state path required}
ROOT=${SLURM_SUBMIT_DIR:-$(pwd)}
SAFE_REL=${REL//\//_}
VASP_EXE=${VASP_EXE:-vasp_ncl}

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "$ROOT/stage_status_${SAFE_REL}.log"; }

converged_outcar() {
  local outcar="$1"
  [ -s "$outcar" ] || return 1
  grep -q 'General timing and accounting' "$outcar" || return 1
  if grep -q 'EDIFF was not reached' "$outcar"; then
    return 1
  fi
  return 0
}

save_outputs() {
  local suffix="$1"
  [ -f output ] && cp -f output "output.$suffix"
  [ -f OUTCAR ] && cp -f OUTCAR "OUTCAR.$suffix"
  [ -f OSZICAR ] && cp -f OSZICAR "OSZICAR.$suffix"
  [ -f CONTCAR ] && cp -f CONTCAR "CONTCAR.$suffix"
  [ -f vasprun.xml ] && cp -f vasprun.xml "vasprun.xml.$suffix"
  [ -f REPORT ] && cp -f REPORT "REPORT.$suffix"
}

clean_runtime_outputs() {
  rm -f OUTCAR OSZICAR vasprun.xml CONTCAR XDATCAR DOSCAR EIGENVAL IBZKPT PCDAT PROCAR REPORT output
}

if [ -f /software/compiler/intel/oneapi/setvars.sh ]; then
  set +u
  source /software/compiler/intel/oneapi/setvars.sh >/dev/null 2>&1 || true
  set -u
fi
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export I_MPI_HYDRA_BOOTSTRAP=${I_MPI_HYDRA_BOOTSTRAP:-slurm}
if [ -e /opt/gridview/slurm/lib/libpmi.so ]; then
  export I_MPI_PMI_LIBRARY=${I_MPI_PMI_LIBRARY:-/opt/gridview/slurm/lib/libpmi.so}
fi

run_one() {
  local dir="$1"
  local suffix="$2"
  cd "$dir"
  if converged_outcar "OUTCAR.$suffix"; then
    log "$REL $suffix already converged, skip"
    return
  fi
  rm -f "output.$suffix" "OUTCAR.$suffix" "OSZICAR.$suffix" "CONTCAR.$suffix" "vasprun.xml.$suffix" "REPORT.$suffix"
  clean_runtime_outputs
  log "Run $suffix in $dir"
  srun "$VASP_EXE" > output 2>&1
  save_outputs "$suffix"
  if ! converged_outcar "OUTCAR.$suffix"; then
    log "ERROR $suffix did not converge normally"
    exit 10
  fi
}

PBE_DIR="$ROOT/pbe/$REL"
FINAL_DIR="$ROOT/$REL"
if [ -d "$PBE_DIR" ]; then
  run_one "$PBE_DIR" pbe_pre
  if [ ! -s "$PBE_DIR/WAVECAR" ] || [ ! -s "$PBE_DIR/CHGCAR" ]; then
    log "ERROR PBE pre finished but WAVECAR/CHGCAR missing"
    exit 11
  fi
  cd "$FINAL_DIR"
  if converged_outcar OUTCAR.hse_no_u; then
    log "$REL hse_no_u already converged, skip"
  else
    rm -f output.hse_no_u OUTCAR.hse_no_u OSZICAR.hse_no_u CONTCAR.hse_no_u vasprun.xml.hse_no_u REPORT.hse_no_u
    clean_runtime_outputs
    cp -f "$PBE_DIR/WAVECAR" WAVECAR
    cp -f "$PBE_DIR/CHGCAR" CHGCAR
    [ -s "$PBE_DIR/CHG" ] && cp -f "$PBE_DIR/CHG" CHG || true
    log "Run hse_no_u in $FINAL_DIR"
    srun "$VASP_EXE" > output 2>&1
    save_outputs hse_no_u
    if ! converged_outcar OUTCAR.hse_no_u; then
      log "ERROR hse_no_u did not converge normally"
      exit 20
    fi
  fi
else
  run_one "$FINAL_DIR" single
fi
log "DONE $REL"
"""
    submit_all = r"""#!/usr/bin/env bash
set -euo pipefail
ROOT=${SLURM_SUBMIT_DIR:-$(pwd)}
cd "$ROOT"
: > submitted_jobs.tsv
tail -n +2 state_jobs.tsv | while IFS=$'\t' read -r relpath job_name description; do
  [ -n "${relpath:-}" ] || continue
  output=$(sbatch --job-name="$job_name" ./run_state.sh "$relpath")
  jobid=$(printf '%s\n' "$output" | awk '{print $NF}')
  printf '%s\t%s\t%s\n' "$relpath" "$job_name" "$jobid" | tee -a submitted_jobs.tsv
done
"""
    postprocess_sh = r"""#!/usr/bin/env bash
set -euo pipefail
ROOT=${SLURM_SUBMIT_DIR:-$(pwd)}
cd "$ROOT"
python3 ./postprocess.py
"""
    (out_root / "run_state.sh").write_text(run_state)
    (out_root / "submit_all.sh").write_text(submit_all)
    (out_root / "postprocess.sh").write_text(postprocess_sh)
    for name in ["run_state.sh", "submit_all.sh", "postprocess.sh"]:
        try:
            (out_root / name).chmod(0o755)
        except OSError:
            pass


POSTPROCESS_PY = r'''#!/usr/bin/env python3
from __future__ import annotations
import json
import math
import os
import re
import sys
from pathlib import Path

FLOAT_RE_TEXT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
ENERGY_RE = re.compile(r"E0=\s*(" + FLOAT_RE_TEXT + r")")
TOTEN_RE = re.compile(r"free\s+energy\s+TOTEN\s+=\s+(" + FLOAT_RE_TEXT + r")")
PENALTY_RE = re.compile(r"\bE_p\s*=\s*(" + FLOAT_RE_TEXT + r")")
LAMBDA_RE = re.compile(r"\blambda\s*=\s*(" + FLOAT_RE_TEXT + r")", re.I)
ION_MW_HEADER_RE = re.compile(r"^\s*ion\s+MW_int\s+M_int\s*$", re.I)
ION_LAMBDA_MW_PERP_HEADER_RE = re.compile(r"^\s*ion\s+lambda\s*\*\s*MW_perp\s*$", re.I)
MAX_PENALTY_EV = float(os.environ.get("FOUR_STATE_MAX_PENALTY_EV", "1e-4"))
MAX_TARGET_ANGLE_DEG = float(os.environ.get("FOUR_STATE_MAX_TARGET_ANGLE_DEG", "5.0"))
STRICT_CONSTRAINTS = os.environ.get("FOUR_STATE_STRICT_CONSTRAINTS", "1").strip().lower() not in {"0", "false", "no"}
NONZERO_MCONSTR_TOL = 1e-8

def parse_float(raw: str) -> float:
    return float(raw.replace("D", "E").replace("d", "e"))

def parse_floats(raw: str) -> list[float]:
    return [parse_float(match) for match in re.findall(FLOAT_RE_TEXT, raw)]

def vector_norm(vec) -> float:
    return math.sqrt(sum(value * value for value in vec))

def vector_dot(a, b) -> float:
    return sum(a[idx] * b[idx] for idx in range(3))

def vector_angle_degrees(a, b) -> float | None:
    na = vector_norm(a)
    nb = vector_norm(b)
    if na <= 1e-12 or nb <= 1e-12:
        return None
    dot = vector_dot(a, b) / (na * nb)
    dot = max(-1.0, min(1.0, dot))
    return math.degrees(math.acos(dot))

def as_vectors(values: list[float]) -> list[tuple[float, float, float]]:
    n = len(values) // 3
    return [tuple(values[3 * idx : 3 * idx + 3]) for idx in range(n)]

def strip_incar_comment(line: str) -> str:
    return re.split(r"[#!]", line, maxsplit=1)[0].strip()

def read_incar_tag(path: Path, tag: str) -> str | None:
    if not path.exists():
        return None
    tag_upper = tag.upper()
    for raw in path.read_text(errors="ignore").splitlines():
        line = strip_incar_comment(raw)
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        if key.strip().upper() == tag_upper:
            return value.strip()
    return None

def read_text_candidates(path: Path, suffix: str) -> list[tuple[Path, str]]:
    candidates = [
        path / f"output.{suffix}",
        path / f"OSZICAR.{suffix}",
        path / f"OUTCAR.{suffix}",
        path / "output",
        path / "OSZICAR",
        path / "OUTCAR",
    ]
    texts = []
    for cand in candidates:
        if cand.exists() and cand.stat().st_size > 0:
            texts.append((cand, cand.read_text(errors="ignore")))
    return texts

def last_regex_float(texts: list[tuple[Path, str]], pattern: re.Pattern[str]) -> float | None:
    for _, text in reversed(texts):
        matches = pattern.findall(text)
        if matches:
            return parse_float(matches[-1])
    return None

def parse_ion_moment_tables(text: str) -> list[tuple[dict[int, tuple[float, float, float]], dict[int, tuple[float, float, float]]]]:
    tables = []
    lines = text.splitlines()
    idx = 0
    while idx < len(lines):
        if not ION_MW_HEADER_RE.match(lines[idx]):
            idx += 1
            continue
        idx += 1
        mw_table: dict[int, tuple[float, float, float]] = {}
        m_table: dict[int, tuple[float, float, float]] = {}
        while idx < len(lines):
            match = re.match(r"^\s*(\d+)\s+(.+?)\s*$", lines[idx])
            if not match:
                break
            vals = parse_floats(match.group(2))
            if len(vals) < 6:
                break
            ion = int(match.group(1))
            mw_table[ion] = tuple(vals[0:3])
            m_table[ion] = tuple(vals[3:6])
            idx += 1
        if mw_table:
            tables.append((mw_table, m_table))
    return tables

def parse_lambda_mw_perp_tables(text: str) -> list[dict[int, tuple[float, float, float]]]:
    tables = []
    lines = text.splitlines()
    idx = 0
    while idx < len(lines):
        if not ION_LAMBDA_MW_PERP_HEADER_RE.match(lines[idx]):
            idx += 1
            continue
        idx += 1
        table: dict[int, tuple[float, float, float]] = {}
        while idx < len(lines):
            match = re.match(r"^\s*(\d+)\s+(.+?)\s*$", lines[idx])
            if not match:
                break
            vals = parse_floats(match.group(2))
            if len(vals) < 3:
                break
            ion = int(match.group(1))
            table[ion] = tuple(vals[0:3])
            idx += 1
        if table:
            tables.append(table)
    return tables

def last_ion_moment_table(texts: list[tuple[Path, str]]) -> tuple[dict[int, tuple[float, float, float]], dict[int, tuple[float, float, float]]]:
    for _, text in reversed(texts):
        tables = parse_ion_moment_tables(text)
        if tables:
            return tables[-1]
    return {}, {}

def last_lambda_mw_perp_table(texts: list[tuple[Path, str]]) -> dict[int, tuple[float, float, float]]:
    for _, text in reversed(texts):
        tables = parse_lambda_mw_perp_tables(text)
        if tables:
            return tables[-1]
    return {}

def constrained_indices_from_mconstr(vectors: list[tuple[float, float, float]]) -> list[int]:
    return [idx + 1 for idx, vec in enumerate(vectors) if vector_norm(vec) > NONZERO_MCONSTR_TOL]

def vectors_for_indices(table: dict[int, tuple[float, float, float]], indices_1based: list[int]) -> list[tuple[float, float, float]]:
    return [table[idx] for idx in indices_1based if idx in table]

def format_vector_list(vectors: list[tuple[float, float, float]]) -> str:
    if not vectors:
        return "NA"
    return ";".join("(" + ",".join(f"{value:.8f}" for value in vec) + ")" for vec in vectors)

def format_index_list(indices: list[int]) -> str:
    return ",".join(str(idx) for idx in indices) if indices else "NA"

def max_norm(vectors: list[tuple[float, float, float]]) -> float | None:
    if not vectors:
        return None
    return max(vector_norm(vec) for vec in vectors)

def extract_energy(path: Path, suffix: str) -> float:
    for cand, text in read_text_candidates(path, suffix):
        matches = ENERGY_RE.findall(text)
        if matches:
            return parse_float(matches[-1])
        matches = TOTEN_RE.findall(text)
        if matches:
            return parse_float(matches[-1])
    raise FileNotFoundError(f"No energy found in {path} for suffix {suffix}")

def extract_state_record(base: Path, rel: str, suffix: str, formula: dict, state_label: str, stage_name: str, metadata: dict) -> dict:
    path = base / rel
    texts = read_text_candidates(path, suffix)
    energy = extract_energy(path, suffix)
    penalty_ep = last_regex_float(texts, PENALTY_RE)
    reported_lambda = last_regex_float(texts, LAMBDA_RE)
    incar_lambda = read_incar_tag(path / "INCAR", "LAMBDA")
    incar_lambda_values = parse_floats(incar_lambda or "")
    lambda_value = reported_lambda if reported_lambda is not None else (max(incar_lambda_values) if incar_lambda_values else None)
    mw_table, m_table = last_ion_moment_table(texts)
    lambda_mw_perp_table = last_lambda_mw_perp_table(texts)
    m_constr_raw = read_incar_tag(path / "INCAR", "M_CONSTR")
    target_vectors = as_vectors(parse_floats(m_constr_raw or ""))
    constrained_indices = constrained_indices_from_mconstr(target_vectors)
    formula_target_indices = [int(idx) for idx in formula.get("target_global_indices_1based", [])]
    target_indices = formula_target_indices or constrained_indices
    constrained_refs = [target_vectors[idx - 1] for idx in constrained_indices if idx <= len(target_vectors)]
    constrained_mw = vectors_for_indices(mw_table, constrained_indices)
    target_mw = vectors_for_indices(mw_table, target_indices)
    constrained_m = vectors_for_indices(m_table, constrained_indices)
    target_m = vectors_for_indices(m_table, target_indices)
    constrained_lambda_mw_perp = vectors_for_indices(lambda_mw_perp_table, constrained_indices)
    target_lambda_mw_perp = vectors_for_indices(lambda_mw_perp_table, target_indices)
    constrained_angles = []
    target_angles = []
    zero_mw_indices = []
    sign_flip_indices = []
    for idx, ref in zip(constrained_indices, constrained_refs):
        measured = mw_table.get(idx)
        if measured is None:
            continue
        angle = vector_angle_degrees(ref, measured)
        if angle is None:
            zero_mw_indices.append(idx)
            continue
        constrained_angles.append(angle)
        if idx in target_indices:
            target_angles.append(angle)
        if vector_dot(ref, measured) < 0.0:
            sign_flip_indices.append(idx)
    max_constrained_angle = max(constrained_angles) if constrained_angles else None
    max_target_angle = max(target_angles) if target_angles else None
    missing_mw_indices = [idx for idx in constrained_indices if idx not in mw_table]
    missing_target_vector_indices = [idx for idx in target_indices if idx > len(target_vectors)]
    unconstrained_target_indices = [
        idx
        for idx in target_indices
        if idx <= len(target_vectors) and vector_norm(target_vectors[idx - 1]) <= NONZERO_MCONSTR_TOL
    ]
    missing_target_mw_indices = [idx for idx in target_indices if idx not in mw_table]
    messages = []
    passed = True

    def fail(message: str) -> None:
        nonlocal passed
        passed = False
        messages.append(message)

    def note(message: str) -> None:
        messages.append(message)

    if penalty_ep is None:
        fail("E_p_missing")
    elif abs(penalty_ep) > MAX_PENALTY_EV:
        fail(f"E_p>{MAX_PENALTY_EV:g}")
    if reported_lambda is None:
        fail("lambda_missing")
    if m_constr_raw is None:
        fail("M_CONSTR_missing")
    elif not target_vectors:
        fail("M_CONSTR_empty_or_unparsed")
    elif not constrained_indices:
        fail("M_CONSTR_has_no_nonzero_vectors")
    if missing_target_vector_indices:
        fail("target_M_CONSTR_missing_ions=" + format_index_list(missing_target_vector_indices))
    if unconstrained_target_indices:
        fail("target_M_CONSTR_zero_ions=" + format_index_list(unconstrained_target_indices))
    if not mw_table:
        fail("MW_int_missing")
    elif missing_mw_indices:
        fail("MW_int_missing_ions=" + format_index_list(missing_mw_indices))
    elif missing_target_mw_indices:
        fail("target_MW_int_missing_ions=" + format_index_list(missing_target_mw_indices))
    if zero_mw_indices:
        fail("MW_int_zero_ions=" + format_index_list(zero_mw_indices))
    if sign_flip_indices:
        fail("spin_sign_flip_ions=" + format_index_list(sign_flip_indices))
    if max_constrained_angle is None:
        if mw_table and constrained_indices and not missing_mw_indices:
            fail("constraint_angle_unavailable")
    elif max_constrained_angle > MAX_TARGET_ANGLE_DEG:
        fail(f"angle>{MAX_TARGET_ANGLE_DEG:g}")
    if not m_table:
        note("M_int_missing")
    if not lambda_mw_perp_table:
        note("lambda_MW_perp_missing")
    if int(metadata.get("constraint_mode", 0) or 0) == 1:
        note("constraint_mode_1_axis_only")
    return {
        "stage": stage_name,
        "formula_label": formula["label"],
        "state_label": state_label,
        "relpath": rel,
        "energy_eV": energy,
        "constraint_mode": metadata.get("constraint_mode", "NA"),
        "penalty_Ep_eV": penalty_ep,
        "max_constrained_angle_deg": max_constrained_angle,
        "max_target_angle_deg": max_target_angle,
        "constrained_global_indices_1based": format_index_list(constrained_indices),
        "target_global_indices_1based": format_index_list(target_indices),
        "constrained_site_MW_int": format_vector_list(constrained_mw),
        "target_site_MW_int": format_vector_list(target_mw),
        "constrained_site_M_int": format_vector_list(constrained_m),
        "target_site_M_int": format_vector_list(target_m),
        "lambda": lambda_value,
        "max_lambda_MW_perp_norm": max_norm(constrained_lambda_mw_perp),
        "lambda_MW_perp": max_norm(constrained_lambda_mw_perp),
        "constrained_site_lambda_MW_perp": format_vector_list(constrained_lambda_mw_perp),
        "target_site_lambda_MW_perp": format_vector_list(target_lambda_mw_perp),
        "constraint_pass": passed,
        "messages": ";".join(messages) if messages else "ok",
    }

def calc_four_state(es, denominator, prefactor):
    return prefactor * (es[0] - es[1] - es[2] + es[3]) * 1000.0 / denominator

def collect_stage(root: Path, metadata: dict, stage: dict) -> tuple[str, list[dict], list[dict]]:
    base = root / stage["base"] if stage["base"] else root
    suffix = stage["suffix"]
    lines = [f"# stage {stage['name']}", "# label quantity E1 E2 E3 E4 value_meV"]
    diagnostics = []
    failures = []
    for formula in metadata["formulas"]:
        state_labels = formula.get("state_labels", [f"E{idx + 1}" for idx in range(len(formula["states"]))])
        state_records = [
            extract_state_record(base, rel, suffix, formula, label, stage["name"], metadata)
            for rel, label in zip(formula["states"], state_labels)
        ]
        diagnostics.extend(state_records)
        failures.extend([record for record in state_records if not record["constraint_pass"]])
        energies = [record["energy_eV"] for record in state_records]
        prefactor = float(formula.get("hamiltonian_prefactor", 1.0))
        if formula["kind"] == "biquadratic":
            j_denominator = float(formula.get("j_energy_denominator", 2.0))
            b_denominator = float(formula.get("b_energy_denominator", 1.0))
            j = prefactor * (energies[0] - energies[1]) * 1000.0 / j_denominator
            b = prefactor * (energies[0] + energies[1] - energies[2] - energies[3]) * 1000.0 / b_denominator
            lines.append(
                f"{formula['label']} J_meV {energies[0]:.12f} {energies[1]:.12f} "
                f"{energies[2]:.12f} {energies[3]:.12f} {j:.8f}"
            )
            lines.append(
                f"{formula['label']} Biquadratic_B_meV {energies[0]:.12f} {energies[1]:.12f} "
                f"{energies[2]:.12f} {energies[3]:.12f} {b:.8f}"
            )
        else:
            denominator = float(formula.get("energy_denominator", 4.0))
            value = calc_four_state(energies, denominator, prefactor)
            lines.append(
                f"{formula['label']} {formula['quantity']} {energies[0]:.12f} {energies[1]:.12f} "
                f"{energies[2]:.12f} {energies[3]:.12f} {value:.8f}"
            )
    return "\n".join(lines) + "\n", diagnostics, failures

def tsv_value(value) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.12g}"
    return str(value).replace("\t", " ").replace("\n", " ")

def write_constraint_diagnostics(path: Path, rows: list[dict]) -> None:
    headers = [
        "stage",
        "formula_label",
        "state_label",
        "relpath",
        "energy_eV",
        "constraint_mode",
        "penalty_Ep_eV",
        "max_constrained_angle_deg",
        "max_target_angle_deg",
        "constrained_global_indices_1based",
        "target_global_indices_1based",
        "constrained_site_MW_int",
        "target_site_MW_int",
        "constrained_site_M_int",
        "target_site_M_int",
        "lambda",
        "max_lambda_MW_perp_norm",
        "lambda_MW_perp",
        "constrained_site_lambda_MW_perp",
        "target_site_lambda_MW_perp",
        "constraint_pass",
        "messages",
    ]
    lines = ["\t".join(headers)]
    for row in rows:
        lines.append("\t".join(tsv_value(row.get(header)) for header in headers))
    path.write_text("\n".join(lines) + "\n")

def main():
    root = Path.cwd()
    metadata = json.loads((root / "metadata.json").read_text())
    results = root / "results"
    results.mkdir(exist_ok=True)
    summary = []
    all_diagnostics = []
    all_failures = []
    for stage in metadata["stages"]:
        text, diagnostics, failures = collect_stage(root, metadata, stage)
        all_diagnostics.extend(diagnostics)
        all_failures.extend(failures)
        out = results / f"{stage['name']}_{metadata['kind']}_energy.dat"
        out.write_text(text)
        summary.append(f"[{stage['name']}]\n{text}")
    write_constraint_diagnostics(results / "constraint_diagnostics.tsv", all_diagnostics)
    if STRICT_CONSTRAINTS and all_failures:
        details = ", ".join(
            f"{row['stage']}:{row['relpath']}({row['messages']})" for row in all_failures[:8]
        )
        raise RuntimeError(
            "Constraint diagnostics failed; refusing final interaction output. "
            f"See {results / 'constraint_diagnostics.tsv'}. First failures: {details}"
        )
    (root / "final_summary.txt").write_text("\n".join(summary))
    print(f"Wrote {root / 'final_summary.txt'}")
    print(f"Wrote {results / 'constraint_diagnostics.tsv'}")

if __name__ == "__main__":
    main()
'''


def write_postprocess_py(out_root: Path) -> None:
    path = out_root / "postprocess.py"
    path.write_text(POSTPROCESS_PY)
    try:
        path.chmod(0o755)
    except OSError:
        pass


def effective_background_axis(kind: str, args: argparse.Namespace) -> str:
    if kind == "sia":
        return "component-specific"
    if args.background_axis:
        return args.background_axis
    if kind == "jiso":
        return "x"
    return "z"


def write_metadata(
    out_root: Path,
    kind: str,
    workflow: str,
    info: PoscarInfo,
    mag: list[int],
    args: argparse.Namespace,
    formulas: list[dict[str, object]],
) -> None:
    stages = STAGE_PBE_HSE if workflow == "pbe-hse" else STAGE_SINGLE
    metadata = {
        "kind": kind,
        "workflow": workflow,
        "poscar": str(Path(args.poscar).resolve()),
        "natoms": info.natoms,
        "elements": info.elements,
        "counts": info.counts,
        "magnetic_elements": args.magnetic_elements,
        "ligand_elements": getattr(args, "ligand_elements", None),
        "magnetic_indices_global_1based": [idx + 1 for idx in mag],
        "moment": args.moment,
        "hamiltonian_sign": getattr(args, "hamiltonian_sign", "plus"),
        "spin_convention": getattr(args, "spin_convention", "unit_vector"),
        "spin_length_S": getattr(args, "spin_length_S", 1.0),
        "bond_distance_tol_A": getattr(args, "bond_distance_tol", 0.02),
        "requested_pair_image_shift": getattr(args, "pair_image_shift", None),
        "constraint_mode": int(getattr(args, "constraint_mode", 4)),
        "vasp_version": getattr(args, "vasp_version", None),
        "max_kitaev_axis_consistency_deg": getattr(args, "max_kitaev_axis_consistency_deg", None),
        "background_axis_effective": effective_background_axis(kind, args),
        "pair_axis": args.pair_axis,
        "index_mode": args.index_mode,
        "stages": stages,
        "formulas": formulas,
    }
    (out_root / "metadata.json").write_text(json.dumps(metadata, indent=2))


def write_readme(out_root: Path, kind: str, workflow: str, args: argparse.Namespace) -> None:
    text = f"""Four-state VASP calculation
kind = {kind}
workflow = {workflow}
moment = {args.moment}
constraint_mode = {int(getattr(args, "constraint_mode", 4))}
vasp_version = {getattr(args, "vasp_version", None)}
background_axis = {effective_background_axis(kind, args)}
pair_axis = {args.pair_axis}
index_mode = {args.index_mode}

Files:
- state_jobs.tsv: relative state directories and suggested Slurm job names.
- run_state.sh: runs one state. Use: sbatch run_state.sh <relpath>
- submit_all.sh: submits every row in state_jobs.tsv.
- postprocess.sh/postprocess.py: collect energies and constraint diagnostics after jobs finish.
- metadata.json: machine-readable formulas and directory map.

Energy formulas:
- Each formula in metadata.json records the actual denominator, Hamiltonian prefactor, spin convention, and pair metadata.
- SIA diagonal differences use denominator 2; off-diagonal SIA and two-site four-state terms use denominator 4 before spin scaling.
- Two-site calculations require a unique explicit POSCAR pair by default. If --allow-periodic-bond-sum was used, the reported value is a summed coupling over periodic translations and is not divided by multiplicity.
- Jiso output is the selected axis component Jaa from --pair-axis, not a tensor trace average.

Indexing:
- input pair/atom labels use the selected index_mode.
- pair_indexing.tsv records global POSCAR indices and magnetic-ion ordinals.
"""
    (out_root / "README.txt").write_text(text)


def prepare(args: argparse.Namespace) -> None:
    info = read_poscar(Path(args.poscar).resolve())
    mag = magnetic_indices(info, args.magnetic_elements)
    if not mag:
        raise ValueError("No magnetic atoms selected")
    out_root = Path(args.out).resolve()
    if out_root.exists() and any(out_root.iterdir()) and not args.force:
        raise FileExistsError(f"Output directory is not empty: {out_root}. Use --force to write into it.")
    validate_cli_saxis(args)
    validate_constraint_mode(args)
    validate_input_templates(info, args)
    out_root.mkdir(parents=True, exist_ok=True)

    pairs: list[PairInfo] = []
    jobs: list[dict[str, str]]
    formulas: list[dict[str, object]]
    if args.kind in {"jani", "jiso", "biqua", "kitaev"}:
        if not args.pair:
            raise ValueError(f"--pair is required for --kind {args.kind}")
        pairs = parse_pairs(args.pair, args.index_mode, mag, info.natoms)
    if args.kind == "jani":
        jobs, formulas = generate_jani(info, mag, pairs, out_root, args)
    elif args.kind == "jiso":
        jobs, formulas = generate_jiso(info, mag, pairs, out_root, args)
    elif args.kind == "kitaev":
        jobs, formulas = generate_kitaev(info, mag, pairs, out_root, args)
    elif args.kind == "sia":
        jobs, formulas = generate_sia(info, mag, out_root, args)
    elif args.kind == "biqua":
        jobs, formulas = generate_biqua(info, mag, pairs, out_root, args)
    else:
        raise ValueError(f"Unknown kind: {args.kind}")

    write_jobs_tsv(out_root, jobs)
    write_pair_indexing(out_root, pairs)
    write_run_scripts(out_root, args.workflow)
    write_postprocess_py(out_root)
    write_metadata(out_root, args.kind, args.workflow, info, mag, args, formulas)
    write_readme(out_root, args.kind, args.workflow, args)
    print(f"[OK] Prepared {len(jobs)} state directories under {out_root}")
    print(f"[OK] Next: cd {out_root} && bash submit_all.sh")


def parse_vaspkit_sequence(raw: str) -> str:
    parts = [part.strip() for part in re.split(r"[,;\s]+", raw) if part.strip()]
    return "\n".join(parts) + "\n"


def extract_rwig_values(potcar: Path) -> list[str]:
    text = potcar.read_text(errors="ignore")
    values = re.findall(r"RWIGS\s*=\s*([-+]?\d+(?:\.\d+)?)", text)
    return values


def update_bootstrap_incar(out_dir: Path, info: PoscarInfo, use_ldau: bool) -> None:
    incar = out_dir / "INCAR"
    tags: dict[str, str] = {}
    potcar = out_dir / "POTCAR"
    if potcar.exists():
        rwigs = extract_rwig_values(potcar)
        if rwigs:
            tags["RWIGS"] = " ".join(rwigs[: len(info.elements)])
    if use_ldau:
        ldaul = []
        ldauu = []
        for elem in info.elements:
            if elem in DEFAULT_U:
                ldaul.append("2")
                ldauu.append(str(DEFAULT_U[elem]))
            else:
                ldaul.append("-1")
                ldauu.append("0")
        tags["LDAU"] = ".TRUE."
        tags["LDAUTYPE"] = "2"
        tags["LDAUL"] = " ".join(ldaul)
        tags["LDAUU"] = " ".join(ldauu)
    if tags:
        set_incar_tags(incar, tags)


def bootstrap(args: argparse.Namespace) -> None:
    poscar = Path(args.poscar).resolve()
    info = read_poscar(poscar)
    out_dir = Path(args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(poscar, out_dir / "POSCAR")
    sequence = parse_vaspkit_sequence(args.vaspkit_sequence)
    print(f"[INFO] Running {args.vaspkit_command} in {out_dir}")
    subprocess.run([args.vaspkit_command], input=sequence, cwd=out_dir, text=True, check=True)
    update_bootstrap_incar(out_dir, info, bool(args.ldau and not args.no_ldau))
    missing = [name for name in REQUIRED_INPUTS if not (out_dir / name).exists()]
    if missing:
        print(f"[WARN] vaspkit finished but missing: {', '.join(missing)}", file=sys.stderr)
    else:
        print(f"[OK] Bootstrapped inputs in {out_dir}")


def collect(args: argparse.Namespace) -> None:
    root = Path(args.root).resolve()
    namespace: dict[str, object] = {"__name__": "__postprocess__"}
    exec(POSTPROCESS_PY, namespace)
    cwd = Path.cwd()
    os.chdir(root)
    try:
        namespace["main"]()
    finally:
        os.chdir(cwd)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Prepare VASP four-state magnetic calculations and analyze magnetic-neighbor geometry.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """\
            Examples:
              python3 four_state_vasp.py neighbors --poscar POSCAR --magnetic-elements Cr --cutoff 10
              python3 four_state_vasp.py prepare --kind jani --poscar POSCAR --input-dir inputs --out pair_14_15 --magnetic-elements Cr --pair 14-15 --moment 6 --workflow pbe-hse
              python3 four_state_vasp.py prepare --kind sia --poscar POSCAR --input-dir inputs --out sia_Cr14 --magnetic-elements Cr --atom 14 --moment 6
              python3 four_state_vasp.py prepare --kind kitaev --poscar POSCAR --input-dir inputs --out kitaev_14_15 --magnetic-elements Cr --ligand-elements I --pair 14-15
              python3 four_state_vasp.py bootstrap --poscar POSCAR --out inputs --vaspkit-sequence "1 102 2 0.04"
            """
        ),
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    neigh = sub.add_parser("neighbors", help="Analyze magnetic-pair neighbor shells and representative pairs")
    neigh.add_argument("--poscar", required=True, help="POSCAR path")
    neigh.add_argument("--magnetic-elements", nargs="+", help="Magnetic element symbols, e.g. Cr")
    neigh.add_argument("--cutoff", type=float, default=10.0, help="Magnetic-pair cutoff in Angstrom")
    neigh.add_argument("--shell-tol", type=float, default=0.08, help="Distance tolerance for shell grouping in Angstrom")
    neigh.add_argument("--center-atom", type=int, help="Preferred center atom for representative shell checks")
    neigh.add_argument("--index-mode", choices=["global", "magnetic"], default="global")
    neigh.add_argument("--boundary-margin", type=float, default=0.5, help="Extra Angstrom margin for expansion heuristic")
    neigh.add_argument("--max-rows", type=int, default=300, help="Maximum contact rows to print to stdout")
    neigh.add_argument("--out", help="Optional directory for neighbor TSV outputs")
    neigh.set_defaults(func=neighbors)

    prep = sub.add_parser("prepare", help="Generate four-state input directories")
    prep.add_argument("--kind", choices=["jani", "jiso", "sia", "biqua", "kitaev"], required=True)
    prep.add_argument("--poscar", required=True, help="POSCAR path")
    prep.add_argument("--input-dir", required=True, help="Directory containing INCAR/KPOINTS/POTCAR/POSCAR templates")
    prep.add_argument("--pbe-input-dir", help="PBE+U template directory for --workflow pbe-hse")
    prep.add_argument("--hse-input-dir", help="HSE no-U template directory for --workflow pbe-hse")
    prep.add_argument("--out", required=True, help="Output calculation root")
    prep.add_argument("--magnetic-elements", nargs="+", help="Magnetic element symbols, e.g. Cr")
    prep.add_argument("--ligand-elements", nargs="+", help="Ligand element symbols for Kitaev octahedral-axis detection, e.g. I")
    prep.add_argument("--pair", action="append", help="Pair I-J. Repeat or comma-separate for multiple pairs.")
    prep.add_argument("--atom", type=int, help="Target atom for SIA")
    prep.add_argument("--index-mode", choices=["global", "magnetic"], default="global")
    prep.add_argument("--moment", type=float, default=6.0)
    prep.add_argument(
        "--constraint-mode",
        type=int,
        choices=[1, 2, 4],
        default=4,
        help="VASP I_CONSTRAINED_M mode. Default 4 constrains direction and sign and requires VASP >= 6.4.0.",
    )
    prep.add_argument("--vasp-version", help="VASP version used for these jobs, e.g. 6.4.3; stored in metadata and checked for mode 4")
    prep.add_argument(
        "--hamiltonian-sign",
        choices=["plus", "minus"],
        default="plus",
        help="Model convention: plus means H = + coefficient term; minus flips extracted coefficient signs.",
    )
    prep.add_argument(
        "--spin-convention",
        choices=["unit_vector", "spin_S"],
        default="unit_vector",
        help="unit_vector reports coefficients for normalized spin directions; spin_S divides by --spin-length-S powers.",
    )
    prep.add_argument("--spin-length-S", type=float, default=1.0, help="Spin length used when --spin-convention spin_S")
    prep.add_argument(
        "--pair-image-shift",
        nargs=3,
        type=int,
        help="Use this image shift only for geometry/reporting; the four-state energy still sums all periodic images of the POSCAR atom pair",
    )
    prep.add_argument(
        "--bond-distance-tol",
        type=float,
        default=0.02,
        help="Distance tolerance in Angstrom for detecting non-unique equal-distance periodic bonds",
    )
    prep.add_argument(
        "--allow-periodic-bond-sum",
        action="store_true",
        help="Allow non-unique POSCAR atom pairs and report the summed coupling over periodic images without multiplicity division",
    )
    prep.add_argument(
        "--background-axis",
        choices=["x", "y", "z"],
        help="Background spin axis. Defaults by kind: jiso=x, jani/biqua=z; SIA uses component defaults.",
    )
    prep.add_argument("--pair-axis", choices=["x", "y", "z"], default="z", help="Pair spin axis for Jiso")
    prep.add_argument("--metal-ligand-cutoff", type=float, default=4.5, help="Metal-ligand cutoff for Kitaev shared-ligand detection")
    prep.add_argument("--kitaev-pair-cutoff", type=float, default=10.0, help="Maximum metal-metal pair distance for Kitaev detection")
    prep.add_argument(
        "--kitaev-no-component-permutation",
        action="store_true",
        help="Use the legacy discrete row/sign match without Cartesian component permutation",
    )
    prep.add_argument(
        "--allow-kitaev-ligand-fallback",
        action="store_true",
        help="Allow nearest-ligand fallback if two shared ligands are not found within --metal-ligand-cutoff",
    )
    prep.add_argument(
        "--max-kitaev-axis-consistency-deg",
        type=float,
        default=25.0,
        help="Reject Kitaev compression if matched local x/y/z axes differ by more than this angle between the two sites",
    )
    prep.add_argument("--workflow", choices=["single", "pbe-hse"], default="single")
    prep.add_argument("--saxis", nargs=3, help="Optional SAXIS values; only the default 0 0 1 is currently supported")
    prep.add_argument("--no-stage-tag-edits", action="store_true", help="Do not add PBE/HSE stage INCAR tag edits")
    prep.add_argument("--force", action="store_true", help="Allow writing into a non-empty output directory")
    prep.set_defaults(func=prepare)

    boot = sub.add_parser("bootstrap", help="Run vaspkit to create a reusable input template directory")
    boot.add_argument("--poscar", required=True)
    boot.add_argument("--out", required=True)
    boot.add_argument("--vaspkit-command", default="vaspkit")
    boot.add_argument("--vaspkit-sequence", default="1 102 2 0.04")
    boot.add_argument("--ldau", action="store_true", help="Apply built-in example LDAU tags after vaspkit; review before production")
    boot.add_argument("--no-ldau", action="store_true", help=argparse.SUPPRESS)
    boot.set_defaults(func=bootstrap)

    coll = sub.add_parser("collect", help="Collect energies using metadata.json")
    coll.add_argument("--root", required=True)
    coll.set_defaults(func=collect)

    krep = sub.add_parser("kitaev-report", help="Detect a Kitaev frame and optionally rotate a Jani tensor into it")
    krep.add_argument("--poscar", required=True)
    krep.add_argument("--magnetic-elements", nargs="+", help="Magnetic element symbols, e.g. Cr")
    krep.add_argument("--ligand-elements", nargs="+", help="Ligand element symbols, e.g. I")
    krep.add_argument("--pair", action="append", required=True, help="Exactly one pair I-J")
    krep.add_argument("--index-mode", choices=["global", "magnetic"], default="global")
    krep.add_argument("--metal-ligand-cutoff", type=float, default=4.5)
    krep.add_argument("--kitaev-pair-cutoff", type=float, default=10.0)
    krep.add_argument(
        "--kitaev-no-component-permutation",
        action="store_true",
        help="Use the legacy discrete row/sign match without Cartesian component permutation",
    )
    krep.add_argument(
        "--allow-kitaev-ligand-fallback",
        action="store_true",
        help="Allow nearest-ligand fallback if two shared ligands are not found within --metal-ligand-cutoff",
    )
    krep.add_argument(
        "--max-kitaev-axis-consistency-deg",
        type=float,
        default=25.0,
        help="Reject common-frame model compression if matched local axes differ too much between the two sites",
    )
    krep.add_argument(
        "--pair-image-shift",
        nargs=3,
        type=int,
        help="Use this image shift only for geometry/reporting; tensor extraction still refers to POSCAR atom indices",
    )
    krep.add_argument("--bond-distance-tol", type=float, default=0.02)
    krep.add_argument(
        "--allow-periodic-bond-sum",
        action="store_true",
        help="Allow non-unique POSCAR atom pairs and label tensors as periodic-image sums",
    )
    krep.add_argument("--jani-root", help="Root of a completed Jani calculation with final_summary.txt")
    krep.add_argument("--stage", help="Stage name in final_summary.txt, e.g. single or hse_no_u")
    krep.add_argument("--out", help="Optional report file")
    krep.set_defaults(func=kitaev_report)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        args.func(args)
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
