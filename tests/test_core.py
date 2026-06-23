import json
import importlib.util
import math
import os
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "four_state_vasp.py"
SPEC = importlib.util.spec_from_file_location("four_state_vasp", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
fsv = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = fsv
SPEC.loader.exec_module(fsv)


def write_mock_vasp_state(
    root: Path,
    rel: str,
    energy: float,
    m_constr: str,
    mw_rows: list[tuple[float, float, float]] | None = None,
    m_rows: list[tuple[float, float, float]] | None = None,
    lambda_mw_perp_rows: list[tuple[float, float, float]] | None = None,
    include_ep: bool = True,
    include_lambda: bool = True,
    include_mw: bool = True,
    include_lambda_mw_perp: bool = True,
) -> None:
    state = root / rel
    state.mkdir(parents=True)
    state.joinpath("INCAR").write_text(f"LAMBDA = 10\nM_CONSTR = {m_constr}\n")
    lines = [f" 1 F= {energy:.8f} E0= {energy:.8f} d E =0\n"]
    if include_ep:
        lines.append("E_p = 0.0\n")
    if include_lambda:
        lines.append("lambda = 10\n")
    rows = mw_rows or []
    if include_mw:
        m_rows = m_rows or rows
        lines.append("ion        MW_int                 M_int\n")
        for ion, (mw, m_int) in enumerate(zip(rows, m_rows), start=1):
            lines.append(
                f"{ion:3d} "
                f"{mw[0]: .8f} {mw[1]: .8f} {mw[2]: .8f} "
                f"{m_int[0]: .8f} {m_int[1]: .8f} {m_int[2]: .8f}\n"
            )
    if include_lambda_mw_perp:
        lambda_rows = lambda_mw_perp_rows or [(0.0, 0.0, 0.0) for _ in rows]
        lines.append("ion             lambda*MW_perp\n")
        for ion, vec in enumerate(lambda_rows, start=1):
            lines.append(f"{ion:3d} {vec[0]: .8f} {vec[1]: .8f} {vec[2]: .8f}\n")
    state.joinpath("OSZICAR.single").write_text("".join(lines))


class CoreFormulaAndGeometryTests(unittest.TestCase):
    def test_sia_diagonal_denominators_are_not_four(self):
        args = SimpleNamespace(spin_convention="unit_vector")
        self.assertEqual(fsv.sia_energy_denominator("Axy"), 4.0)
        self.assertEqual(fsv.sia_energy_denominator("Ayy_minus_Axx"), 2.0)
        self.assertEqual(fsv.energy_denominator(args, fsv.sia_energy_denominator("Azz_minus_Axx"), 2), 2.0)

    def test_energy_denominator_scales_spin_but_not_bond_multiplicity(self):
        args = SimpleNamespace(spin_convention="spin_S", spin_length_S=1.5)
        self.assertEqual(fsv.energy_denominator(args, 4.0, 2), 9.0)
        self.assertAlmostEqual(fsv.energy_denominator(args, 1.0, 4), 5.0625)
        with self.assertRaisesRegex(ValueError, "multiplicity"):
            fsv.energy_denominator(args, 4.0, 2, 3)

    def test_validate_template_incar_rejects_unsafe_tags(self):
        info = fsv.PoscarInfo(
            path=Path("POSCAR"),
            elements=["Cr", "I"],
            counts=[1, 2],
            atom_symbols=["Cr", "I", "I"],
            lattice=[(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)],
            frac_coords=[(0.0, 0.0, 0.0), (0.5, 0.0, 0.0), (0.0, 0.5, 0.0)],
        )
        with tempfile.TemporaryDirectory() as tmp:
            incar = Path(tmp) / "INCAR"
            incar.write_text("LAMBDA = 10\nRWIGS = 1.2 1.8\nLSORBIT = .TRUE.\n")
            fsv.validate_template_incar(incar, info, "jani", "single")

            incar.write_text("LAMBDA = 10\nRWIGS = 1.2 1.8\nLSORBIT = .FALSE.\n")
            with self.assertRaisesRegex(ValueError, "requires LSORBIT"):
                fsv.validate_template_incar(incar, info, "kitaev", "single")

            incar.write_text("LAMBDA = 10\nRWIGS = 1.2 1.8\nLSORBIT = .TRUE.\n")
            with self.assertRaisesRegex(ValueError, "without LSORBIT"):
                fsv.validate_template_incar(incar, info, "biqua", "single")

            incar.write_text("RWIGS = 1.2\nISPIN = 2\n")
            with self.assertRaisesRegex(ValueError, "missing positive LAMBDA"):
                fsv.validate_template_incar(incar, info, "jiso", "single")

            incar.write_text("LAMBDA = 10\nRWIGS = 1.2 1.8\nLSORBIT = .TRUE.\nSAXIS = 1 0 0\n")
            with self.assertRaisesRegex(ValueError, "non-default SAXIS"):
                fsv.validate_template_incar(incar, info, "jani", "single")

        with self.assertRaisesRegex(ValueError, "Only the default SAXIS"):
            fsv.validate_cli_saxis(SimpleNamespace(saxis=["0", "1", "0"]))

    def test_nearest_image_delta_matches_bruteforce_for_skew_cell(self):
        lattice = [(2.0, 0.0, 0.0), (1.7, 0.45, 0.0), (0.2, 0.1, 2.5)]
        frac_i = (0.02, 0.04, 0.1)
        frac_j = (0.71, 0.77, 0.14)
        disp_frac, shift = fsv.nearest_image_delta(frac_i, frac_j, lattice)
        dist = fsv.vec_norm(fsv.frac_to_cart(disp_frac, lattice))

        brute = []
        raw = fsv.vec_sub(frac_j, frac_i)
        for tx in range(-4, 5):
            for ty in range(-4, 5):
                for tz in range(-4, 5):
                    trial = (raw[0] + tx, raw[1] + ty, raw[2] + tz)
                    brute.append((fsv.vec_norm(fsv.frac_to_cart(trial, lattice)), (tx, ty, tz)))
        best_dist, best_shift = min(brute, key=lambda item: item[0])
        self.assertAlmostEqual(dist, best_dist, places=10)
        self.assertEqual(shift, best_shift)

    def test_nearest_image_delta_finds_audited_counterexample_and_random_cells(self):
        lattice = [
            (1.11338991, 0.0, 0.0),
            (2.68612083, 2.27728193, 0.0),
            (2.09824066, -3.98315157, 2.33785853),
        ]
        frac_i = (0.0, 0.0, 0.0)
        frac_j = (0.22154003, -0.27123778, 0.44527070)
        disp_frac, shift = fsv.nearest_image_delta(frac_i, frac_j, lattice)
        self.assertEqual(shift, (-3, 1, 0))
        self.assertAlmostEqual(fsv.vec_norm(fsv.frac_to_cart(disp_frac, lattice)), 1.06644576, places=7)

        rng = random.Random(7)
        for _ in range(20):
            a = (rng.uniform(1.0, 3.0), 0.0, 0.0)
            b = (rng.uniform(-2.0, 2.0), rng.uniform(1.0, 3.0), 0.0)
            c = (rng.uniform(-2.0, 2.0), rng.uniform(-2.0, 2.0), rng.uniform(1.0, 3.0))
            lattice = [a, b, c]
            frac_i = (rng.random(), rng.random(), rng.random())
            frac_j = (rng.random(), rng.random(), rng.random())
            disp_frac, shift = fsv.nearest_image_delta(frac_i, frac_j, lattice)
            dist = fsv.vec_norm(fsv.frac_to_cart(disp_frac, lattice))
            raw = fsv.vec_sub(frac_j, frac_i)
            brute = []
            for tx in range(-6, 7):
                for ty in range(-6, 7):
                    for tz in range(-6, 7):
                        trial = (raw[0] + tx, raw[1] + ty, raw[2] + tz)
                        brute.append((fsv.vec_norm(fsv.frac_to_cart(trial, lattice)), (tx, ty, tz)))
            best_dist, best_shift = min(
                brute,
                key=lambda item: (
                    item[0],
                    abs(item[1][0]) + abs(item[1][1]) + abs(item[1][2]),
                    abs(item[1][0]),
                    abs(item[1][1]),
                    abs(item[1][2]),
                ),
            )
            self.assertAlmostEqual(dist, best_dist, places=10)
            self.assertEqual(shift, best_shift)

    def test_non_unique_periodic_pair_is_rejected_by_default(self):
        info = fsv.PoscarInfo(
            path=Path("POSCAR"),
            elements=["Cr"],
            counts=[2],
            atom_symbols=["Cr", "Cr"],
            lattice=[(1.0, 0.0, 0.0), (0.0, 4.0, 0.0), (0.0, 0.0, 4.0)],
            frac_coords=[(0.0, 0.0, 0.0), (0.5, 0.0, 0.0)],
        )
        pair = fsv.PairInfo(1, 2, 0, 1, 1, 2)
        args = SimpleNamespace(pair_image_shift=None, bond_distance_tol=1e-8, allow_periodic_bond_sum=False)
        bond = fsv.pair_bond_context(info, pair, args)
        self.assertEqual(bond.multiplicity, 2)
        with self.assertRaisesRegex(ValueError, "not a unique explicit pair"):
            fsv.validate_pair_bond_context("jani", pair, bond, args)

    def test_non_nearest_pair_image_shift_is_rejected_by_default(self):
        info = fsv.PoscarInfo(
            path=Path("POSCAR"),
            elements=["Cr"],
            counts=[2],
            atom_symbols=["Cr", "Cr"],
            lattice=[(10.0, 0.0, 0.0), (0.0, 10.0, 0.0), (0.0, 0.0, 10.0)],
            frac_coords=[(0.0, 0.0, 0.0), (0.2, 0.0, 0.0)],
        )
        pair = fsv.PairInfo(1, 2, 0, 1, 1, 2)
        args = SimpleNamespace(pair_image_shift=(1, 0, 0), bond_distance_tol=1e-8, allow_periodic_bond_sum=False)
        with self.assertRaisesRegex(ValueError, "not the nearest image"):
            fsv.pair_bond_context(info, pair, args)

    def test_constraint_mode_version_guard(self):
        self.assertEqual(fsv.version_tuple("vasp.6.4.1"), (6, 4, 1))
        fsv.validate_constraint_mode(SimpleNamespace(constraint_mode=4, vasp_version="6.4.0"))
        fsv.validate_constraint_mode(SimpleNamespace(constraint_mode=1, vasp_version="6.3.2"))
        with self.assertRaisesRegex(ValueError, "requires VASP >= 6.4.0"):
            fsv.validate_constraint_mode(SimpleNamespace(constraint_mode=4, vasp_version="6.3.2"))

    def test_postprocess_parses_vasp_tables_and_checks_background_moments(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            metadata = {
                "kind": "jiso",
                "constraint_mode": 4,
                "stages": [{"name": "single", "base": "", "suffix": "single"}],
                "formulas": [
                    {
                        "kind": "jiso",
                        "label": "Jzz",
                        "quantity": "Jzz_meV",
                        "states": ["s1", "s2", "s3", "s4"],
                        "state_labels": ["upup", "updn", "dnup", "dndn"],
                        "energy_denominator": 4.0,
                        "hamiltonian_prefactor": 1.0,
                        "target_global_indices_1based": [1, 2],
                    }
                ],
            }
            (root / "metadata.json").write_text(json.dumps(metadata))
            (root / "postprocess.py").write_text(fsv.POSTPROCESS_PY)
            energies = [-10.0, -11.0, -12.0, -13.0]
            good_mw = [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
            bad_background_mw = [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 1.0, 0.0)]
            for idx, energy in enumerate(energies, start=1):
                write_mock_vasp_state(
                    root,
                    f"s{idx}",
                    energy,
                    "1 0 0 0 1 0 0 0 1",
                    bad_background_mw if idx == 4 else good_mw,
                )

            failed = subprocess.run(
                [sys.executable, "postprocess.py"],
                cwd=root,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertNotEqual(failed.returncode, 0)
            diagnostics = root / "results" / "constraint_diagnostics.tsv"
            self.assertTrue(diagnostics.exists())
            diag_text = diagnostics.read_text()
            self.assertIn("1,2,3", diag_text)
            self.assertIn("angle>5", diag_text)
            self.assertIn("False", diag_text)
            self.assertFalse((root / "final_summary.txt").exists())

            env = os.environ.copy()
            env["FOUR_STATE_STRICT_CONSTRAINTS"] = "0"
            relaxed = subprocess.run(
                [sys.executable, "postprocess.py"],
                cwd=root,
                env=env,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(relaxed.returncode, 0, relaxed.stderr)
            self.assertIn("Jzz_meV", (root / "final_summary.txt").read_text())

    def test_postprocess_fails_on_signed_flip(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            metadata = {
                "kind": "jiso",
                "constraint_mode": 1,
                "stages": [{"name": "single", "base": "", "suffix": "single"}],
                "formulas": [
                    {
                        "kind": "jiso",
                        "label": "Jzz",
                        "quantity": "Jzz_meV",
                        "states": ["s1", "s2", "s3", "s4"],
                        "state_labels": ["upup", "updn", "dnup", "dndn"],
                        "energy_denominator": 4.0,
                        "hamiltonian_prefactor": 1.0,
                        "target_global_indices_1based": [1],
                    }
                ],
            }
            (root / "metadata.json").write_text(json.dumps(metadata))
            (root / "postprocess.py").write_text(fsv.POSTPROCESS_PY)
            for idx, energy in enumerate([-1.0, -2.0, -3.0, -4.0], start=1):
                mw = [(0.0, 0.0, -1.0)] if idx == 1 else [(0.0, 0.0, 1.0)]
                write_mock_vasp_state(root, f"s{idx}", energy, "0 0 1", mw)

            failed = subprocess.run(
                [sys.executable, "postprocess.py"],
                cwd=root,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertNotEqual(failed.returncode, 0)
            diag_text = (root / "results" / "constraint_diagnostics.tsv").read_text()
            self.assertIn("spin_sign_flip_ions=1", diag_text)
            self.assertIn("180", diag_text)
            self.assertFalse((root / "final_summary.txt").exists())

    def test_postprocess_fails_when_formula_target_is_not_constrained(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            metadata = {
                "kind": "jiso",
                "constraint_mode": 4,
                "stages": [{"name": "single", "base": "", "suffix": "single"}],
                "formulas": [
                    {
                        "kind": "jiso",
                        "label": "Jzz",
                        "quantity": "Jzz_meV",
                        "states": ["s1", "s2", "s3", "s4"],
                        "state_labels": ["upup", "updn", "dnup", "dndn"],
                        "energy_denominator": 4.0,
                        "hamiltonian_prefactor": 1.0,
                        "target_global_indices_1based": [2],
                    }
                ],
            }
            (root / "metadata.json").write_text(json.dumps(metadata))
            (root / "postprocess.py").write_text(fsv.POSTPROCESS_PY)
            rows = [(0.0, 0.0, 1.0), (0.0, 1.0, 0.0)]
            for idx, energy in enumerate([-1.0, -2.0, -3.0, -4.0], start=1):
                write_mock_vasp_state(root, f"s{idx}", energy, "0 0 1 0 0 0", rows)

            failed = subprocess.run(
                [sys.executable, "postprocess.py"],
                cwd=root,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertNotEqual(failed.returncode, 0)
            diag_text = (root / "results" / "constraint_diagnostics.tsv").read_text()
            self.assertIn("target_M_CONSTR_zero_ions=2", diag_text)
            self.assertFalse((root / "final_summary.txt").exists())

    def test_postprocess_fails_closed_when_constraint_data_are_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            metadata = {
                "kind": "jiso",
                "constraint_mode": 4,
                "stages": [{"name": "single", "base": "", "suffix": "single"}],
                "formulas": [
                    {
                        "kind": "jiso",
                        "label": "Jzz",
                        "quantity": "Jzz_meV",
                        "states": ["s1", "s2", "s3", "s4"],
                        "state_labels": ["upup", "updn", "dnup", "dndn"],
                        "energy_denominator": 4.0,
                        "hamiltonian_prefactor": 1.0,
                        "target_global_indices_1based": [1],
                    }
                ],
            }
            (root / "metadata.json").write_text(json.dumps(metadata))
            (root / "postprocess.py").write_text(fsv.POSTPROCESS_PY)
            for idx, energy in enumerate([-1.0, -2.0, -3.0, -4.0], start=1):
                write_mock_vasp_state(
                    root,
                    f"s{idx}",
                    energy,
                    "0 0 1",
                    include_ep=False,
                    include_lambda=False,
                    include_mw=False,
                    include_lambda_mw_perp=False,
                )

            failed = subprocess.run(
                [sys.executable, "postprocess.py"],
                cwd=root,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertNotEqual(failed.returncode, 0)
            diag_text = (root / "results" / "constraint_diagnostics.tsv").read_text()
            self.assertIn("E_p_missing", diag_text)
            self.assertIn("lambda_missing", diag_text)
            self.assertIn("MW_int_missing", diag_text)
            self.assertFalse((root / "final_summary.txt").exists())

    def test_physical_kitaev_decomposition_ignores_local_gauge_rotation(self):
        theta = math.radians(37.0)
        local_x_j = (math.cos(theta), math.sin(theta), 0.0)
        local_y_j = (-math.sin(theta), math.cos(theta), 0.0)
        frame = fsv.KitaevFrame(
            pair=fsv.PairInfo(1, 2, 0, 1, 1, 2),
            gamma_axis=(0.0, 0.0, 1.0),
            local_x=(1.0, 0.0, 0.0),
            local_y=(0.0, 1.0, 0.0),
            gamma_axis_j=(0.0, 0.0, 1.0),
            local_x_j=local_x_j,
            local_y_j=local_y_j,
            bond_axis=(1.0, 0.0, 0.0),
            pair_shift=(0, 0, 0),
            shared_ligands=[],
            gamma_label="z",
            alpha_label="x",
            beta_label="y",
            kitaev_axes=[("x", (1.0, 0.0, 0.0)), ("y", (0.0, 1.0, 0.0)), ("z", (0.0, 0.0, 1.0))],
            kitaev_axes_j=[("x", local_x_j), ("y", local_y_j), ("z", (0.0, 0.0, 1.0))],
            reference_basis=[],
            reference_basis_j=[],
            axis_overlaps=(1.0, 1.0, 1.0),
            axis_overlaps_j=(1.0, 1.0, 1.0),
            axis_bond_dots=(1.0, 0.0, 0.0),
            axis_consistency_degrees=(0.0, 0.0, 0.0),
            ligand_method="test",
            axis_match="test",
            axis_match_j="test",
            octahedral_ligands=[],
            octahedral_ligands_j=[],
            octahedral_angle_error=None,
            octahedral_angle_error_j=None,
        )
        j0 = 12.3
        global_matrix = [[j0, 0.0, 0.0], [0.0, j0, 0.0], [0.0, 0.0, j0]]
        common_matrix, local_gauge_matrix, _ = fsv.kitaev_exchange_matrices(global_matrix, frame)
        physical = fsv.decompose_common_exchange(common_matrix)
        self.assertAlmostEqual(physical["K_gamma_minus_alpha_beta_avg_meV"], 0.0, places=10)
        self.assertAlmostEqual(physical["Gamma_alpha_beta_meV"], 0.0, places=10)
        self.assertAlmostEqual(physical["Gamma_prime_avg_meV"], 0.0, places=10)
        self.assertAlmostEqual(physical["DMI_alpha_meV"], 0.0, places=10)
        self.assertAlmostEqual(physical["DMI_beta_meV"], 0.0, places=10)
        self.assertAlmostEqual(physical["DMI_gamma_meV"], 0.0, places=10)
        self.assertNotAlmostEqual(fsv.antisymmetric_dmi(local_gauge_matrix)[2], 0.0, places=6)

    def test_kitaev_axes_rotate_with_reference_basis(self):
        reference = [(0.0, 1.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 0.0, 1.0)]
        basis, _, match = fsv.align_ideal_kitaev_basis(reference)
        expected = [
            fsv.basis_linear_combination(row, reference, f"expected {label}")
            for label, row in zip(fsv.KITAEV_AXIS_LABELS, fsv.IDEAL_KITAEV_BASIS)
        ]
        self.assertEqual(match, "continuous_reference_projection")
        for got, want in zip(basis, expected):
            self.assertAlmostEqual(abs(fsv.vec_dot(got, want)), 1.0, places=10)

    def test_parse_jani_summary_requires_stage_when_multiple_are_present(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            lines = []
            for stage, offset in (("pbe_pre", 0.0), ("hse_no_u", 100.0)):
                lines.append(f"[{stage}]")
                lines.append(f"# stage {stage}")
                lines.append("# label quantity E1 E2 E3 E4 value_meV")
                for idx, comp in enumerate(fsv.JANI_COMPONENTS, start=1):
                    lines.append(f"{comp} {comp}_meV 0 0 0 0 {offset + idx:.8f}")
            (root / "final_summary.txt").write_text("\n".join(lines) + "\n")

            with self.assertRaisesRegex(ValueError, "multiple stages"):
                fsv.parse_jani_summary(root, "pair_1_2", None)
            matrix, stage = fsv.parse_jani_summary(root, "pair_1_2", "hse_no_u")
            self.assertEqual(stage, "hse_no_u")
            self.assertEqual(matrix[0][0], 101.0)


if __name__ == "__main__":
    unittest.main()
