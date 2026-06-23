import importlib.util
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


class CoreFormulaAndGeometryTests(unittest.TestCase):
    def test_sia_diagonal_denominators_are_not_four(self):
        args = SimpleNamespace(spin_convention="unit_vector")
        self.assertEqual(fsv.sia_energy_denominator("Axy"), 4.0)
        self.assertEqual(fsv.sia_energy_denominator("Ayy_minus_Axx"), 2.0)
        self.assertEqual(fsv.energy_denominator(args, fsv.sia_energy_denominator("Azz_minus_Axx"), 2), 2.0)

    def test_energy_denominator_includes_spin_and_bond_multiplicity(self):
        args = SimpleNamespace(spin_convention="spin_S", spin_length_S=1.5)
        self.assertEqual(fsv.energy_denominator(args, 4.0, 2, 3), 27.0)
        self.assertAlmostEqual(fsv.energy_denominator(args, 1.0, 4, 2), 10.125)

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
