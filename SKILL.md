---
name: four-state-vasp
description: Analyze POSCAR magnetic-neighbor shells and generate, submit, and postprocess VASP four-state magnetic-interaction calculations for Jani exchange tensors, Jiso isotropic exchange, single-ion anisotropy (SIA), biquadratic, and Kitaev interactions. Use when the user provides a POSCAR/input template and asks Codex to identify magnetic near-neighbor pairs within a cutoff, choose representative pairs, warn about supercell/boundary issues, calculate interactions between selected atoms, detect octahedral Kitaev axes from shared ligands, rotate Jani tensors into local Kitaev frames, or prepare VASP jobs using vaspkit-generated inputs, PBE+U preconvergence, HSE06 no-U, MAGMOM/M_CONSTR state generation, or energy extraction.
---

# Four-State VASP

## Core Workflow

1. Identify the requested interaction:
   - `jani`: anisotropic exchange tensor for a pair.
   - `jiso`: isotropic four-state exchange for one or more pairs.
   - `sia`: single-ion anisotropy for one atom.
   - `biqua`: biquadratic interaction for a pair.
   - `kitaev`: octahedral-axis Kitaev projection for a pair, or rotate finished Jani into a Kitaev local frame.
2. If the user gives a POSCAR but no pair, run `scripts/four_state_vasp.py neighbors` first. Report all magnetic-pair contacts within the cutoff, representative pairs by shell, and supercell/boundary warnings before preparing calculations.
3. Locate the POSCAR and input templates. If only POSCAR is provided, create an input template directory with `scripts/four_state_vasp.py bootstrap`, which runs vaspkit and updates RWIGS/LDAU defaults.
4. Confirm or infer magnetic element(s), moment, and indexing. Default to POSCAR global 1-based atom indices unless the user explicitly says magnetic-ion indices.
5. Generate the calculation tree with `scripts/four_state_vasp.py prepare`.
6. Submit jobs from the output root with `bash submit_all.sh`, or submit selected rows from `state_jobs.tsv`.
7. After VASP jobs finish, run `bash postprocess.sh` in the output root, or run `scripts/four_state_vasp.py collect --root <output-root>`.

Read `references/methods.md` when you need formulas, state definitions, or the provenance of the older reference scripts.

## Input Generation

Use vaspkit only when the input template directory does not already contain `INCAR`, `KPOINTS`, `POSCAR`, and `POTCAR`.

Linux/HPC example:

```bash
python3 /path/to/four-state-vasp/scripts/four_state_vasp.py bootstrap \
  --poscar POSCAR \
  --out inputs \
  --vaspkit-sequence "1 102 2 0.04"
```

Desktop Windows example:

```powershell
py -X utf8 C:\Users\admin\.codex\skills\four-state-vasp\scripts\four_state_vasp.py bootstrap `
  --poscar X:\path\POSCAR `
  --out X:\path\inputs `
  --vaspkit-sequence "1 102 2 0.04"
```

If vaspkit is unavailable, ask the user for a directory containing the four required VASP input files.

## Neighbor Analysis

When a user asks which pair to calculate, inspect magnetic-neighbor shells before generating VASP jobs:

```bash
python3 scripts/four_state_vasp.py neighbors \
  --poscar POSCAR \
  --magnetic-elements Cr \
  --cutoff 10 \
  --out neighbor_analysis
```

Read the stdout and generated files:

- `neighbor_contacts.tsv`: every magnetic-magnetic contact within the cutoff, including image shifts.
- `neighbor_representatives.tsv`: one representative pair per shell.
- `supercell_warnings.txt`: finite-size and boundary warnings.

Use `image_shift_j != 0,0,0` as a boundary warning. If the representative shell around the center atom crosses the periodic boundary, tell the user that four-state pair calculations may require a larger supercell or a more central equivalent pair. The expansion heuristic reports minimum multiples from `2*r + margin` compared with the cell heights.

## Generation Examples

Jani for pair 14-15 using POSCAR global indices:

```bash
python3 scripts/four_state_vasp.py prepare \
  --kind jani \
  --poscar POSCAR \
  --input-dir inputs \
  --out pair_14_15_4state \
  --magnetic-elements Cr \
  --pair 14-15 \
  --moment 6 \
  --workflow pbe-hse \
  --saxis 0 0 1
```

Jiso for several pairs:

```bash
python3 scripts/four_state_vasp.py prepare \
  --kind jiso \
  --poscar POSCAR \
  --input-dir inputs \
  --out jiso_pairs \
  --magnetic-elements Cr \
  --pair 14-23,14-20,14-24,14-22 \
  --moment 6 \
  --background-axis x \
  --pair-axis z \
  --workflow pbe-hse
```

SIA for atom 14:

```bash
python3 scripts/four_state_vasp.py prepare \
  --kind sia \
  --poscar POSCAR \
  --input-dir inputs \
  --out sia_Cr14 \
  --magnetic-elements Cr \
  --atom 14 \
  --moment 6 \
  --workflow pbe-hse
```

Biquadratic for pair 14-15:

```bash
python3 scripts/four_state_vasp.py prepare \
  --kind biqua \
  --poscar POSCAR \
  --input-dir inputs \
  --out biqua_pair_14_15 \
  --magnetic-elements Cr \
  --pair 14-15 \
  --moment 6 \
  --workflow pbe-hse
```

Kitaev projection for pair 14-15 using shared iodine ligands:

```bash
python3 scripts/four_state_vasp.py prepare \
  --kind kitaev \
  --poscar POSCAR \
  --input-dir inputs \
  --out kitaev_pair_14_15 \
  --magnetic-elements Cr \
  --ligand-elements I \
  --pair 14-15 \
  --moment 6 \
  --workflow pbe-hse
```

Rotate a completed Jani tensor into the octahedral Kitaev frame:

```bash
python3 scripts/four_state_vasp.py kitaev-report \
  --poscar POSCAR \
  --magnetic-elements Cr \
  --ligand-elements I \
  --pair 14-15 \
  --jani-root pair_14_15_4state \
  --stage hse_no_u
```

## Indexing Rules

- Use `--index-mode global` by default: atom labels are POSCAR 1-based atom indices.
- Use `--index-mode magnetic` only when the user is following the older scripts where pair labels are 1-based positions inside the magnetic-ion list.
- Always inspect `pair_indexing.tsv` or `sia_target.tsv` after generation. These files map user labels to POSCAR global indices and magnetic-ion ordinals.

The older Jiso/Jani/biquadratic scripts often used magnetic-ion indices. This skill defaults to global POSCAR indices to avoid silent mistakes when magnetic atoms are not the first POSCAR species.

## Output Layout

The generator writes:

- `state_jobs.tsv`: one state per row, with relative path and suggested Slurm job name.
- `run_state.sh`: runs one state; it supports both single-stage and PBE+U -> HSE no-U workflows.
- `submit_all.sh`: submits all states.
- `postprocess.sh` and `postprocess.py`: collect final energies, write `results/constraint_diagnostics.tsv`, and write `final_summary.txt` only when strict constraint diagnostics pass.
- `metadata.json`: directory map, formulas, stages, and POSCAR/magnetic-index metadata.
- For `kitaev`, `kitaev_frames.tsv` records the two-site matched `x/y/z` Kitaev axes, selected `gamma` axis, auxiliary `alpha/beta` axes, pair image shift, shared ligands, octahedral ligand references, axis overlaps, axis-consistency angles, and bond-axis dot products.

For `--workflow pbe-hse`, each state has both `pbe/<relpath>` and `<relpath>`. `run_state.sh` runs PBE+U first, then copies `WAVECAR/CHGCAR` into the HSE no-U directory.

## Defaults To Preserve

- Magnetic moment default: `6.0`.
- Jani and biquadratic background axis default: `z`.
- Jiso default: non-pair magnetic atoms along `x`, pair spins along `z`.
- SIA uses component-specific orthogonal backgrounds following the `single_fix/4state.py` logic.
- Kitaev axis detection defaults to nonmagnetic atoms as ligands unless `--ligand-elements` is provided; prefer explicit ligand elements for mixed-anion structures. The default metal-ligand cutoff is `4.5` A for Kitaev detection.
- Kitaev axis detection follows the improved `generated_materials` workflow but is now rotation-covariant: select two cutoff-shared edge ligands, build separate reference bases for the two magnetic sites, use nearby octahedral ligand bonds to stabilize handedness, continuously project the ideal `x/y/z` Kitaev basis into each local basis, then choose `gamma` as the site-i axis most perpendicular to the metal-metal bond.
- If two cutoff-shared ligands are not found, stop and ask the user to check `--ligand-elements`/`--metal-ligand-cutoff`; use `--allow-kitaev-ligand-fallback` only for exploratory geometry checks.
- Use `--kitaev-no-component-permutation` when reproducing older `generated_materials` runs that used a discrete row/sign match without Cartesian component permutation.
- Direct `prepare --kind kitaev` gives the `J_gamma_gamma` projection using site-i and site-j local `gamma` axes. For the physically useful Kitaev anisotropy, run full `jani` and then `kitaev-report`; inspect `physical_K_gamma_minus_alpha_beta_avg_meV`, `physical_traceless_gamma_anisotropy_meV`, `physical_Gamma_alpha_beta_meV`, `physical_Gamma_prime_avg_meV`, and `physical_DMI_*`. Treat `local_gauge_J_meV` and `kitaev_gauge_J_meV` as diagnostics only because their left/right bases differ.
- `prepare` validates templates before writing state trees: positive `LAMBDA`, enough positive `RWIGS`, no `ISPIN=2`, and default `SAXIS = 0 0 1` only; `jani`/`sia`/`kitaev` require `LSORBIT=.TRUE.`, while `biqua` rejects `LSORBIT=.TRUE.`.
- `prepare --kind jiso` reports the selected axis component `Jaa`, not the trace-average isotropic exchange.
- Formula denominators are stored per entry in `metadata.json`; SIA diagonal differences use base denominator 2, off-diagonal SIA and two-site four-state terms use base denominator 4, then spin convention is applied. Never divide by periodic bond multiplicity: `jani`/`jiso`/`kitaev`/`biqua` reject non-unique POSCAR pair images by default, and `--allow-periodic-bond-sum` reports the summed coupling over periodic translations.
- `bootstrap` only fills POTCAR-derived `RWIGS` by default; add `--ldau` only when the built-in example U values are intentionally desired and will be reviewed.
- `kitaev-report --jani-root` must be given `--stage` when `final_summary.txt` contains more than one stage.
- Energy unit in postprocessing: meV. Constraint diagnostics are written to `results/constraint_diagnostics.tsv`; strict thresholds can be adjusted with `FOUR_STATE_MAX_PENALTY_EV`, `FOUR_STATE_MAX_TARGET_ANGLE_DEG`, and `FOUR_STATE_STRICT_CONSTRAINTS=0`.
