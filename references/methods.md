# Four-State Method Notes

## Reference Provenance

The skill consolidates logic from these user reference locations:

- SIA: `/public/home/gaoshiyan/cri2/single_fix/4state.py`
- Jiso: `/public/home/gaoshiyan/cri2hse_4state/iso_4state`
- Jani: `/public/home/gaoshiyan/cri2/notsame`, especially `do.sh`, `4state.py`, `get-E.sh`
- Biquadratic: `/public/home/gaoshiyan/cri2/biqua`, especially `biquadratic_4state.py`
- Kitaev axes: `/public/home/gaoshiyan/generated_materials`, especially `4state_easy/axis1.py` and the material-specific `axis1.py` / `axis.txt` examples.

In the Codex desktop workspace these paths usually correspond to `X:\cri2\...` and `X:\cri2hse_4state\...`.

## Suggested Citations

If results from this workflow are used in a manuscript, cite the method papers that match the calculation. BibTeX entries are collected in `bibliography.bib`.

- Four-state / spin-lattice energy mapping: Xiang et al., Phys. Rev. B 84, 224429 (2011), DOI `10.1103/PhysRevB.84.224429`.
- Energy-mapping analysis review: Xiang et al., Dalton Trans. 42, 823-853 (2013), DOI `10.1039/C2DT31662E`.
- Generic four-state exchange mapping: Sabani et al., Phys. Rev. B 102, 014457 (2020), DOI `10.1103/PhysRevB.102.014457`.
- Kitaev local-axis model context: Kitaev, Ann. Phys. 321, 2-111 (2006), DOI `10.1016/j.aop.2005.10.005`.

## Conventions

The default mapping uses a positive coefficient convention. For the bilinear two-site tensor this is:

```text
H = + S_i^T J S_j
```

Use `--hamiltonian-sign minus` for the opposite coefficient sign convention. The generated metadata records the chosen convention, the spin convention, and the actual energy denominator used for every component.

By default spin directions are unit vectors. With `--spin-convention spin_S --spin-length-S S`, bilinear terms are divided by `S^2` and biquadratic terms by `S^4`. If the selected periodic pair represents multiple equivalent bonds within `--bond-distance-tol`, the denominator also includes that bond multiplicity. This matters in small supercells and for boundary-crossing pairs.

## Jani

Components:

```text
Jxx Jxy Jyy Jyz Jzz Jzx Jxz Jyx Jzy
```

For component `Jab`, the pair states are:

```text
E1: +a, +b
E2: +a, -b
E3: -a, +b
E4: -a, -b
```

The output directories are `<component>/1` through `<component>/4`.

Formula:

```text
Jab_meV = prefactor * (E1 - E2 - E3 + E4) * 1000 / denominator
```

For the default unit-vector convention, `denominator = 4 * bond_multiplicity`. For `--hamiltonian-sign minus`, `prefactor = -1`; otherwise `prefactor = +1`.

## Jiso

Pair states:

```text
upup: +z, +z
updn: +z, -z
dnup: -z, +z
dndn: -z, -z
```

The old `iso_4state/gen_files.py` used a non-pair magnetic background along `+x` and the selected pair along `+/-z`. In the consolidated script this mode computes the selected axis component `Jaa`, not the true trace average `(Jxx + Jyy + Jzz) / 3`. The default output key for `--pair-axis z` is therefore `Jzz_meV`.

Formula:

```text
Jaa_meV = prefactor * (E_upup - E_updn - E_dnup + E_dndn) * 1000 / denominator
```

For the default unit-vector convention, `denominator = 4 * bond_multiplicity`.

## SIA

Components:

```text
Axy Axz Ayz Ayy_minus_Axx Azz_minus_Axx
```

Offdiagonal components use the ordinary four-state denominator:

```text
Aab_meV = prefactor * (E1 - E2 - E3 + E4) * 1000 / denominator
```

For the default unit-vector convention, offdiagonal `denominator = 4`. Diagonal differences such as `Ayy_minus_Axx` and `Azz_minus_Axx` use `denominator = 2`, because the energy difference maps directly to the difference of two squared direction cosines. SIA is a single-site term and does not include pair bond multiplicity.

State definitions:

```text
Axy:
  E1 (+x,+y), E2 (+x,-y), E3 (-x,+y), E4 (-x,-y), background +z
Axz:
  E1 (+x,+z), E2 (+x,-z), E3 (-x,+z), E4 (-x,-z), background +y
Ayz:
  E1 (+y,+z), E2 (+y,-z), E3 (-y,+z), E4 (-y,-z), background +x
Ayy_minus_Axx:
  E1 +y, E2 +x, E3 -x, E4 -y, background +z
Azz_minus_Axx:
  E1 +z, E2 +x, E3 -x, E4 -z, background +y
```

The old `single_fix/4state.py` required an even number of special indices and used the first half as the measured atoms. The consolidated script directly targets one atom, with all other magnetic atoms as background.

## Biquadratic

This mode is intended for scalar non-SOC energy mapping. The template validator rejects `LSORBIT = .TRUE.` for `--kind biqua` to avoid mixing anisotropic SOC terms into the scalar fit.

Pair states:

```text
E1: S1=(1,0,0), S2=(1,0,0)
E2: S1=(1,0,0), S2=(-1,0,0)
E3: S1=(1,0,0), S2=(1/sqrt2,1/sqrt2,0)
E4: S1=(1,0,0), S2=(-1/sqrt2,-1/sqrt2,0)
```

Formulas:

```text
J_meV = prefactor * (E1 - E2) * 1000 / J_denominator
Biquadratic_B_meV = prefactor * (E1 + E2 - E3 - E4) * 1000 / B_denominator
```

For the default unit-vector convention, `J_denominator = 2 * bond_multiplicity` and `B_denominator = bond_multiplicity`.

## Indexing Caution

The old Jani/Jiso/biquadratic scripts often accepted 1-based magnetic-ion ordinals and then padded zeros for nonmagnetic atoms. The consolidated script defaults to 1-based POSCAR global atom indices. Use `--index-mode magnetic` only when reproducing an old magnetic-ion-index calculation.

## Magnetic Neighbor Shells

Run `four_state_vasp.py neighbors` before choosing a pair when the user provides only a POSCAR or asks which neighbor shell to calculate.

The script searches magnetic-magnetic contacts under periodic boundary conditions up to the requested cutoff, clusters distances into shells, and writes:

```text
neighbor_contacts.tsv
neighbor_representatives.tsv
supercell_warnings.txt
```

`image_shift_j = 0,0,0` means the representative contact is inside the current POSCAR cell. Nonzero shifts mean the nearest image crosses a periodic boundary. For four-state pair calculations, warn the user if a shell around the chosen center atom crosses the boundary: skipping those equivalent contacts or choosing a boundary pair can give a finite-size artifact. The expansion suggestion is a heuristic based on `2*r + margin` compared with the three cell heights.

## Kitaev Axis and Rotation

For edge-sharing octahedra, the script detects all three local Kitaev `x/y/z` axes before selecting the pair's `gamma` direction:

1. Find the nearest periodic image of the selected magnetic pair.
2. Find ligand images bonded to both magnetic atoms within `--metal-ligand-cutoff`; the Kitaev default is `4.5` A.
3. If fewer than two cutoff-shared ligands are found, the default behavior is to stop with an error. Use `--allow-kitaev-ligand-fallback` only for exploratory screening; it restores the nearest-six distance-sum fallback inspired by `generated_materials/4state_easy/axis1.py`.
4. Use the two shared edge ligands to build a pair reference basis for site `i` and a separate periodic-image-aware reference basis for site `j`.
5. Find local octahedral references around both magnetic sites by choosing three nearly orthogonal ligand directions. This stabilizes handedness and catches distorted octahedra.
6. Align the ideal Kitaev `x/y/z` basis by continuous projection into the detected reference basis. This is rotation-covariant, so a globally rotated POSCAR rotates the detected Kitaev axes instead of changing their labels through a Cartesian component permutation.
7. Choose `gamma` as the site-`i` matched Kitaev axis with the smallest absolute dot product against the metal-metal bond, then use the same label on the site-`j` local basis. The remaining labels are reported as `alpha` and `beta`; the report includes site-`i`/site-`j` axis consistency angles.

Use `--kitaev-no-component-permutation` only when deliberately comparing against older material-specific `axis1.py` behavior that used a discrete matching path.

`prepare --kind kitaev` creates four states along the detected `gamma` axis:

```text
pp: +gamma_i, +gamma_j
pm: +gamma_i, -gamma_j
mp: -gamma_i, +gamma_j
mm: -gamma_i, -gamma_j
```

The direct result is:

```text
J_gamma_gamma_meV = prefactor * (E_pp - E_pm - E_mp + E_mm) * 1000 / denominator
```

For the default unit-vector convention, `denominator = 4 * bond_multiplicity`.

To estimate the Kitaev anisotropy from a full anisotropic exchange tensor, first run `jani`, then run `kitaev-report --jani-root <jani-root>`. If `final_summary.txt` contains multiple bracketed stage blocks, pass `--stage` explicitly so the report does not silently select the wrong PBE/HSE stage. The report rotates the global J matrix into a two-site Kitaev frame, with rows in the site-`i` basis and columns in the site-`j` basis. It reports:

```text
kitaev_J_meV rows_i=(x,y,z) cols_j=(x,y,z)
local_J_meV rows_i=(alpha,beta,gamma) cols_j=(alpha,beta,gamma)
J_trace_iso_meV = trace(J_local) / 3
K_gamma_minus_alpha_beta_avg_meV = J_gamma_gamma - (J_alpha_alpha + J_beta_beta) / 2
traceless_gamma_anisotropy_meV = J_gamma_gamma - J_trace_iso
Gamma_alpha_beta_meV
Gamma_prime_avg_meV
DMI_alpha_meV, DMI_beta_meV, DMI_gamma_meV
```

Use `K_gamma_minus_alpha_beta_avg_meV` when the convention is to compare the gamma component against the average of the other two local diagonal components. The DMI values are extracted from the antisymmetric part of the two-site local tensor.
