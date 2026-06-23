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
Jab_meV = (E1 - E2 - E3 + E4) * 1000 / 4
```

## Jiso

Pair states:

```text
upup: +z, +z
updn: +z, -z
dnup: -z, +z
dndn: -z, -z
```

The old `iso_4state/gen_files.py` used a non-pair magnetic background along `+x` and the selected pair along `+/-z`.

Formula:

```text
Jiso_meV = (E_upup - E_updn - E_dnup + E_dndn) * 1000 / 4
```

## SIA

Components:

```text
Axy Axz Ayz Ayy_minus_Axx Azz_minus_Axx
```

All use:

```text
value_meV = (E1 - E2 - E3 + E4) * 1000 / 4
```

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

Pair states:

```text
E1: S1=(1,0,0), S2=(1,0,0)
E2: S1=(1,0,0), S2=(-1,0,0)
E3: S1=(1,0,0), S2=(1/sqrt2,1/sqrt2,0)
E4: S1=(1,0,0), S2=(-1/sqrt2,-1/sqrt2,0)
```

Formulas:

```text
J_meV = (E1 - E2) * 1000 / 2
Biquadratic_B_meV = (E1 + E2 - E3 - E4) * 1000
```

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
3. If fewer than two cutoff-shared ligands are found, retry with the nearest six ligand images around each magnetic atom and choose the best shared candidates by distance sum. This follows the robust fallback style used in `generated_materials/4state_easy/axis1.py`.
4. Use the two shared edge ligands to build a pair reference basis.
5. Find a local octahedral reference from nearby metal-ligand bond vectors by choosing three nearly orthogonal ligand directions. This stabilizes handedness and catches distorted octahedra.
6. Align the fixed ideal Kitaev `x/y/z` basis to the pair reference. By default the matcher allows row permutations, row signs, Cartesian component permutations, and component signs, which reproduces the improved generated-materials recognition for cases such as Co/Br and Ir/Te examples.
7. Choose `gamma` as the matched Kitaev axis with the smallest absolute dot product against the metal-metal bond. The remaining two axes are reported as `alpha` and `beta`.

Use `--kitaev-no-component-permutation` to reproduce older material-specific `axis1.py` outputs that only matched ideal rows by row order and sign.

`prepare --kind kitaev` creates four states along the detected `gamma` axis:

```text
pp: +gamma, +gamma
pm: +gamma, -gamma
mp: -gamma, +gamma
mm: -gamma, -gamma
```

The direct result is:

```text
J_gamma_gamma_meV = (E_pp - E_pm - E_mp + E_mm) * 1000 / 4
```

To estimate the Kitaev anisotropy from a full anisotropic exchange tensor, first run `jani`, then run `kitaev-report --jani-root <jani-root>`. It rotates the global J matrix into both the matched Kitaev `x/y/z` frame and the selected `(alpha,beta,gamma)` frame, then reports:

```text
kitaev_J_meV rows=(x,y,z)
local_J_meV rows=(alpha,beta,gamma)
Jiso_trace = trace(J_local) / 3
K_gamma_minus_trace_iso = J_gamma_gamma - Jiso_trace
K_gamma_minus_alpha_beta_avg = J_gamma_gamma - (J_alpha_alpha + J_beta_beta) / 2
```

Use the last quantity when the convention is to compare the gamma component against the average of the other two local diagonal components.
