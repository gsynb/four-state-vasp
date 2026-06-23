# four-state-vasp

一个用于 VASP 四态法磁相互作用计算的 Codex / Claude Code skill。它可以从 POSCAR 分析磁性原子近邻关系，生成 Jani、Jiso、SIA、biquadratic 和 Kitaev 相互作用的四态法计算目录，并在作业结束后提取能量与相互作用参数。

This is a Codex / Claude Code skill for VASP four-state magnetic-interaction workflows. It analyzes magnetic-neighbor shells from POSCAR files, generates VASP calculation trees for Jani, Jiso, SIA, biquadratic, and Kitaev interactions, and postprocesses finished jobs into meV-scale interaction parameters.

## Contributors / 贡献者

- [@gsynb](https://github.com/gsynb) - project owner; provided the scientific workflow requirements, reference-script provenance, and validation targets.

## 功能概览

- `neighbors`: 给定 POSCAR 后，列出 10 A 或指定 cutoff 内的所有磁性原子近邻，按壳层给出代表 pair，并提示周期边界和扩包风险。
- `bootstrap`: 在没有完整 VASP 输入模板时调用 `vaspkit` 生成 `INCAR`、`KPOINTS`、`POTCAR`、`POSCAR` 模板，并补充常用 LDAU/RWIGS 设置。
- `prepare --kind jani`: 为指定 pair 生成各向异性交换张量 `Jab` 的四态法目录。
- `prepare --kind jiso`: 为一个或多个 pair 生成各向同性交换 `Jiso` 的四态法目录。
- `prepare --kind sia`: 为指定磁性原子生成单离子各向异性 `SIA` 的四态法目录。
- `prepare --kind biqua`: 为指定 pair 生成 biquadratic 相互作用目录。
- `prepare --kind kitaev`: 根据共享配体和八面体方向识别 Kitaev `x/y/z` 轴，选择当前键的 `gamma` 轴并生成四态目录。
- `kitaev-report`: 将已完成的 Jani 张量旋转到 Kitaev 局域坐标，输出 `K_gamma_minus_alpha_beta_avg_meV` 等指标。
- `collect`: 收集 VASP 完成后的能量并写出 `final_summary.txt`。

## Feature Summary

- `neighbors`: list all magnetic-neighbor contacts within a cutoff, group them into shells, choose representative pairs, and warn about boundary/supercell issues.
- `bootstrap`: call `vaspkit` to create reusable VASP input templates when `INCAR`, `KPOINTS`, `POTCAR`, and `POSCAR` are not already available.
- `prepare --kind jani`: generate four-state anisotropic exchange tensor calculations.
- `prepare --kind jiso`: generate four-state isotropic exchange calculations.
- `prepare --kind sia`: generate single-ion anisotropy calculations.
- `prepare --kind biqua`: generate biquadratic interaction calculations.
- `prepare --kind kitaev`: detect octahedral Kitaev axes from shared ligands and generate the selected `gamma` projection.
- `kitaev-report`: rotate a completed Jani tensor into the matched Kitaev local frame.
- `collect`: extract final VASP energies and write `final_summary.txt`.

## 目录结构

```text
four-state-vasp/
├── SKILL.md
├── CONTRIBUTORS.md
├── agents/
│   └── openai.yaml
├── references/
│   ├── methods.md
│   └── bibliography.bib
└── scripts/
    └── four_state_vasp.py
```

`SKILL.md` 是 agent 入口说明；`scripts/four_state_vasp.py` 是实际生成和后处理脚本；`references/methods.md` 记录公式、状态定义和旧脚本来源；`references/bibliography.bib` 提供可复制的文献引用。

## 安装到 Codex

### Linux / macOS / HPC

```bash
mkdir -p ~/.codex/skills
git clone https://github.com/gsynb/four-state-vasp.git ~/.codex/skills/four-state-vasp
```

### Windows PowerShell

```powershell
New-Item -ItemType Directory -Force "$env:USERPROFILE\.codex\skills"
git clone https://github.com/gsynb/four-state-vasp.git "$env:USERPROFILE\.codex\skills\four-state-vasp"
```

安装后开启一个新的 Codex 会话。如果 skill 没有立刻出现在可用列表中，重启 Codex。

使用时可以直接说：

```text
用 four-state-vasp 分析这个 POSCAR 中 Cr 原子 10 A 内所有近邻，并建议代表 pair。
```

或者：

```text
用 four-state-vasp 给 POSCAR 里的 Cr14-Cr15 生成 Jani 四态法计算，workflow 用 pbe-hse。
```

## Install in Codex

Clone this repository into your personal Codex skills directory:

```bash
mkdir -p ~/.codex/skills
git clone https://github.com/gsynb/four-state-vasp.git ~/.codex/skills/four-state-vasp
```

Then start a new Codex session and ask for the skill by name, for example:

```text
Use four-state-vasp to analyze magnetic neighbor shells in this POSCAR and prepare a Jani calculation for pair 14-15.
```

## 安装到 Claude Code

Claude Code skills 是包含 `SKILL.md` 的目录。个人 skill 可放在 `~/.claude/skills/<skill-name>/SKILL.md`，项目 skill 可放在仓库的 `.claude/skills/<skill-name>/SKILL.md`。Claude Code 可以根据描述自动调用 skill，也可以用 `/skill-name` 手动调用。官方说明见 [Claude Code skills documentation](https://docs.anthropic.com/en/docs/claude-code/skills)。

### 个人安装，全项目可用

```bash
mkdir -p ~/.claude/skills
git clone https://github.com/gsynb/four-state-vasp.git ~/.claude/skills/four-state-vasp
```

进入任意项目后启动 Claude Code：

```bash
claude
```

然后可以手动调用：

```text
/four-state-vasp analyze POSCAR neighbor shells and prepare Jiso for pair 14-15
```

### 项目内安装，只对当前项目可用

```bash
mkdir -p .claude/skills
git clone https://github.com/gsynb/four-state-vasp.git .claude/skills/four-state-vasp
```

如果 Claude Code 已经在运行，新建顶层 skills 目录后建议重启一次 Claude Code；之后修改 `SKILL.md` 通常会被自动检测。

## Install in Claude Code

Claude Code skills are directories containing a `SKILL.md` file. Personal skills live at `~/.claude/skills/<skill-name>/SKILL.md` and are available across projects. Project skills live at `.claude/skills/<skill-name>/SKILL.md` and apply only to that repository.

Personal install:

```bash
mkdir -p ~/.claude/skills
git clone https://github.com/gsynb/four-state-vasp.git ~/.claude/skills/four-state-vasp
```

Project install:

```bash
mkdir -p .claude/skills
git clone https://github.com/gsynb/four-state-vasp.git .claude/skills/four-state-vasp
```

Invoke it directly with:

```text
/four-state-vasp prepare Kitaev local-frame analysis for pair 14-15
```

or ask naturally and let Claude Code load the skill when relevant.

## 基本依赖

- Python 3.10 或更新版本。
- VASP 输入模板目录，至少包含 `INCAR`、`KPOINTS`、`POSCAR`、`POTCAR`。
- `vaspkit` 可选，仅在需要 `bootstrap` 自动生成输入模板时使用。
- `bash`，用于运行生成的 `submit_all.sh`、`run_state.sh`、`postprocess.sh`。
- VASP 与集群提交环境由用户本地系统提供。

注意：不要把商业/授权限制文件如 `POTCAR` 上传到公开仓库。本仓库的 `.gitignore` 默认忽略常见 VASP 大输出和 `POTCAR`。

## Requirements

- Python 3.10 or newer.
- A VASP input template directory containing `INCAR`, `KPOINTS`, `POSCAR`, and `POTCAR`.
- Optional `vaspkit` for `bootstrap`.
- `bash` for generated submission and postprocessing scripts.
- VASP and the cluster scheduler are provided by your local/HPC environment.

Do not commit licensed or large VASP files such as `POTCAR`, `WAVECAR`, or `CHGCAR` to a public repository.

## 命令行使用

下面示例假设已经在仓库根目录：

```bash
python3 scripts/four_state_vasp.py --help
```

### 1. 分析磁性近邻

```bash
python3 scripts/four_state_vasp.py neighbors \
  --poscar POSCAR \
  --magnetic-elements Cr \
  --cutoff 10 \
  --out neighbor_analysis
```

输出：

- `neighbor_contacts.tsv`: cutoff 内所有磁性 pair，包括周期 image shift。
- `neighbor_representatives.tsv`: 每个近邻壳层的代表 pair。
- `supercell_warnings.txt`: 扩包、中心原子和周期边界提醒。

如果代表 shell 绕中心原子一圈时有等效 pair 穿过周期边界，四态法计算可能需要更大超胞，或者选择更中心的等效 pair。

### 2. 生成 VASP 输入模板

```bash
python3 scripts/four_state_vasp.py bootstrap \
  --poscar POSCAR \
  --out inputs \
  --vaspkit-sequence "1 102 2 0.04"
```

如果已有 `INCAR/KPOINTS/POSCAR/POTCAR` 模板目录，可以跳过这一步。

### 3. Jani

```bash
python3 scripts/four_state_vasp.py prepare \
  --kind jani \
  --poscar POSCAR \
  --input-dir inputs \
  --out pair_14_15_jani \
  --magnetic-elements Cr \
  --pair 14-15 \
  --moment 6 \
  --workflow pbe-hse \
  --saxis 0 0 1
```

结果为 `Jxx Jxy Jxz Jyx Jyy Jyz Jzx Jzy Jzz` 四态目录。公式见 `references/methods.md`。

### 4. Jiso

```bash
python3 scripts/four_state_vasp.py prepare \
  --kind jiso \
  --poscar POSCAR \
  --input-dir inputs \
  --out jiso_pairs \
  --magnetic-elements Cr \
  --pair 14-23,14-20,14-24 \
  --moment 6 \
  --background-axis x \
  --pair-axis z \
  --workflow pbe-hse
```

默认非 pair 磁性原子沿 `x`，选中 pair 沿 `z` 做 `upup/updn/dnup/dndn`。

### 5. SIA

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

输出 `Axy Axz Ayz Ayy_minus_Axx Azz_minus_Axx`。

### 6. Biquadratic

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

### 7. Kitaev

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

Kitaev 轴识别默认使用 `4.5 A` 金属-配体 cutoff。脚本会：

1. 找到选中金属 pair 的最近周期 image。
2. 搜索共享配体；不足两个时用 nearest-six 和距离和 fallback。
3. 用共享边配体建立 pair reference basis。
4. 用附近八面体配体方向稳定手性。
5. 将理想 Kitaev `x/y/z` 轴通过行/符号和笛卡尔分量置换匹配到局域结构。
6. 选择与金属-金属键最垂直的轴作为 `gamma`。

如需复现旧的 generated_materials material-specific `axis1.py` 行/符号匹配方式：

```bash
python3 scripts/four_state_vasp.py kitaev-report \
  --poscar POSCAR \
  --magnetic-elements Cr \
  --ligand-elements Cl \
  --pair 7-26 \
  --kitaev-no-component-permutation
```

### 8. 提交和后处理

进入输出目录后：

```bash
bash submit_all.sh
```

等 VASP 作业结束后：

```bash
bash postprocess.sh
```

或者：

```bash
python3 scripts/four_state_vasp.py collect --root pair_14_15_jani
```

最终结果写入输出目录下的 `final_summary.txt`，能量单位为 meV。

## CLI Quick Start

Analyze neighbors:

```bash
python3 scripts/four_state_vasp.py neighbors \
  --poscar POSCAR \
  --magnetic-elements Cr \
  --cutoff 10 \
  --out neighbor_analysis
```

Generate a Jani calculation:

```bash
python3 scripts/four_state_vasp.py prepare \
  --kind jani \
  --poscar POSCAR \
  --input-dir inputs \
  --out pair_14_15_jani \
  --magnetic-elements Cr \
  --pair 14-15 \
  --moment 6 \
  --workflow pbe-hse
```

Generate a Kitaev projection:

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

Rotate a completed Jani tensor into the Kitaev frame:

```bash
python3 scripts/four_state_vasp.py kitaev-report \
  --poscar POSCAR \
  --magnetic-elements Cr \
  --ligand-elements I \
  --pair 14-15 \
  --jani-root pair_14_15_jani \
  --stage hse_no_u
```

Collect results:

```bash
python3 scripts/four_state_vasp.py collect --root pair_14_15_jani
```

## 索引规则

默认使用 POSCAR 全局 1-based 原子编号：

```bash
--index-mode global
```

如果你要复现旧脚本中“磁性原子列表内部编号”的 pair 写法，使用：

```bash
--index-mode magnetic
```

每次生成后都建议检查：

- `pair_indexing.tsv`
- `sia_target.tsv`
- `metadata.json`

这些文件会明确记录用户输入编号、POSCAR 全局编号和磁性原子序号之间的映射。

## Indexing Rules

The default index mode is POSCAR global 1-based atom indexing:

```bash
--index-mode global
```

Use magnetic-ion ordinal indexing only when reproducing older scripts:

```bash
--index-mode magnetic
```

Always inspect `pair_indexing.tsv`, `sia_target.tsv`, and `metadata.json` after generation.

## 常见工作流

### 用户只给 POSCAR，还没选 pair

1. 运行 `neighbors --cutoff 10`。
2. 查看所有近邻和代表 pair。
3. 如果某个 shell 穿过周期边界，先提醒扩包倍数或选择更中心的等效 pair。
4. 再运行 `prepare`。

### 已经知道 pair

1. 确认编号是 POSCAR 全局编号还是磁性原子内部编号。
2. 确认 `--magnetic-elements` 和 `--moment`。
3. 有配体依赖时显式传入 `--ligand-elements`。
4. 运行对应 `prepare --kind ...`。
5. 提交，完成后 `postprocess.sh` 或 `collect`。

### 想算真正 Kitaev 强度

1. 对 pair 先做完整 `jani`。
2. VASP 完成后 `collect`。
3. 运行 `kitaev-report --jani-root <jani-root>`。
4. 优先查看 `K_gamma_minus_alpha_beta_avg_meV`；也可同时比较 `K_gamma_minus_trace_iso_meV`。

## Notes

- Generated energies are reported in meV.
- `--workflow pbe-hse` creates both PBE+U preconvergence and HSE no-U directories.
- `prepare --kind kitaev` directly computes the selected `J_gamma_gamma` projection. For a more complete Kitaev anisotropy estimate, run full `jani` first and then `kitaev-report`.
- The methods and state formulas are documented in `references/methods.md`.

## Literature / 文献引用

If you publish results generated with this workflow, cite the papers that match the calculation. BibTeX entries are provided in `references/bibliography.bib`.

- Four-state / spin-lattice energy mapping: H. Xiang, E. Kan, S.-H. Wei, M.-H. Whangbo, and X. G. Gong, Phys. Rev. B 84, 224429 (2011), DOI: [10.1103/PhysRevB.84.224429](https://doi.org/10.1103/PhysRevB.84.224429).
- Energy-mapping review: H. Xiang, C. Lee, H.-J. Koo, X. G. Gong, and M.-H. Whangbo, Dalton Trans. 42, 823-853 (2013), DOI: [10.1039/C2DT31662E](https://doi.org/10.1039/C2DT31662E).
- Generic four-state mapping reference: D. Sabani, C. Bacaksiz, and M. V. Milosevic, Phys. Rev. B 102, 014457 (2020), DOI: [10.1103/PhysRevB.102.014457](https://doi.org/10.1103/PhysRevB.102.014457).
- Kitaev model/local-axis context: A. Kitaev, Ann. Phys. 321, 2-111 (2006), DOI: [10.1016/j.aop.2005.10.005](https://doi.org/10.1016/j.aop.2005.10.005).
