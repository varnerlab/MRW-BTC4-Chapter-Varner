# Mathematical Models in Biotechnology — Chapter Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Produce a build-ready review chapter ("Mathematical Models in Biotechnology," Comprehensive Biotechnology 4e) with six sections and one reproducible Julia worked example per pillar.

**Architecture:** Prose lives in `chapter/` (plain `article` class, `Chapter.tex` + `chapter/sections/*.tex`, built by a `Makefile`). Code lives in `code/` behind a single `Include.jl` entry point with shared `src/` utilities and per-pillar subfolders. Each example script regenerates figures into `code/figs/`, which are copied into `chapter/figures/` for `\includegraphics`. Source code is vendored and adapted from four Varner-lab repos.

**Tech Stack:** Julia 1.10+ (DifferentialEquations.jl, JuMP + HiGHS for FBA LP, Flux.jl for deep learning, CairoMakie for figures, CSV/DataFrames/JSON). LaTeX (`article`, `natbib`/`amsmath`/`graphicx`), BibTeX.

## Global Constraints

- Language: **Julia throughout**; no Python in the chapter pipeline.
- LaTeX: plain **`article`** document class (publisher template swapped later).
- Build: **`Makefile`** in `chapter/`: `pdflatex → bibtex → pdflatex → pdflatex`.
- Math: **full derivations inline** (main text or `appendix.tex`).
- Author: **solo — Jeffrey D. Varner**, Robert Frederick Smith School of Chemical and Biomolecular Engineering, Cornell University, Ithaca, NY.
- Deep learning: **Flux.jl** (LSTM built-in; S4/HiPPO-LegS as a custom layer). CHO training data is **synthetic**, from the §2 ODE.
- Transformers: **brief conceptual subsection only**.
- Reproducibility: vendor source into `code/`; **commit `code/Manifest.toml`** (relax the `.gitignore` rule).
- Figures: generated into `code/figs/`, copied into `chapter/figures/`.
- Branch: `chapter-draft`. Deadline: **2026-07-24**.
- Source repos (clone fresh if `/tmp` copies are gone):
  - 5820 (local): `/Users/jdv27/Desktop/julia_work/CHEME-5820-Instances/Spring-2026/`
  - FBA: `github.com/varnerlab/Lecture-5430-FluxBalanceAnalysis`
  - Kompala (cybernetic): `github.com/varnerlab/Kompala-Model-LP-Paper`
  - §4 network motifs are written fresh (no repo). Gluconate model (`/Users/jdv27/Desktop/papers/Gluconate-Sensor-Model-Paper`) and Adhikari 2020 are **cited only**, not vendored.

---

## Phase 1 — Scaffold (Week 1, now–Jun 29)

### Task 1.1: Julia environment and entry point

**Files:**
- Create: `code/Project.toml`
- Create: `code/Include.jl`
- Create: `code/src/Runtime.jl`
- Modify: `.gitignore` (relax `Manifest*.toml` so `code/Manifest.toml` is tracked)

**Produces:** `Include.jl` defines globals `_ROOT`, `_PATH_TO_SRC`, `_PATH_TO_DATA`, `_PATH_TO_FIGS`; activates the env; includes `src/Runtime.jl`. `Runtime.jl` exports `figpath(name)::String` (returns `joinpath(_PATH_TO_FIGS, name)`) and `datapath(name)::String`.

- [ ] **Step 1: Create the project file**

```toml
# code/Project.toml
name = "BiotechChapter"
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
Flux = "587475ba-b771-5e3f-ad9e-33799f191a9c"
HiGHS = "87dc4568-4c63-4d18-b0c0-bb2238e4078b"
JSON = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
JuMP = "4076af6c-e467-56ae-b986-b466b2749572"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
```

- [ ] **Step 2: Create the entry point**

```julia
# code/Include.jl
const _ROOT = @__DIR__
const _PATH_TO_SRC  = joinpath(_ROOT, "src")
const _PATH_TO_DATA = joinpath(_ROOT, "data")
const _PATH_TO_FIGS = joinpath(_ROOT, "figs")

using Pkg
Pkg.activate(_ROOT)

using CSV, DataFrames, JSON, Statistics, LinearAlgebra, Random
using DifferentialEquations, JuMP, HiGHS, Flux, CairoMakie

include(joinpath(_PATH_TO_SRC, "Runtime.jl"))
```

- [ ] **Step 3: Create the runtime helper**

```julia
# code/src/Runtime.jl
figpath(name::AbstractString)  = joinpath(_PATH_TO_FIGS, name)
datapath(name::AbstractString) = joinpath(_PATH_TO_DATA, name)
```

- [ ] **Step 4: Instantiate and verify the environment resolves**

```bash
cd code && mkdir -p src data figs && \
julia --project=. -e 'using Pkg; Pkg.instantiate(); include("Include.jl"); println("OK: ", figpath("x.pdf"))'
```
Expected: ends with `OK: .../code/figs/x.pdf` and no resolver errors. (First run downloads packages; allow several minutes.)

- [ ] **Step 5: Relax the Manifest ignore rule**

In `.gitignore`, replace the `Manifest*.toml` line with `!code/Manifest.toml` exception (keep ignoring Manifests elsewhere):
```
Manifest*.toml
!code/Manifest.toml
```

- [ ] **Step 6: Commit**

```bash
git add code/Project.toml code/Manifest.toml code/Include.jl code/src/Runtime.jl .gitignore
git commit -m "Scaffold Julia environment and Include.jl entry point"
```

---

### Task 1.2: LaTeX skeleton and Makefile

**Files:**
- Create: `chapter/Chapter.tex`
- Create: `chapter/sections/{introduction,kinetics,fba,geneexpression,deeplearning,hybrid,appendix}.tex` (stubs)
- Create: `chapter/References.bib` (empty placeholder with one entry)
- Create: `chapter/Makefile`
- Create: `chapter/figures/.gitkeep`

**Produces:** `make` in `chapter/` builds `Chapter.pdf`. Section files are included via `\input{sections/<name>}`.

- [ ] **Step 1: Create the main file**

```latex
% chapter/Chapter.tex
\documentclass[11pt]{article}
\usepackage[margin=1in]{geometry}
\usepackage{amsmath,amssymb,graphicx,booktabs,siunitx}
\usepackage[numbers,sort&compress]{natbib}
\usepackage[hidelinks]{hyperref}
\graphicspath{{figures/}}
\title{Mathematical Models in Biotechnology}
\author{Jeffrey D. Varner\\
  Robert Frederick Smith School of Chemical and Biomolecular Engineering,\\
  Cornell University, Ithaca, NY, USA}
\date{}
\begin{document}
\maketitle
\begin{abstract}
Placeholder abstract.
\end{abstract}
\input{sections/introduction}
\input{sections/kinetics}
\input{sections/fba}
\input{sections/geneexpression}
\input{sections/deeplearning}
\input{sections/hybrid}
\appendix
\input{sections/appendix}
\bibliographystyle{unsrtnat}
\bibliography{References}
\end{document}
```

- [ ] **Step 2: Create section stubs**

Each `chapter/sections/<name>.tex` gets a single line so the build resolves, e.g. `introduction.tex`:
```latex
\section{Introduction}\label{sec:intro}
Placeholder.
```
Repeat with appropriate `\section{...}\label{...}` for `kinetics` (`\section{Kinetic and Mechanistic Models}\label{sec:kinetics}`), `fba` (`\section{Flux Balance Analysis}\label{sec:fba}`), `geneexpression` (`\section{Models of Gene Expression}\label{sec:gene}`), `deeplearning` (`\section{Deep Time-Series Models}\label{sec:dl}`), `hybrid` (`\section{Hybrid Models and Outlook}\label{sec:hybrid}`), `appendix` (`\section{Derivations}\label{sec:appendix}`).

- [ ] **Step 3: Seed the bibliography**

```bibtex
% chapter/References.bib
@article{Orth2010,
  author = {Orth, Jeffrey D. and Thiele, Ines and Palsson, Bernhard {\O}.},
  title = {What is flux balance analysis?},
  journal = {Nature Biotechnology}, volume = {28}, number = {3},
  pages = {245--248}, year = {2010}, doi = {10.1038/nbt.1614}
}
```

- [ ] **Step 4: Create the Makefile**

```makefile
# chapter/Makefile
DOC = Chapter
.PHONY: all clean
all: $(DOC).pdf
$(DOC).pdf: $(DOC).tex sections/*.tex References.bib
	pdflatex -interaction=nonstopmode $(DOC).tex
	-bibtex $(DOC)
	pdflatex -interaction=nonstopmode $(DOC).tex
	pdflatex -interaction=nonstopmode $(DOC).tex
clean:
	rm -f $(DOC).aux $(DOC).bbl $(DOC).blg $(DOC).log $(DOC).out $(DOC).pdf
```

- [ ] **Step 5: Build and verify the PDF is produced**

```bash
cd chapter && touch figures/.gitkeep && make 2>&1 | tail -5 && ls -la Chapter.pdf
```
Expected: `Chapter.pdf` exists; no fatal LaTeX errors (undefined-citation warnings are fine at this stage).

- [ ] **Step 6: Commit**

```bash
git add chapter/Chapter.tex chapter/sections chapter/References.bib chapter/Makefile chapter/figures/.gitkeep
git commit -m "Scaffold LaTeX skeleton, section stubs, and Makefile"
```

---

### Task 1.3: Vendor and smoke-test the kinetics example (CHO fed-batch ODE)

**Files:**
- Create: `code/kinetics/cho_model.jl` (the mechanistic ODE; adapted from 5820 Wk12 `src/{Kinetics,MassBalances,Parameters,Types}.jl` and `fedbatch-cho-model.md`)
- Create: `code/kinetics/run_cho.jl` (simulate + save figure + save synthetic dataset)
- Create: `code/data/.gitkeep`

**Source:** `/Users/jdv27/Desktop/julia_work/CHEME-5820-Instances/Spring-2026/CHEME-5820-Labs-Spring-2026/labs/week-12/L12d/src/` (Kinetics.jl, MassBalances.jl, Parameters.jl, Types.jl) and the model notes at `.../lectures/week-12/L12c/fedbatch-cho-model.md`.

**Produces:** `cho_rhs!(du,u,p,t)` for a 7-state system `u = [V, X, Glc, Gln, Lac, Amm, mAb]`; `default_cho_params()::NamedTuple`; `simulate_cho(; tspan, saveat)::DataFrame` with columns `t,V,X,Glc,Gln,Lac,Amm,mAb`. `run_cho.jl` writes `code/data/cho_trajectories.csv` (consumed by Phase 4) and `code/figs/cho_kinetics.pdf`.

- [ ] **Step 1: Vendor and adapt the model into `code/kinetics/cho_model.jl`**

Copy the rate laws (Monod growth with lactate/ammonia inhibition, death rate, Luedeking–Piret product formation, substrate uptake) and the fed-batch mass balances from the L12d `src/` files into a single `cho_model.jl`, replacing their module headers with plain functions that match the **Produces** signatures above. Implement the glucose-triggered feed as a callback or feed term `F(t)` per `fedbatch-cho-model.md`.

- [ ] **Step 2: Write the run script**

```julia
# code/kinetics/run_cho.jl
include(joinpath(@__DIR__, "..", "Include.jl"))
include(joinpath(@__DIR__, "cho_model.jl"))
df = simulate_cho(tspan=(0.0, 240.0), saveat=1.0)   # hours
CSV.write(datapath("cho_trajectories.csv"), df)
let fig = Figure()
    ax = Axis(fig[1,1], xlabel="time (h)", ylabel="conc.")
    lines!(ax, df.t, df.X;   label="biomass")
    lines!(ax, df.t, df.mAb; label="mAb")
    axislegend(ax)
    save(figpath("cho_kinetics.pdf"), fig)
end
println("rows=", nrow(df), " final_mAb=", df.mAb[end])
```

- [ ] **Step 3: Run it and verify outputs exist and are physical**

```bash
cd code && julia --project=. kinetics/run_cho.jl && \
test -f figs/cho_kinetics.pdf && test -f data/cho_trajectories.csv && echo FIGURE_OK
```
Expected: prints `rows=241 final_mAb=<positive number>`, then `FIGURE_OK`. Sanity invariants to confirm in the printout: biomass and mAb are non-negative and mAb is monotonically non-decreasing.

- [ ] **Step 4: Commit**

```bash
git add code/kinetics code/data/.gitkeep
git commit -m "Vendor and smoke-test CHO fed-batch mechanistic model"
```

---

### Task 1.4: Vendor and smoke-test the FBA example (urea cycle)

**Files:**
- Create: `code/fba/urea_cycle.jl` (stoichiometry + bounds; adapted from the 5430 repo urea-cycle example and `src/Stoichiometric.jl`)
- Create: `code/fba/run_fba.jl` (build LP with JuMP/HiGHS, solve, save flux figure)

**Source:** clone `github.com/varnerlab/Lecture-5430-FluxBalanceAnalysis`; use `CHEME-5450-Example-Solution-UreaCycle-S2026.ipynb` and `src/{Stoichiometric,Network,Types}.jl`.

**Produces:** `urea_cycle_model()::NamedTuple` with fields `S::Matrix{Float64}`, `reactions::Vector{String}`, `metabolites::Vector{String}`, `lb::Vector{Float64}`, `ub::Vector{Float64}`, `c::Vector{Float64}` (objective). `solve_fba(model)::DataFrame` with columns `reaction, flux`.

- [ ] **Step 1: Adapt the urea-cycle stoichiometry into `urea_cycle.jl`**

Translate the notebook's stoichiometric matrix, flux bounds, and objective into the `urea_cycle_model()` NamedTuple above (pure data; no notebook dependencies).

- [ ] **Step 2: Write the LP solve + run script**

```julia
# code/fba/run_fba.jl
include(joinpath(@__DIR__, "..", "Include.jl"))
include(joinpath(@__DIR__, "urea_cycle.jl"))
m = urea_cycle_model()
function solve_fba(m)
    model = Model(HiGHS.Optimizer); set_silent(model)
    n = length(m.reactions)
    @variable(model, m.lb[i] <= v[i=1:n] <= m.ub[i])
    @constraint(model, m.S * v .== 0)
    @objective(model, Max, sum(m.c[i]*v[i] for i in 1:n))
    optimize!(model)
    DataFrame(reaction=m.reactions, flux=value.(v))
end
res = solve_fba(m)
CSV.write(datapath("urea_fba_solution.csv"), res)
let fig = Figure()
    ax = Axis(fig[1,1], xticks=(1:nrow(res), res.reaction), ylabel="flux",
              xticklabelrotation=pi/4)
    barplot!(ax, 1:nrow(res), res.flux)
    save(figpath("urea_fba.pdf"), fig)
end
println("objective_flux=", res.flux[argmax(m.c)])
```

- [ ] **Step 3: Run and verify a feasible, mass-balanced solution**

```bash
cd code && julia --project=. fba/run_fba.jl && test -f figs/urea_fba.pdf && echo FBA_OK
```
Expected: prints a finite `objective_flux=...`, `FBA_OK`. Add a Julia assertion in the script-end that `maximum(abs.(m.S * res.flux)) < 1e-6` (steady-state mass balance holds).

- [ ] **Step 4: Commit**

```bash
git add code/fba
git commit -m "Vendor and smoke-test urea-cycle FBA example"
```

---

### Task 1.5: Vendor and smoke-test the cybernetic (Kompala) model

**Files:**
- Create: `code/geneexpression/cybernetic.jl` (Kompala diauxic-growth / LP-choice model; adapted from the Kompala repo)
- Create: `code/geneexpression/run_cybernetic.jl` (simulate mixed-substrate diauxie, save figure)

**Source:** clone `github.com/varnerlab/Kompala-Model-LP-Paper`; adapt `code/choice-simulation` and `code/discrete-cellmass-simulation`. (This is the only §4 model that needs vendoring; smoke-testing it in Phase 1 de-risks the one remote dependency early.)

**Produces:** `simulate_diauxie(; tspan)::DataFrame` with columns `t, biomass, S1, S2` (two substrates, sequential uptake). Writes `code/figs/cybernetic_diauxie.pdf`.

- [ ] **Step 1: Adapt the Kompala model into `cybernetic.jl`**

Translate the cybernetic variables / LP allocation of enzyme synthesis across the two substrates from the Kompala repo into `simulate_diauxie` matching the **Produces** signature (pure Julia; no notebook dependencies).

- [ ] **Step 2: Write the run script**

```julia
# code/geneexpression/run_cybernetic.jl
include(joinpath(@__DIR__, "..", "Include.jl"))
include(joinpath(@__DIR__, "cybernetic.jl"))
df = simulate_diauxie(tspan=(0.0, 20.0))
CSV.write(datapath("cybernetic_diauxie.csv"), df)
let fig = Figure()
    ax = Axis(fig[1,1], xlabel="time (h)", ylabel="conc.")
    lines!(ax, df.t, df.biomass; label="biomass")
    lines!(ax, df.t, df.S1;      label="substrate 1")
    lines!(ax, df.t, df.S2;      label="substrate 2")
    axislegend(ax)
    save(figpath("cybernetic_diauxie.pdf"), fig)
end
println("S1_end=", df.S1[end], " S2_end=", df.S2[end])
```

- [ ] **Step 3: Run and verify sequential substrate consumption (diauxie)**

```bash
cd code && julia --project=. geneexpression/run_cybernetic.jl && \
test -f figs/cybernetic_diauxie.pdf && echo CYB_OK
```
Expected: `CYB_OK`; the printout/figure shows S1 is depleted before S2 begins declining (diauxic lag); both substrates end near zero.

- [ ] **Step 4: Commit**

```bash
git add code/geneexpression/cybernetic.jl code/geneexpression/run_cybernetic.jl
git commit -m "Vendor and smoke-test Kompala cybernetic model"
```

---

### Task 1.6: Draft the Introduction (§1)

**Files:**
- Modify: `chapter/sections/introduction.tex`

**Produces:** Full Introduction prose establishing the mechanistic→data-driven spectrum, model purposes (prediction/design/control/understanding), cross-cutting themes (identifiability, calibration vs. validation, extrapolation), the two recurring systems, and a roadmap referencing `\ref{sec:kinetics}`…`\ref{sec:hybrid}`.

- [ ] **Step 1: Write the Introduction section** (replace the stub with ~1.5–2 pages of prose; introduce the CHO bioreactor and E. coli gene-regulatory-circuit running examples; forward-reference each section by label).

- [ ] **Step 2: Build and verify references resolve**

```bash
cd chapter && make 2>&1 | grep -iE "undefined|error" | head; ls -la Chapter.pdf
```
Expected: no "Reference ... undefined" for `sec:*` labels; `Chapter.pdf` regenerated.

- [ ] **Step 3: Commit**

```bash
git add chapter/sections/introduction.tex && git commit -m "Draft Introduction (section 1)"
```

---

### Task 1.7: Seed the bibliography

**Files:**
- Modify: `chapter/References.bib`

**Produces:** BibTeX entries for all anchor references named in the spec so later sections can `\cite` them immediately.

- [ ] **Step 1: Add entries** for: Michaelis & Menten 1913; Monod 1949; Luedeking & Piret 1959; Thiele & Palsson 2013; Heirendt 2019; **Alon 2006 (*An Introduction to Systems Biology*); Shen-Orr et al. 2002; Mangan & Alon 2003; Goentoro et al. 2009; Shoval & Alon 2010 (SnapShot: Network Motifs, *Cell* 143:326)**; Kompala 1986; Adhikari 2020 (`10.3389/fbioe.2020.539081`, cited as biophysical next level); Hochreiter & Schmidhuber 1997; Gu, Goel & Ré 2022; Vaswani 2017. Pull bibliographic data from the PDFs in the 5430 repo `docs/` and the Alon SnapShot reference list where available.

- [ ] **Step 2: Verify BibTeX parses**

```bash
cd chapter && make 2>&1 | grep -iE "bibtex|warning--" | head; echo done
```
Expected: no BibTeX syntax errors.

- [ ] **Step 3: Commit**

```bash
git add chapter/References.bib && git commit -m "Seed bibliography with anchor references"
```

---

## Phase 2 — Kinetics and FBA (Week 2, Jun 30–Jul 6)

### Task 2.1: §2 Kinetics prose, derivations, and figure

**Files:**
- Modify: `chapter/sections/kinetics.tex`
- Modify: `chapter/sections/appendix.tex` (MM quasi-steady-state derivation)
- Possibly modify: `code/kinetics/run_cho.jl` (final publication figure styling)

**Interfaces — Consumes:** `code/figs/cho_kinetics.pdf` (Task 1.3).

- [ ] **Step 1: Regenerate the publication CHO figure**

```bash
cd code && julia --project=. kinetics/run_cho.jl && cp figs/cho_kinetics.pdf ../chapter/figures/
```
Expected: `chapter/figures/cho_kinetics.pdf` updated.

- [ ] **Step 2: Write §2 prose** — enzyme kinetics, Monod growth, mass-action, Luedeking–Piret; build to the CHO fed-batch model; include the governing ODE system and a `\ref{fig:cho}` to the figure; add the Hockin–Mann sidebar. Put the full Michaelis–Menten quasi-steady-state derivation in `appendix.tex` and reference it.

- [ ] **Step 3: Build and verify figure + derivation render**

```bash
cd chapter && make 2>&1 | grep -iE "undefined|error|Warning: File" | head; ls figures/cho_kinetics.pdf
```
Expected: figure included; no undefined `fig:cho` reference.

- [ ] **Step 4: Commit**

```bash
git add chapter/sections/kinetics.tex chapter/sections/appendix.tex chapter/figures/cho_kinetics.pdf code/kinetics/run_cho.jl
git commit -m "Write section 2 (kinetics) with CHO model figure and MM derivation"
```

---

### Task 2.2: §3 FBA prose, derivation, and figure

**Files:**
- Modify: `chapter/sections/fba.tex`
- Modify: `chapter/sections/appendix.tex` (derivation of `S·v = 0`)

**Interfaces — Consumes:** `code/figs/urea_fba.pdf` and `code/data/urea_fba_solution.csv` (Task 1.4).

- [ ] **Step 1: Regenerate the FBA figure into the chapter**

```bash
cd code && julia --project=. fba/run_fba.jl && cp figs/urea_fba.pdf ../chapter/figures/
```

- [ ] **Step 2: Write §3 prose** — stoichiometric matrix, the FBA LP, exchange reactions, FVA, genome-scale/SVD note; present the urea-cycle result with `\ref{fig:fba}`. Put the open-species mole-balance derivation of `S·v = 0` (steady-state, constant volume) in `appendix.tex`.

- [ ] **Step 3: Build and verify**

```bash
cd chapter && make 2>&1 | grep -iE "undefined|error" | head; ls figures/urea_fba.pdf
```

- [ ] **Step 4: Commit**

```bash
git add chapter/sections/fba.tex chapter/sections/appendix.tex chapter/figures/urea_fba.pdf
git commit -m "Write section 3 (FBA) with urea-cycle example and S*v=0 derivation"
```

---

## Phase 3 — Gene expression: network motifs + cybernetic coda (Week 3, Jul 7–13)

### Task 3.1: Network-motif models and figure (the worked example)

**Files:**
- Create: `code/geneexpression/motifs.jl` (Hill input functions + four motif ODEs)
- Create: `code/geneexpression/run_motifs.jl` (simulate all four motifs, save one multi-panel figure)

**Produces:** `hill_act(x; K, n)`, `hill_rep(x; K, n)`; `simulate_nar(; tspan)`, `simulate_c1ffl(; pulse, tspan)`, `simulate_i1ffl(; tspan)`, `simulate_oscillator(; tspan)`, each returning a `DataFrame` keyed on `t`. `run_motifs.jl` writes `code/figs/motifs.pdf` (4 panels: NAR vs. simple regulation; C1-FFL persistence detector; I1-FFL pulse; oscillator).

- [ ] **Step 1: Implement Hill functions and the four motif ODEs**

```julia
# code/geneexpression/motifs.jl
hill_act(x; K=1.0, n=2) = x^n / (K^n + x^n)
hill_rep(x; K=1.0, n=2) = K^n / (K^n + x^n)

# Negative autoregulation vs. simple regulation (β chosen so both share steady state)
function simulate_nar(; tspan=(0.0, 10.0), β=10.0, α=1.0, K=1.0, n=2)
    f!(du,u,p,t) = begin
        du[1] = β - α*u[1]                       # simple regulation
        du[2] = β*hill_rep(u[2]; K=K, n=n) - α*u[2]  # NAR (note: rescale β to match Yst if desired)
    end
    sol = solve(ODEProblem(f!, [0.0,0.0], tspan), Tsit5(), saveat=0.05)
    DataFrame(t=sol.t, simple=getindex.(sol.u,1), nar=getindex.(sol.u,2))
end

# Coherent type-1 FFL with AND gate: X->Y, X->Z, Y->Z (persistence detector)
function simulate_c1ffl(; tspan=(0.0, 16.0), βy=1.0, βz=1.0, αy=1.0, αz=1.0,
                          Kxy=0.5, Kxz=0.5, Kyz=0.5, n=2, pulse=(2.0,4.0), step_on=8.0)
    Xsig(t) = (pulse[1] ≤ t ≤ pulse[2]) || (t ≥ step_on) ? 1.0 : 0.0   # transient pulse, then sustained step
    f!(du,u,p,t) = begin
        X = Xsig(t)
        du[1] = βy*hill_act(X; K=Kxy, n=n) - αy*u[1]                       # Y
        du[2] = βz*hill_act(X; K=Kxz, n=n)*hill_act(u[1]; K=Kyz, n=n) - αz*u[2]  # Z (AND)
    end
    sol = solve(ODEProblem(f!, [0.0,0.0], tspan), Tsit5(), saveat=0.05)
    DataFrame(t=sol.t, X=Xsig.(sol.t), Y=getindex.(sol.u,1), Z=getindex.(sol.u,2))
end

# Incoherent type-1 FFL: X->Y, X->Z, Y-|Z (pulse generator)
function simulate_i1ffl(; tspan=(0.0, 12.0), βy=1.0, βz=1.0, αy=1.0, αz=1.0,
                          Kxy=0.5, Kxz=0.5, Kyz=0.5, n=2, step_on=1.0)
    Xsig(t) = t ≥ step_on ? 1.0 : 0.0
    f!(du,u,p,t) = begin
        X = Xsig(t)
        du[1] = βy*hill_act(X; K=Kxy, n=n) - αy*u[1]                       # Y
        du[2] = βz*hill_act(X; K=Kxz, n=n)*hill_rep(u[1]; K=Kyz, n=n) - αz*u[2]  # Z (X AND not-Y)
    end
    sol = solve(ODEProblem(f!, [0.0,0.0], tspan), Tsit5(), saveat=0.05)
    DataFrame(t=sol.t, X=Xsig.(sol.t), Y=getindex.(sol.u,1), Z=getindex.(sol.u,2))
end

# Goodwin-style negative-feedback oscillator (3 variables, high Hill coefficient)
function simulate_oscillator(; tspan=(0.0, 80.0), a=1.0, b=0.1, c=1.0, d=0.1, e=1.0, f=0.1, K=1.0, n=9)
    g!(du,u,p,t) = begin
        du[1] = a*hill_rep(u[3]; K=K, n=n) - b*u[1]
        du[2] = c*u[1] - d*u[2]
        du[3] = e*u[2] - f*u[3]
    end
    sol = solve(ODEProblem(g!, [0.1,0.1,0.1], tspan), Tsit5(), saveat=0.1)
    DataFrame(t=sol.t, X=getindex.(sol.u,1), Y=getindex.(sol.u,2), Z=getindex.(sol.u,3))
end
```
(`solve`, `ODEProblem`, `Tsit5` come from DifferentialEquations.jl, already loaded by `Include.jl`.)

- [ ] **Step 2: Write the run script (multi-panel figure)**

```julia
# code/geneexpression/run_motifs.jl
include(joinpath(@__DIR__, "..", "Include.jl"))
include(joinpath(@__DIR__, "motifs.jl"))
nar  = simulate_nar();  c1 = simulate_c1ffl();  i1 = simulate_i1ffl();  osc = simulate_oscillator()
fig = Figure(size=(900,650))
let ax = Axis(fig[1,1], title="Negative autoregulation", xlabel="time", ylabel="Y")
    lines!(ax, nar.t, nar.simple; label="simple"); lines!(ax, nar.t, nar.nar; label="NAR"); axislegend(ax)
end
let ax = Axis(fig[1,2], title="C1-FFL (persistence detector)", xlabel="time")
    lines!(ax, c1.t, c1.X; label="X"); lines!(ax, c1.t, c1.Z; label="Z"); axislegend(ax)
end
let ax = Axis(fig[2,1], title="I1-FFL (pulse generator)", xlabel="time")
    lines!(ax, i1.t, i1.X; label="X"); lines!(ax, i1.t, i1.Z; label="Z"); axislegend(ax)
end
let ax = Axis(fig[2,2], title="Negative-feedback oscillator", xlabel="time")
    lines!(ax, osc.t, osc.X; label="X"); axislegend(ax)
end
save(figpath("motifs.pdf"), fig)
println("nar_halftime_simple<nar? ", nothing)  # visual check; see expectations below
```

- [ ] **Step 3: Run and verify the signature dynamics**

```bash
cd code && julia --project=. geneexpression/run_motifs.jl && test -f figs/motifs.pdf && echo MOTIFS_OK
```
Expected: `MOTIFS_OK`. Visual/assertion checks: (a) NAR reaches half its steady state faster than simple regulation; (b) C1-FFL Z stays ~0 during the brief X pulse but rises under the sustained X step (transient filtered); (c) I1-FFL Z rises then falls back (a pulse); (d) the oscillator's X shows ≥2 peaks over the window. Add `@assert` statements for (c) `maximum(i1.Z) > 1.5*i1.Z[end]` and (d) a simple peak count.

- [ ] **Step 4: Commit**

```bash
git add code/geneexpression/motifs.jl code/geneexpression/run_motifs.jl
git commit -m "Implement network-motif worked example (NAR, C1/I1-FFL, oscillator)"
```

---

### Task 3.2: §4 prose, derivations, and figure assembly

**Files:**
- Modify: `chapter/sections/geneexpression.tex`
- Modify: `chapter/sections/appendix.tex` (Hill input-function derivation; C1-FFL sign-sensitive-delay and I1-FFL pulse analysis; cybernetic allocation objective)

**Interfaces — Consumes:** `code/figs/motifs.pdf` (Task 3.1) and `code/figs/cybernetic_diauxie.pdf` (Task 1.5).

- [ ] **Step 1: Refresh both figures into the chapter**

```bash
cd code && julia --project=. geneexpression/run_motifs.jl && julia --project=. geneexpression/run_cybernetic.jl && \
cp figs/motifs.pdf figs/cybernetic_diauxie.pdf ../chapter/figures/
```

- [ ] **Step 2: Write §4 prose** — (4.1) regulatory input functions (Hill activation/repression, cooperativity); (4.2) network motifs as design principles, presenting NAR, C1-FFL, I1-FFL, and the oscillator with their governing ODEs and `\ref{fig:motifs}`; note that more granular biophysical promoter models (Adhikari 2020; the structured-promoter gluconate biosensor) are the next level of detail, with citations; (4.3) the cybernetic/Kompala coda for whole-cell resource allocation with `\ref{fig:cybernetic}`. Put the Hill-function derivation, the C1-FFL sign-sensitive-delay argument, the I1-FFL pulse analysis, and the cybernetic objective in `appendix.tex`.

- [ ] **Step 3: Build and verify both figures and cross-refs resolve**

```bash
cd chapter && make 2>&1 | grep -iE "undefined|error|Warning: File" | head
```
Expected: no undefined `fig:motifs` / `fig:cybernetic`; both PDFs included.

- [ ] **Step 4: Commit**

```bash
git add chapter/sections/geneexpression.tex chapter/sections/appendix.tex chapter/figures/motifs.pdf chapter/figures/cybernetic_diauxie.pdf
git commit -m "Write section 4 (network motifs + cybernetic coda) with figures and derivations"
```

---

## Phase 4 — Deep time-series and hybrid (Week 4, Jul 14–20)

### Task 4.1: LSTM forecaster on CHO data (Flux)

**Files:**
- Create: `code/deeplearning/data_prep.jl` (load `cho_trajectories.csv`, window into sequences, normalize, train/test split)
- Create: `code/deeplearning/lstm.jl` (Flux LSTM model + train loop)
- Create: `code/deeplearning/run_lstm.jl` (train, forecast, save figure + metrics)

**Interfaces — Consumes:** `code/data/cho_trajectories.csv` (Task 1.3).
**Produces:** `make_windows(df; lookback, horizon)::Tuple` returning `(Xtrain, Ytrain, Xtest, Ytest, scaler)`; `train_lstm(Xtrain, Ytrain; epochs)::Chain`; `forecast(model, Xtest)::Matrix`. Writes `code/figs/lstm_cho.pdf` and `code/data/lstm_metrics.csv` (columns `state, rmse`).

- [ ] **Step 1: Implement windowing/normalization** in `data_prep.jl` (fixed `lookback`, multi-step `horizon`; return a reversible `scaler`).

- [ ] **Step 2: Implement the Flux LSTM** in `lstm.jl` (`Chain(LSTM(in=>hidden), Dense(hidden=>out))`, MSE loss, ADAM).

- [ ] **Step 3: Write `run_lstm.jl`** to train, forecast the held-out horizon, compute per-state RMSE, and save the overlay figure (truth vs. forecast).

- [ ] **Step 4: Run and verify training converges and a forecast is produced**

```bash
cd code && julia --project=. deeplearning/run_lstm.jl && \
test -f figs/lstm_cho.pdf && test -f data/lstm_metrics.csv && echo LSTM_OK
```
Expected: loss decreases across epochs (printed), `LSTM_OK`, finite RMSEs in the CSV.

- [ ] **Step 5: Commit**

```bash
git add code/deeplearning/data_prep.jl code/deeplearning/lstm.jl code/deeplearning/run_lstm.jl
git commit -m "Implement LSTM forecaster on CHO data"
```

---

### Task 4.2: S4 / structured state-space forecaster (custom Flux layer)

**Files:**
- Create: `code/deeplearning/s4.jl` (HiPPO-LegS SSM layer; bilinear discretization; trainable readout)
- Create: `code/deeplearning/run_s4.jl` (train, forecast, save figure + metrics)

**Source:** 5820 Wk14 `L14d/src/Compute.jl` (`build_legS_matrices_mimo`, discretization) and the L14a HiPPO example.
**Interfaces — Consumes:** `make_windows` (Task 4.1), `code/data/cho_trajectories.csv`.
**Produces:** `S4Layer` (Flux-compatible, with `Flux.@functor`); `train_s4(...)::Chain`; writes `code/figs/s4_cho.pdf` and `code/data/s4_metrics.csv` (columns `state, rmse`).

- [ ] **Step 1: Implement the HiPPO-LegS matrices and bilinear discretization** in `s4.jl` (adapt `build_legS_matrices_mimo` and the discretization from Wk14 `Compute.jl`).

- [ ] **Step 2: Wrap the SSM as a Flux layer** (`struct S4Layer; A; B; C; D; end` with `Flux.@functor S4Layer` and a callable applying the discrete recurrence over a sequence), followed by a `Dense` readout.

- [ ] **Step 3: Write `run_s4.jl`** mirroring `run_lstm.jl` (same windows/scaler) so results are directly comparable; save figure and per-state RMSE.

- [ ] **Step 4: Run and verify**

```bash
cd code && julia --project=. deeplearning/run_s4.jl && \
test -f figs/s4_cho.pdf && test -f data/s4_metrics.csv && echo S4_OK
```
Expected: loss decreases; `S4_OK`; finite RMSEs.

- [ ] **Step 5: Commit**

```bash
git add code/deeplearning/s4.jl code/deeplearning/run_s4.jl
git commit -m "Implement S4/HiPPO-LegS forecaster on CHO data"
```

---

### Task 4.3: S4-vs-LSTM comparison figure

**Files:**
- Create: `code/deeplearning/run_comparison.jl`

**Interfaces — Consumes:** `code/data/{lstm_metrics,s4_metrics}.csv`; the trained models' forecasts.
**Produces:** `code/figs/s4_vs_lstm.pdf` (per-state RMSE bars and a representative forecast overlay).

- [ ] **Step 1: Write the comparison script** loading both metric CSVs and plotting grouped RMSE bars + one overlay panel (truth, LSTM, S4).

- [ ] **Step 2: Run and verify**

```bash
cd code && julia --project=. deeplearning/run_comparison.jl && test -f figs/s4_vs_lstm.pdf && echo CMP_OK
```

- [ ] **Step 3: Commit**

```bash
git add code/deeplearning/run_comparison.jl && git commit -m "Add S4-vs-LSTM comparison figure"
```

---

### Task 4.4: §5 prose and derivations (deep time-series)

**Files:**
- Modify: `chapter/sections/deeplearning.tex`
- Modify: `chapter/sections/appendix.tex` (BPTT; LSTM gate equations; S4 continuous→discrete bilinear recurrence + HiPPO)

**Interfaces — Consumes:** `code/figs/{lstm_cho,s4_cho,s4_vs_lstm}.pdf`.

- [ ] **Step 1: Copy figures into the chapter**

```bash
cd code && for f in lstm_cho s4_cho s4_vs_lstm; do cp figs/$f.pdf ../chapter/figures/; done
```

- [ ] **Step 2: Write §5 prose** — RNN→LSTM→S4 spine, data needs and extrapolation limits, the CHO S4-vs-LSTM result with figure refs, and the brief Transformer/attention conceptual subsection. Put BPTT, LSTM gate, and S4 discretization derivations in `appendix.tex`.

- [ ] **Step 3: Build and verify**

```bash
cd chapter && make 2>&1 | grep -iE "undefined|error|Warning: File" | head
```

- [ ] **Step 4: Commit**

```bash
git add chapter/sections/deeplearning.tex chapter/sections/appendix.tex chapter/figures/{lstm_cho,s4_cho,s4_vs_lstm}.pdf
git commit -m "Write section 5 (deep time-series) with S4-vs-LSTM example"
```

---

### Task 4.5: §6 Hybrid example and prose

**Files:**
- Create: `code/deeplearning/run_hybrid.jl` (mechanistic+ML: e.g., ML residual on top of the §2 ODE, or mechanism-informed features)
- Modify: `chapter/sections/hybrid.tex`

**Interfaces — Consumes:** `cho_model.jl`, `make_windows`, the trained forecasters.
**Produces:** `code/figs/hybrid_cho.pdf` comparing pure-mechanistic, pure-data, and hybrid forecasts; the §6 prose and outlook.

- [ ] **Step 1: Implement the hybrid** in `run_hybrid.jl` — train an ML model to predict the residual between the mechanistic CHO prediction and the (synthetic, noised) data, then add it back; plot all three forecasts.

- [ ] **Step 2: Run and verify the hybrid improves on at least one pure approach**

```bash
cd code && julia --project=. deeplearning/run_hybrid.jl && test -f figs/hybrid_cho.pdf && echo HYB_OK
```
Expected: printout shows hybrid RMSE ≤ pure-data RMSE on the held-out horizon; `HYB_OK`.

- [ ] **Step 3: Copy figure and write §6 prose + outlook**

```bash
cd code && cp figs/hybrid_cho.pdf ../chapter/figures/
```
Then write `hybrid.tex` (gray-box/SciML framing, the three-way comparison, outlook on identifiability, data scarcity, digital twins).

- [ ] **Step 4: Build and commit**

```bash
cd chapter && make 2>&1 | grep -iE "undefined|error" | head
cd .. && git add code/deeplearning/run_hybrid.jl chapter/sections/hybrid.tex chapter/figures/hybrid_cho.pdf
git commit -m "Write section 6 (hybrid + outlook) with mechanistic+ML example"
```

---

## Phase 5 — Integration and polish (Week 5, Jul 21–24)

### Task 5.1: Abstract, full reference pass, and end-to-end build

**Files:**
- Modify: `chapter/Chapter.tex` (real abstract)
- Modify: `chapter/References.bib` (fill any `\cite` gaps to ~60–100 entries)

- [ ] **Step 1: Write the abstract** (replace placeholder) summarizing the arc and contributions.

- [ ] **Step 2: Find undefined citations and fix them**

```bash
cd chapter && make 2>&1 | grep -i "Citation.*undefined" | sort -u
```
Expected after fixes: no undefined citations.

- [ ] **Step 3: Full clean build**

```bash
cd chapter && make clean && make 2>&1 | tail -3 && ls -la Chapter.pdf
```
Expected: `Chapter.pdf` builds clean; no undefined references/citations.

- [ ] **Step 4: Commit**

```bash
git add chapter/Chapter.tex chapter/References.bib
git commit -m "Add abstract and complete reference pass"
```

---

### Task 5.2: Reproducibility check and README

**Files:**
- Create: `README.md` (root; how to reproduce figures and build the PDF)

- [ ] **Step 1: Regenerate every figure from scratch** to confirm reproducibility

```bash
cd code && for s in kinetics/run_cho fba/run_fba geneexpression/run_motifs \
  geneexpression/run_cybernetic \
  deeplearning/run_lstm deeplearning/run_s4 deeplearning/run_comparison \
  deeplearning/run_hybrid; do echo "== $s =="; julia --project=. $s.jl || exit 1; done
```
Expected: all scripts exit 0; every referenced figure exists in `code/figs/`.

- [ ] **Step 2: Copy refreshed figures into the chapter and rebuild**

```bash
cd code && cp figs/*.pdf ../chapter/figures/ && cd ../chapter && make 2>&1 | tail -3
```

- [ ] **Step 3: Write the root README** (project overview; `julia --project=code code/<pillar>/run_*.jl` to regenerate figures; `cd chapter && make` to build the PDF; pointer to the spec).

- [ ] **Step 4: Commit**

```bash
git add README.md chapter/figures/*.pdf && git commit -m "Add README and verify end-to-end reproducibility"
```

---

### Task 5.3: Final read-through and submission build

- [ ] **Step 1: Self-review pass** — read `Chapter.pdf` for flow, consistent notation across sections, and that every figure/derivation is referenced in text.

- [ ] **Step 2: Fix issues** found (prose/notation/cross-refs) in the relevant `sections/*.tex`.

- [ ] **Step 3: Final clean build and confirm page count is reasonable**

```bash
cd chapter && make clean && make >/dev/null 2>&1 && \
julia -e 'println("pages≈", read(`pdfinfo Chapter.pdf`, String))' 2>/dev/null || echo "Chapter.pdf built"
```

- [ ] **Step 4: Commit and tag the submission**

```bash
git add -A && git commit -m "Final read-through fixes; submission build"
git tag chapter-submission-2026-07-24
```

---

## Self-Review (against the spec)

**Spec coverage:** §1 Intro → Task 1.6; §2 Kinetics → 1.3, 2.1; §3 FBA → 1.4, 2.2; §4 gene expression (Hill input functions + network motifs + cybernetic coda) → 1.5 (Kompala vendoring), 3.1 (motif models), 3.2 (prose + derivations); §5 deep TS (RNN/LSTM/S4 + Transformer note) → 4.1, 4.2, 4.3, 4.4; §6 hybrid/outlook → 4.5. Repo layout → 1.1, 1.2. Build/Makefile → 1.2. References → 1.7, 5.1. Reproducibility/Manifest → 1.1, 5.2. Two recurring systems → CHO (1.3→4.x), E. coli gene circuits (3.1, 3.2). Derivations inline → appendix updated in 2.1, 2.2, 3.2, 4.4. All spec sections map to tasks.

**De-risking cuts (from spec §8), if behind by Jul 18:** make the Kompala cybernetic coda (Tasks 1.5/3.2) review-only (cite, do not re-run); drop the negative-feedback-oscillator motif in Task 3.1 (keep NAR + C1-FFL + I1-FFL); lighten Task 4.5 (hybrid) to a conceptual figure. None block the build.

**Placeholder scan:** code steps show real Julia; LaTeX steps specify exact files and content; verification commands have expected output. Vendoring tasks name exact source files to copy and the exact target interface to expose — not "implement later."

**Type consistency:** `simulate_cho`→`cho_trajectories.csv`→`make_windows` (4.1) reused by `run_s4.jl` (4.2) and `run_hybrid.jl` (4.5); `solve_fba` returns `DataFrame(reaction,flux)` used in 1.4/2.2; metric CSVs `{lstm,s4}_metrics.csv` (cols `state,rmse`) consumed by `run_comparison.jl` (4.3). Names consistent across tasks.
