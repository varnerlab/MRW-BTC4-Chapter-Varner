# FBA-Integrative-Framework Chapter Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rewrite the *Comprehensive Biotechnology 4e* chapter around flux balance analysis as an integrative framework ("bounds are the gateway; keep `μ` everywhere"), with two reproducible steady-state Julia examples (urea-cycle metabolic FBA; a `BSTModelKit.jl` branched-feedback toy whose regulation is recovered by opening the FBA bound) and two reviewed capstones.

**Architecture:** Julia code in `code/` produces figures + numbers first; LaTeX prose in `chapter/` is written against those artifacts. Example 1 reuses existing `code/fba` (JuMP+HiGHS, `S·v=0`). Example 2 depends on the registered `BSTModelKit.jl` package to generate synthetic "truth" data on a branched-feedback network, then adds a JuMP+HiGHS FBA layer on the same stoichiometry and shows the regulation gateway recovering the truth. Deep-learning / motif / hybrid content and the CHO ODE code are deleted; Michaelis–Menten and gene-expression material are re-homed into the bounds-gateway section.

**Tech Stack:** LaTeX (`article` class, `Makefile`, `pdflatex`→`bibtex`→`pdflatex`×2), `natbib`; Julia 1.10+ with `JuMP`, `HiGHS`, `BSTModelKit`, `CairoMakie`, `CSV`, `DataFrames`, `LinearAlgebra`.

## Global Constraints

- Deadline 2026-07-24; solo author Jeffrey D. Varner (Cornell).
- Official title stays `Mathematical Models in Biotechnology`.
- Thesis (must recur): flux bounds are the integrative gateway.
- Intellectual signature (must appear and stay consistent): keep the growth-dilution term `μ` — the honest metabolite balance is `S·v̂ = μx`, and `S·v̂ = 0` is a labeled special case; the expression balance keeps the `(θ+μ)` denominators.
- References policy: minimal, cite only what is built on; no survey-padding.
- Solver: `JuMP` + `HiGHS` everywhere (never GLPK).
- Style: no em-dashes in prose (use commas/parentheses/colons); no subsection headings inside Results/Discussion-style example sections; no section self-references ("as this section shows"); every prose file uses LaTeX line-numberable paragraphs. Run the `style-check` and `review-section` skills on each prose file before its commit.
- Reproducibility: each example runs via `cd code && julia --project=. <example>/run_*.jl` and regenerates its figure into `code/figs/`, which is copied to `chapter/figures/`. Commit `code/Manifest.toml`.
- Every prose task ends with a clean `cd chapter && make` (no undefined references, no missing citations in the `.log`).

---

## File Structure

**Code (`code/`):**
- `code/Project.toml` — modify: drop `Flux`, `DifferentialEquations`, the `OrdinaryDiffEqCore` pin; add `BSTModelKit`.
- `code/Include.jl` — modify: drop `DifferentialEquations`, `Flux`; add `using BSTModelKit`.
- `code/Manifest.toml` — regenerate (committed).
- `code/fba/urea_cycle.jl`, `code/fba/run_fba.jl` — keep; verify unchanged behavior (Example 1).
- `code/feedback/Branched-Feedback.toml` — create (copied from BSTModelKit paper example).
- `code/feedback/run_feedback.jl` — create: BST synthetic-data generation + FBA layer + figure.
- Delete: `code/kinetics/`, `code/deeplearning/`, `code/geneexpression/`.

**Prose (`chapter/`):**
- `chapter/Chapter.tex` — modify: new abstract, new `\input` list.
- Create: `sections/derivation.tex` (§2), `sections/linearprogram.tex` (§3), `sections/gateways.tex` (§4), `sections/example_urea.tex` (§5), `sections/example_feedback.tex` (§6), `sections/capstones.tex` (§7), `sections/outlook.tex` (§8).
- Rewrite: `sections/introduction.tex`, `sections/appendix.tex`.
- Delete: `sections/kinetics.tex`, `sections/fba.tex`, `sections/geneexpression.tex`, `sections/deeplearning.tex`, `sections/hybrid.tex`.
- `chapter/References.bib` — modify: trim to the minimal cited set below.
- `chapter/figures/` — receives `urea_fba.pdf`, `feedback_gateway.pdf`.

**Minimal reference set (`References.bib`) — these and no others unless a task adds one:**
Orth-Palsson-Thiele 2010 (What is FBA); Palsson/Bordbar 2014 (constraint-based review, for exchange-reaction figure lineage); Savageau 1976 + Savageau-Voit-Irvine 1987 (BST/S-systems); Vadhin-Varner 2026 (BSTModelKit.jl, arXiv:2603.19115); Vilkhovoy 2018 (sequence-specific CFPS); Adhikari 2020 (effective biophysical / partition-function TX-TL); Vilkhovoy 2023 (integrated dynamic cell-free, bioRxiv 2023.02.10.528035); Wayman 2019 (designer glycans); Allen-Palsson 2003 (sequence-based demands); eQuilibrator (Beber 2022); BRENDA (Chang 2021). Drop all deep-learning, motif, cybernetic, and CHO citations from the old bib.

---

## Phase 0 — Teardown and environment

### Task 1: Delete cut content and reduce `Chapter.tex` to the new skeleton

**Files:**
- Delete: `chapter/sections/kinetics.tex`, `chapter/sections/fba.tex`, `chapter/sections/geneexpression.tex`, `chapter/sections/deeplearning.tex`, `chapter/sections/hybrid.tex`
- Delete: `code/kinetics/`, `code/deeplearning/`, `code/geneexpression/`
- Delete figures: `code/figs/cho_kinetics.pdf`, `cybernetic_diauxie.pdf`, `hybrid_cho.pdf`, `lstm_cho.pdf`, `motifs.pdf`, `s4_cho.pdf`, `s4_vs_lstm.pdf` and the same names under `chapter/figures/`
- Modify: `chapter/Chapter.tex`

**Interfaces:**
- Produces: a `Chapter.tex` whose body is a new abstract placeholder plus `\input` lines for the eight new section files and the appendix; empty stub section files so the document compiles.

- [ ] **Step 1: Create empty stub section files** so the skeleton compiles

```bash
cd chapter/sections
for s in derivation linearprogram gateways example_urea example_feedback capstones outlook; do
  printf '%%%% %s.tex (stub)\n' "$s" > "$s.tex"
done
```

- [ ] **Step 2: Delete cut prose, code, and figures**

```bash
cd /Users/jdv27/Desktop/papers/MRW-BTC4-Chapter-Varner
git rm chapter/sections/kinetics.tex chapter/sections/fba.tex chapter/sections/geneexpression.tex chapter/sections/deeplearning.tex chapter/sections/hybrid.tex
git rm -r code/kinetics code/deeplearning code/geneexpression
git rm code/figs/cho_kinetics.pdf code/figs/cybernetic_diauxie.pdf code/figs/hybrid_cho.pdf code/figs/lstm_cho.pdf code/figs/motifs.pdf code/figs/s4_cho.pdf code/figs/s4_vs_lstm.pdf
git rm chapter/figures/cho_kinetics.pdf chapter/figures/cybernetic_diauxie.pdf chapter/figures/hybrid_cho.pdf chapter/figures/lstm_cho.pdf chapter/figures/motifs.pdf chapter/figures/s4_cho.pdf chapter/figures/s4_vs_lstm.pdf
```

- [ ] **Step 3: Replace the `\input` list and abstract in `Chapter.tex`**

In `chapter/Chapter.tex`, replace the entire `\begin{abstract}...\end{abstract}` block with a one-line placeholder `\begin{abstract}\input{sections/abstract_placeholder}\end{abstract}` OR a temporary sentence (the real abstract is Task 15), and replace the six old `\input` lines plus appendix with:

```latex
\input{sections/introduction}
\input{sections/derivation}
\input{sections/linearprogram}
\input{sections/gateways}
\input{sections/example_urea}
\input{sections/example_feedback}
\input{sections/capstones}
\input{sections/outlook}
\appendix
\input{sections/appendix}
```

- [ ] **Step 4: Blank the introduction and appendix to stubs** (real content in later tasks)

```bash
printf '%%%% introduction.tex (stub)\n' > chapter/sections/introduction.tex
printf '%%%% appendix.tex (stub)\n\\section{Derivations}\\label{sec:appendix}\n' > chapter/sections/appendix.tex
```

- [ ] **Step 5: Build to verify the skeleton compiles**

Run: `cd chapter && make`
Expected: `Chapter.pdf` builds; `Chapter.log` has no "File ... not found" and no undefined-input errors. (Undefined *references* are fine at this stage.)

- [ ] **Step 6: Commit**

```bash
git add -A
git commit -m "Tear down old chapter: delete DL/motif/hybrid/kinetics content, stub new spine"
```

### Task 2: Rewire the Julia environment to the new dependency set

**Files:**
- Modify: `code/Project.toml`, `code/Include.jl`
- Regenerate: `code/Manifest.toml`

**Interfaces:**
- Produces: an environment where `include("code/Include.jl")` loads `CSV, DataFrames, LinearAlgebra, Statistics, Random, JuMP, HiGHS, BSTModelKit, CairoMakie` and the `src/Runtime.jl` helpers (`datapath`, `figpath`). No `Flux`, no `DifferentialEquations`.

- [ ] **Step 1: Edit `code/Include.jl`** — change the `using` line to drop DiffEq/Flux and add BSTModelKit

```julia
using CSV, DataFrames, JSON, Statistics, LinearAlgebra, Random
using JuMP, HiGHS, BSTModelKit, CairoMakie
```

- [ ] **Step 2: Edit `code/Project.toml`** — remove `[deps]` entries `DifferentialEquations`, `Flux`, and the `[compat]`/`[extras]`/`[targets]` `OrdinaryDiffEqCore` blocks; leave the rest.

- [ ] **Step 3: Add BSTModelKit and resolve**

Run:
```bash
cd code && julia --project=. -e 'using Pkg; Pkg.add("BSTModelKit"); Pkg.rm("Flux"); Pkg.rm("DifferentialEquations"); Pkg.resolve()'
```
Expected: resolves cleanly; `BSTModelKit` appears in `Project.toml [deps]`.
(If `BSTModelKit` is unregistered, instead run `Pkg.add(url="https://github.com/varnerlab/BSTModelKit.jl.git")`.)

- [ ] **Step 4: Verify the environment loads**

Run: `cd code && julia --project=. -e 'include("Include.jl"); println("env ok")'`
Expected: prints `env ok` with no load errors.

- [ ] **Step 5: Commit**

```bash
git add code/Project.toml code/Manifest.toml code/Include.jl
git commit -m "Rewire Julia env: drop Flux/DiffEq, add BSTModelKit"
```

---

## Phase 1 — Build the reproducible examples (figures + numbers first)

### Task 3: Verify and regenerate Example 1 (urea-cycle FBA)

**Files:**
- Verify: `code/fba/urea_cycle.jl`, `code/fba/run_fba.jl`
- Produces: `code/figs/urea_fba.pdf`, `code/data/urea_fba_solution.csv`

**Interfaces:**
- Consumes: `urea_cycle_model()` returning NamedTuple `(S, reactions, metabolites, lb, ub, c)`.
- Produces: optimal urea-export flux magnitude `0.0328` mmol/gDW/h (the argininosuccinate-lyase `v2` capacity), used verbatim by §5 prose.

- [ ] **Step 1: Run the example**

Run: `cd code && julia --project=. fba/run_fba.jl`
Expected: prints `objective_flux=` with magnitude `0.0328` (to 3 s.f.); the `@assert maximum(abs.(m.S*res.flux)) < 1e-6` passes; `figs/urea_fba.pdf` regenerates.

- [ ] **Step 2: If the run errors** (e.g., a dropped dependency), read the error, restore only the needed `using` in `Include.jl`, re-run. Do not re-introduce Flux/DiffEq.

- [ ] **Step 3: Copy the figure into the chapter**

Run: `cp code/figs/urea_fba.pdf chapter/figures/urea_fba.pdf`

- [ ] **Step 4: Commit**

```bash
git add code/figs/urea_fba.pdf chapter/figures/urea_fba.pdf code/data/urea_fba_solution.csv
git commit -m "Regenerate urea-cycle FBA figure under new env"
```

### Task 4: Example 2a — synthetic "truth" from the BSTModelKit branched-feedback model

**Files:**
- Create: `code/feedback/Branched-Feedback.toml` (copy the model below)
- Create: `code/feedback/run_feedback.jl` (data-generation half; FBA half added in Task 5)

**Interfaces:**
- Produces: `feedback_truth()` returning `(species, Xss, reactions, vss)` where `species = ["A","B","C","D","E"]`, `reactions = ["r1","r2","r3","r4","r0","r5","r6"]`, `Xss` the steady-state concentrations, `vss` the per-reaction steady-state fluxes. Consumed by Task 5.

- [ ] **Step 1: Create `code/feedback/Branched-Feedback.toml`** with exactly:

```toml
[metadata]
author = "jdv27@cornell.edu"
version = "0.1"
description = "Branched pathway with product feedback inhibition (chapter Example 2)"

[model]
list_of_static_species = ["E1", "E2", "E3", "E4"]
list_of_dynamic_species = ["A", "B", "C", "D", "E"]
list_of_connection_records = [
    "r1::A --> B",
    "r2::B --> C",
    "r3::C --> D",
    "r4::C --> E",
    "r0::{} --> A",
    "r5::D --> {}",
    "r6::E --> {}"
]
list_of_kinetics_records = [
    "r0::{}",
    "r1::{A,E,E1}",
    "r2::{B,D,E2}",
    "r3::{C,E3}",
    "r4::{C,E4}",
    "r5::{D}",
    "r6::{E}"
]
list_of_stoichiometry_records = []
```

- [ ] **Step 2: Read the BSTModelKit rate API** so the flux extraction is correct

Run: `cd code && julia --project=. -e 'using BSTModelKit; println(pathof(BSTModelKit))'`
Then read the package's `src/Kinetics.jl` and `src/Balances.jl` at that path to find the function that returns per-reaction rates given a model and a state vector. Note its exact name and signature for Step 4.

- [ ] **Step 3: Write the data-generation function** in `code/feedback/run_feedback.jl`

```julia
include(joinpath(@__DIR__, "..", "Include.jl"))

# Build the branched-feedback S-system and return steady-state concentrations + fluxes.
function feedback_truth()
    path = joinpath(@__DIR__, "Branched-Feedback.toml")
    model = build(path)

    # reaction order (connection records): [r1, r2, r3, r4, r0, r5, r6]
    reactions = ["r1","r2","r3","r4","r0","r5","r6"]
    species   = ["A","B","C","D","E"]

    model.initial_condition_array = [1.0, 0.1, 0.1, 0.1, 0.1]   # A,B,C,D,E
    model.static_factors_array    = [1.0, 1.0, 1.0, 1.0]        # E1,E2,E3,E4
    model.α = [5.0, 5.0, 5.0, 5.0, 2.0, 1.0, 1.0]               # per-reaction rate constants

    # feedback via negative kinetic orders (as in the BSTModelKit paper example)
    G = model.G
    G[5, 1] = -0.5   # E (species 5) inhibits r1 (reaction 1)
    G[4, 2] = -0.5   # D (species 4) inhibits r2 (reaction 2)
    model.G = G

    Xss = steadystate(model; tspan=(0.0, 200.0))
    vss = reaction_fluxes(model, Xss)     # <-- replace with the actual API found in Step 2
    return (species=species, Xss=Xss, reactions=reactions, vss=vss)
end
```

- [ ] **Step 4: Replace `reaction_fluxes(model, Xss)`** with the actual BSTModelKit call found in Step 2. If the package exposes only species-level balances, compute per-reaction fluxes directly as `v_j = α_j * prod(Xss[i]^G[i,j] for i in dynamic species) * prod(static_factors[k]^G_static...)` using the kinetic records; encode this as a small local `reaction_fluxes` helper in the same file. Verify it satisfies mass balance: `r0 ≈ r1 ≈ r2`, `r3 ≈ r5`, `r4 ≈ r6`, `r3 + r4 ≈ r2`.

- [ ] **Step 5: Add a self-check at the bottom of the data-gen path**

```julia
let t = feedback_truth()
    A,B,C,D,E = t.Xss
    v = Dict(t.reactions .=> t.vss)
    @assert isapprox(v["r0"], v["r1"]; rtol=1e-3)
    @assert isapprox(v["r1"], v["r2"]; rtol=1e-3)
    @assert isapprox(v["r3"] + v["r4"], v["r2"]; rtol=1e-3)
    @assert isapprox(v["r3"], v["r5"]; rtol=1e-3)
    @assert isapprox(v["r4"], v["r6"]; rtol=1e-3)
    println("BST truth OK: T=", v["r1"], " split D/E=", v["r3"], "/", v["r4"])
end
```

- [ ] **Step 6: Run and verify**

Run: `cd code && julia --project=. feedback/run_feedback.jl`
Expected: prints `BST truth OK:` with a positive throughput `T` and a nonzero split between `r3` and `r4` (feedback keeps both branches active). All asserts pass.

- [ ] **Step 7: Commit**

```bash
git add code/feedback/Branched-Feedback.toml code/feedback/run_feedback.jl
git commit -m "Example 2a: BSTModelKit branched-feedback synthetic-data generator"
```

### Task 5: Example 2b — FBA layer and the gateway-closed-vs-open comparison

**Files:**
- Modify: `code/feedback/run_feedback.jl` (append the FBA layer + figure)
- Produces: `code/figs/feedback_gateway.pdf`, `code/data/feedback_fba.csv`

**Interfaces:**
- Consumes: `feedback_truth()` from Task 4.
- Produces: `feedback_fba(truth; gateway::Bool)` returning a flux vector aligned to `truth.reactions`; the figure and CSV used by §6 prose.

- [ ] **Step 1: Append the stoichiometry + FBA solver** to `code/feedback/run_feedback.jl`

```julia
# Stoichiometric matrix S (species x reactions), rows [A,B,C,D,E],
# cols [r1,r2,r3,r4,r0,r5,r6]:
const S_FEEDBACK = [
    -1  0  0  0  1  0  0    # A
     1 -1  0  0  0  0  0    # B
     0  1 -1 -1  0  0  0    # C
     0  0  1  0  0 -1  0    # D
     0  0  0  1  0  0 -1    # E
]

# Allosteric control function theta(P) in [0,1]: two-microstate partition
# function (free-active vs inhibitor-bound-inactive), reduces to a Hill term.
theta(P; K, n) = 1 / (1 + (P/K)^n)

function feedback_fba(truth; gateway::Bool)
    rxn = truth.reactions
    idx = Dict(rxn .=> eachindex(rxn))
    A,B,C,D,E = truth.Xss

    Vmax = 10.0                 # generous common capacity
    lb = fill(0.0, length(rxn)) # all irreversible here
    ub = fill(Vmax, length(rxn))

    if gateway
        # open the regulation gateway: cap the committed steps r1,r2 using the
        # measured end-product concentrations E,D through the control function,
        # and cap the branch reactions r3,r4 at their true capacities.
        ub[idx["r1"]] = Vmax * theta(E; K=0.5, n=2)
        ub[idx["r2"]] = Vmax * theta(D; K=0.5, n=2)
        ub[idx["r3"]] = truth.vss[idx["r3"]]
        ub[idx["r4"]] = truth.vss[idx["r4"]]
    end

    model = Model(HiGHS.Optimizer); set_silent(model)
    @variable(model, lb[i] <= v[i=1:length(rxn)] <= ub[i])
    @constraint(model, S_FEEDBACK * v .== 0)
    @objective(model, Max, v[idx["r5"]])      # maximise D-branch export
    optimize!(model)
    value.(v)
end
```

- [ ] **Step 2: Add the comparison + figure + assertions** at the bottom

```julia
let t = feedback_truth()
    v_closed = feedback_fba(t; gateway=false)
    v_open   = feedback_fba(t; gateway=true)

    df = DataFrame(reaction=t.reactions, truth=t.vss,
                   gateway_closed=v_closed, gateway_open=v_open)
    CSV.write(datapath("feedback_fba.csv"), df)

    # gateway-open must match the BST truth; gateway-closed must not
    @assert maximum(abs.(v_open .- t.vss)) < 1e-2 "gateway-open should match truth"
    @assert maximum(abs.(v_closed .- t.vss)) > 1e-1 "gateway-closed should diverge"

    fig = Figure(size=(720,320))
    ax  = Axis(fig[1,1], xticks=(1:length(t.reactions), t.reactions),
               ylabel="flux (AU)", title="Feedback recovered through the bound")
    w = 0.25
    barplot!(ax, (1:length(t.reactions)) .- w, t.vss;       width=w, label="BST truth")
    barplot!(ax, (1:length(t.reactions)),      v_closed;    width=w, label="FBA, gateway closed")
    barplot!(ax, (1:length(t.reactions)) .+ w, v_open;      width=w, label="FBA, gateway open")
    axislegend(ax; position=:lt)
    save(figpath("feedback_gateway.pdf"), fig)
    println("feedback FBA OK")
end
```

- [ ] **Step 3: Run and verify**

Run: `cd code && julia --project=. feedback/run_feedback.jl`
Expected: prints `feedback FBA OK`; both asserts pass; `figs/feedback_gateway.pdf` shows gateway-open bars matching truth while gateway-closed diverges (drives D-branch high, starves the E-branch).
If the asserts fail, tune `theta`'s `K,n` and the `Vmax` so the closed case over-drives `r5` and the open case matches; the `θ(E)`, `θ(D)` caps must reproduce `truth.vss` for `r1,r2`.

- [ ] **Step 4: Copy figure to chapter and commit**

```bash
cp code/figs/feedback_gateway.pdf chapter/figures/feedback_gateway.pdf
git add code/feedback/run_feedback.jl code/figs/feedback_gateway.pdf chapter/figures/feedback_gateway.pdf code/data/feedback_fba.csv
git commit -m "Example 2b: FBA layer recovers BST feedback through the bound"
```

---

## Phase 2 — Prose (written against the artifacts)

Each prose task: write the section to its content spec, `cd chapter && make`, confirm no undefined refs / missing citations in `Chapter.log`, run `style-check` + `review-section`, commit. Equations below are the required load-bearing content and must appear verbatim (LaTeX).

### Task 6: §2 Derivation — from mole balances to the constraint (`sections/derivation.tex`)

**Content spec:** `\section{From Mole Balances to the Flux Constraint}\label{sec:derivation}`. Draw from CHEME-5430 advanced-derivation notebook and L5c/L6a notes. Must contain, in order:
1. Open species mole balance → concentration balance: `\sum_{s}d_s C_{i,s}\dot V_s + \sum_j \sigma_{ij}\hat v_j V = \frac{d}{dt}(C_i V)`.
2. Specific units `V=B\bar V`; expand accumulation; steady intracellular pools + constant culture volume + no transport to reach `C_i\mu = \sum_j \sigma_{ij}\hat v_j`, i.e. `\mathbf{S}\hat{\mathbf v}=\mu\mathbf{x}` (number this equation; it is *the* honest constraint).
3. State plainly that dropping `\mu\mathbf{x}` gives the Palsson special case `\mathbf{S}\hat{\mathbf v}=\mathbf 0`, valid when growth is slow vs. metabolism and pools are dilute, and **failing** for fed-batch (`d\bar V/dt\neq0`) and cell-free (`\mu=0` but `dx/dt\neq0`). Forward-reference §7 capstones.
4. Exchange reactions: how hypothetical boundary reactions make a closed steady-state statement represent an open system. Cite Orth 2010, Bordbar 2014.
Full algebra goes in the appendix (Task 14); §2 shows the result and the interpretation. ~1.5 pp. Cite: Orth 2010, Bordbar 2014.

- [ ] Step 1: Write `sections/derivation.tex` to the spec above.
- [ ] Step 2: `cd chapter && make`; confirm build + label `sec:derivation` resolves.
- [ ] Step 3: Run `style-check` and `review-section` on the file; fix flags.
- [ ] Step 4: `git add chapter/sections/derivation.tex && git commit -m "Write section 2: keep-mu derivation"`

### Task 7: §3 The linear program and its geometry (`sections/linearprogram.tex`)

**Content spec:** `\section{The Linear Program and Its Geometry}\label{sec:lp}`. Must contain:
1. The FBA LP: `\max_{\hat v} \mathbf c^\top\hat{\mathbf v}` s.t. `\mathbf S\hat{\mathbf v}=\mathbf 0` (note: or `\mu\mathbf x`), `\boldsymbol{\ell}\le\hat{\mathbf v}\le\mathbf u`. Objective choice: biomass vs. product export.
2. Underdetermination: `r>m`; the null space `\mathcal N(\mathbf S)`.
3. SVD `\mathbf S=\mathbf U\boldsymbol\Sigma\mathbf V^\top`: right null space (zero-singular-value columns of `\mathbf V`) spans feasible flux modes; left null space spans conserved moieties. Genome-scale: `\mathrm{rank}(\mathbf S)\ll r`.
4. FVA in two sentences (min/max each `v_j` at fixed optimum) to characterize what constraints leave undetermined.
~1.5 pp. Cite: Orth 2010; Palsson/Bordbar 2014.

- [ ] Step 1: Write the file. Step 2: build + refs. Step 3: style-check + review-section. Step 4: commit `-m "Write section 3: LP and geometry"`.

### Task 8: §4 Bounds as gateways (`sections/gateways.tex`) — the heart

**Content spec:** `\section{Bounds as Gateways}\label{sec:gateways}`. Draw from L5c/L6a/L6c. Must contain, in order:
1. The general bound (number it): `-\delta_j\big[V^{\circ}_{\max,j}(e/e^{\circ})\theta_j(\cdot)f_j(\cdot)\big]\le\hat v_j\le V^{\circ}_{\max,j}(e/e^{\circ})\theta_j(\cdot)f_j(\cdot)`, and the framing sentence: each factor is a gateway for a distinct class of information.
2. Thermodynamic gateway `\delta_j\in\{0,1\}` from `\mathrm{sign}(\Delta G^{\circ}-\Delta G^{*})` or a `K_{eq}` cutoff (eQuilibrator). 
3. Kinetic gateway: `V^{\circ}_{\max,j}=k^{\circ}_{\mathrm{cat},j}e^{\circ}` (BRENDA); substrate saturation `f_j` (Michaelis–Menten form — this is where the old kinetics material is re-homed; include the MM rate law and its quasi-steady-state derivation pointer to the appendix).
4. Expression gateway `(e/e^{\circ})`: the mRNA/protein balances `\dot m_j=r_{X,j}u_j-(\theta_{m,j}+\mu)m_j+\lambda_j`, `\dot p_j=r_{L,j}w_j-(\theta_{p,j}+\mu)p_j`; steady state `m^{\star}_j=\frac{r_{X,j}u_j+\lambda_j}{\theta_{m,j}+\mu}`, `p^{\star}_j=\frac{r_{L,j}w_j}{\theta_{p,j}+\mu}`; `(e/e^{\circ})=p^{\star}_j/e^{\circ}` via GPR. Emphasize the kept `\mu` in the denominators (the signature).
5. Regulation gateway `\theta_j`: the partition-function control function `p_s=\frac{f_s e^{-\beta\epsilon_s}}{Z}`, `Z=\sum_s f_s e^{-\beta\epsilon_s}`, `\bar u=\sum_{s\in\mathcal A}p_s+\sum_{s\in\mathcal B}p_s=u+u^{\dagger}`, `\lambda\equiv r_X u^{\dagger}`; `f_s` Hill-type in effector concentration. State that this is the Adhikari effective biophysical (= Boltzmann) model — one object.
6. Simplified bounds model `-\delta_j V^{\circ}_{\max,j}\le\hat v_j\le V^{\circ}_{\max,j}` as the tractable baseline (all correction factors `\sim 1`).
7. One paragraph pointer: sequence-specific transcription/translation (Allen-Palsson demands; Vilkhovoy 2018 CFPS) as the mechanism behind the §7 capstone; not developed here.
~5 pp (the longest section). Cite: eQuilibrator (Beber 2022), BRENDA (Chang 2021), Adhikari 2020, Allen-Palsson 2003, Vilkhovoy 2018.

- [ ] Step 1: Write the file to the full spec. Step 2: build + refs. Step 3: style-check + review-section. Step 4: commit `-m "Write section 4: bounds as gateways"`.

### Task 9: §5 Example 1 — urea-cycle metabolism (`sections/example_urea.tex`)

**Content spec:** `\section{Worked Example: Urea-Cycle Metabolism}\label{sec:urea}`. This exercises the thermodynamic + kinetic gateways from databases. Reuse the technical content from the (now-deleted) old `fba.tex` urea material as raw source, but reframed: it is the *metabolism half* demonstrating database-informed gateways, not a standalone FBA intro. Must state: the network (5 enzymatic + 14 exchange, `\mathbf S\in\mathbb R^{18\times19}`); reversibility from eQuilibrator (`\Delta G^{\circ}` table), `V_{\max}` from BRENDA `k_{\mathrm{cat}}`; the objective (maximize urea export, `c_{b4}=-1`); the result (v1–v4 saturate at `0.0328` mmol/gDW/h, the argininosuccinate-lyase bottleneck; NOS branch carries zero flux). Reference `Figure~\ref{fig:urea}` (`\includegraphics{urea_fba.pdf}`). The `0.0328` value must match Task 3 output. ~2.5 pp. Cite: eQuilibrator, BRENDA, Orth 2010.

- [ ] Step 1: Write the file (recover equations/numbers from git history `git show HEAD~N:chapter/sections/fba.tex` if useful). Step 2: build; confirm `fig:urea` + `0.0328` present. Step 3: style-check + review-section. Step 4: commit `-m "Write section 5: urea-cycle example"`.

### Task 10: §6 Example 2 — branched-feedback toy (`sections/example_feedback.tex`)

**Content spec:** `\section{Worked Example: Regulation Through the Bound}\label{sec:feedback}`. This exercises the regulation gateway with synthetic data. Must contain:
1. The branched network `r0:\emptyset\to A`, `r1:A\to B`, `r2:B\to C`, `r3:C\to D`, `r4:C\to E`, `r5:D\to\emptyset`, `r6:E\to\emptyset`, with product feedback (E inhibits r1, D inhibits r2).
2. The two representations (the duality): the S-system power-law form `\frac{dX_i}{dt}=\alpha_i\prod_j X_j^{g_{ij}}-\beta_i\prod_j X_j^{h_{ij}}` encodes feedback as a negative kinetic order; FBA encodes the same fact as the `\theta` control function in the bound. Cite BSTModelKit (Vadhin-Varner 2026), Savageau.
3. The experiment: `BSTModelKit.jl` generates steady-state concentrations and fluxes (the synthetic truth); naive FBA (gateway closed) maximizing D-export overshoots and starves the E-branch; opening the gateway (bounds on `r1,r2` from `\theta` evaluated at the measured `E,D`; branch capacities on `r3,r4`) recovers the truth. Reference `Figure~\ref{fig:feedback}` (`feedback_gateway.pdf`).
4. Punchline sentence tying back to the thesis: regulation is real and, in the constraint-based world, it enters through the bound.
~3 pp. Cite: Vadhin-Varner 2026, Savageau 1976, Adhikari 2020.

- [ ] Step 1: Write the file; numbers/behavior must match Task 5 (`feedback_fba.csv`). Step 2: build; confirm `fig:feedback`. Step 3: style-check + review-section. Step 4: commit `-m "Write section 6: feedback toy example"`.

### Task 11: §7 Capstones — reviewed (`sections/capstones.tex`)

**Content spec:** `\section{Integration in Practice}\label{sec:capstones}`. Two reviewed subsections (prose only, one figure each optional — reuse published figures only if license permits; otherwise describe):
1. **Vilkhovoy 2023** — the dynamic integrated cell-free model: §5 (metabolism) and §6 (regulation) coupled and made dynamic, driven by real time-resolved intracellular data (63 metabolites + mRNA/protein/enzyme activity). Make the point: cell-free is where the kept term is honestly measurable because the intracellular compartment *is* the reactor. This is where `S\hat v=0` is not assumed.
2. **Wayman 2019** — bounds as engineering design levers: model-guided metabolic engineering for designer glycans; changing bounds = changing the achievable phenotype.
~2.5 pp. Cite: Vilkhovoy 2023, Wayman 2019, Adhikari 2020, Vilkhovoy 2018.

- [ ] Step 1: Write the file. Step 2: build + refs. Step 3: style-check + review-section. Step 4: commit `-m "Write section 7: reviewed capstones"`.

### Task 12: §1 Introduction (`sections/introduction.tex`)

**Content spec:** `\section{Introduction}\label{sec:intro}`. Written last-but-one so it can promise exactly what the chapter delivers. Must: frame flux as the phenotype and its estimation as the problem; state the two load-bearing ideas (bounds are the gateway; keep `\mu`); preview the spine (derivation → LP/geometry → gateways → two worked examples → capstones); set the tone (a rigorous, self-contained account, not a survey). No section self-references. ~1.5 pp. Cite: Orth 2010 only.

- [ ] Step 1: Write the file. Step 2: build. Step 3: style-check + review-section. Step 4: commit `-m "Write section 1: introduction"`.

### Task 13: §8 Outlook (`sections/outlook.tex`)

**Content spec:** `\section{Outlook}\label{sec:outlook}`. Brief (~0.75 pp). Where the framework goes: learning the control functions and `k_{\mathrm{cat}}` from data; scaling the gateway idea to genome scale; the dynamic/data-rich regime (cell-free) as the honest setting for the kept term. Non-crap, no hype. Cite nothing new.

- [ ] Step 1: Write. Step 2: build. Step 3: style-check + review-section. Step 4: commit `-m "Write section 8: outlook"`.

### Task 14: Appendix — full derivations (`sections/appendix.tex`)

**Content spec:** `\section{Derivations}\label{sec:appendix}`. Three worked derivations referenced from the body: (A) the full open-system mole balance → `\mathbf S\hat{\mathbf v}=\mu\mathbf x` (every algebra step from L5c/L6a); (B) SVD null-space structure and conserved moieties; (C) steady-state expression `m^{\star},p^{\star}` from the mRNA/protein balances. ~3 pp. Cite: as in body, no new refs.

- [ ] Step 1: Write. Step 2: build; confirm appendix labels referenced by §2/§3/§4 resolve. Step 3: style-check. Step 4: commit `-m "Write appendix: full derivations"`.

---

## Phase 3 — Front matter, references, integration

### Task 15: Abstract, title block, and References.bib

**Files:** Modify `chapter/Chapter.tex` (abstract), `chapter/References.bib`.

- [ ] **Step 1: Write the new abstract** into `Chapter.tex` (replace the Task-1 placeholder). ~180 words: FBA as an integrative framework; bounds as gateways for thermodynamic/kinetic/expression/regulatory information; the kept-`\mu` constraint `S\hat v=\mu x`; two worked examples (urea-cycle metabolism; a branched-feedback toy where opening the bound recovers regulation); reviewed cell-free and metabolic-engineering capstones. No em-dashes.
- [ ] **Step 2: Trim `References.bib`** to exactly the minimal set listed in the File Structure section; delete all deep-learning/motif/cybernetic/CHO entries; add missing entries (BSTModelKit arXiv:2603.19115; Vilkhovoy 2018; Adhikari 2020 PMID 33324619; Vilkhovoy 2023 bioRxiv 2023.02.10.528035; Wayman 2019; Allen-Palsson 2003; Savageau 1976 + 1987).
- [ ] **Step 3: Build and check citations**

Run: `cd chapter && make`
Expected: no "Citation undefined" and no "There were undefined references" in `Chapter.log`; every `\cite` resolves; no uncited bib bloat (optional: check with `checkcites` if available).

- [ ] **Step 4: Commit** `git add chapter/Chapter.tex chapter/References.bib && git commit -m "Add abstract; trim references to minimal cited set"`

### Task 16: Whole-chapter integration pass

**Files:** any prose file needing a cross-reference fix.

- [ ] **Step 1: Clean build from scratch**

Run: `cd chapter && make clean && make`
Expected: `Chapter.pdf` builds with zero warnings for undefined references/citations in `Chapter.log`.

- [ ] **Step 2: Cross-reference audit** — grep every `\ref`/`\eqref`/`\label` and confirm each `\ref` has a target and each figure (`fig:urea`, `fig:feedback`) is referenced in text.

Run: `cd chapter && grep -rn '\\ref{' sections/ && grep -rn '\\label{' sections/`

- [ ] **Step 3: Reproducibility check** — regenerate both figures end-to-end

Run: `cd code && julia --project=. fba/run_fba.jl && julia --project=. feedback/run_feedback.jl`
Expected: both print their OK/values; figures regenerate; copy any changed figure into `chapter/figures/` and rebuild.

- [ ] **Step 4: Thesis-consistency read** — confirm "bounds are the gateway" and the kept-`\mu` statement (`S\hat v=\mu x`) each appear in the intro, §2, and §4, and that no section reintroduces `S\hat v=0` as *the* constraint without the special-case caveat.

- [ ] **Step 5: Run the `audit-magic-numbers` skill** on §5 and §6 to confirm every numeric result traces to a figure/table/CSV.

- [ ] **Step 6: Final commit**

```bash
git add -A
git commit -m "Integration pass: clean build, cross-refs, reproducibility verified"
```

---

## Self-Review (completed against the spec)

- **Spec coverage:** §1 Task 12; §2 Task 6; §3 Task 7; §4 Task 8; §5 Task 9; §6 Task 10; §7 Task 11; §8 Task 13; appendix Task 14; abstract/refs Task 15; both built examples Tasks 3–5; teardown/re-home Tasks 1–2; keep-`μ` present in Tasks 6/8/14/16; partition-function regulation in Task 8; BSTModelKit reuse in Tasks 2/4/5. No spec section is unmapped.
- **Placeholder scan:** the only deliberately deferred item is the BSTModelKit rate-extraction API (Task 4 Steps 2/4), which is a read-the-source instruction with a concrete fallback formula, not a "TODO".
- **Type consistency:** `feedback_truth()` fields (`species, Xss, reactions, vss`) and `feedback_fba(truth; gateway)` signature are used identically in Tasks 4 and 5; reaction order `[r1,r2,r3,r4,r0,r5,r6]` and `S_FEEDBACK` columns are aligned; `0.0328` urea value shared by Task 3 and Task 9.
