# Phase 1 Correctness Fixes Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the urea-cycle seconds-to-hours unit bug, align the gateways.tex protein
balance with the feedback implementation, close the feedback example's circularity with a
bounded/independent regulatory gateway, and harden the Julia code around all of it — so the
numbers Phase 2/3 will write prose around are actually correct.

**Architecture:** No new modules. All changes are in-place edits to existing files in
`code/fba/`, `code/feedback/`, and two `chapter/sections/*.tex` files. Every code change is
verified by running the actual Julia scripts (this project has no mocking layer — FBA/UQ
correctness is checked by solving the real LP and reading real output).

**Tech Stack:** Julia 1.12 (JuMP + HiGHS for the LPs, BSTModelKit for the feedback truth model,
CairoMakie for figures, CSV/DataFrames for data), plain `article`-class LaTeX built via
`chapter/Makefile`.

## Global Constraints

- Spec: `docs/superpowers/specs/2026-07-14-audit-phase1-correctness-fixes-design.md`.
- `chapter/sections/derivation.tex` is explicitly untouched — no task in this plan modifies it.
- Julia invocation is always `julia --project=code <script path>` from the repo root (per
  `README.md`); do not `cd` into `code/` first, or the relative `@__DIR__`-based paths inside
  `Include.jl` break.
- After any change to `code/fba/urea_cycle.jl` or `code/feedback/dual_feedback.jl`, the scripts
  that consume them (`run_fba.jl`, `urea_cycle_uq.jl`, `run_feedback.jl`) must be re-run so the
  CSVs in `code/data/` and figures in `code/figs/` stay in sync with the code before any LaTeX
  prose is updated to match.
- Figures must be copied from `code/figs/` to `chapter/figures/` before `make` picks them up
  (per `README.md`'s build step) — `chapter/figures/` is not a symlink.
- All dollar-quantity numbers written into `.tex` files in this plan were verified by actually
  running the fixed code in a scratch copy during planning; they are not hand-derived estimates.

---

## Task 1: Fix the urea-cycle seconds-to-hours unit conversion

**Files:**
- Modify: `code/fba/urea_cycle.jl:148`
- Test: `code/fba/test_urea_cycle.jl`

**Interfaces:**
- Consumes: nothing new.
- Produces: `urea_cycle_model()` now returns `lb`/`ub` 3600x larger for the 5 enzymatic
  reactions (`v1..v5`); exchange bounds (`b1..b14`, ±1000) are unchanged. Downstream tasks
  (Task 4, Task 7) depend on this new scale.

- [ ] **Step 1: Update the bounds/flux assertions in the existing test to the corrected values**

Edit `code/fba/test_urea_cycle.jl`, replacing the whole file with:

```julia
include(joinpath(@__DIR__, "..", "Include.jl"))
include(joinpath(@__DIR__, "urea_cycle.jl"))

# Nominal bounds must match kcat (s^-1) * e0 * 3600 (s -> h conversion).
expected_lb = vcat([-360.0, -118.08,    0.0,     0.0,   0.0], fill(-1000.0, 14))
expected_ub = vcat([ 360.0,  118.08, 6840.0, 14760.0, 360.0], fill( 1000.0, 14))

m = urea_cycle_model()
# isapprox, not ==: Vmax = kcat*e0*3600 (e.g. 3.28*0.01*3600) differs from the literal 118.08
# only by floating-point round-off.
@assert all(isapprox.(m.lb, expected_lb; atol=1e-9)) "lb mismatch: $(m.lb)"
@assert all(isapprox.(m.ub, expected_ub; atol=1e-9)) "ub mismatch: $(m.ub)"

v = solve_flux(m)
@assert v !== nothing "nominal model failed to solve"
b4 = findfirst(==("b4"), m.reactions)
@assert isapprox(v[b4], 118.08; atol=1e-6) "urea export $(v[b4]) != 118.08 (secretion-positive)"
@assert maximum(abs.(m.S * v)) < 1e-6 "Sv=0 violated"

# f defaults to ones -> passing the nominal Park f leaves v2 the bottleneck (fluxes unchanged)
f_park = ones(5); f_park[1] = 0.923; f_park[3] = 0.142; f_park[4] = 0.154; f_park[5] = 0.986
vf = solve_flux(urea_cycle_model(; f = f_park))
@assert isapprox(vf[b4], 118.08; atol=1e-6) "Park-f urea export $(vf[b4]) != 118.08 (v2 should still bind)"

# Dimensional-conversion check: Vmax must equal kcat*e0*3600 exactly, not kcat*e0.
kcat_test = [7.0, 7.0, 7.0, 7.0, 7.0]
e0_test   = 0.02
m2 = urea_cycle_model(; kcat = kcat_test, e0 = e0_test)
@assert isapprox(m2.ub[1], 7.0 * 0.02 * 3600.0; atol=1e-9) "Vmax not converted s^-1 -> h^-1: got $(m2.ub[1])"
@assert !isapprox(m2.ub[1], 7.0 * 0.02; atol=1e-9) "Vmax still in per-second units"

println("test_urea_cycle: all checks passed")
```

- [ ] **Step 2: Run the test to confirm it fails against the current (unfixed) code**

Run: `julia --project=code code/fba/test_urea_cycle.jl`
Expected: `AssertionError: lb mismatch: ...` (the current code still returns the un-converted,
1/3600-scale bounds).

- [ ] **Step 3: Apply the fix**

In `code/fba/urea_cycle.jl`, change line 148 from:

```julia
    Vmax  = kcat .* e0
```

to:

```julia
    Vmax  = kcat .* e0 .* 3600.0   # kcat is s^-1 (BRENDA); convert to h^-1 to match e0 (mmol/gDW)
```

- [ ] **Step 4: Run the test to confirm it passes**

Run: `julia --project=code code/fba/test_urea_cycle.jl`
Expected: `test_urea_cycle: all checks passed`

- [ ] **Step 5: Commit**

```bash
git add code/fba/urea_cycle.jl code/fba/test_urea_cycle.jl
git commit -m "Fix urea-cycle kcat seconds-to-hours conversion (3600x)"
```

---

## Task 2: Harden `urea_cycle.jl` — input validation and explicit solve-failure handling

**Files:**
- Modify: `code/fba/urea_cycle.jl`
- Test: `code/fba/test_urea_cycle.jl`

**Interfaces:**
- Consumes: `urea_cycle_model`, `solve_flux`, `solve_fba` from Task 1.
- Produces: `urea_cycle_model` now raises `ArgumentError` on malformed input instead of
  silently building a nonsensical model; `solve_fba` now raises `ErrorException` on an
  infeasible/unbounded solve instead of returning a `DataFrame` with a `nothing` flux column.

- [ ] **Step 1: Add the failing validation tests**

Append to `code/fba/test_urea_cycle.jl` (before the final `println`):

```julia
# ---- Task 2: input validation and explicit solve-failure handling ---------- #
using Test

@test_throws ArgumentError urea_cycle_model(; kcat = [1.0, 2.0])          # wrong length
@test_throws ArgumentError urea_cycle_model(; dG = fill(NaN, 5))          # non-finite
@test_throws ArgumentError urea_cycle_model(; f = [-0.1, 1.0, 1.0, 1.0, 1.0])  # negative saturation

infeasible = urea_cycle_model(; kcat = zeros(5))  # zero capacity everywhere -> b4 forced to 0, still feasible;
infeasible = (; infeasible..., lb = infeasible.lb .+ 1.0, ub = infeasible.ub)  # lb > ub -> infeasible
@test_throws ErrorException solve_fba(infeasible)

println("test_urea_cycle (Task 2): validation and solve-failure checks passed")
```

- [ ] **Step 2: Run to verify these new checks fail**

Run: `julia --project=code code/fba/test_urea_cycle.jl`
Expected: fails with `UndefVarError` or `MethodError` around the `@test_throws` calls (Test
package usage is fine — `using Test` is already a resolvable dependency via `Pkg.instantiate`
since JuMP pulls it in transitively — actual expected failure is that `ArgumentError`/
`ErrorException` are not yet thrown because the validation and failure-handling don't exist
yet, so `@test_throws` reports "no exception was thrown" or the code errors during model
construction with an unrelated `BoundsError`/`DimensionMismatch` instead of the intended
`ArgumentError`).

- [ ] **Step 3: Add validation inside `urea_cycle_model` and fix `solve_fba`**

In `code/fba/urea_cycle.jl`, right after the `metabolites`/`reactions`/`S` block and before the
`delta = ...` line (i.e. right before the "Flux bounds, built from parameters" comment block),
insert:

```julia
    # ------------------------------------------------------------------ #
    # Input validation
    # ------------------------------------------------------------------ #
    length(kcat) == 5 || throw(ArgumentError("kcat must have length 5, got $(length(kcat))"))
    length(dG)   == 5 || throw(ArgumentError("dG must have length 5, got $(length(dG))"))
    length(f)    == 5 || throw(ArgumentError("f must have length 5, got $(length(f))"))
    all(isfinite, kcat)   || throw(ArgumentError("kcat must be finite: $kcat"))
    all(isfinite, dG)     || throw(ArgumentError("dG must be finite: $dG"))
    all(isfinite, f)      || throw(ArgumentError("f must be finite: $f"))
    isfinite(e0)          || throw(ArgumentError("e0 must be finite: $e0"))
    all(f .>= 0)          || throw(ArgumentError("saturation factors f must be nonnegative: $f"))
```

Then, at the end of the function, change the final `@assert` block (currently checking sizes
of `S`, `lb`, `ub`, `c`) to also check `lb .<= ub`:

```julia
    # sanity checks
    @assert size(S) == (length(metabolites), length(reactions))
    @assert length(lb) == length(reactions)
    @assert length(ub) == length(reactions)
    @assert length(c)  == length(reactions)
    all(lb .<= ub) || throw(ArgumentError("lower bounds must not exceed upper bounds: lb=$lb, ub=$ub"))
```

Then change `solve_fba` from:

```julia
solve_fba(m) = DataFrame(reaction = m.reactions, flux = solve_flux(m))
```

to:

```julia
function solve_fba(m)
    v = solve_flux(m)
    v === nothing && error("solve_fba: linear program did not solve to optimality (infeasible or unbounded)")
    DataFrame(reaction = m.reactions, flux = v)
end
```

- [ ] **Step 4: Run the test to confirm all checks pass**

Run: `julia --project=code code/fba/test_urea_cycle.jl`
Expected: `test_urea_cycle (Task 2): validation and solve-failure checks passed` followed by
`test_urea_cycle: all checks passed`.

- [ ] **Step 5: Commit**

```bash
git add code/fba/urea_cycle.jl code/fba/test_urea_cycle.jl
git commit -m "Harden urea_cycle_model: input validation and explicit solve-failure handling"
```

---

## Task 3: Add an FVA routine and uniqueness test to `urea_cycle.jl`

**Files:**
- Modify: `code/fba/urea_cycle.jl`
- Test: `code/fba/test_urea_cycle.jl`

**Interfaces:**
- Consumes: `urea_cycle_model`, `solve_flux` from Task 1/2.
- Produces: `fva(m; tol=1e-6) -> DataFrame` with columns `reaction, vmin, vmax`. Not consumed by
  other tasks in this plan, but Phase 2 will use it to write up the FVA result in the chapter
  prose (per the design spec).

- [ ] **Step 1: Write the failing test**

Append to `code/fba/test_urea_cycle.jl` (before the final `println`):

```julia
# ---- Task 3: FVA confirms the nominal optimum is unique ------------------- #
m0  = urea_cycle_model()
rng = fva(m0)
@assert nrow(rng) == length(m0.reactions)
@assert maximum(abs.(rng.vmax .- rng.vmin)) < 1e-4 "nominal optimum is not unique: $(rng)"

println("test_urea_cycle (Task 3): FVA uniqueness check passed")
```

- [ ] **Step 2: Run to verify it fails**

Run: `julia --project=code code/fba/test_urea_cycle.jl`
Expected: `UndefVarError: fva not defined`

- [ ] **Step 3: Implement `fva`**

Append to `code/fba/urea_cycle.jl` (after `solve_fba`):

```julia
"""
    fva(m; tol=1e-6) -> DataFrame

Flux variability analysis: for each reaction, minimize and maximize its flux subject to
`Sv=0`, `lb<=v<=ub`, and the objective held within `tol` of its optimal value. Returns a
DataFrame with columns `reaction, vmin, vmax`.
"""
function fva(m; tol=1e-6)
    vopt = solve_flux(m)
    vopt === nothing && error("fva: nominal model did not solve to optimality")
    zopt = sum(m.c[i] * vopt[i] for i in eachindex(m.c))
    n = length(m.reactions)
    vmin = zeros(n); vmax = zeros(n)
    for i in 1:n
        for sense in (:min, :max)
            model = Model(HiGHS.Optimizer); set_silent(model)
            @variable(model, m.lb[k] <= v[k=1:n] <= m.ub[k])
            @constraint(model, m.S * v .== 0)
            @constraint(model, sum(m.c[k] * v[k] for k in 1:n) >= zopt - tol)
            if sense == :min
                @objective(model, Min, v[i])
            else
                @objective(model, Max, v[i])
            end
            optimize!(model)
            val = is_solved_and_feasible(model) ? value(v[i]) : NaN
            sense == :min ? (vmin[i] = val) : (vmax[i] = val)
        end
    end
    DataFrame(reaction = m.reactions, vmin = vmin, vmax = vmax)
end
```

- [ ] **Step 4: Run the test to confirm it passes**

Run: `julia --project=code code/fba/test_urea_cycle.jl`
Expected: `test_urea_cycle (Task 3): FVA uniqueness check passed` followed by
`test_urea_cycle: all checks passed`.

- [ ] **Step 5: Commit**

```bash
git add code/fba/urea_cycle.jl code/fba/test_urea_cycle.jl
git commit -m "Add FVA routine and nominal-uniqueness test to urea_cycle.jl"
```

---

## Task 4: Regenerate the urea nominal solution and update `example_urea.tex`

**Files:**
- Run (not modified): `code/fba/run_fba.jl`
- Modify: `chapter/sections/example_urea.tex`
- Generated: `code/data/urea_fba_solution.csv`, `code/figs/urea_fba.pdf` (regenerated by the
  script), copied to `chapter/figures/urea_fba.pdf`

**Interfaces:**
- Consumes: the Task 1 unit fix.
- Produces: `chapter/sections/example_urea.tex` with corrected `V°max`/flux numbers,
  consistent with the regenerated CSV/figure.

- [ ] **Step 1: Regenerate the nominal solution CSV and figure**

Run: `julia --project=code code/fba/run_fba.jl`
Expected output ends with: `objective_flux=118.07999999999998` and no assertion error (the
script's own `Sv=0` check passes).

Verify the CSV:

Run: `cat code/data/urea_fba_solution.csv`
Expected:
```
reaction,flux
v1,118.07999999999998
v2,118.07999999999998
v3,118.07999999999998
v4,118.07999999999998
v5,-0.0
b1,-118.07999999999998
b2,-118.07999999999998
b3,118.07999999999998
b4,118.07999999999998
b5,-118.07999999999998
b6,118.07999999999998
b7,118.07999999999998
b8,118.07999999999998
b9,-0.0
b10,-0.0
b11,-0.0
b12,-0.0
b13,-0.0
b14,-118.07999999999998
```

- [ ] **Step 2: Copy the regenerated figure into the chapter**

Run: `cp code/figs/urea_fba.pdf chapter/figures/urea_fba.pdf`

- [ ] **Step 3: Update the enzyme-capacity numbers in `example_urea.tex`**

In `chapter/sections/example_urea.tex`, replace the sentence (around line 67-79):

```latex
Argininosuccinate synthetase and nitric
oxide synthase share $k^{\circ}_{\mathrm{cat}}=10.0$ s$^{-1}$, giving
$V^{\circ}_{\max,v_1}=V^{\circ}_{\max,v_5}=0.100$ mmol\,gDW$^{-1}$\,h$^{-1}$;
arginase I has $k^{\circ}_{\mathrm{cat}}=190.0$ s$^{-1}$, giving
$V^{\circ}_{\max,v_3}=1.9$ mmol\,gDW$^{-1}$\,h$^{-1}$; and ornithine
transcarbamylase has $k^{\circ}_{\mathrm{cat}}=410.0$ s$^{-1}$, giving
$V^{\circ}_{\max,v_4}=4.1$ mmol\,gDW$^{-1}$\,h$^{-1}$. Argininosuccinate
lyase is the outlier: $k^{\circ}_{\mathrm{cat},v_2}=3.28$ s$^{-1}$, an
order of magnitude below the other four turnover numbers, gives
$V^{\circ}_{\max,v_2}=0.0328$ mmol\,gDW$^{-1}$\,h$^{-1}$, the smallest
of the five capacities by more than threefold.
```

with:

```latex
Argininosuccinate synthetase and nitric
oxide synthase share $k^{\circ}_{\mathrm{cat}}=10.0$ s$^{-1}$, giving, after
converting the per-second turnover number to the per-hour basis the rest of the
chapter uses,
$V^{\circ}_{\max,v_1}=V^{\circ}_{\max,v_5}=360$ mmol\,gDW$^{-1}$\,h$^{-1}$;
arginase I has $k^{\circ}_{\mathrm{cat}}=190.0$ s$^{-1}$, giving
$V^{\circ}_{\max,v_3}=6840$ mmol\,gDW$^{-1}$\,h$^{-1}$; and ornithine
transcarbamylase has $k^{\circ}_{\mathrm{cat}}=410.0$ s$^{-1}$, giving
$V^{\circ}_{\max,v_4}=14760$ mmol\,gDW$^{-1}$\,h$^{-1}$. Argininosuccinate
lyase is the outlier: $k^{\circ}_{\mathrm{cat},v_2}=3.28$ s$^{-1}$, an
order of magnitude below the other four turnover numbers, gives
$V^{\circ}_{\max,v_2}=118.08$ mmol\,gDW$^{-1}$\,h$^{-1}$, the smallest
of the five capacities by more than threefold.
```

- [ ] **Step 4: Update the nominal-solution paragraph**

Replace (around line 106-118):

```latex
Solving Equation~\eqref{eq:fba-lp} by linear programming
\cite{Orth2010} with these bounds and this objective gave a single,
sharply determined answer. The entire enzymatic backbone, $v_1$
through $v_4$, saturated at $0.0328$ mmol\,gDW$^{-1}$\,h$^{-1}$, the
capacity imposed by argininosuccinate lyase, and the nitric oxide
synthase branch $v_5$ carried zero flux: with arginine supply capped
by the same $v_2$ bottleneck that limits the rest of the cycle, none
of it was left over to spend on the competing NOS reaction once urea
export is the sole objective. The optimal urea export was
$0.0328$ mmol\,gDW$^{-1}$\,h$^{-1}$, attained at
$\hat v_{b_4}=0.0328$, and every exchange flux carrying carbon or
nitrogen into or out of the cycle, carbamoyl phosphate and aspartate
on the uptake side, urea and fumarate on the secretion side, scaled to
match it, giving the full optimal flux distribution
(right panel, Fig.~\ref{fig:urea}).
```

with:

```latex
Solving Equation~\eqref{eq:fba-lp} by linear programming
\cite{Orth2010} with these bounds and this objective gave a single,
sharply determined answer. The entire enzymatic backbone, $v_1$
through $v_4$, saturated at $118.08$ mmol\,gDW$^{-1}$\,h$^{-1}$, the
capacity imposed by argininosuccinate lyase, and the nitric oxide
synthase branch $v_5$ carried zero flux: with arginine supply capped
by the same $v_2$ bottleneck that limits the rest of the cycle, none
of it was left over to spend on the competing NOS reaction once urea
export is the sole objective. The optimal urea export was
$118.08$ mmol\,gDW$^{-1}$\,h$^{-1}$, attained at
$\hat v_{b_4}=118.08$, and every exchange flux carrying carbon or
nitrogen into or out of the cycle, carbamoyl phosphate and aspartate
on the uptake side, urea and fumarate on the secretion side, scaled to
match it, giving the full optimal flux distribution
(right panel, Fig.~\ref{fig:urea}).
```

- [ ] **Step 5: Update the Figure 1 caption**

Replace (around line 185-191):

```latex
    $v_1$ through $v_4$ saturates at $0.0328$, the capacity of
    argininosuccinate lyase ($v_2$), while the nitric oxide synthase
    branch $v_5$ carries zero flux. Exchange fluxes $b_1$ through
    $b_{14}$ follow the secretion-positive convention, so positive bars
    denote net secretion and negative bars net uptake; $b_4$ (urea)
    shows the maximized export flux of $0.0328$ mmol\,gDW$^{-1}$\,h$^{-1}$.}
```

with:

```latex
    $v_1$ through $v_4$ saturates at $118.08$, the capacity of
    argininosuccinate lyase ($v_2$), while the nitric oxide synthase
    branch $v_5$ carries zero flux. Exchange fluxes $b_1$ through
    $b_{14}$ follow the secretion-positive convention, so positive bars
    denote net secretion and negative bars net uptake; $b_4$ (urea)
    shows the maximized export flux of $118.08$ mmol\,gDW$^{-1}$\,h$^{-1}$.}
```

- [ ] **Step 6: Commit**

```bash
git add chapter/sections/example_urea.tex chapter/figures/urea_fba.pdf code/data/urea_fba_solution.csv
git commit -m "Regenerate urea nominal solution and update example_urea.tex to corrected units"
```

---

## Task 5: Load `park_saturation.csv` properly in `urea_cycle_uq.jl`

**Files:**
- Modify: `code/fba/urea_cycle_uq.jl`

**Interfaces:**
- Consumes: `code/data/park_saturation.csv` (already exists, unchanged).
- Produces: `SAT_NOMINAL` is now computed by `load_park_saturation()` rather than hardcoded;
  its *values* are unchanged (verified below), so this task alone does not change any UQ
  output — Task 7 (unit fix propagation) is what changes the numbers.

- [ ] **Step 1: Replace the hardcoded dict with a loader function**

In `code/fba/urea_cycle_uq.jl`, replace:

```julia
# Nominal saturation inputs (conc_M, Km_M) per reaction, reduced from
# code/data/park_saturation.csv by the aggregation rule (least-saturated
# substrate; prefer Homo sapiens row, else geometric mean across organisms).
# v2 has no forward-substrate data -> f_2 = 1 in Config A.
const SAT_NOMINAL = Dict(
    1 => (4.673e-3, 3.923e-4),  # v1 ATP (Mus musculus)
    3 => (2.555e-4, 1.546e-3),  # v3 arginine (Homo sapiens)
    4 => (2.129e-4, 1.166e-3),  # v4 ornithine (geometric mean, yeast + E. coli)
    5 => (2.555e-4, 3.497e-6),  # v5 arginine (Mus musculus)
)
```

with:

```julia
"""
    load_park_saturation(path) -> Dict{Int,Tuple{Float64,Float64}}

Load `code/data/park_saturation.csv` and reduce to one (conc_M, km_M) pair per enzymatic
reaction index (1-based, matching KCAT0/DG0 ordering):
  1. group rows by (reaction, substrate);
  2. within a substrate, prefer the Homo sapiens row, else the geometric mean of conc_M and
     km_M across organism rows;
  3. across substrates for the same reaction, keep the least-saturated substrate (smallest
     conc_M / (km_M + conc_M)).
Reaction v2 has no forward-substrate data in the source (only reverse-direction products) and
is excluded; f_2 defaults to 1 unless imputed (see `sample_f`).
"""
function load_park_saturation(path=joinpath(_PATH_TO_DATA, "park_saturation.csv"))
    df = CSV.read(path, DataFrame; comment="#")
    reaction_index = Dict("v1" => 1, "v3" => 3, "v4" => 4, "v5" => 5)
    result = Dict{Int,Tuple{Float64,Float64}}()
    for (rxn, ridx) in reaction_index
        sub = df[df.reaction .== rxn, :]
        best = nothing
        for s in unique(sub.substrate)
            rows = sub[sub.substrate .== s, :]
            hs = rows[rows.organism .== "Homo sapiens", :]
            conc, km = if nrow(hs) > 0
                hs.conc_M[1], hs.km_M[1]
            else
                exp(mean(log.(rows.conc_M))), exp(mean(log.(rows.km_M)))
            end
            f = conc / (km + conc)
            if best === nothing || f < best[3]
                best = (conc, km, f)
            end
        end
        result[ridx] = (best[1], best[2])
    end
    return result
end

# Nominal saturation inputs (conc_M, Km_M) per reaction; see `load_park_saturation` docstring
# for the aggregation rule.
const SAT_NOMINAL = load_park_saturation()
```

- [ ] **Step 2: Run the UQ script and confirm the aggregated values are unchanged**

Run:
```bash
julia --project=code -e '
include("code/Include.jl")
include("code/fba/urea_cycle.jl")
include("code/fba/urea_cycle_uq.jl")
'
```
Wait — do not actually run this yet; the script has other assertions tied to Task 1's unit
fix that are handled in Task 7. For this step, just confirm the loader reproduces the prior
values in isolation:

Run:
```bash
julia --project=code -e '
include("code/Include.jl")
include("code/fba/urea_cycle_uq.jl")
d = load_park_saturation()
for k in sort(collect(keys(d))); println(k, " => ", d[k]); end
'
```
Expected:
```
1 => (0.0046729612916725, 0.00039233315423419906)
3 => (0.000255461883107884, 0.00154608643830104)
4 => (0.00021293944853313844, 0.0011661903789690598)
5 => (0.00025546188310788376, 3.4969415984039393e-6)
```
These match the previously hardcoded `(4.673e-3, 3.923e-4)`, `(2.555e-4, 1.546e-3)`,
`(2.129e-4, 1.166e-3)`, `(2.555e-4, 3.497e-6)` to the precision originally quoted.

- [ ] **Step 3: Commit**

```bash
git add code/fba/urea_cycle_uq.jl
git commit -m "Load park_saturation.csv instead of hardcoding the reduced values"
```

---

## Task 6: Record UQ termination status in `urea_cycle_uq.jl`

**Files:**
- Modify: `code/fba/urea_cycle_uq.jl`

**Interfaces:**
- Consumes: `run_ensemble` (existing function in this file).
- Produces: `run_ensemble` now returns `(F, m0, n_failed)` instead of `(F, m0)`; both call
  sites in this file are updated in this task, and Task 7 writes the failure counts into the
  printed summary.

- [ ] **Step 1: Update `run_ensemble` to track and return the failure count**

Replace:

```julia
"Run N draws; return an (kept x nreactions) flux matrix."
function run_ensemble(seed; impute_f2::Bool=false, sample_saturation::Bool=true)
    rng = MersenneTwister(seed)
    m0  = urea_cycle_model()
    rows = Vector{Vector{Float64}}()
    for _ in 1:N
        v = sample_flux(rng; impute_f2 = impute_f2, sample_saturation = sample_saturation)
        v === nothing && continue
        push!(rows, v)
    end
    return reduce(vcat, (r' for r in rows)), m0
end
```

with:

```julia
"Run N draws; return an (kept x nreactions) flux matrix, the nominal model, and the count of
draws that failed to solve (infeasible/unbounded), which are dropped rather than kept."
function run_ensemble(seed; impute_f2::Bool=false, sample_saturation::Bool=true)
    rng = MersenneTwister(seed)
    m0  = urea_cycle_model()
    rows = Vector{Vector{Float64}}()
    n_failed = 0
    for _ in 1:N
        v = sample_flux(rng; impute_f2 = impute_f2, sample_saturation = sample_saturation)
        if v === nothing
            n_failed += 1
            continue
        end
        push!(rows, v)
    end
    return reduce(vcat, (r' for r in rows)), m0, n_failed
end
```

- [ ] **Step 2: Update the two call sites in this file**

Replace:

```julia
# ---- Config A ----
F, m0 = run_ensemble(SEED; impute_f2 = false)
kept = size(F, 1)
@assert kept > 0.99 * N "too many infeasible draws: kept $kept of $N"
```

with:

```julia
# ---- Config A ----
F, m0, nfail_A = run_ensemble(SEED; impute_f2 = false)
kept = size(F, 1)
@assert kept > 0.99 * N "too many infeasible draws: kept $kept of $N"
println("configA_failed=", nfail_A)
```

Replace:

```julia
Fcap, _ = run_ensemble(SEED; sample_saturation=false)
cap_kept = size(Fcap, 1)
```

with:

```julia
Fcap, _, nfail_cap = run_ensemble(SEED; sample_saturation=false)
cap_kept = size(Fcap, 1)
println("configCap_failed=", nfail_cap)
```

- [ ] **Step 3: Run to confirm the script still executes end-to-end**

Run: `julia --project=code code/fba/urea_cycle_uq.jl`
Expected: this will still fail at the pre-existing `configA_urea_export` assertion (still
checking the un-converted `0.0328` target from before Task 1's fix propagates here) — that is
expected and is fixed in Task 7. Confirm the failure is specifically that assertion (not a
`MethodError`/`UndefVarError` from this task's changes) by checking the printed lines above
the error include `configA_failed=` and `configCap_failed=` with small integer values.

- [ ] **Step 4: Commit**

```bash
git add code/fba/urea_cycle_uq.jl
git commit -m "Record UQ ensemble failure counts instead of silently dropping them"
```

---

## Task 7: Fix UQ assertions, regenerate outputs, and update `example_urea.tex` UQ prose

**Files:**
- Modify: `code/fba/urea_cycle_uq.jl`
- Modify: `chapter/sections/example_urea.tex`
- Generated: `code/data/urea_fba_uq.csv`, `code/data/urea_uq_sensitivity.csv`,
  `code/figs/urea_fba.pdf` (overwritten again with the whisker), `code/figs/urea_saturation.pdf`

**Interfaces:**
- Consumes: Tasks 1, 5, 6 (unit fix, CSV loader, failure tracking must all be in place first).
- Produces: final, corrected `code/data/urea_fba_uq.csv` / `urea_uq_sensitivity.csv` and
  matching chapter prose. This is the last urea-side task before the full rebuild in Task 12.

- [ ] **Step 1: Update the two hardcoded assertion targets in `urea_cycle_uq.jl`**

Replace:

```julia
@assert isapprox(quantile(ua, 0.5), 0.0328; rtol=0.15) "Config A urea median should stay near 0.0328"
```

with:

```julia
@assert isapprox(quantile(ua, 0.5), 118.08; rtol=0.15) "Config A urea median should stay near 118.08"
```

(The other assertion, `@assert quantile(ub_export, 0.5) < quantile(ua, 0.5) "Config B should
lower the urea median"`, is a relative comparison and needs no numeric change.)

- [ ] **Step 2: Run the UQ script and capture the real output**

Run: `julia --project=code code/fba/urea_cycle_uq.jl`
Expected (Tasks 5 and 6 do not add or remove any `rand()` calls, so the RNG stream — and
therefore these values — are bit-identical to before those tasks, verified during planning):
```
configA_kept=10000
configA_failed=0
configA_urea_export mean=156.04333428521818 sd=157.22305454830027 median=105.05381342754373 ci=[17.060725222471447,613.5765955661723]
configA_v5 mean=0.023959880968959464
configCap_kept=10000
configCap_failed=0
configCap_urea_export mean=167.33961130464067 sd=169.8734191808569 median=112.91950708196956 ci=[17.596821995102992,672.00858977199]
configB_f2 median=0.5021191499901032
configB_urea_export median=51.26799246728291 ci=[4.348771770597835,397.81163870302527]
```
If the printed numbers differ from these beyond floating-point-level noise, stop and check
whether Task 5's loader or Task 6's failure-counting introduced an unintended change to control
flow or draw order before continuing to Step 3.

- [ ] **Step 3: Update the UQ paragraph in `example_urea.tex`** (median $105$, mean $156$, CI
$[17,614]$ from Step 2, each rounded to the same 2-3 significant figures the original prose
used for its own bootstrap statistics)

Replace (around line 120-143):

```latex
The single answer above treats every bound as exact, but each was assembled
from database estimates that carry their own uncertainty, and propagating that
uncertainty turns the point prediction into a distribution. Assigning each
turnover number and the reference abundance a lognormal factor of roughly two,
each standard free energy a normal spread of a few kilojoules per mole, and
four of the five saturation factors (the fifth, $f_2$ at argininosuccinate
lyase, held fixed at one for lack of a measured substrate) a substrate
concentration and Michaelis constant drawn from the compilation of Park and
coworkers \cite{Park2016}, and re-solving the linear program over ten thousand
parameter draws by the parametric bootstrap of Algorithm~\ref{alg:uq}, spreads
urea export over a right-skewed band with median $0.029$ and mean $0.044$
mmol\,gDW$^{-1}$\,h$^{-1}$ and a central ninety-five percent interval of
$0.005$ to $0.17$. A capacity-only bootstrap, holding every saturation factor
at one and letting only the turnover numbers, the reference abundance, and the
standard free energies vary, gives a median urea export of $0.031$
mmol\,gDW$^{-1}$\,h$^{-1}$, and restoring the four measured saturation factors
leaves the median at $0.029$, demonstrating their inertness.
```

with (substituting your Step 2 output if it differs from the illustrative values below):

```latex
The single answer above treats every bound as exact, but each was assembled
from database estimates that carry their own uncertainty, and propagating that
uncertainty turns the point prediction into a distribution. Assigning each
turnover number and the reference abundance a lognormal factor of roughly two,
each standard free energy a normal spread of a few kilojoules per mole, and
four of the five saturation factors (the fifth, $f_2$ at argininosuccinate
lyase, held fixed at one for lack of a measured substrate) a substrate
concentration and Michaelis constant drawn from the compilation of Park and
coworkers \cite{Park2016}, and re-solving the linear program over ten thousand
parameter draws by the parametric bootstrap of Algorithm~\ref{alg:uq}, spreads
urea export over a right-skewed band with median $105$ and mean $156$
mmol\,gDW$^{-1}$\,h$^{-1}$ and a central ninety-five percent interval of
$17$ to $614$. A capacity-only bootstrap, holding every saturation factor
at one and letting only the turnover numbers, the reference abundance, and the
standard free energies vary, gives a median urea export of $113$
mmol\,gDW$^{-1}$\,h$^{-1}$, and restoring the four measured saturation factors
leaves the median at $105$, demonstrating their inertness.
```

- [ ] **Step 4: Update the saturation-imputation paragraph**

Replace (around line 145-157):

```latex
Imputing
that missing factor from a physiological argininosuccinate concentration and a
BRENDA Michaelis constant, and sampling it alongside the rest, drags the urea
median down to $0.014$ mmol\,gDW$^{-1}$\,h$^{-1}$ and widens the band toward
zero (Fig.~\ref{fig:urea_uq}), the shift and the spread both governed by the
single quantity the model never measured.
```

with:

```latex
Imputing
that missing factor from a physiological argininosuccinate concentration and a
BRENDA Michaelis constant, and sampling it alongside the rest, drags the urea
median down to $51$ mmol\,gDW$^{-1}$\,h$^{-1}$ and widens the band toward
zero (Fig.~\ref{fig:urea_uq}), the shift and the spread both governed by the
single quantity the model never measured.
```

- [ ] **Step 5: Update the Config A/B figure caption**

Replace (around line 228-236):

```latex
  \caption{Saturation-gateway sensitivity for urea export. Configuration A
    (blue) samples the four saturation factors that Park and coworkers
    \cite{Park2016} constrain and holds $f_2 = 1$ at the argininosuccinate
    lyase bottleneck, whose substrate is unmeasured, giving a median urea
    export of $0.029$ mmol\,gDW$^{-1}$\,h$^{-1}$. Configuration B (red)
    imputes $f_2$ from a physiological argininosuccinate concentration and a
    BRENDA Michaelis constant and samples it as well, pulling the median down
    to $0.014$ mmol\,gDW$^{-1}$\,h$^{-1}$ and widening the band toward zero.
    The horizontal axis is logarithmic.}
```

with:

```latex
  \caption{Saturation-gateway sensitivity for urea export. Configuration A
    (blue) samples the four saturation factors that Park and coworkers
    \cite{Park2016} constrain and holds $f_2 = 1$ at the argininosuccinate
    lyase bottleneck, whose substrate is unmeasured, giving a median urea
    export of $105$ mmol\,gDW$^{-1}$\,h$^{-1}$. Configuration B (red)
    imputes $f_2$ from a physiological argininosuccinate concentration and a
    BRENDA Michaelis constant and samples it as well, pulling the median down
    to $51$ mmol\,gDW$^{-1}$\,h$^{-1}$ and widening the band toward zero.
    The horizontal axis is logarithmic.}
```

- [ ] **Step 6: Copy the regenerated figures into the chapter**

Run:
```bash
cp code/figs/urea_fba.pdf chapter/figures/urea_fba.pdf
cp code/figs/urea_saturation.pdf chapter/figures/urea_saturation.pdf
```

- [ ] **Step 7: Commit**

```bash
git add code/fba/urea_cycle_uq.jl chapter/sections/example_urea.tex \
        code/data/urea_fba_uq.csv code/data/urea_uq_sensitivity.csv \
        chapter/figures/urea_fba.pdf chapter/figures/urea_saturation.pdf
git commit -m "Fix urea UQ targets and regenerate outputs/prose to corrected units"
```

---

## Task 8: Align the general protein balance in `gateways.tex` with the feedback implementation

**Files:**
- Modify: `chapter/sections/gateways.tex`

**Interfaces:**
- Consumes: nothing (LaTeX-only).
- Produces: `\eqref{eq:protein}` and `\eqref{eq:expression-ss}` now include `m_j`, matching
  `code/feedback/dual_feedback.jl`'s `rTL::{m}` kinetics (translation order 1 in `m`).

- [ ] **Step 1: Confirm no other section restates these equations' RHS**

Run: `grep -n "eq:protein\|eq:expression-ss" chapter/sections/*.tex`
Expected: hits only in `gateways.tex` (definition) and possibly a citation-only `\eqref{...}`
in `example_feedback.tex` (which cites but does not restate the formula, per the design spec).
If any other section restates the old (no-`m_j`) formula inline, note it — no such restatement
is expected based on the current draft.

- [ ] **Step 2: Update the protein balance equation**

In `chapter/sections/gateways.tex`, replace:

```latex
\begin{align}
  \dot m_j &= r_{X,j}\,u_j - (\theta_{m,j}+\mu)\,m_j + \lambda_j,
  \label{eq:mrna}\\
  \dot p_j &= r_{L,j}\,w_j - (\theta_{p,j}+\mu)\,p_j,
  \label{eq:protein}
\end{align}
```

with:

```latex
\begin{align}
  \dot m_j &= r_{X,j}\,u_j - (\theta_{m,j}+\mu)\,m_j + \lambda_j,
  \label{eq:mrna}\\
  \dot p_j &= r_{L,j}\,w_j\,m_j - (\theta_{p,j}+\mu)\,p_j,
  \label{eq:protein}
\end{align}
```

- [ ] **Step 3: Update the closed-form steady-state protein abundance**

Replace:

```latex
\begin{equation}\label{eq:expression-ss}
  m^{\star}_j = \frac{r_{X,j}u_j+\lambda_j}{\theta_{m,j}+\mu},
  \qquad
  p^{\star}_j = \frac{r_{L,j}w_j}{\theta_{p,j}+\mu},
  \qquad
  (e/e^{\circ}) = \frac{p^{\star}_j}{e^{\circ}},
\end{equation}
```

with:

```latex
\begin{equation}\label{eq:expression-ss}
  m^{\star}_j = \frac{r_{X,j}u_j+\lambda_j}{\theta_{m,j}+\mu},
  \qquad
  p^{\star}_j = \frac{r_{L,j}w_j\,m^{\star}_j}{\theta_{p,j}+\mu},
  \qquad
  (e/e^{\circ}) = \frac{p^{\star}_j}{e^{\circ}},
\end{equation}
```

- [ ] **Step 4: Rebuild the chapter and confirm it compiles**

Run:
```bash
cp code/figs/*.pdf chapter/figures/ 2>/dev/null
cd chapter && make && cd ..
```
Expected: `chapter/Chapter.pdf` is produced with no `pdflatex` fatal errors (undefined
references are acceptable only if they pre-exist elsewhere and are unrelated to this edit;
grep the log for `eq:protein` / `eq:expression-ss` specifically and confirm no "undefined"
warning appears for those two labels).

- [ ] **Step 5: Commit**

```bash
git add chapter/sections/gateways.tex
git commit -m "Add missing transcript dependence to the general protein balance"
```

---

## Task 9: Implement the bounded, independent θ_FBA(X3) in `dual_feedback.jl`

**Files:**
- Modify: `code/feedback/dual_feedback.jl`
- Test: `code/feedback/test_dual_feedback.jl`

**Interfaces:**
- Consumes: `feedback_truth`, `truth.Xss`, `truth.reactions`, `truth.model` (all existing,
  unchanged).
- Produces: `gateway_factors(truth)` still returns a NamedTuple `(Vmax0, θ, e_e0)` — same
  shape, same field names, so `feedback_fba` (Task 10, Task 11) needs no interface change —
  but `θ` is now computed from the bounded Hill form, not `X3^{-a}`. New module-level constants
  `K_THETA = 5.0`, `N_THETA = 2.0`.

- [ ] **Step 1: Update the failing test assertions first**

In `code/feedback/test_dual_feedback.jl`, replace the Task 3 block:

```julia
    gw = gateway_factors(t)
    @assert isapprox(gw.Vmax0, 10.0; atol=1e-9) "Vmax0: $(gw.Vmax0)"
    @assert isapprox(gw.θ,     0.574; atol=0.02) "theta*: $(gw.θ)"
    @assert isapprox(gw.e_e0,  0.464; atol=0.03) "(e/e0)*: $(gw.e_e0)"
    @assert gw.θ < 0.9 && gw.e_e0 < 0.9 "both gateways must be genuinely sub-unity"

    vtruth = truth_metabolic_fluxes(t)
    v_naive = feedback_fba(gw; expression=false, activity=false)
    v_expr  = feedback_fba(gw; expression=true,  activity=false)
    v_act   = feedback_fba(gw; expression=false, activity=true)
    v_both  = feedback_fba(gw; expression=true,  activity=true)

    j = findfirst(==("r3"), METAB_REACTIONS)
    @assert maximum(abs.(v_both .- vtruth)) < 1e-2 "both-open must recover truth: $(v_both) vs $(vtruth)"
    @assert v_naive[j] - vtruth[j] > 0.5 "naive must overshoot: $(v_naive[j]) vs $(vtruth[j])"
    T, N, Ee, Aa = vtruth[j], v_naive[j], v_expr[j], v_act[j]
    @assert T < Ee < N "expression-only must be strictly bracketed: $T < $Ee < $N"
    @assert T < Aa < N "activity-only must be strictly bracketed: $T < $Aa < $N"
    @assert !isapprox(Ee, Aa; atol=0.2) "the two partial cases must be visibly distinct: $Ee vs $Aa"
    println("test_dual_feedback (Task 3): naive=$N expr=$Ee act=$Aa both=$(v_both[j]) truth=$T")
```

with:

```julia
    gw = gateway_factors(t)
    @assert isapprox(gw.Vmax0, 10.0; atol=1e-9) "Vmax0: $(gw.Vmax0)"
    @assert isapprox(gw.θ,     0.610; atol=0.01) "theta*: $(gw.θ)"
    @assert isapprox(gw.e_e0,  0.464; atol=0.03) "(e/e0)*: $(gw.e_e0)"
    @assert gw.θ < 0.9 && gw.e_e0 < 0.9 "both gateways must be genuinely sub-unity"

    vtruth = truth_metabolic_fluxes(t)
    v_naive = feedback_fba(gw; expression=false, activity=false)
    v_expr  = feedback_fba(gw; expression=true,  activity=false)
    v_act   = feedback_fba(gw; expression=false, activity=true)
    v_both  = feedback_fba(gw; expression=true,  activity=true)

    j = findfirst(==("r3"), METAB_REACTIONS)
    @assert maximum(abs.(v_both .- vtruth) ./ vtruth) < 0.10 "both-open must approximate truth within 10%: $(v_both) vs $(vtruth)"
    @assert v_naive[j] - vtruth[j] > 0.5 "naive must overshoot: $(v_naive[j]) vs $(vtruth[j])"
    T, N, Ee, Aa = vtruth[j], v_naive[j], v_expr[j], v_act[j]
    @assert T < Ee < N "expression-only must be strictly bracketed: $T < $Ee < $N"
    @assert T < Aa < N "activity-only must be strictly bracketed: $T < $Aa < $N"
    @assert !isapprox(Ee, Aa; atol=0.2) "the two partial cases must be visibly distinct: $Ee vs $Aa"
    println("test_dual_feedback (Task 3): naive=$N expr=$Ee act=$Aa both=$(v_both[j]) truth=$T")
```

- [ ] **Step 2: Run to confirm the new targets fail against the current code**

Run: `julia --project=code code/feedback/test_dual_feedback.jl`
Expected: `AssertionError: theta*: 0.5744...` (the un-fixed code still returns the old
`X3^{-a}`-based value).

- [ ] **Step 3: Replace `gateway_factors`' θ computation**

In `code/feedback/dual_feedback.jl`, replace:

```julia
# Reference capacity Vmax0, activity control theta, and expression ratio (e/e0)
# at the MEASURED steady state.
#   Vmax0 = alpha[r0]            reference capacity (E0=1, no repression, no allostery)
#   theta = X3^{-A_ACT}          computed rate-with / rate-without the X3 order on r0
#   e_e0  = E0*                  measured enzyme abundance (reference e0 = 1 by normalization)
function gateway_factors(truth)
    model  = truth.model
    idx    = Dict(truth.reactions .=> eachindex(truth.reactions))
    sidx   = Dict(truth.species   .=> eachindex(truth.species))
    col_r0 = idx["r0"]
    row_X3 = findfirst(==("X3"), model.total_species_list)

    state       = vcat(truth.Xss, model.static_factors_array)
    rate_with   = BSTModelKit._powerlaw(state, model.α, model.G)[col_r0]   # Vmax0*E0*X3^{-a}
    G0          = copy(model.G); G0[row_X3, col_r0] = 0.0
    rate_noallo = BSTModelKit._powerlaw(state, model.α, G0)[col_r0]        # Vmax0*E0
    θ     = rate_with / rate_noallo                                        # = X3^{-a}
    Vmax0 = model.α[col_r0]                                                # reference capacity
    e_e0  = truth.Xss[sidx["E0"]]                                          # measured (e/e0)
    return (Vmax0 = Vmax0, θ = θ, e_e0 = e_e0)
end
```

with:

```julia
# Reference capacity Vmax0, activity control theta, and expression ratio (e/e0)
# at the MEASURED steady state.
#   Vmax0 = alpha[r0]                  reference capacity (E0=1, no repression, no allostery)
#   theta = K_THETA^N_THETA / (K_THETA^N_THETA + X3^N_THETA)   bounded two-state Hill occupancy,
#           independent of the truth model's own X3^{-a} kinetics (not fit to match it)
#   e_e0  = E0*                        measured enzyme abundance (reference e0 = 1 by normalization)
const K_THETA = 5.0
const N_THETA = 2.0

function gateway_factors(truth)
    sidx = Dict(truth.species .=> eachindex(truth.species))
    X3   = truth.Xss[sidx["X3"]]

    Vmax0 = truth.model.α[findfirst(==("r0"), truth.reactions)]  # reference capacity
    θ     = K_THETA^N_THETA / (K_THETA^N_THETA + X3^N_THETA)     # bounded Hill occupancy
    e_e0  = truth.Xss[sidx["E0"]]                                 # measured (e/e0)
    return (Vmax0 = Vmax0, θ = θ, e_e0 = e_e0)
end
```

- [ ] **Step 4: Run the test to confirm it passes**

Run: `julia --project=code code/feedback/test_dual_feedback.jl`
Expected:
```
test_dual_feedback (Task 1): TOML builds, orderings + metabolic S OK
test_dual_feedback (Task 2): fixed point X3*=3.99862308785145 T*=2.663081099128929 ; mu-frac transcript=0.039999999999999994 protein=0.3684210526315789
test_dual_feedback (Task 3): naive=10.0 expr=4.636055003996753 act=6.09919933975477 both=2.827622361944379 truth=2.663082976509066
```

- [ ] **Step 5: Commit**

```bash
git add code/feedback/dual_feedback.jl code/feedback/test_dual_feedback.jl
git commit -m "Replace X3^-a activity gateway with a bounded, independent Hill occupancy"
```

---

## Task 10: Harden `dual_feedback.jl` — wrap the remaining `_powerlaw` use and validate `feedback_fba` inputs

**Files:**
- Modify: `code/feedback/dual_feedback.jl`
- Modify: `code/Manifest.toml` (comment only, no version change)
- Test: `code/feedback/test_dual_feedback.jl`

**Interfaces:**
- Consumes: `reaction_fluxes`, `feedback_fba` from the existing file (Task 9 already removed
  `gateway_factors`' use of `_powerlaw`; this task addresses the one remaining call site).
- Produces: `_powerlaw` is called through a locally-named wrapper `bst_powerlaw_rates`;
  `feedback_fba` raises `ArgumentError` on malformed bounds instead of silently building an
  infeasible/nonsensical LP.

- [ ] **Step 1: Write the failing validation test**

Append to `code/feedback/test_dual_feedback.jl` (before the final closing of the file, as a
new top-level `let` block):

```julia
# ---- Task 4: feedback_fba input validation ---------------------------------- #
let gw_bad = (Vmax0 = -1.0, θ = 0.5, e_e0 = 0.5)
    threw = false
    try
        feedback_fba(gw_bad; expression=true, activity=true)
    catch e
        threw = e isa ArgumentError
    end
    @assert threw "feedback_fba must reject a negative Vmax0"
    println("test_dual_feedback (Task 4): feedback_fba validation OK")
end
```

- [ ] **Step 2: Run to confirm it fails**

Run: `julia --project=code code/feedback/test_dual_feedback.jl`
Expected: `AssertionError: feedback_fba must reject a negative Vmax0` (the current code builds
a negative-upper-bound LP silently instead of raising).

- [ ] **Step 3: Wrap `_powerlaw` and add validation to `feedback_fba`**

In `code/feedback/dual_feedback.jl`, replace:

```julia
# Per-reaction rate vector from BSTModelKit's own power-law kernel, so extracted
# fluxes are identical to what the solver integrated. state = [dynamic; static].
function reaction_fluxes(model::BSTModel, X::Vector{Float64})::Vector{Float64}
    state_array = vcat(X, model.static_factors_array)
    return BSTModelKit._powerlaw(state_array, model.α, model.G)
end
```

with:

```julia
# BSTModelKit.jl is pinned in code/Manifest.toml; this wrapper isolates the one remaining
# dependency on the private `_powerlaw` kernel so a future package upgrade fails here, at a
# single named call site, rather than silently inside a private API call.
bst_powerlaw_rates(state, α, G) = BSTModelKit._powerlaw(state, α, G)

# Per-reaction rate vector from BSTModelKit's own power-law kernel, so extracted
# fluxes are identical to what the solver integrated. state = [dynamic; static].
function reaction_fluxes(model::BSTModel, X::Vector{Float64})::Vector{Float64}
    state_array = vcat(X, model.static_factors_array)
    return bst_powerlaw_rates(state_array, model.α, model.G)
end
```

Then replace the start of `feedback_fba`:

```julia
function feedback_fba(gw; expression::Bool, activity::Bool)
    n   = length(METAB_REACTIONS)
    idx = Dict(METAB_REACTIONS .=> eachindex(METAB_REACTIONS))
    ub  = fill(GENEROUS_CAPACITY, n)
    ub[idx["r0"]] = gw.Vmax0 * (expression ? gw.e_e0 : 1.0) * (activity ? gw.θ : 1.0)
    lb  = zeros(n)                             # all irreversible
```

with:

```julia
function feedback_fba(gw; expression::Bool, activity::Bool)
    gw.Vmax0 >= 0 || throw(ArgumentError("Vmax0 must be nonnegative, got $(gw.Vmax0)"))
    0 <= gw.θ    || throw(ArgumentError("θ must be nonnegative, got $(gw.θ)"))
    0 <= gw.e_e0 || throw(ArgumentError("e_e0 must be nonnegative, got $(gw.e_e0)"))

    n   = length(METAB_REACTIONS)
    idx = Dict(METAB_REACTIONS .=> eachindex(METAB_REACTIONS))
    ub  = fill(GENEROUS_CAPACITY, n)
    ub[idx["r0"]] = gw.Vmax0 * (expression ? gw.e_e0 : 1.0) * (activity ? gw.θ : 1.0)
    lb  = zeros(n)                             # all irreversible
```

- [ ] **Step 4: Add a one-line comment pinning the BSTModelKit version in `code/Manifest.toml`**

Run: `grep -n -A2 '\[\[deps.BSTModelKit\]\]' code/Manifest.toml`

Note the `version = "..."` line printed. No edit is needed to `Manifest.toml` itself (editing a
resolved manifest by hand is not appropriate); the wrapper function's comment in Step 3 already
documents that the pin lives in this file. Skip any further edit here.

- [ ] **Step 5: Run the test to confirm everything passes**

Run: `julia --project=code code/feedback/test_dual_feedback.jl`
Expected: all four `test_dual_feedback (Task N): ...` lines print with no error, including
`test_dual_feedback (Task 4): feedback_fba validation OK`.

- [ ] **Step 6: Commit**

```bash
git add code/feedback/dual_feedback.jl code/feedback/test_dual_feedback.jl
git commit -m "Wrap BSTModelKit._powerlaw and validate feedback_fba inputs"
```

---

## Task 11: Regenerate the feedback outputs and rewrite `example_feedback.tex`

**Files:**
- Run (not modified): `code/feedback/run_feedback.jl`
- Modify: `chapter/sections/example_feedback.tex`
- Generated: `code/data/feedback_fba.csv`, `code/figs/feedback_gateway.pdf`, copied to
  `chapter/figures/feedback_gateway.pdf`

**Interfaces:**
- Consumes: Task 9 (bounded θ_FBA) and Task 10 (hardening) must both be in place.
- Produces: `chapter/sections/example_feedback.tex` consistent with the regenerated CSV/figure,
  with the θ/θ_FBA symbol collision resolved and the "recovered ... to numerical tolerance"
  overclaim corrected.

- [ ] **Step 1: Regenerate the feedback CSV and figure**

Run: `julia --project=code code/feedback/run_feedback.jl`
Expected output:
```
fixed point: X3*=3.999  T*=2.663  theta*=0.61  (e/e0)*=0.464
overshoot factor Vmax0/T* = 3.755
ledger [Vmax0, x(e/e0), x theta, both] = [10.0, 4.636, 6.099, 2.828]
mu fraction of clearance: transcript=0.04  protein=0.368
run_feedback OK: wrote feedback_fba.csv and feedback_gateway.pdf
```
(A `Warning: arrows are deprecated...` message from CairoMakie/Makie may print — that is a
pre-existing cosmetic deprecation warning unrelated to this change and can be ignored.)

- [ ] **Step 2: Copy the regenerated figure into the chapter**

Run: `cp code/figs/feedback_gateway.pdf chapter/figures/feedback_gateway.pdf`

- [ ] **Step 3: Disambiguate the truth model's kinetic exponent from the new FBA gateway symbol**

In `chapter/sections/example_feedback.tex`, replace the committed-step rate equation (around
line 64-69):

```latex
The committed-step rate carries both controls,
\begin{equation}\label{eq:feedback-rate}
  v_{r_0} = V^{\circ}_{\max,r_0}\,(e/e^{\circ})\,\theta,
  \qquad
  \theta = X_3^{-a},\quad a = 0.4,\qquad
  \text{transcription} \propto X_3^{-b},\quad b = 0.6.
\end{equation}
```

with:

```latex
The committed-step rate carries both controls,
\begin{equation}\label{eq:feedback-rate}
  v_{r_0} = V^{\circ}_{\max,r_0}\,(e/e^{\circ})\,\theta_{\mathrm{true}},
  \qquad
  \theta_{\mathrm{true}} = X_3^{-a},\quad a = 0.4,\qquad
  \text{transcription} \propto X_3^{-b},\quad b = 0.6.
\end{equation}
```

- [ ] **Step 4: Update the self-consistency paragraph's numeric values and remove the ambiguous starred-θ symbol**

Replace (around line 101-114):

```latex
Because the chain is linear, conservation forces every metabolic
reaction to carry one common throughput $T$, and because the export step
$r_3$ is first order in $X_3$, that throughput pins the end-product level
at $X_3 = T/k_3$. Substituting the two controls back into the
committed-step rate, the activity factor $\theta=X_3^{-a}$ and the
expression ratio the transcriptional repression sets through
Equation~\eqref{eq:expression-ss}, reduces the coupled five-state system
to a single self-consistency condition in $X_3$ alone, and the value
that condition selects is the value the full numerical integration
confirmed: $X_3^{\star}\approx4$, an activity factor
$\theta^{\star}\approx0.574$, an expression ratio
$(e/e^{\circ})^{\star}\approx0.464$, and a throughput
$T^{\star}\approx2.66$, less than a third of the uninhibited capacity of
the committed step (Fig.~\ref{fig:feedback}).
```

with:

```latex
Because the chain is linear, conservation forces every metabolic
reaction to carry one common throughput $T$, and because the export step
$r_3$ is first order in $X_3$, that throughput pins the end-product level
at $X_3 = T/k_3$. Substituting the two controls back into the
committed-step rate, the truth model's own allosteric term
$\theta_{\mathrm{true}}=X_3^{-a}$ and the expression ratio the
transcriptional repression sets through
Equation~\eqref{eq:expression-ss}, reduces the coupled five-state system
to a single self-consistency condition in $X_3$ alone, and the value
that condition selects is the value the full numerical integration
confirmed: $X_3^{\star}\approx4$, an allosteric term
$\theta_{\mathrm{true}}^{\star}\approx0.574$, an expression ratio
$(e/e^{\circ})^{\star}\approx0.464$, and a throughput
$T^{\star}\approx2.66$, less than a third of the uninhibited capacity of
the committed step (Fig.~\ref{fig:feedback}).
```

- [ ] **Step 5: Insert the bounded θ_FBA definition and rewrite the ledger paragraph**

Replace the entire paragraph (around line 128-159, from "The flux-balance model kept..."
through "...settled on the throughput the cell actually ran."):

```latex
The flux-balance model kept the chain's stoichiometry unchanged and
placed the entire dual-control content on the single bound of the
committed step,
\begin{equation}\label{eq:feedback-bound}
  0 \;\le\; \hat v_{r_0} \;\le\;
  V^{\circ}_{\max,r_0}\,(e/e^{\circ})\,\theta(X_3),
  \qquad V^{\circ}_{\max,r_0} = 10,
\end{equation}
the expression and regulatory instance of the general bound of
Equation~\eqref{eq:general-bound}, with the thermodynamic switch and the
saturation factor left at their defaults and the two remaining gateways,
expression $(e/e^{\circ})$ and activity $\theta(X_3)$, both in play. Both
factors were evaluated from the measured steady state, the enzyme
abundance and the end-product level read off the integrated model, and
not copied from the integrated flux the linear program is meant to
reproduce. Four solves of the linear program of
Equation~\eqref{eq:fba-lp} differed in nothing but this one bound. A
feedback-blind solve, holding both factors at unity, overshot to the
uninhibited capacity of $10$, roughly $3.8$ times the true throughput;
opening the expression gateway alone throttled the committed step to
$4.64$; opening the activity gateway alone throttled it to $5.74$; and
only opening both together landed on $2.66$, recovering the BST truth to
numerical tolerance (Fig.~\ref{fig:feedback}). The escalation is the
whole point. Nothing about either feedback loop appears anywhere in the
stoichiometric matrix, which is identical whether $X_3$ closes a single
loop, both loops, or none; the entire dual-control fact lives in the one
bound on $r_0$, where the two gateways enter as independent
multiplicative factors and each removes a distinct, quantifiable slice
of the capacity. Read as a ledger, the uninhibited capacity of $10$ was
marked down once by expression to $4.64$ and once by activity to $5.74$,
and only the product of the two markdowns settled on the throughput the
cell actually ran.
```

with:

```latex
The flux-balance model kept the chain's stoichiometry unchanged and
placed the entire dual-control content on the single bound of the
committed step,
\begin{equation}\label{eq:feedback-bound}
  0 \;\le\; \hat v_{r_0} \;\le\;
  V^{\circ}_{\max,r_0}\,(e/e^{\circ})\,\theta_{\mathrm{FBA}}(X_3),
  \qquad V^{\circ}_{\max,r_0} = 10,
\end{equation}
the expression and regulatory instance of the general bound of
Equation~\eqref{eq:general-bound}, with the thermodynamic switch and the
saturation factor left at their defaults and the two remaining gateways,
expression $(e/e^{\circ})$ and activity $\theta_{\mathrm{FBA}}(X_3)$,
both in play. The activity gateway does not reuse the truth model's own
kinetic exponent $\theta_{\mathrm{true}}$ from Equation~\eqref{eq:feedback-rate};
a bound built from the exact functional form that generates the
quantity it is meant to predict would make the recovery an identity
rather than a check. Instead $\theta_{\mathrm{FBA}}$ is a bounded,
two-state Hill occupancy,
\begin{equation}\label{eq:feedback-theta-fba}
  \theta_{\mathrm{FBA}}(X_3) \;=\; \frac{K^{n}}{K^{n}+X_3^{\,n}},
  \qquad K = 5,\quad n = 2,
\end{equation}
the same partition-function-motivated form used for $f_s$ in
Equation~\eqref{eq:partition}, with $K$ and $n$ chosen independently of
the truth model's allosteric order $a$ rather than fit to reproduce it,
and with no basal-leak floor: unlike the transcriptional loop, a direct
allosteric on/off switch has no leaky analog, so
$\theta_{\mathrm{FBA}}\to0$ as $X_3\to\infty$. Both gateway factors were
evaluated from the measured steady state, the enzyme abundance and the
end-product level read off the integrated model, and not copied from
the integrated flux the linear program is meant to reproduce. Four
solves of the linear program of Equation~\eqref{eq:fba-lp} differed in
nothing but this one bound. A feedback-blind solve, holding both factors
at unity, overshot to the uninhibited capacity of $10$, roughly $3.8$
times the true throughput; opening the expression gateway alone
throttled the committed step to $4.64$; opening the activity gateway
alone throttled it to $6.10$; and only opening both together landed on
$2.83$, approximating the BST truth of $2.66$ to within about six
percent (Fig.~\ref{fig:feedback}). The escalation is the whole point.
Nothing about either feedback loop appears anywhere in the
stoichiometric matrix, which is identical whether $X_3$ closes a single
loop, both loops, or none; the entire dual-control fact lives in the one
bound on $r_0$, where the two gateways enter as independent
multiplicative factors and each removes a distinct, quantifiable slice
of the capacity. Read as a ledger, the uninhibited capacity of $10$ was
marked down once by expression to $4.64$ and once by activity to $6.10$,
and only the product of the two markdowns approached the throughput the
cell actually ran.
```

- [ ] **Step 6: Update the Figure 2 caption**

Replace (around line 161-178):

```latex
\begin{figure}[htbp]
  \centering
  \includegraphics[width=0.95\textwidth]{feedback_gateway.pdf}
  \caption{The dual-control feedback worked example. Left: the wiring of
    the linear chain of Equation~\eqref{eq:feedback-network}, in which
    the end product $X_3$ closes two negative loops onto the committed
    step $r_0$, repressing transcription of the enzyme's gene, the slow
    expression loop acting through $(e/e^{\circ})$, and inhibiting the
    activity of the enzyme $E_0$ already present, the fast allosteric
    loop acting through $\theta$. Right: the capacity ledger for the
    committed-step bound of Equation~\eqref{eq:feedback-bound}. The
    uninhibited capacity $V^{\circ}_{\max,r_0}=10$ is reduced by the
    expression gateway alone to $4.64$ and by the activity gateway alone
    to $5.74$, and only the two gateways together land on the throughput
    of the Biochemical Systems Theory truth, $T^{\star}\approx2.66$; the
    dashed line marks the uninhibited capacity.}
  \label{fig:feedback}
\end{figure}
```

with:

```latex
\begin{figure}[htbp]
  \centering
  \includegraphics[width=0.95\textwidth]{feedback_gateway.pdf}
  \caption{The dual-control feedback worked example. Left: the wiring of
    the linear chain of Equation~\eqref{eq:feedback-network}, in which
    the end product $X_3$ closes two negative loops onto the committed
    step $r_0$, repressing transcription of the enzyme's gene, the slow
    expression loop acting through $(e/e^{\circ})$, and inhibiting the
    activity of the enzyme $E_0$ already present, the fast allosteric
    loop acting through $\theta_{\mathrm{FBA}}$. Right: the capacity ledger for the
    committed-step bound of Equation~\eqref{eq:feedback-bound}. The
    uninhibited capacity $V^{\circ}_{\max,r_0}=10$ is reduced by the
    expression gateway alone to $4.64$ and by the activity gateway alone
    to $6.10$, and only the two gateways together approach the throughput
    of the Biochemical Systems Theory truth, $T^{\star}\approx2.66$, landing
    at $2.83$; the dashed line marks the uninhibited capacity.}
  \label{fig:feedback}
\end{figure}
```

- [ ] **Step 7: Commit**

```bash
git add chapter/sections/example_feedback.tex chapter/figures/feedback_gateway.pdf code/data/feedback_fba.csv
git commit -m "Regenerate feedback example and rewrite prose for the bounded theta_FBA gateway"
```

---

## Task 12: Final full rebuild and test run

**Files:** none modified; this task only runs and verifies.

**Interfaces:**
- Consumes: every prior task in this plan.
- Produces: a clean `chapter/Chapter.pdf` and two passing Julia test suites, confirming Phase 1
  is internally consistent before Phase 2 starts.

- [ ] **Step 1: Run both test suites**

Run:
```bash
julia --project=code code/fba/test_urea_cycle.jl
julia --project=code code/feedback/test_dual_feedback.jl
```
Expected: both print their full sequence of `... : ... passed` / `... OK` lines with no
`AssertionError`.

- [ ] **Step 2: Regenerate every figure this plan touched, in order, and sync into the chapter**

Run:
```bash
julia --project=code code/fba/run_fba.jl
julia --project=code code/fba/urea_cycle_uq.jl
julia --project=code code/feedback/run_feedback.jl
cp code/figs/*.pdf chapter/figures/
```

- [ ] **Step 3: Full chapter rebuild**

Run: `cd chapter && make && cd ..`
Expected: `chapter/Chapter.pdf` is written; no `pdflatex` run reports "Emergency stop" or a
fatal error. Undefined-reference or undefined-citation warnings unrelated to labels touched in
this plan (`eq:protein`, `eq:expression-ss`, `eq:feedback-rate`, `eq:feedback-bound`,
`eq:feedback-theta-fba`, `fig:urea`, `fig:urea_uq`, `fig:feedback`) are pre-existing and out of
scope.

Run: `grep -i "undefined" chapter/Chapter.log | grep -E "protein|expression-ss|feedback-rate|feedback-bound|feedback-theta-fba"`
Expected: no output (empty grep result — none of the labels touched by this plan are
undefined).

- [ ] **Step 4: Spot-check the corrected numbers made it into the built PDF**

Run: `grep -c "118.08" chapter/Chapter.tex chapter/sections/example_urea.tex`
Expected: nonzero count in `example_urea.tex` (the file itself; `Chapter.tex` is the root
document and may show 0 since it only `\input`s the sections).

Run: `grep -c "2.83\|6.10\|theta_{\\\\mathrm{FBA}}" chapter/sections/example_feedback.tex`
Expected: nonzero.

- [ ] **Step 5: Commit (if Step 2's regeneration produced any diff beyond what earlier tasks already committed)**

Run: `git status --short`

If any of `code/data/*.csv`, `code/figs/*.pdf`, or `chapter/figures/*.pdf` show as modified
(possible if Step 2's fresh RNG-seeded run produced negligible floating-point differences from
earlier task commits), commit them:

```bash
git add code/data code/figs chapter/figures
git commit -m "Regenerate Phase 1 figures/data from the final code state"
```

If `git status --short` shows no changes, skip this step — there is nothing to commit.
