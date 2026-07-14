# Dual-Control Feedback Example Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rebuild the linear-feedback worked example so one committed enzyme is controlled by its own end product through two gateways at once — fast allosteric activity regulation (`theta`) and slow transcriptional expression control (`(e/e0)`, an Alon negative-feedback loop) — and both are recovered through the product bound `Vmax0 * (e/e0) * theta`.

**Architecture:** Extend the BSTModelKit S-system to a five-species cascade (`X1,X2,X3,m,E0`) integrated to a self-consistent dual-feedback steady state (the "truth"). A constraint-based layer recovers that truth through the committed-step bound in four escalating cases (naive → expression-only → activity-only → both). Expression capacities come from sequence-derived quantities; `(e/e0)` and `theta` are computed from the measured state, never copied from the integrated fluxes. Code splits into `dual_feedback.jl` (functions), `test_dual_feedback.jl` (assert gates), `run_feedback.jl` (driver: CSV + figure).

**Tech Stack:** Julia; BSTModelKit (GMA power-law S-systems), JuMP + HiGHS (LP), CairoMakie (figure), CSV/DataFrames. LaTeX (`article` class, `chapter/Makefile`).

## Global Constraints

- Entry point: every script's first line is `include(joinpath(@__DIR__, "..", "Include.jl"))`, which `Pkg.activate`s `code/` and brings in `CSV, DataFrames, Statistics, LinearAlgebra, Random, JuMP, HiGHS, BSTModelKit, CairoMakie` plus `figpath`/`datapath` (`code/src/Runtime.jl`).
- Run commands from the repo root: `julia --project=code code/feedback/<script>.jl`.
- BSTModelKit facts (from the existing `run_feedback.jl`): `build(path)::BSTModel`; `steadystate(model; tspan)` returns the state vector; `BSTModelKit._powerlaw(state, model.α, model.G)` returns the per-reaction rate vector where `state = vcat(dynamic_state, model.static_factors_array)`. `model.G` rows index `model.total_species_list = vcat(dynamic, static)`; `model.G` cols index `model.list_of_reactions`. Kinetics records default each listed species to order `+1`; set feedback orders by name after `build`.
- Never read fluxes back from `truth.vss`: `Vmax0`, `theta`, `(e/e0)` are all computed from the model at the measured steady state.
- Feedback strengths and scales (verbatim): `A_ACT = 0.4`, `B_EXPR = 0.6`, `PHI = 0.05`, `VMAX0 = 10.0`, `K3 = 0.666`, `FAST = 100.0`. Fixed point targets: `X3* ≈ 4.0`, `T* ≈ 2.66`, `theta* ≈ 0.574`, `(e/e0)* ≈ 0.464`, overshoot `≈ 3.8x`. Ledger: `Vmax0=10 → x(e/e0)=4.64 → x theta=5.74 → both=2.66`.
- Style (from `docs/superpowers/specs/2026-07-14-dual-control-feedback-example-design.md` and the chapter conventions): no em dashes; no subsection headings in the example; no self-references ("this section..."); long complete paragraphs; `mu` is one culture-wide growth rate identical in both expression denominators — never imply a species-dependent `mu`; numeric results in prose tie to the figure or are structural setup numbers.
- Commit message trailer on every commit: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.

---

### Task 1: Dual-Feedback BST model spec

**Files:**
- Create: `code/feedback/Dual-Feedback.toml`
- Delete: `code/feedback/Linear-Feedback.toml`
- Test: `code/feedback/test_dual_feedback.jl`

**Interfaces:**
- Produces: a BSTModelKit TOML that `build` parses into a model with `list_of_dynamic_species == ["X1","X2","X3","m","E0"]`, `list_of_static_species == ["E"]`, and `list_of_reactions == ["rTX","rTXb","rMdeg","rTL","rEdeg","r0","r1","r2","r3"]` (connection-record order). The metabolic submatrix `model.S[rows X1,X2,X3][cols r0,r1,r2,r3]` equals `[1 -1 0 0; 0 1 -1 0; 0 0 1 -1]`.

- [ ] **Step 1: Write the failing test**

Create `code/feedback/test_dual_feedback.jl`:

```julia
include(joinpath(@__DIR__, "..", "Include.jl"))

# ---- Task 1: the Dual-Feedback TOML builds with the expected orderings ----- #
let model = build(joinpath(@__DIR__, "Dual-Feedback.toml"))
    @assert model.list_of_dynamic_species == ["X1","X2","X3","m","E0"] "dyn species: $(model.list_of_dynamic_species)"
    @assert model.list_of_reactions == ["rTX","rTXb","rMdeg","rTL","rEdeg","r0","r1","r2","r3"] "reactions: $(model.list_of_reactions)"
    # metabolic submatrix (rows X1,X2,X3 ; cols r0..r3) is the linear chain
    srow(s) = findfirst(==(s), model.list_of_dynamic_species)
    scol(r) = findfirst(==(r), model.list_of_reactions)
    Smet = model.S[[srow("X1"),srow("X2"),srow("X3")], [scol("r0"),scol("r1"),scol("r2"),scol("r3")]]
    @assert Smet == [1.0 -1.0 0.0 0.0; 0.0 1.0 -1.0 0.0; 0.0 0.0 1.0 -1.0] "metabolic S mismatch:\n$Smet"
    println("test_dual_feedback (Task 1): TOML builds, orderings + metabolic S OK")
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=code code/feedback/test_dual_feedback.jl`
Expected: FAIL — `SystemError`/`ArgumentError` opening `Dual-Feedback.toml` (file does not exist yet).

- [ ] **Step 3: Create the TOML**

Create `code/feedback/Dual-Feedback.toml`:

```toml
[metadata]
author = "jdv27@cornell.edu"
version = "0.1"
date = "2026-07-14"
description = "Linear pathway under dual end-product control: expression (e/e0) + activity (theta) on one committed enzyme (chapter Example 2)"

[model]
list_of_static_species = ["E"]
list_of_dynamic_species = ["X1", "X2", "X3", "m", "E0"]
list_of_connection_records = [
    "rTX::{} --> m",
    "rTXb::{} --> m",
    "rMdeg::m --> {}",
    "rTL::{} --> E0",
    "rEdeg::E0 --> {}",
    "r0::{} --> X1",
    "r1::X1 --> X2",
    "r2::X2 --> X3",
    "r3::X3 --> {}"
]
list_of_kinetics_records = [
    "rTX::{X3}",
    "rTXb::{E}",
    "rMdeg::{m}",
    "rTL::{m}",
    "rEdeg::{E0}",
    "r0::{E0, X3}",
    "r1::{X1}",
    "r2::{X2}",
    "r3::{X3}"
]
list_of_stoichiometry_records = []
```

- [ ] **Step 4: Delete the obsolete TOML**

Run: `git rm code/feedback/Linear-Feedback.toml`

- [ ] **Step 5: Run test to verify it passes**

Run: `julia --project=code code/feedback/test_dual_feedback.jl`
Expected: PASS — `test_dual_feedback (Task 1): TOML builds, orderings + metabolic S OK`.

**Fallback (only if Step 5 fails on `build` because a static species may not appear in a kinetics record):** remove the `"rTXb::{} --> m"` connection record and the `"rTXb::{E}"` kinetics record, and drop `"rTXb"` from the expected `list_of_reactions` in the test. This sets the basal leak to zero (`PHI = 0`); Task 2 then uses `PHI = 0.0`, `K3 = 0.625`, and the alternate targets `T* ≈ 2.5`, `(e/e0)* ≈ 0.435` (adjust the Task 2 assert atols accordingly). The leak stays a prose concept only. Record which path was taken in the commit message.

- [ ] **Step 6: Commit**

```bash
git add code/feedback/Dual-Feedback.toml code/feedback/test_dual_feedback.jl
git commit -m "$(cat <<'EOF'
Add dual-feedback BST model spec (5-species cascade)

Replace the linear-feedback TOML with a five-species cascade
(X1,X2,X3,m,E0): transcription (repressed by X3), translation, mRNA and
enzyme turnover, plus the metabolic chain whose committed step r0 is
catalyzed by E0 and allosterically inhibited by X3.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 2: Truth model — sequence params, build, integrate

**Files:**
- Create: `code/feedback/dual_feedback.jl`
- Modify: `code/feedback/test_dual_feedback.jl` (append Task 2 gates)

**Interfaces:**
- Consumes: `Dual-Feedback.toml` (Task 1); `build`, `steadystate`, `BSTModelKit._powerlaw`.
- Produces:
  - `const A_ACT=0.4, B_EXPR=0.6, PHI=0.05, VMAX0=10.0, K3=0.666, FAST=100.0`.
  - `struct SeqParams` with fields `L_prot, L_gene, tau_X, tau_L, theta_m, theta_p, mu` (all `Float64`).
  - `sequence_params()::SeqParams`.
  - `reaction_fluxes(model::BSTModel, X::Vector{Float64})::Vector{Float64}`.
  - `feedback_truth()::NamedTuple` with fields `species::Vector{String}`, `Xss::Vector{Float64}` (order `[X1,X2,X3,m,E0]`), `reactions::Vector{String}`, `vss::Vector{Float64}` (9-vector in reaction order), `model::BSTModel`, `seq::SeqParams`.

- [ ] **Step 1: Write the failing test (append to `test_dual_feedback.jl`)**

Append:

```julia
include(joinpath(@__DIR__, "dual_feedback.jl"))

# ---- Task 2: the integrated truth hits the dual-feedback fixed point ------- #
let t = feedback_truth()
    X1, X2, X3, m, E0 = t.Xss
    idx = Dict(t.reactions .=> eachindex(t.reactions))
    vmet = [t.vss[idx[r]] for r in ["r0","r1","r2","r3"]]
    @assert all(isapprox.(vmet, vmet[1]; rtol=1e-3)) "linear chain not balanced: $vmet"
    T = vmet[1]
    @assert isapprox(X3, 4.0; atol=0.15) "X3* != 4: $X3"
    @assert isapprox(T, 2.66; atol=0.12) "T* != 2.66: $T"
    @assert X3 > 1.0 "end product must accumulate so both gateways are sub-unity"
    @assert isapprox(E0, 0.464; atol=0.03) "E0* (=(e/e0)) off target: $E0"

    sp = t.seq
    frac_m = sp.mu / (sp.theta_m + sp.mu)
    frac_p = sp.mu / (sp.theta_p + sp.mu)
    @assert frac_m < 0.10 "mu should be a small fraction of transcript clearance: $frac_m"
    @assert frac_p > 0.30 "mu should be a dominant fraction of enzyme clearance: $frac_p"
    println("test_dual_feedback (Task 2): fixed point X3*=$X3 T*=$T ; mu-frac transcript=$frac_m protein=$frac_p")
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=code code/feedback/test_dual_feedback.jl`
Expected: FAIL — `UndefVarError: feedback_truth not defined`.

- [ ] **Step 3: Create `dual_feedback.jl` (truth half)**

Create `code/feedback/dual_feedback.jl`:

```julia
# --------------------------------------------------------------------------- #
# Chapter Example 2: dual end-product control of ONE committed enzyme.
#   activity   (theta) : X3 allosterically inhibits E0             (fast)
#   expression (e/e0)  : X3 transcriptionally represses gene(E0)   (slow)
#
# The "truth" is an extended BSTModelKit S-system over [X1,X2,X3,m,E0]. The FBA
# layer (below) recovers it through the product bound Vmax0*(e/e0)*theta on r0.
# Nothing is copied from the integrated fluxes: (e/e0) is the MEASURED enzyme
# abundance (reference e0 = 1 by the nondimensional normalization), and theta is
# computed rate-with / rate-without on the allosteric order.
#
# Expression states are nondimensional ratios normalized to their de-repressed
# reference (mhat = m/m_ref, ehat = E0/e0, both = 1 when repression is off), so
# the metabolic chain stays in the arbitrary units of the existing example
# (reference capacity Vmax0 = 10).
# --------------------------------------------------------------------------- #

const A_ACT  = 0.4     # allosteric (activity) feedback order on X3 in r0
const B_EXPR = 0.6     # transcriptional (expression) repression order on X3 in rTX
const PHI    = 0.05    # basal transcription (leak) fraction; lambda = (theta_m+mu)*PHI
const VMAX0  = 10.0    # reference committed-step capacity (= kcat0 * e0, AU)
const K3     = 0.666   # export rate constant; places the fixed point at X3* = 4
const FAST   = 100.0   # r1, r2 internal steps: fast enough never to limit

# --------------------------------------------------------------------------- #
# Sequence-specific expression parameters (realistic-but-fake). Capacities and
# timescales are computed from gene + protein length and elongation rates; NO
# precursors enter S (that is the cell-free capstone's construction). These set
# the residence times that carry mu and calibrate the reference abundance; the
# (e/e0) response to X3 is carried by B_EXPR.
# --------------------------------------------------------------------------- #
struct SeqParams
    L_prot::Float64; L_gene::Float64
    tau_X::Float64;  tau_L::Float64
    theta_m::Float64; theta_p::Float64; mu::Float64
end

function sequence_params()::SeqParams
    L_prot = 330.0                 # aa, fake committed-step enzyme
    L_gene = 3 * L_prot            # 990 nt (3 nt / codon)
    e_X = 60.0                     # nt/s transcription elongation (E. coli 40-70)
    e_L = 16.0                     # aa/s translation elongation   (E. coli 12-20)
    tau_X = L_gene / e_X           # 16.5 s transcript elongation delay
    tau_L = L_prot / e_L           # 20.6 s protein elongation delay
    t_half_m = 2.5                 # min, mRNA half-life
    t_half_p = 35.0                # min, stable-enzyme half-life
    T_double = 60.0                # min, doubling time
    theta_m = log(2) / t_half_m    # ~0.277 /min
    theta_p = log(2) / t_half_p    # ~0.0198 /min
    mu      = log(2) / T_double    # ~0.0116 /min  (ONE growth rate, both denominators)
    return SeqParams(L_prot, L_gene, tau_X, tau_L, theta_m, theta_p, mu)
end

# Per-reaction rate vector from BSTModelKit's own power-law kernel, so extracted
# fluxes are identical to what the solver integrated. state = [dynamic; static].
function reaction_fluxes(model::BSTModel, X::Vector{Float64})::Vector{Float64}
    state_array = vcat(X, model.static_factors_array)
    return BSTModelKit._powerlaw(state_array, model.α, model.G)
end

# Build the 5-species dual-feedback S-system and integrate to steady state.
function feedback_truth()
    sp    = sequence_params()
    model = build(joinpath(@__DIR__, "Dual-Feedback.toml"))

    species   = ["X1","X2","X3","m","E0"]
    reactions = ["rTX","rTXb","rMdeg","rTL","rEdeg","r0","r1","r2","r3"]
    @assert model.list_of_dynamic_species == species   "dyn-species order changed: $(model.list_of_dynamic_species)"
    @assert model.list_of_reactions == reactions       "reaction order changed: $(model.list_of_reactions)"

    model.initial_condition_array = ones(length(species))   # all > 0 (X3>0 keeps X3^-a finite)
    model.static_factors_array    = [1.0]                    # E: constant driver for the leak

    θm, θp, μ = sp.theta_m, sp.theta_p, sp.mu
    idx = Dict(reactions .=> eachindex(reactions))

    # Rate constants. Expression coefficients are tied to (theta+mu) so the
    # de-repressed reference steady state is mhat = ehat = 1; the repression /
    # leak split sets the operating fraction.
    α = zeros(length(reactions))
    α[idx["rTX"]]   = (θm + μ) * (1 - PHI)   # repressible transcription
    α[idx["rTXb"]]  = (θm + μ) * PHI         # basal leak (rate = coeff * E, E=1)
    α[idx["rMdeg"]] = (θm + μ)               # mRNA turnover + dilution
    α[idx["rTL"]]   = (θp + μ)               # translation (m catalytic)
    α[idx["rEdeg"]] = (θp + μ)               # enzyme turnover + dilution
    α[idx["r0"]]    = VMAX0                   # reference committed-step capacity
    α[idx["r1"]]    = FAST
    α[idx["r2"]]    = FAST
    α[idx["r3"]]    = K3
    model.α = α

    # Feedback kinetic orders, set by name; everything else keeps the default +1
    # the kinetics records declared. Assert those defaults so a package change
    # fails loudly instead of silently mistargeting a feedback.
    trow(s) = findfirst(==(s), model.total_species_list)
    G = model.G
    G[trow("X3"), idx["rTX"]] = -B_EXPR      # X3 represses transcription
    G[trow("X3"), idx["r0"]]  = -A_ACT       # X3 allosterically inhibits r0
    model.G = G
    @assert model.G[trow("E0"), idx["r0"]]   == 1.0 "E0 order in r0 not default +1"
    @assert model.G[trow("m"),  idx["rTL"]]  == 1.0 "m order in rTL not default +1"
    @assert model.G[trow("E"),  idx["rTXb"]] == 1.0 "E order in rTXb not default +1"

    # Stiff system (metabolic alpha ~100 vs expression ~0.03): integrate long
    # enough for the slow protein mode (tau ~ 1/(theta_p+mu) ~ 31 min).
    Xss = steadystate(model; tspan=(0.0, 1000.0))
    vss = reaction_fluxes(model, Xss)
    return (species=species, Xss=Xss, reactions=reactions, vss=vss, model=model, seq=sp)
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=code code/feedback/test_dual_feedback.jl`
Expected: PASS — prints the Task 2 line with `X3*≈4.0`, `T*≈2.66`, `mu-frac transcript≈0.04 protein≈0.37`.

If `X3` lands outside `4.0 ± 0.15`, the integration has not converged: raise `tspan` to `(0.0, 3000.0)`. Do not retune `K3` unless the analytic fixed point itself is wrong (it is verified in the spec).

- [ ] **Step 5: Commit**

```bash
git add code/feedback/dual_feedback.jl code/feedback/test_dual_feedback.jl
git commit -m "$(cat <<'EOF'
Integrate the dual-feedback truth (sequence-specific expression cascade)

Build the 5-species S-system and integrate to the dual-feedback steady
state. Sequence-derived capacities set the (theta+mu) residence times;
the same mu is ~4% of fast transcript clearance but ~37% of stable
enzyme clearance, the visible signature of growth-aware bookkeeping.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 3: Gateway factors and the four-case FBA recovery

**Files:**
- Modify: `code/feedback/dual_feedback.jl` (append FBA half)
- Modify: `code/feedback/test_dual_feedback.jl` (append Task 3 gates)

**Interfaces:**
- Consumes: `feedback_truth()` output (Task 2); `JuMP`, `HiGHS`.
- Produces:
  - `const S_FEEDBACK::Matrix{Float64}` (3x4), `const GENEROUS_CAPACITY=100.0`, `const METAB_REACTIONS=["r0","r1","r2","r3"]`.
  - `gateway_factors(truth)::NamedTuple` with fields `Vmax0::Float64`, `θ::Float64`, `e_e0::Float64`.
  - `feedback_fba(gw; expression::Bool, activity::Bool)::Vector{Float64}` (length 4, order `METAB_REACTIONS`).
  - `truth_metabolic_fluxes(truth)::Vector{Float64}` (length 4, order `METAB_REACTIONS`).

- [ ] **Step 1: Write the failing test (append to `test_dual_feedback.jl`)**

Append:

```julia
# ---- Task 3: FBA recovers the truth only when BOTH gateways are open ------- #
let t = feedback_truth()
    @assert S_FEEDBACK == [1.0 -1.0 0.0 0.0; 0.0 1.0 -1.0 0.0; 0.0 0.0 1.0 -1.0] "S_FEEDBACK wrong"
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
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=code code/feedback/test_dual_feedback.jl`
Expected: FAIL — `UndefVarError: S_FEEDBACK not defined`.

- [ ] **Step 3: Append the FBA half to `dual_feedback.jl`**

Append:

```julia
# =========================================================================== #
# Constraint-based recovery. The SAME linear chain as a flux-balance model; the
# feedback enters ONLY through the committed-step bound, as the product
#   ub[r0] = Vmax0 * (e/e0) * theta
# with each gateway switchable so its marginal effect is visible.
# =========================================================================== #

# Metabolic stoichiometry (rows X1,X2,X3 ; cols r0,r1,r2,r3). Asserted equal to
# the metabolic submatrix of the BST model's S in the test.
const S_FEEDBACK = [
     1.0 -1.0  0.0  0.0    # X1: r0 in, r1 out
     0.0  1.0 -1.0  0.0    # X2: r1 in, r2 out
     0.0  0.0  1.0 -1.0    # X3: r2 in, r3 out
]
const GENEROUS_CAPACITY = 100.0
const METAB_REACTIONS   = ["r0","r1","r2","r3"]

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

# Metabolic throughput fluxes [r0,r1,r2,r3] read out of the BST truth.
function truth_metabolic_fluxes(truth)
    idx = Dict(truth.reactions .=> eachindex(truth.reactions))
    return [truth.vss[idx[r]] for r in METAB_REACTIONS]
end

# Solve the constraint-based problem; `expression`/`activity` gate each factor on
# the committed-step bound. Every other bound is identical between runs.
function feedback_fba(gw; expression::Bool, activity::Bool)
    n   = length(METAB_REACTIONS)
    idx = Dict(METAB_REACTIONS .=> eachindex(METAB_REACTIONS))
    ub  = fill(GENEROUS_CAPACITY, n)
    ub[idx["r0"]] = gw.Vmax0 * (expression ? gw.e_e0 : 1.0) * (activity ? gw.θ : 1.0)
    lb  = zeros(n)                             # all irreversible

    m = Model(HiGHS.Optimizer); set_silent(m)
    @variable(m, lb[i] <= v[i=1:n] <= ub[i])
    @constraint(m, S_FEEDBACK * v .== 0)
    @objective(m, Max, v[idx["r3"]])           # maximize export = throughput
    optimize!(m)
    return value.(v)
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `julia --project=code code/feedback/test_dual_feedback.jl`
Expected: PASS — Task 3 line prints `naive≈10 expr≈4.64 act≈5.74 both≈2.66 truth≈2.66`.

- [ ] **Step 5: Commit**

```bash
git add code/feedback/dual_feedback.jl code/feedback/test_dual_feedback.jl
git commit -m "$(cat <<'EOF'
Recover the dual-feedback truth through the product bound Vmax*(e/e0)*theta

Add the constraint-based layer: gateway_factors computes Vmax0, theta,
and (e/e0) from the measured state (never from vss), and feedback_fba
runs the four escalating cases. Only both gateways open reaches the
truth; each alone leaves a strictly bracketed, distinct gap.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 4: Driver — CSV, figure, prose numbers

**Files:**
- Rewrite: `code/feedback/run_feedback.jl`
- Output: `code/data/feedback_fba.csv`, `code/figs/feedback_gateway.pdf`
- Copy: `chapter/figures/feedback_gateway.pdf`

**Interfaces:**
- Consumes: all of `dual_feedback.jl`; `datapath`, `figpath`, `CSV`, `DataFrame`, `CairoMakie`.
- Produces: `feedback_fba.csv` with columns `reaction, truth, naive, expression_only, activity_only, both_open`; the two-panel `feedback_gateway.pdf`; stdout prose numbers.

- [ ] **Step 1: Rewrite `run_feedback.jl`**

Replace the entire contents of `code/feedback/run_feedback.jl`:

```julia
include(joinpath(@__DIR__, "..", "Include.jl"))
include(joinpath(@__DIR__, "dual_feedback.jl"))

# --------------------------------------------------------------------------- #
# Example 2 driver: integrate the dual-feedback truth, recover it through the
# committed-step bound in four escalating cases, and emit the CSV + figure.
# --------------------------------------------------------------------------- #
t  = feedback_truth()
gw = gateway_factors(t)

vtruth  = truth_metabolic_fluxes(t)
v_naive = feedback_fba(gw; expression=false, activity=false)
v_expr  = feedback_fba(gw; expression=true,  activity=false)
v_act   = feedback_fba(gw; expression=false, activity=true)
v_both  = feedback_fba(gw; expression=true,  activity=true)

df = DataFrame(reaction=METAB_REACTIONS, truth=vtruth, naive=v_naive,
               expression_only=v_expr, activity_only=v_act, both_open=v_both)
CSV.write(datapath("feedback_fba.csv"), df)

# ---- prose numbers -------------------------------------------------------- #
sp = t.seq
X3 = t.Xss[findfirst(==("X3"), t.species)]
T  = vtruth[findfirst(==("r3"), METAB_REACTIONS)]
ledger = [gw.Vmax0, gw.Vmax0*gw.e_e0, gw.Vmax0*gw.θ, gw.Vmax0*gw.e_e0*gw.θ]
println("fixed point: X3*=", round(X3,digits=3), "  T*=", round(T,digits=3),
        "  theta*=", round(gw.θ,digits=3), "  (e/e0)*=", round(gw.e_e0,digits=3))
println("overshoot factor Vmax0/T* = ", round(gw.Vmax0/T,digits=3))
println("ledger [Vmax0, x(e/e0), x theta, both] = ", round.(ledger,digits=3))
println("mu fraction of clearance: transcript=", round(sp.mu/(sp.theta_m+sp.mu),digits=3),
        "  protein=", round(sp.mu/(sp.theta_p+sp.mu),digits=3))

# ---- figure: [wiring schematic | capacity ledger] ------------------------- #
fig = Figure(size=(940, 380))

# panel (a): dual-feedback wiring schematic
axS = Axis(fig[1,1], title="Dual end-product control of the committed enzyme")
hidedecorations!(axS); hidespines!(axS); xlims!(axS, 0, 10); ylims!(axS, 0, 6)
nodes = Dict("X1"=>(2.0,1.2), "X2"=>(4.0,1.2), "X3"=>(6.0,1.2),
             "gene"=>(6.0,4.6), "m"=>(4.0,4.6), "E0"=>(2.0,4.6))
for (lbl,(x,y)) in nodes
    scatter!(axS, [x],[y]; markersize=34, color=(:steelblue,0.20), strokecolor=:steelblue, strokewidth=1.5)
    text!(axS, x, y; text=lbl, align=(:center,:center), fontsize=13)
end
# metabolic chain (solid) and expression cascade (solid, rightward gene->m->E0)
arrows!(axS, [0.5,2.3,4.3,6.3], [1.2,1.2,1.2,1.2], [1.1,1.4,1.4,1.1], [0,0,0,0]; color=:black)  # r0 in, X1->X2, X2->X3, r3 out
arrows!(axS, [5.6,3.6],[4.6,4.6],[-1.2,-1.2],[0,0]; color=:seagreen)                             # gene->m, m->E0
arrows!(axS, [2.0],[4.1],[0.0],[-2.2]; color=:seagreen)                                          # E0 catalyzes r0 (down)
text!(axS, 1.0,1.55; text="r0", fontsize=11); text!(axS, 7.1,1.55; text="r3", fontsize=11)
text!(axS, 2.35,3.0; text="catalysis", fontsize=10, color=:seagreen, align=(:left,:center))
# two repression arrows from X3 (dashed, crimson): slow expression + fast activity
lines!(axS, [6.0,6.0],[1.6,4.2]; color=:crimson, linestyle=:dash)                               # X3 -| gene (expression)
lines!(axS, [5.7,2.3],[1.5,4.2]; color=:crimson, linestyle=:dash)                               # X3 -| E0   (activity)
text!(axS, 6.15,3.0; text="(e/e0)\nexpression, slow", color=:crimson, fontsize=10, align=(:left,:center))
text!(axS, 3.7,2.4; text="theta\nactivity, fast", color=:crimson, fontsize=10, align=(:left,:center))

# panel (b): capacity ledger
axL = Axis(fig[1,2], xticks=(1:4, ["Vmax0","x(e/e0)","x theta","both"]),
           ylabel="committed-step capacity (AU)",
           title="Capacity ledger: the two gateways multiply")
barplot!(axL, 1:4, ledger; color=[:gray70,:steelblue,:orange,:seagreen])
hlines!(axL, [gw.Vmax0]; color=:gray, linestyle=:dash, label="un-inhibited capacity")
hlines!(axL, [T];        color=:crimson, linestyle=:dot, label="BST truth T*")
axislegend(axL; position=:rt, framevisible=false)

save(figpath("feedback_gateway.pdf"), fig)
println("run_feedback OK: wrote feedback_fba.csv and feedback_gateway.pdf")
```

**Fallback (only if the panel-(a) schematic renders poorly):** delete the `axS` block and place the ledger in `fig[1,1]`; then in Task 5 add a TikZ wiring diagram in the `.tex` instead. The figure remains `feedback_gateway.pdf`.

- [ ] **Step 2: Run the driver**

Run: `julia --project=code code/feedback/run_feedback.jl`
Expected: prints the fixed point (`X3*≈4.0 T*≈2.66 theta*≈0.574 (e/e0)*≈0.464`), overshoot `≈3.8`, ledger `≈[10, 4.64, 5.74, 2.66]`, mu fractions `≈[0.04, 0.37]`, and `run_feedback OK`.

- [ ] **Step 3: Verify the artifacts**

Run: `julia --project=code -e 'using CSV,DataFrames; show(CSV.read("code/data/feedback_fba.csv", DataFrame)); println(); @assert isfile("code/figs/feedback_gateway.pdf")'`
Expected: a 4-row table (`r0,r1,r2,r3`) with `both_open ≈ truth ≈ 2.66`, `naive ≈ 10`; and the assert passes (PDF exists).

- [ ] **Step 4: Copy the figure into the chapter**

Run: `cp code/figs/feedback_gateway.pdf chapter/figures/feedback_gateway.pdf`

- [ ] **Step 5: Commit**

```bash
git add code/feedback/run_feedback.jl code/data/feedback_fba.csv code/figs/feedback_gateway.pdf chapter/figures/feedback_gateway.pdf
git commit -m "$(cat <<'EOF'
Emit dual-feedback CSV and the two-panel capacity-ledger figure

Driver integrates the truth, runs the four FBA cases, writes
feedback_fba.csv, and renders feedback_gateway.pdf: the dual-feedback
wiring schematic and the capacity ledger showing Vmax0 whittled by
(e/e0) then theta down onto the BST truth.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 5: Rewrite the section prose

**Files:**
- Rewrite: `chapter/sections/example_feedback.tex`
- Verify: `chapter/Chapter.tex` builds via `chapter/Makefile`

**Interfaces:**
- Consumes: `feedback_gateway.pdf` (Task 4); the labels `eq:general-bound`, `eq:fba-lp`, `eq:mrna`, `eq:protein`, `eq:expression-ss`, `eq:partition`, `eq:control-split`, `eq:mm` from `gateways.tex`; keeps `\label{sec:feedback}` and `\label{fig:feedback}`.
- Produces: the rewritten section; no new preamble packages required.

- [ ] **Step 1: Confirm the numbers to quote**

Re-read the Task 4 stdout. The prose must quote only: exponents `a=0.4`/`b=0.6`; sequence inputs (330 aa, 990 nt, 60 nt/s, 16 aa/s, half-lives 2.5/35 min, doubling 60 min); the derived `theta_m+mu`, `theta_p+mu` and the two mu-fractions (~4%, ~37%); the fixed point (`X3*≈4`, `T*≈2.66`, `theta*≈0.574`, `(e/e0)*≈0.464`); the ledger `10 → 4.64 → 5.74 → 2.66`; overshoot `≈3.8x`. Every capacity/flux number ties to Fig.~\ref{fig:feedback}; the exponents, lengths, rates, and half-lives are structural setup numbers.

- [ ] **Step 2: Rewrite `example_feedback.tex`**

Replace the file with five long paragraphs following this arc. Write real, continuous prose (no bullet lists, no subsection headings, no em dashes, no self-references). Include the display equations shown below verbatim.

Paragraph 1 — dual control motivation. The urea example opened the thermodynamic and kinetic gateways; the regulatory factor `theta_j` was exercised by the first feedback pass but the expression ratio `(e/e0)` was left at unity. A committed biosynthetic step is controlled by its own end product at two levels and on two timescales: fast allosteric inhibition of the enzyme's activity, and slow transcriptional repression of the enzyme's gene, the transcriptional loop being an Alon negative-feedback motif. Regulated as the trp operon is. Keep the abstract chain and open both gateways at once on the single committed step `r0`. Introduce the chain with the existing display:

```latex
\begin{equation}\label{eq:feedback-network}
  \emptyset \xrightarrow{\,r_0\,} X_1 \xrightarrow{\,r_1\,} X_2
  \xrightarrow{\,r_2\,} X_3 \xrightarrow{\,r_3\,} \emptyset,
\end{equation}
```

Paragraph 2 — the extended BST truth. The committed enzyme is now modeled explicitly: its transcript `m` and the enzyme `E0` are dynamic species alongside `X1,X2,X3` in one BSTModelKit S-system, integrated to a self-consistent steady state with BSTModelKit.jl \cite{Vadhin2026}. `X3` enters twice: as a negative kinetic order `-a` in the rate of `r0` (allosteric activity) and as a negative kinetic order `-b` in the rate of transcription (expression repression). State that power-law repression is the concrete operating-point instance of the partition-function control of Eqs.~\eqref{eq:partition}-\eqref{eq:control-split}, the same reduction the chapter already makes for `theta`. State the modeling convention: the metabolic balances keep `S v = 0`, while the `+mu` term lives in the expression balances of Eqs.~\eqref{eq:mrna}-\eqref{eq:protein}. Give the committed-step rate:

```latex
\begin{equation}\label{eq:feedback-rate}
  v_{r_0} = V^{\circ}_{\max,r_0}\,(e/e^{\circ})\,\theta,
  \qquad
  \theta = X_3^{-a},\quad a = 0.4,\qquad
  \text{transcription} \propto X_3^{-b},\quad b = 0.6.
\end{equation}
```

Paragraph 3 — sequence-specific capacities and the mu signature. The transcription and translation capacities `r_X`, `r_L` follow from the gene and protein the enzyme is made from, not from a lookup: a 330-residue enzyme, a 990-nucleotide gene, transcription and translation elongating at 60 nt/s and 16 aa/s. Together with mRNA and protein half-lives (2.5 and 35 min) and a 60-min doubling time these set the residence times in the steady states of Eqs.~\eqref{eq:mrna}-\eqref{eq:protein}, giving `(e/e0)` through Eq.~\eqref{eq:expression-ss}. The mu payoff, stated so it cannot be misread: mu is one number, the culture's specific growth rate, identical in both denominators; only each species' own turnover differs. The same mu is about four percent of the fast-degrading transcript's clearance `theta_m + mu` but more than a third of the stable enzyme's clearance `theta_p + mu`, so for the enzyme growth dilution is a dominant clearance term. This is why the chapter carries mu in every intracellular balance.

Paragraph 4 — coupled steady state and the control scheme. Because the chain is linear every metabolic flux equals the throughput `T`, and `r3` first-order pins `X3 = T/k3`. Substituting the two controls into the committed-step rate gives a single self-consistency condition in `X3`, whose solution the integration confirms: `X3* ≈ 4`, `theta* ≈ 0.574`, `(e/e0)* ≈ 0.464`, throughput `T* ≈ 2.66`, less than a third of the un-inhibited capacity. Name the control scheme: a repressor corepressed by `X3` throttles transcription down to a basal leak floor (the `lambda`/`u_dagger` of Eqs.~\eqref{eq:mrna} and~\eqref{eq:control-split}); the strengths `a=0.4`, `b=0.6` were chosen so each gateway is individually load-bearing.

Paragraph 5 — FBA escalation and the ledger payoff. The flux-balance model keeps the chain's stoichiometry and puts all regulatory content on one bound:

```latex
\begin{equation}\label{eq:feedback-bound}
  0 \;\le\; \hat v_{r_0} \;\le\;
  V^{\circ}_{\max,r_0}\,(e/e^{\circ})\,\theta(X_3),
  \qquad V^{\circ}_{\max,r_0} = 10,
\end{equation}
```

the expression and regulatory instance of Eq.~\eqref{eq:general-bound} with the thermodynamic and saturation gateways at their defaults. Both factors are computed from the measured steady state, not copied from the integrated flux. Four solves of Eq.~\eqref{eq:fba-lp} differing only in this bound: feedback-blind (both factors one) overshoots to the un-inhibited capacity of 10; expression alone throttles it to 4.64, activity alone to 5.74, and only both together reach 2.66, matching the BST truth. Nothing about either loop appears in the stoichiometric matrix; the entire dual-control fact lives in the one bound, where the two gateways multiply. Close on the ledger reading of Fig.~\ref{fig:feedback}.

Update the figure caption to describe the two panels: (left) the dual-feedback wiring, `X3` repressing both transcription of `gene(E0)` (slow, `(e/e0)`) and `E0` activity (fast, `theta`); (right) the capacity ledger, `Vmax0 = 10` reduced by `(e/e0)` to 4.64 and by `theta` to 5.74, both together landing on the BST truth `T* ≈ 2.66`, with the un-inhibited capacity as a dashed line.

- [ ] **Step 3: Build the chapter**

Run: `cd chapter && make`
Expected: `Chapter.pdf` regenerates; no `Undefined control sequence`, no unresolved `??` for `fig:feedback`/`sec:feedback`. Check the log:
Run: `grep -Ei "undefined|multiply defined|LaTeX Warning: Reference" chapter/Chapter.log || echo "no reference/undefined warnings"`
Expected: `no reference/undefined warnings` (or only pre-existing unrelated warnings).

- [ ] **Step 4: Style + magic-number pass**

Run the `style-check` skill on `chapter/sections/example_feedback.tex` (em dashes, subsection headings, self-references, line numbers) and the `audit-magic-numbers` skill on the section. Fix any flagged item inline. Optionally run `review-section` for a flow pass and merge any choppy paragraphs (Varner prefers long complete paragraphs).

- [ ] **Step 5: Commit**

```bash
git add chapter/sections/example_feedback.tex
git commit -m "$(cat <<'EOF'
Rewrite feedback example for dual expression + activity control

The committed enzyme is now controlled by its own end product through
two gateways: fast allosteric activity (theta) and slow transcriptional
expression ((e/e0)) computed from sequence-specific capacities. The bound
Vmax*(e/e0)*theta recovers the extended BST truth; each gateway alone
leaves a distinct gap. Foreground the mu signature: the same growth
dilution dominates stable-enzyme clearance but not transcript clearance.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Self-Review

**Spec coverage:**
- Dual control on one enzyme, product bound `Vmax0*(e/e0)*theta` — Tasks 2-3, Eq.~\eqref{eq:feedback-bound}. ✓
- Extended 5-species BST truth, power-law repression as the Hill/partition instance — Task 1 (TOML), Task 2 (integrate). ✓
- Nondimensional expression states, `S v = 0` metabolic / `+mu` expression convention — Task 2 (alpha coefficients tied to `theta+mu`), Paragraph 2. ✓
- Sequence-specific capacities, no precursors in S — Task 2 `sequence_params`, Paragraph 3. ✓
- The mu signature (~4% vs ~37%, same mu) — Task 2 assert + Paragraph 3 with the explicit "one mu" framing. ✓
- Control scheme: repressor + basal leak `lambda` (`PHI`) — Task 2 `rTXb`, Paragraph 4. ✓
- Four escalating FBA cases, computed-not-copied factors — Task 3 `gateway_factors`/`feedback_fba`. ✓
- Two-panel figure (schematic + ledger) — Task 4, with TikZ fallback. ✓
- Five long paragraphs, style conventions — Task 5. ✓
- Out of scope honored: no precursor columns in S; `w_j` (translation control) fixed at 1; steady-state only. ✓

**Placeholder scan:** No TBD/TODO; every code step shows complete code; every prose paragraph names its content and equations. Fallbacks are concrete (exact edits), not "handle errors."

**Type consistency:** `feedback_truth()` returns `(species, Xss, reactions, vss, model, seq)`, consumed with those exact names in `gateway_factors`, `truth_metabolic_fluxes`, and the driver. `gateway_factors` returns `(Vmax0, θ, e_e0)`, consumed by `feedback_fba(gw; expression, activity)` and the ledger. `METAB_REACTIONS`, `S_FEEDBACK`, `SeqParams` fields (`theta_m, theta_p, mu, ...`) match across tasks. Reaction order `["rTX","rTXb","rMdeg","rTL","rEdeg","r0","r1","r2","r3"]` and species order `["X1","X2","X3","m","E0"]` are asserted identically in Task 1 test and Task 2 `feedback_truth`.
