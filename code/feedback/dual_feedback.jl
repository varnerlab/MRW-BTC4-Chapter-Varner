# --------------------------------------------------------------------------- #
# Chapter Example 2: dual end-product control of ONE committed enzyme.
#   activity   (theta) : X3 allosterically inhibits E0             (fast)
#   expression (e/e0)  : X3 transcriptionally represses gene(E0)   (slow)
#
# The "truth" is an extended BSTModelKit S-system over [X1,X2,X3,m,E0]. The FBA
# layer (below) approximates it, to within about 6%, through the product bound
# Vmax0*(e/e0)*theta on r0. Nothing is copied from the integrated fluxes: (e/e0)
# is the MEASURED enzyme abundance (reference e0 = 1 by the nondimensional
# normalization), and theta is a bounded, independent Hill occupancy
# K^n/(K^n+X3^n) evaluated at the measured X3 -- not the truth model's own
# X3^-a allosteric kinetics.
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
    @assert model.G[trow("m"),  idx["rMdeg"]] == 1.0 "m order in rMdeg not default +1"
    @assert model.G[trow("E0"), idx["rEdeg"]] == 1.0 "E0 order in rEdeg not default +1"
    @assert model.G[trow("X1"), idx["r1"]]    == 1.0 "X1 order in r1 not default +1"
    @assert model.G[trow("X2"), idx["r2"]]    == 1.0 "X2 order in r2 not default +1"
    @assert model.G[trow("X3"), idx["r3"]]    == 1.0 "X3 order in r3 not default +1"

    # Stiff system (metabolic alpha ~100 vs expression ~0.03): integrate long
    # enough for the slow protein mode (tau ~ 1/(theta_p+mu) ~ 31 min).
    Xss = steadystate(model; tspan=(0.0, 1000.0))
    vss = reaction_fluxes(model, Xss)
    return (species=species, Xss=Xss, reactions=reactions, vss=vss, model=model, seq=sp)
end

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

# Metabolic throughput fluxes [r0,r1,r2,r3] read out of the BST truth.
function truth_metabolic_fluxes(truth)
    idx = Dict(truth.reactions .=> eachindex(truth.reactions))
    return [truth.vss[idx[r]] for r in METAB_REACTIONS]
end

# Solve the constraint-based problem; `expression`/`activity` gate each factor on
# the committed-step bound. Every other bound is identical between runs.
function feedback_fba(gw; expression::Bool, activity::Bool)
    gw.Vmax0 >= 0 || throw(ArgumentError("Vmax0 must be nonnegative, got $(gw.Vmax0)"))
    0 <= gw.θ    || throw(ArgumentError("θ must be nonnegative, got $(gw.θ)"))
    0 <= gw.e_e0 || throw(ArgumentError("e_e0 must be nonnegative, got $(gw.e_e0)"))

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
