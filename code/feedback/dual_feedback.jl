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
