include(joinpath(@__DIR__, "..", "Include.jl"))

# --------------------------------------------------------------------------- #
# Example 2a: synthetic "truth" from a BSTModelKit LINEAR feedback S-system.
#
# A linear power-law (generalized-mass-action) chain
#
#     r0 : {} -> X1     uptake / input, INHIBITED by the end product X3
#     r1 : X1 -> X2     first-order in X1
#     r2 : X2 -> X3     first-order in X2
#     r3 : X3 -> {}     export, first-order in X3
#
# is driven to steady state. Because the chain is linear, at steady state every
# flux is equal to the throughput T (r0=r1=r2=r3=T) — there is no branch split.
# The single control point is the input reaction r0, whose rate carries a
# NEGATIVE kinetic order on the end product X3: the more X3 accumulates, the
# harder the cell throttles uptake. That is textbook end-product inhibition.
#
# We return both the steady-state concentrations (Xss) and the per-reaction
# steady-state fluxes (vss). Those fluxes are the "truth" the FBA layer
# (Example 2b) must reproduce THROUGH THE BOUND, never by copying vss.
#
# Reaction order (connection-record order): [r0, r1, r2, r3]
# Species order:                            [X1, X2, X3]
# --------------------------------------------------------------------------- #

# Per-reaction steady-state fluxes.
#
# BSTModelKit is generalized-mass-action under the hood: the ODE right-hand
# side is S * v, where the per-reaction rate vector v is computed by the
# package's own power-law kernel `_powerlaw(state, α, G)` with
#   v[j] = α[j] * prod_i state[i]^G[i,j].
# The kinetic state the balances see is the dynamic species stacked on top of
# the static factors, i.e. state = vcat(x, static_factors_array) (see
# BSTModelKit/src/Balances.jl). We reuse that exact kernel here so the extracted
# fluxes are, by construction, identical to what the solver integrated.
function reaction_fluxes(model::BSTModel, X::Vector{Float64})::Vector{Float64}
    state_array = vcat(X, model.static_factors_array)   # [X1,X2,X3, E(static)]
    return BSTModelKit._powerlaw(state_array, model.α, model.G)
end

# Build the linear feedback S-system and return steady-state concentrations + fluxes.
function feedback_truth()
    path  = joinpath(@__DIR__, "Linear-Feedback.toml")
    model = build(path)

    # reaction order (connection records): [r0, r1, r2, r3]
    reactions = ["r0","r1","r2","r3"]
    species   = ["X1","X2","X3"]

    model.initial_condition_array = [1.0, 1.0, 1.0]   # X1,X2,X3 (X3>0 keeps X3^-1/2 finite)
    model.static_factors_array    = [1.0]             # single unused static enzyme E

    # Per-reaction rate constants.
    #   α[r0]=10  : un-inhibited uptake capacity Vmax.
    #   α[r1]=α[r2]=100 : fast internal steps so they never limit throughput.
    #   α[r3]=1   : export. With export order 1 the self-consistent steady state
    #               is X3 = T and T = 10·X3^(-1/2), i.e. T = 10^(2/3) ≈ 4.64,
    #               X3 ≈ 4.64, θ = X3^(-1/2) ≈ 0.46, capacity Vmax = 10. The gap
    #               between Vmax (10) and the true throughput (4.64) is exactly
    #               what lets a feedback-blind FBA overshoot in Example 2b.
    model.α = [10.0, 100.0, 100.0, 1.0]

    # Feedback via a NEGATIVE kinetic order. G rows index model.total_species_list
    # (= vcat(dynamic, static) = [X1,X2,X3]); G cols index model.list_of_reactions
    # (connection-record order = [r0,r1,r2,r3]). BSTModelKit's Factory.jl builds G
    # as `_build_exponent_matrix(total_species_list, sorted_rate_dict_array)` and
    # exposes both orderings on the model object, so we look up the feedback entry
    # BY NAME (findfirst) rather than hardcoding indices, and assert the orderings
    # match what the rest of this file assumes so a package change fails loudly
    # instead of silently mistargeting the feedback.
    @assert model.list_of_dynamic_species == species "BSTModelKit dynamic-species order changed: expected $(species), got $(model.list_of_dynamic_species)"
    @assert model.list_of_reactions == reactions "BSTModelKit reaction order changed: expected $(reactions), got $(model.list_of_reactions)"

    row_X3 = findfirst(==("X3"), model.total_species_list)
    col_r0 = findfirst(==("r0"), model.list_of_reactions)

    G = model.G
    G[row_X3, col_r0] = -0.5   # end product X3 inhibits the input r0
    model.G = G

    println("feedback: X3 -| r0 via G[$(row_X3),$(col_r0)] = ", model.G[row_X3, col_r0])

    Xss = steadystate(model; tspan=(0.0, 200.0))
    vss = reaction_fluxes(model, Xss)
    # The built model (α, G, static factors, orderings) is returned so the FBA
    # layer below can evaluate the bound models at the measured Xss with the exact
    # same power-law kernel — never reading fluxes back from vss.
    return (species=species, Xss=Xss, reactions=reactions, vss=vss, model=model)
end

# --------------------------------------------------------------------------- #
# Self-check (TDD gate): the extracted fluxes must satisfy the linear chain's
# mass balance at steady state (all fluxes equal the throughput T), and the end
# product must sit above one so the feedback factor θ = X3^(-1/2) is genuinely
# sub-unity (not an inert regulation term).
# --------------------------------------------------------------------------- #
let t = feedback_truth()
    X1,X2,X3 = t.Xss
    v = Dict(t.reactions .=> t.vss)
    @assert isapprox(v["r0"], v["r1"]; rtol=1e-3)
    @assert isapprox(v["r1"], v["r2"]; rtol=1e-3)
    @assert isapprox(v["r2"], v["r3"]; rtol=1e-3)
    @assert X3 > 1.0   # end product accumulates ⇒ θ = X3^(-1/2) < 1 (feedback live)
    T = v["r3"]
    θ = X3^(-0.5)
    println("BST truth OK: T = ", T, "  X3 = ", X3, "  θ = X3^(-1/2) = ", θ)
end

# =========================================================================== #
# Example 2b: the SAME linear pathway as a flux-balance (constraint-based)
# model. The point of the section: product feedback is real, and in
# constraint-based modeling it does NOT live in the objective — it enters
# through the flux BOUND on the regulated input.
#
#   * gateway CLOSED (feedback OFF): the input is bounded at its UN-inhibited
#     capacity Vmax (θ forced to 1). Regulation is invisible, so throughput
#     overshoots up to the uptake capacity (≈10).
#   * gateway OPEN  (regulation ON): the input bound is Vmax·θ(X3_measured),
#     the capacity throttled by the measured end product. The FBA optimum then
#     reproduces the BST "truth" of Example 2a (≈4.64).
#
# The ONLY difference between the two scenarios is the feedback factor θ on the
# input bound (closed uses θ=1, open uses θ measured from Xss). Nothing is copied
# from truth.vss: both Vmax and θ are COMPUTED from Xss with the same power-law
# kernel Example 2a used. That makes the demonstration load-bearing — turning θ
# back to 1 IS the closed case, and it provably changes the answer.
#
# Stoichiometric matrix S (species x reactions). Rows [X1,X2,X3];
# columns [r0,r1,r2,r3] (same reaction order as feedback_truth()).
# =========================================================================== #
const S_FEEDBACK = [
     1 -1  0  0    # X1:  r0 in,  r1 out
     0  1 -1  0    # X2:  r1 in,  r2 out
     0  0  1 -1    # X3:  r2 in,  r3 out
]

# Generous constant capacity for the NON-regulated reactions. The linear chain
# forces every flux to equal r0's, so any bound comfortably above the uptake
# capacity (10) leaves the input reaction as the sole throughput determinant.
const GENEROUS_CAPACITY = 100.0

# --------------------------------------------------------------------------- #
# Bound model — how the measured metabolic state sets the input flux bound.
#
# Capacity / regulation decomposition for the feedback-controlled input step:
#
#     ub[r0] = Vmax * θ(X3)
#
#   Vmax : the UN-inhibited uptake capacity — r0's power-law rate with its
#          feedback kinetic-order zeroed, evaluated at the measured state (⇒ α[r0]).
#   θ    : the regulation control factor ∈ (0,1], = (rate WITH feedback) /
#          (rate WITHOUT feedback) at the measured end product (⇒ X3^(-1/2)).
#
# Both factors come from the SAME `_powerlaw` kernel Example 2a used, evaluated
# at the MEASURED Xss — never read back from truth.vss.
# --------------------------------------------------------------------------- #

# (Vmax, θ) for reaction column `col` whose feedback kinetic order lives at
# G[fb_row, col].
function capacity_and_control(model::BSTModel, Xss::Vector{Float64},
                              col::Int, fb_row::Int)
    state     = vcat(Xss, model.static_factors_array)
    rate_with = BSTModelKit._powerlaw(state, model.α, model.G)[col]   # feedback active
    G0 = copy(model.G); G0[fb_row, col] = 0.0                         # remove regulation
    Vmax = BSTModelKit._powerlaw(state, model.α, G0)[col]             # un-inhibited capacity
    return (Vmax = Vmax, θ = rate_with / Vmax)
end

# --------------------------------------------------------------------------- #
# feedback_fba(truth; gateway) — solve the constraint-based problem and return a
# flux vector aligned to truth.reactions. The input bound is Vmax·θ; `gateway`
# selects θ:
#   gateway=false (CLOSED): θ = 1        → input at un-inhibited capacity Vmax.
#   gateway=true  (OPEN):   θ = θ(Xss)   → input throttled by measured end product.
# Every other bound is identical between the two runs, so the two scenarios
# differ ONLY in the feedback factor.
# --------------------------------------------------------------------------- #
function feedback_fba(truth; gateway::Bool)
    model = truth.model
    rxn   = truth.reactions
    idx   = Dict(rxn .=> eachindex(rxn))
    n     = length(rxn)

    # feedback location looked up BY NAME (robust to package reordering)
    row_X3 = findfirst(==("X3"), model.total_species_list)
    g0 = capacity_and_control(model, truth.Xss, idx["r0"], row_X3)

    # closed → θ=1 (regulation blind); open → measured θ (regulation on)
    θ = gateway ? g0.θ : 1.0

    lb = zeros(n)                          # all reactions irreversible here
    ub = fill(GENEROUS_CAPACITY, n)        # non-regulated steps: generous
    ub[idx["r0"]] = g0.Vmax * θ            # the SOLE difference between scenarios

    m = Model(HiGHS.Optimizer); set_silent(m)
    @variable(m, lb[i] <= v[i=1:n] <= ub[i])
    @constraint(m, S_FEEDBACK * v .== 0)
    @objective(m, Max, v[idx["r3"]])       # maximise export = throughput
    optimize!(m)
    return value.(v)
end

# --------------------------------------------------------------------------- #
# Comparison, artifacts, and the TDD gate: gateway-open must recover the BST
# truth; gateway-closed must overshoot it; and θ must be genuinely sub-unity
# (guards against an inert-feedback regression).
# --------------------------------------------------------------------------- #
let t = feedback_truth()
    # sanity: the hand-written FBA stoichiometry matches the package's S.
    @assert S_FEEDBACK == t.model.S "S_FEEDBACK disagrees with BSTModelKit S"

    v_closed = feedback_fba(t; gateway=false)
    v_open   = feedback_fba(t; gateway=true)

    idx  = Dict(t.reactions .=> eachindex(t.reactions))
    row_X3 = findfirst(==("X3"), t.model.total_species_list)
    g0   = capacity_and_control(t.model, t.Xss, idx["r0"], row_X3)
    θ    = g0.θ

    df = DataFrame(reaction=t.reactions, truth=t.vss,
                   gateway_closed=v_closed, gateway_open=v_open)
    CSV.write(datapath("feedback_fba.csv"), df)

    err_open  = maximum(abs.(v_open .- t.vss))
    overshoot = v_closed[idx["r3"]] - t.vss[idx["r3"]]
    println("gateway-open   vs truth: max |Δflux| = ", err_open)
    println("gateway-closed overshoot on r3      = ", overshoot,
            "  (Vmax=", g0.Vmax, ", θ=", θ, ")")

    # ACCEPTANCE ASSERTS (the test)
    @assert err_open  < 1e-2 "gateway-open should recover the BST truth (got $err_open)"
    @assert overshoot > 0.5  "gateway-closed should overshoot throughput (got $overshoot)"
    @assert θ         < 0.9  "feedback factor must be genuinely sub-unity (got θ=$θ)"

    fig = Figure(size=(760,320))
    ax  = Axis(fig[1,1], xticks=(1:length(t.reactions), t.reactions),
               ylabel="flux (AU)", title="End-product feedback recovered through the input bound")
    w = 0.25
    barplot!(ax, (1:length(t.reactions)) .- w, t.vss;    width=w, label="BST truth")
    barplot!(ax, (1:length(t.reactions)),      v_closed; width=w, label="FBA, gateway closed (θ=1)")
    barplot!(ax, (1:length(t.reactions)) .+ w, v_open;   width=w, label="FBA, gateway open (θ=X3^-1/2)")
    hlines!(ax, [g0.Vmax]; color=:gray, linestyle=:dash, label="un-inhibited capacity Vmax")
    Legend(fig[1,2], ax)
    save(figpath("feedback_gateway.pdf"), fig)
    println("feedback FBA OK")
end
