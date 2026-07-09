include(joinpath(@__DIR__, "..", "Include.jl"))

# --------------------------------------------------------------------------- #
# Example 2a: synthetic "truth" from a BSTModelKit branched-feedback S-system.
#
# A branched power-law (generalized-mass-action) pathway with product-feedback
# inhibition is driven to steady state. We return both the steady-state
# concentrations (Xss) and the per-reaction steady-state fluxes (vss). Those
# fluxes are the "truth" the FBA layer (Task 5) must reproduce.
#
# Reaction order (connection-record order): [r1, r2, r3, r4, r0, r5, r6]
# Species order:                            [A, B, C, D, E]
# --------------------------------------------------------------------------- #

# Per-reaction steady-state fluxes.
#
# BSTModelKit is generalized-mass-action under the hood: the ODE right-hand
# side is S * v, where the per-reaction rate vector v is computed by the
# package's own power-law kernel `_powerlaw(state, α, G)` with
#   v[j] = α[j] * prod_i state[i]^G[i,j].
# The kinetic state that the balances see is the dynamic species stacked on
# top of the static factors, i.e. state = vcat(x, static_factors_array)
# (see BSTModelKit/src/Balances.jl). We reuse that exact kernel here so the
# extracted fluxes are, by construction, identical to what the solver used.
function reaction_fluxes(model::BSTModel, X::Vector{Float64})::Vector{Float64}
    state_array = vcat(X, model.static_factors_array)   # [A,B,C,D,E, E1,E2,E3,E4]
    return BSTModelKit._powerlaw(state_array, model.α, model.G)
end

# Build the branched-feedback S-system and return steady-state concentrations + fluxes.
function feedback_truth()
    path  = joinpath(@__DIR__, "Branched-Feedback.toml")
    model = build(path)

    # reaction order (connection records): [r1, r2, r3, r4, r0, r5, r6]
    reactions = ["r1","r2","r3","r4","r0","r5","r6"]
    species   = ["A","B","C","D","E"]

    model.initial_condition_array = [1.0, 0.1, 0.1, 0.1, 0.1]   # A,B,C,D,E
    model.static_factors_array    = [1.0, 1.0, 1.0, 1.0]        # E1,E2,E3,E4

    # Per-reaction rate constants. The export constants for r5,r6 are set below
    # unity so the feedback end-products D,E accumulate to a steady state ABOVE
    # one (here D=E=2). That is deliberate: the feedback control factor is
    # θ=P^{-1/2}, so with D=E=1 (unit export) it would collapse to θ=1 and the
    # regulation would be numerically invisible — the un-inhibited capacity would
    # already equal the true flux. With D=E=2 the control factor is a genuine
    # sub-unity term (θ=1/√2≈0.71) and the un-inhibited committed-step capacity
    # (≈2.83) sits strictly above the inhibited flux (2.0). That gap is exactly
    # what lets the gateway-closed FBA in Example 2b overshoot when it ignores
    # the feedback. Throughput T=2 and the branch split 1/1 are unchanged.
    model.α = [5.0, 5.0, 5.0, 5.0, 2.0, 0.5, 0.5]               # per-reaction rate constants (r5,r6 export < 1)

    # Feedback via negative kinetic orders. G rows index model.total_species_list
    # (= vcat(dynamic, static) = [A,B,C,D,E,E1,E2,E3,E4]); G cols index
    # model.list_of_reactions (connection-record order = [r1,r2,r3,r4,r0,r5,r6]).
    # BSTModelKit's Factory.jl builds G as
    # `_build_exponent_matrix(total_species_list, sorted_rate_dict_array)` and
    # exposes both orderings on the model object itself, so we look up the
    # feedback indices BY NAME rather than hardcoding row/column numbers, and we
    # assert those orderings match what the rest of this file assumes (species,
    # reactions) so a package change would fail loudly instead of silently
    # mistargeting the feedback.
    @assert model.list_of_dynamic_species == species "BSTModelKit dynamic-species order changed: expected $(species), got $(model.list_of_dynamic_species)"
    @assert model.list_of_reactions == reactions "BSTModelKit reaction order changed: expected $(reactions), got $(model.list_of_reactions)"

    row_E  = findfirst(==("E"),  model.total_species_list)
    col_r1 = findfirst(==("r1"), model.list_of_reactions)
    row_D  = findfirst(==("D"),  model.total_species_list)
    col_r2 = findfirst(==("r2"), model.list_of_reactions)

    G = model.G
    G[row_E, col_r1] = -0.5   # E inhibits r1
    G[row_D, col_r2] = -0.5   # D inhibits r2
    model.G = G

    println("feedback: E->r1, D->r2 via G[$(row_E),$(col_r1)], G[$(row_D),$(col_r2)]")

    Xss = steadystate(model; tspan=(0.0, 200.0))
    vss = reaction_fluxes(model, Xss)
    # The built model (α, G, static factors, orderings) is returned so the FBA
    # layer below can evaluate the bound models at the measured Xss with the
    # exact same power-law kernel — never reading fluxes back from vss.
    return (species=species, Xss=Xss, reactions=reactions, vss=vss, model=model)
end

# --------------------------------------------------------------------------- #
# Self-check (TDD gate): the extracted fluxes must satisfy the pathway's
# mass balance at steady state, and the branch point must stay split.
# --------------------------------------------------------------------------- #
let t = feedback_truth()
    A,B,C,D,E = t.Xss
    v = Dict(t.reactions .=> t.vss)
    @assert isapprox(v["r0"], v["r1"]; rtol=1e-3)
    @assert isapprox(v["r1"], v["r2"]; rtol=1e-3)
    @assert isapprox(v["r3"] + v["r4"], v["r2"]; rtol=1e-3)
    @assert isapprox(v["r3"], v["r5"]; rtol=1e-3)
    @assert isapprox(v["r4"], v["r6"]; rtol=1e-3)
    @assert v["r3"] > 0 && v["r4"] > 0   # feedback keeps both branches active
    println("BST truth OK: T=", v["r1"], " split D/E=", v["r3"], "/", v["r4"])
end

# =========================================================================== #
# Example 2b: the SAME branched pathway as a flux-balance (constraint-based)
# model. The point of the section: product feedback is real, and in
# constraint-based modeling it does NOT live in the objective — it enters
# through the flux BOUND.
#
#   * gateway CLOSED: a naive FBA that ignores the metabolic state and puts a
#     generous constant capacity on every internal reaction. Feedback and
#     kinetic saturation are both invisible, so throughput overshoots.
#   * gateway OPEN:   the committed-step bounds are informed by the MEASURED
#     steady-state concentrations (truth.Xss) through the bound models below,
#     and the FBA optimum then reproduces the BST "truth" of Example 2a.
#
# Nothing is copied from truth.vss: every open bound is COMPUTED from Xss with
# the same power-law kernel Example 2a used.
#
# Stoichiometric matrix S (species x reactions). Rows [A,B,C,D,E];
# columns [r1,r2,r3,r4,r0,r5,r6] (same reaction order as feedback_truth()).
# =========================================================================== #
const S_FEEDBACK = [
    -1  0  0  0  1  0  0    # A:  r0 in,  r1 out
     1 -1  0  0  0  0  0    # B:  r1 in,  r2 out
     0  1 -1 -1  0  0  0    # C:  r2 in,  r3+r4 out
     0  0  1  0  0 -1  0    # D:  r3 in,  r5 out
     0  0  0  1  0  0 -1    # E:  r4 in,  r6 out
]

# Generous constant capacity for the gateway-CLOSED bounds: every internal
# reaction may run this fast (all regulation/saturation correction factors = 1).
# It sits far above the true throughput (T=2), so the closed FBA overshoots.
const CLOSED_CAPACITY = 10.0

# --------------------------------------------------------------------------- #
# Bound models — how the measured metabolic state sets an FBA flux bound.
#
# Capacity / regulation decomposition for a feedback-controlled committed step:
#
#     ub_j = Vmax_j * θ_j(inhibitor)
#
#   Vmax_j : the UN-inhibited capacity — the reaction's power-law rate with its
#            feedback kinetic-order zeroed, evaluated at the measured state.
#   θ_j    : the regulation control factor ∈ (0,1], = (rate WITH feedback) /
#            (rate WITHOUT feedback) at the measured inhibitor concentration.
#            This is the partition-function / Hill feedback term of §4; here it
#            evaluates to E^{-1/2} for r1 and D^{-1/2} for r2 (≈0.71 at Xss).
#
# Both factors are produced by the SAME `_powerlaw` kernel Example 2a used,
# evaluated at the MEASURED Xss — never read back from truth.vss.
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

# Kinetic gateway for an UN-regulated reaction: its power-law rate at Xss. For
# the branch steps r3,r4 this caps each branch at the measured C-rate, which
# pins the C-split so the recovered branch fluxes are unique.
kinetic_gateway(model::BSTModel, Xss::Vector{Float64}, col::Int) =
    reaction_fluxes(model, Xss)[col]

# --------------------------------------------------------------------------- #
# feedback_fba(truth; gateway) — solve the constraint-based problem and return a
# flux vector aligned to truth.reactions. `gateway=false` uses the naive
# (closed) constant bounds; `gateway=true` opens the regulation gateway on the
# committed steps r1,r2 and the kinetic gateway on the branches r3,r4.
# --------------------------------------------------------------------------- #
function feedback_fba(truth; gateway::Bool)
    model = truth.model
    rxn   = truth.reactions
    idx   = Dict(rxn .=> eachindex(rxn))
    n     = length(rxn)

    lb = zeros(n)                     # all reactions irreversible here
    ub = fill(CLOSED_CAPACITY, n)     # closed default: generous constant on every reaction

    if gateway
        # feedback locations looked up BY NAME (robust to package reordering)
        row_E = findfirst(==("E"), model.total_species_list)
        row_D = findfirst(==("D"), model.total_species_list)

        # committed steps r1, r2: bound = capacity * control, from measured E, D
        g1 = capacity_and_control(model, truth.Xss, idx["r1"], row_E)
        g2 = capacity_and_control(model, truth.Xss, idx["r2"], row_D)
        ub[idx["r1"]] = g1.Vmax * g1.θ
        ub[idx["r2"]] = g2.Vmax * g2.θ

        # branch steps r3, r4: kinetic gateway at measured C pins the split
        ub[idx["r3"]] = kinetic_gateway(model, truth.Xss, idx["r3"])
        ub[idx["r4"]] = kinetic_gateway(model, truth.Xss, idx["r4"])

        # boundary/export r0, r5, r6 stay generously bounded — mass balance
        # (S v = 0) fixes them from the internal fluxes.
    end

    m = Model(HiGHS.Optimizer); set_silent(m)
    @variable(m, lb[i] <= v[i=1:n] <= ub[i])
    @constraint(m, S_FEEDBACK * v .== 0)
    @objective(m, Max, v[idx["r5"]] + v[idx["r6"]])   # maximise TOTAL export
    optimize!(m)
    return value.(v)
end

# --------------------------------------------------------------------------- #
# Comparison, artifacts, and the TDD gate: gateway-open must recover the BST
# truth; gateway-closed must overshoot it.
# --------------------------------------------------------------------------- #
let t = feedback_truth()
    v_closed = feedback_fba(t; gateway=false)
    v_open   = feedback_fba(t; gateway=true)

    df = DataFrame(reaction=t.reactions, truth=t.vss,
                   gateway_closed=v_closed, gateway_open=v_open)
    CSV.write(datapath("feedback_fba.csv"), df)

    err_open   = maximum(abs.(v_open   .- t.vss))
    err_closed = maximum(abs.(v_closed .- t.vss))
    println("gateway-open   vs truth: max |Δflux| = ", err_open)
    println("gateway-closed vs truth: max |Δflux| = ", err_closed)

    # gateway-open recovers the BST truth; gateway-closed diverges (overshoots)
    @assert err_open   < 1e-2 "gateway-open should recover the BST truth (got $err_open)"
    @assert err_closed > 1e-1 "gateway-closed should overshoot the truth (got $err_closed)"

    fig = Figure(size=(760,320))
    ax  = Axis(fig[1,1], xticks=(1:length(t.reactions), t.reactions),
               ylabel="flux (AU)", title="Feedback recovered through the bound")
    w = 0.25
    barplot!(ax, (1:length(t.reactions)) .- w, t.vss;    width=w, label="BST truth")
    barplot!(ax, (1:length(t.reactions)),      v_closed; width=w, label="FBA, gateway closed")
    barplot!(ax, (1:length(t.reactions)) .+ w, v_open;   width=w, label="FBA, gateway open")
    Legend(fig[1,2], ax)
    save(figpath("feedback_gateway.pdf"), fig)
    println("feedback FBA OK")
end
