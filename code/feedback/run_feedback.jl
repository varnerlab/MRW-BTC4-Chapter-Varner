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
    model.α = [5.0, 5.0, 5.0, 5.0, 2.0, 1.0, 1.0]               # per-reaction rate constants

    # Feedback via negative kinetic orders. G rows index total_species_list
    # [A,B,C,D,E,E1,E2,E3,E4]; G cols index reactions [r1,r2,r3,r4,r0,r5,r6].
    G = model.G
    G[5, 1] = -0.5   # E (species 5) inhibits r1 (reaction 1)
    G[4, 2] = -0.5   # D (species 4) inhibits r2 (reaction 2)
    model.G = G

    Xss = steadystate(model; tspan=(0.0, 200.0))
    vss = reaction_fluxes(model, Xss)
    return (species=species, Xss=Xss, reactions=reactions, vss=vss)
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
