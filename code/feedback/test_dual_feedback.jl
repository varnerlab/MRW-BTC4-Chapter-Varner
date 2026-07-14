include(joinpath(@__DIR__, "..", "Include.jl"))

# ---- Task 1: the Dual-Feedback TOML builds with the expected orderings ----- #
let model = build(joinpath(@__DIR__, "Dual-Feedback.toml"))
    @assert model.list_of_dynamic_species == ["X1","X2","X3","m","E0"] "dyn species: $(model.list_of_dynamic_species)"
    @assert model.list_of_static_species == ["E"] "static species: $(model.list_of_static_species)"
    @assert model.list_of_reactions == ["rTX","rTXb","rMdeg","rTL","rEdeg","r0","r1","r2","r3"] "reactions: $(model.list_of_reactions)"
    # metabolic submatrix (rows X1,X2,X3 ; cols r0..r3) is the linear chain
    srow(s) = findfirst(==(s), model.list_of_dynamic_species)
    scol(r) = findfirst(==(r), model.list_of_reactions)
    Smet = model.S[[srow("X1"),srow("X2"),srow("X3")], [scol("r0"),scol("r1"),scol("r2"),scol("r3")]]
    @assert Smet == [1.0 -1.0 0.0 0.0; 0.0 1.0 -1.0 0.0; 0.0 0.0 1.0 -1.0] "metabolic S mismatch:\n$Smet"
    println("test_dual_feedback (Task 1): TOML builds, orderings + metabolic S OK")
end

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

# ---- Task 3: FBA recovers the truth only when BOTH gateways are open ------- #
let t = feedback_truth()
    @assert S_FEEDBACK == [1.0 -1.0 0.0 0.0; 0.0 1.0 -1.0 0.0; 0.0 0.0 1.0 -1.0] "S_FEEDBACK wrong"
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
end

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
