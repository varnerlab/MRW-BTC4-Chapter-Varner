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

# ---- Task 2: the integrated reference reaches the expected fixed point ----- #
let reference = feedback_reference()
    X1, X2, X3, m, E0 = reference.Xss
    idx = Dict(reference.reactions .=> eachindex(reference.reactions))
    vmet = [reference.vss[idx[r]] for r in ["r0","r1","r2","r3"]]
    @assert all(isapprox.(vmet, vmet[1]; rtol=1e-3)) "linear chain not balanced: $vmet"
    T = vmet[1]
    @assert isapprox(X3, 4.0; atol=0.15) "X3* != 4: $X3"
    @assert isapprox(T, 2.66; atol=0.12) "T* != 2.66: $T"
    @assert X3 > 1.0 "end product must accumulate so both controls are sub-unity"
    @assert isapprox(E0, 0.464; atol=0.03) "E0* (=(e/e0)) off target: $E0"

    sp = reference.seq
    @assert isapprox(sp.transcription_time_s, 16.5; atol=1e-9)
    @assert isapprox(sp.translation_time_s, 20.625; atol=1e-9)
    frac_m = sp.mu / (sp.theta_m + sp.mu)
    frac_p = sp.mu / (sp.theta_p + sp.mu)
    @assert frac_m < 0.10 "mu should be a small fraction of transcript clearance: $frac_m"
    @assert frac_p > 0.30 "mu should be a dominant fraction of enzyme clearance: $frac_p"
    println("test_dual_feedback (Task 2): fixed point X3*=$X3 T*=$T ; mu-frac transcript=$frac_m protein=$frac_p")
end

# ---- Task 3: compare the reference with the four FBA configurations -------- #
let reference = feedback_reference()
    @assert S_FEEDBACK == [1.0 -1.0 0.0 0.0; 0.0 1.0 -1.0 0.0; 0.0 0.0 1.0 -1.0] "S_FEEDBACK wrong"
    controls = control_factors(reference)
    @assert isapprox(controls.Vmax0, 10.0; atol=1e-9) "Vmax0: $(controls.Vmax0)"
    @assert isapprox(controls.θ,     0.610; atol=0.01) "theta*: $(controls.θ)"
    @assert isapprox(controls.e_e0,  0.464; atol=0.03) "(e/e0)*: $(controls.e_e0)"
    @assert controls.θ < 0.9 && controls.e_e0 < 0.9 "both controls must be sub-unity"

    v_reference = reference_metabolic_fluxes(reference)
    v_naive = feedback_fba(controls; expression=false, activity=false)
    v_expr  = feedback_fba(controls; expression=true,  activity=false)
    v_act   = feedback_fba(controls; expression=false, activity=true)
    v_both  = feedback_fba(controls; expression=true,  activity=true)

    j = findfirst(==("r3"), METAB_REACTIONS)
    @assert maximum(abs.(v_both .- v_reference) ./ v_reference) < 0.10 "dual-control result must be within 10% of reference: $(v_both) vs $(v_reference)"
    @assert v_naive[j] - v_reference[j] > 0.5 "naive result must exceed reference: $(v_naive[j]) vs $(v_reference[j])"
    T, N, Ee, Aa = v_reference[j], v_naive[j], v_expr[j], v_act[j]
    @assert T < Ee < N "expression-only must be strictly bracketed: $T < $Ee < $N"
    @assert T < Aa < N "activity-only must be strictly bracketed: $T < $Aa < $N"
    @assert !isapprox(Ee, Aa; atol=0.2) "the two partial cases must be visibly distinct: $Ee vs $Aa"
    println("test_dual_feedback (Task 3): naive=$N expr=$Ee act=$Aa both=$(v_both[j]) reference=$T")
end

# ---- Task 5: dual-control range includes the reference across shapes ------- #
let reference = feedback_reference()
    controls = control_factors(reference)
    T = reference_metabolic_fluxes(reference)[findfirst(==("r3"), METAB_REACTIONS)]
    sw = theta_shape_sweep(reference, controls)
    Tlo, Thi = extrema(sw.Tboth)
    @assert Tlo < T < Thi "dual-control sweep must include reference: [$Tlo, $Thi] vs $T"
    @assert minimum(sw.Tact) > T "activity-only must stay above reference: $(minimum(sw.Tact)) vs $T"
    @assert sw.Texpr > T "expression-only must stay above reference: $(sw.Texpr) vs $T"
    T_chosen = controls.Vmax0 * controls.e_e0 * controls.θ
    @assert Tlo <= T_chosen <= Thi "chosen throughput must lie in sweep range: $T_chosen not in [$Tlo, $Thi]"
    println("test_dual_feedback (Task 5): T_both in [$(round(Tlo,digits=3)), $(round(Thi,digits=3))] includes reference $(round(T,digits=3)); single-control min act=$(round(minimum(sw.Tact),digits=3)) expr=$(round(sw.Texpr,digits=3))")
end

# ---- Task 4: feedback_fba input validation ---------------------------------- #
let controls_bad = (Vmax0 = -1.0, θ = 0.5, e_e0 = 0.5)
    threw = false
    try
        feedback_fba(controls_bad; expression=true, activity=true)
    catch e
        threw = e isa ArgumentError
    end
    @assert threw "feedback_fba must reject a negative Vmax0"
    println("test_dual_feedback (Task 4): feedback_fba validation OK")
end
