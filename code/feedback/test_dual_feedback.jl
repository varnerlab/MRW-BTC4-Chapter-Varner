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
