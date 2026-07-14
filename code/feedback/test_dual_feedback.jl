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
