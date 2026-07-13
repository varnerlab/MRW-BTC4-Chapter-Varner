include(joinpath(@__DIR__, "..", "Include.jl"))
include(joinpath(@__DIR__, "urea_cycle.jl"))

# Nominal bounds must match the pre-refactor hardcoded values.
expected_lb = vcat([-0.1, -0.0328, 0.0, 0.0, 0.0], fill(-1000.0, 14))
expected_ub = vcat([ 0.1,  0.0328, 1.9, 4.1, 0.1], fill( 1000.0, 14))

m = urea_cycle_model()
# isapprox, not ==: Vmax = kcat*e0 (e.g. 3.28*0.01) differs from the literal 0.0328
# only by floating-point round-off (~1e-18).
@assert all(isapprox.(m.lb, expected_lb; atol=1e-12)) "lb mismatch: $(m.lb)"
@assert all(isapprox.(m.ub, expected_ub; atol=1e-12)) "ub mismatch: $(m.ub)"

v = solve_flux(m)
@assert v !== nothing "nominal model failed to solve"
b4 = findfirst(==("b4"), m.reactions)
@assert isapprox(v[b4], -0.0328; atol=1e-6) "urea export $(v[b4]) != -0.0328"
@assert maximum(abs.(m.S * v)) < 1e-6 "Sv=0 violated"

# f defaults to ones -> passing the nominal Park f leaves v2 the bottleneck (fluxes unchanged)
f_park = ones(5); f_park[1] = 0.923; f_park[3] = 0.142; f_park[4] = 0.154; f_park[5] = 0.986
vf = solve_flux(urea_cycle_model(; f = f_park))
@assert isapprox(vf[b4], -0.0328; atol=1e-6) "Park-f urea export $(vf[b4]) != -0.0328 (v2 should still bind)"

println("test_urea_cycle: all checks passed")
