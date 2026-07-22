include(joinpath(@__DIR__, "..", "Include.jl"))
include(joinpath(@__DIR__, "urea_cycle.jl"))

# Nominal bounds must match kcat (s^-1) * e0 * 3600 (s -> h conversion).
expected_lb = vcat([-360.0, -118.08,    0.0,     0.0,   0.0], fill(-1000.0, 14))
expected_ub = vcat([ 360.0,  118.08, 6840.0, 14760.0, 360.0], fill( 1000.0, 14))

m = urea_cycle_model()
# isapprox, not ==: Vmax = kcat*e0*3600 (e.g. 3.28*0.01*3600) differs from the literal 118.08
# only by floating-point round-off.
@assert all(isapprox.(m.lb, expected_lb; atol=1e-9)) "lb mismatch: $(m.lb)"
@assert all(isapprox.(m.ub, expected_ub; atol=1e-9)) "ub mismatch: $(m.ub)"

nominal_result = solve_flux_with_status(m)
@assert nominal_result.status == "OPTIMAL" "unexpected solver status: $(nominal_result.status)"
v = nominal_result.flux
@assert v !== nothing "nominal model failed to solve"
b4 = findfirst(==("b4"), m.reactions)
@assert isapprox(v[b4], 118.08; atol=1e-6) "urea export $(v[b4]) != 118.08 (secretion-positive)"
@assert maximum(abs.(m.S * v)) < 1e-6 "Sv=0 violated"

# f defaults to ones -> passing the nominal Park f leaves v2 the bottleneck (fluxes unchanged)
f_park = ones(5); f_park[1] = 0.923; f_park[3] = 0.142; f_park[4] = 0.154; f_park[5] = 0.986
vf = solve_flux(urea_cycle_model(; f = f_park))
@assert isapprox(vf[b4], 118.08; atol=1e-6) "Park-f urea export $(vf[b4]) != 118.08 (v2 should still bind)"

# Dimensional-conversion check: Vmax must equal kcat*e0*3600 exactly, not kcat*e0.
kcat_test = [7.0, 7.0, 7.0, 7.0, 7.0]
e0_test   = 0.02
m2 = urea_cycle_model(; kcat = kcat_test, e0 = e0_test)
@assert isapprox(m2.ub[1], 7.0 * 0.02 * 3600.0; atol=1e-9) "Vmax not converted s^-1 -> h^-1: got $(m2.ub[1])"
@assert !isapprox(m2.ub[1], 7.0 * 0.02; atol=1e-9) "Vmax still in per-second units"

# ---- Task 2: input validation and explicit solve-failure handling ---------- #
using Test

@test_throws ArgumentError urea_cycle_model(; kcat = [1.0, 2.0])          # wrong length
@test_throws ArgumentError urea_cycle_model(; dG = fill(NaN, 5))          # non-finite
@test_throws ArgumentError urea_cycle_model(; f = [-0.1, 1.0, 1.0, 1.0, 1.0])  # negative saturation
@test_throws ArgumentError urea_cycle_model(; e0 = -0.01)                 # lb > ub (negative e0)

infeasible = urea_cycle_model(; kcat = zeros(5))  # zero capacity everywhere -> b4 forced to 0, still feasible;
infeasible = (; infeasible..., lb = infeasible.lb .+ 1.0, ub = infeasible.ub)  # lb > ub -> infeasible
failed_result = solve_flux_with_status(infeasible)
@assert failed_result.flux === nothing
@assert failed_result.status != "OPTIMAL"
@test_throws ErrorException solve_fba(infeasible)

println("test_urea_cycle (Task 2): validation and solve-failure checks passed")

# ---- Task 3: FVA confirms the nominal optimum is unique ------------------- #
m0  = urea_cycle_model()
rng = fva(m0)
@assert nrow(rng) == length(m0.reactions)
@assert maximum(abs.(rng.vmax .- rng.vmin)) < 1e-4 "nominal optimum is not unique: $(rng)"

println("test_urea_cycle (Task 3): FVA uniqueness check passed")

println("test_urea_cycle: all checks passed")
