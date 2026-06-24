include(joinpath(@__DIR__, "..", "Include.jl"))
include(joinpath(@__DIR__, "urea_cycle.jl"))

m = urea_cycle_model()

function solve_fba(m)
    model = Model(HiGHS.Optimizer); set_silent(model)
    n = length(m.reactions)
    @variable(model, m.lb[i] <= v[i=1:n] <= m.ub[i])
    @constraint(model, m.S * v .== 0)
    @objective(model, Max, sum(m.c[i]*v[i] for i in 1:n))
    optimize!(model)
    DataFrame(reaction=m.reactions, flux=value.(v))
end

res = solve_fba(m)
CSV.write(datapath("urea_fba_solution.csv"), res)

let fig = Figure()
    ax = Axis(fig[1,1], xticks=(1:nrow(res), res.reaction), ylabel="flux",
              xticklabelrotation=pi/4)
    barplot!(ax, 1:nrow(res), res.flux)
    save(figpath("urea_fba.pdf"), fig)
end

println("objective_flux=", res.flux[argmax(abs.(m.c))])

# Steady-state mass balance: Sv = 0 must hold at the optimal solution
@assert maximum(abs.(m.S * res.flux)) < 1e-6 "Mass balance Sv=0 failed at optimal solution"
