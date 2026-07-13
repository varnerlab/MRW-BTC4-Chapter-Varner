include(joinpath(@__DIR__, "..", "Include.jl"))
include(joinpath(@__DIR__, "urea_cycle.jl"))

m = urea_cycle_model()

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
