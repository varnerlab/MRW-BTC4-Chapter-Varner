# run_cho.jl
# Simulate the CHO fed-batch ODE, write synthetic dataset + kinetics figure.
# Usage (from code/):  julia --project=. kinetics/run_cho.jl

include(joinpath(@__DIR__, "..", "Include.jl"))
include(joinpath(@__DIR__, "cho_model.jl"))

df = simulate_cho(tspan=(0.0, 240.0), saveat=1.0)   # hours

CSV.write(datapath("cho_trajectories.csv"), df)

let fig = Figure()
    ax = Axis(fig[1,1], xlabel="time (h)", ylabel="conc.")
    lines!(ax, df.t, df.X;   label="biomass")
    lines!(ax, df.t, df.mAb; label="mAb")
    axislegend(ax)
    save(figpath("cho_kinetics.pdf"), fig)
end

println("rows=", nrow(df), " final_mAb=", df.mAb[end])
