# code/geneexpression/run_cybernetic.jl
include(joinpath(@__DIR__, "..", "Include.jl"))
include(joinpath(@__DIR__, "cybernetic.jl"))
df = simulate_diauxie(tspan=(0.0, 20.0))
CSV.write(datapath("cybernetic_diauxie.csv"), df)
let fig = Figure()
    ax = Axis(fig[1,1], xlabel="time (h)", ylabel="conc.")
    lines!(ax, df.t, df.biomass; label="biomass")
    lines!(ax, df.t, df.S1;      label="substrate 1")
    lines!(ax, df.t, df.S2;      label="substrate 2")
    axislegend(ax)
    save(figpath("cybernetic_diauxie.pdf"), fig)
end
println("S1_end=", df.S1[end], " S2_end=", df.S2[end])
