# Parametric-bootstrap UQ for the urea-cycle FBA example.
# Config A: kcat, e0, dG, and saturation f_j (from Park et al. data) are sampled;
#           f_2 = 1 (argininosuccinate not measured). See design spec 2026-07-13.
include(joinpath(@__DIR__, "..", "Include.jl"))
include(joinpath(@__DIR__, "urea_cycle.jl"))

const N        = 10_000
const SEED     = 20260713
const SIGMA_LN = 0.69     # lognormal spread for kcat, e0, conc, Km
const DG_SIGMA = 2.0      # kJ/mol, normal spread for dG

# Nominal saturation inputs (conc_M, Km_M) per reaction, reduced from
# code/data/park_saturation.csv by the aggregation rule (least-saturated
# substrate; prefer Homo sapiens row, else geometric mean across organisms).
# v2 has no forward-substrate data -> f_2 = 1 in Config A.
const SAT_NOMINAL = Dict(
    1 => (4.673e-3, 3.923e-4),  # v1 ATP (Mus musculus)
    3 => (2.555e-4, 1.546e-3),  # v3 arginine (Homo sapiens)
    4 => (2.129e-4, 1.166e-3),  # v4 ornithine (geometric mean, yeast + E. coli)
    5 => (2.555e-4, 3.497e-6),  # v5 arginine (Mus musculus)
)

lognrand(rng, median, sln) = median * exp(sln * randn(rng))

"Draw a saturation vector f (length 5). If impute_f2, also sample f_2."
function sample_f(rng; impute_f2::Bool=false)
    f = ones(5)
    for (j, (c0, k0)) in SAT_NOMINAL
        c = lognrand(rng, c0, SIGMA_LN)
        k = lognrand(rng, k0, SIGMA_LN)
        f[j] = c / (k + c)
    end
    if impute_f2
        asa = lognrand(rng, 5e-6, 1.0)   # [argininosuccinate], M
        km2 = lognrand(rng, 5e-6, 0.7)   # ASL Km, M
        f[2] = asa / (km2 + asa)
    end
    return f
end

"""
Draw one parameter set and return the optimal flux vector (or nothing).
If `sample_saturation` is false, the saturation vector is held at `ones(5)`
(no substrate/Km draws at all) so kcat, e0, and dG alone set the capacity;
this isolates the capacity-only reference ensemble from the saturation-gateway
effect tested by Config A/B.
"""
function sample_flux(rng; impute_f2::Bool=false, sample_saturation::Bool=true)
    kcat = KCAT0 .* exp.(SIGMA_LN .* randn(rng, 5))
    e0   = lognrand(rng, E0, SIGMA_LN)          # one shared draw per sample
    dG   = DG0 .+ DG_SIGMA .* randn(rng, 5)
    f    = sample_saturation ? sample_f(rng; impute_f2 = impute_f2) : ones(5)
    solve_flux(urea_cycle_model(; kcat = kcat, e0 = e0, dG = dG, f = f))
end

"Run N draws; return an (kept x nreactions) flux matrix."
function run_ensemble(seed; impute_f2::Bool=false, sample_saturation::Bool=true)
    rng = MersenneTwister(seed)
    m0  = urea_cycle_model()
    rows = Vector{Vector{Float64}}()
    for _ in 1:N
        v = sample_flux(rng; impute_f2 = impute_f2, sample_saturation = sample_saturation)
        v === nothing && continue
        push!(rows, v)
    end
    return reduce(vcat, (r' for r in rows)), m0
end

# ---- Config A ----
F, m0 = run_ensemble(SEED; impute_f2 = false)
kept = size(F, 1)
@assert kept > 0.99 * N "too many infeasible draws: kept $kept of $N"
vnom = solve_flux(m0)

stats = DataFrame(reaction=String[], nominal=Float64[], mean=Float64[],
                  std=Float64[], cv=Float64[], q025=Float64[], q50=Float64[], q975=Float64[])
for (j, r) in enumerate(m0.reactions)
    col = F[:, j]
    mu, sd = mean(col), std(col)
    push!(stats, (r, vnom[j], mu, sd, sd / (abs(mu) + eps()),
                  quantile(col, 0.025), quantile(col, 0.5), quantile(col, 0.975)))
end
CSV.write(datapath("urea_fba_uq.csv"), stats)

# Figure: nominal flux bars. Because the cycle is linear, every nonzero flux
# equals the single throughput up to sign, so all bootstrap bands are identical;
# drawing the 2.5-97.5 percentile whisker on b4 (urea export, the objective)
# alone shows that one uncertainty without cloning it across 13 correlated fluxes.
let fig = Figure()
    ax = Axis(fig[1,1], xticks=(1:nrow(stats), stats.reaction), ylabel="flux",
              xticklabelrotation=pi/4)
    barplot!(ax, 1:nrow(stats), stats.nominal)
    b4row = findfirst(==("b4"), stats.reaction)
    rangebars!(ax, [b4row], [stats.q025[b4row]], [stats.q975[b4row]];
               color=:black, whiskerwidth=8)
    save(figpath("urea_fba.pdf"), fig)
end

b4 = findfirst(==("b4"), m0.reactions)
ua = F[:, b4]  # urea export (secretion-positive: b4 > 0 is export)
println("configA_kept=", kept)
println("configA_urea_export mean=", mean(ua), " sd=", std(ua),
        " median=", quantile(ua, 0.5), " ci=[", quantile(ua, 0.025), ",", quantile(ua, 0.975), "]")
println("configA_v5 mean=", mean(F[:, findfirst(==("v5"), m0.reactions)]))
@assert isapprox(quantile(ua, 0.5), 0.0328; rtol=0.15) "Config A urea median should stay near 0.0328"

# ---- Capacity-only reference ensemble (no saturation sampling at all) ----
# Same kcat/e0/dG sampling as Config A, but f held at ones(5) throughout, so
# this isolates the "min of sampled capacities" effect from the saturation-
# gateway effect that Config A/B probe. Not plotted; reported in prose only.
Fcap, _ = run_ensemble(SEED; sample_saturation=false)
cap_kept = size(Fcap, 1)
uacap = Fcap[:, b4]  # urea export (secretion-positive), capacity-only
println("configCap_kept=", cap_kept)
println("configCap_urea_export mean=", mean(uacap), " sd=", std(uacap),
        " median=", quantile(uacap, 0.5), " ci=[", quantile(uacap, 0.025), ",", quantile(uacap, 0.975), "]")

# ---- Config B: impute and sample f_2 (argininosuccinate unmeasured) ----
# Re-run the ensemble capturing f_2 and urea export per draw.
let rng = MersenneTwister(SEED + 1)
    b4idx = findfirst(==("b4"), m0.reactions)
    f2s = Float64[]; ub_export = Float64[]
    for _ in 1:N
        # replicate sample_flux but keep f_2 for reporting
        kcat = KCAT0 .* exp.(SIGMA_LN .* randn(rng, 5))
        e0   = lognrand(rng, E0, SIGMA_LN)
        dG   = DG0 .+ DG_SIGMA .* randn(rng, 5)
        f    = sample_f(rng; impute_f2 = true)
        v    = solve_flux(urea_cycle_model(; kcat=kcat, e0=e0, dG=dG, f=f))
        v === nothing && continue
        push!(f2s, f[2]); push!(ub_export, v[b4idx])
    end

    # Sensitivity CSV: Config A (f2=1) and Config B (sampled f2) urea export
    sens = DataFrame(config=String[], f2=Float64[], urea_export=Float64[])
    for x in ua;        push!(sens, ("A", 1.0, x)); end
    for (g, x) in zip(f2s, ub_export); push!(sens, ("B", g, x)); end
    CSV.write(datapath("urea_uq_sensitivity.csv"), sens)

    println("configB_f2 median=", quantile(f2s, 0.5))
    println("configB_urea_export median=", quantile(ub_export, 0.5),
            " ci=[", quantile(ub_export, 0.025), ",", quantile(ub_export, 0.975), "]")
    @assert quantile(ub_export, 0.5) < quantile(ua, 0.5) "Config B should lower the urea median"

    # Figure: overlaid urea-export densities, Config A vs Config B (log-x for skew)
    fig = Figure()
    ax = Axis(fig[1,1],
              xlabel = rich("urea export (mmol gDW", superscript("-1"), " h", superscript("-1"), ")"),
              ylabel="density", xscale=log10)
    # Explicit positive `boundary` keeps each KDE support off the automatic
    # padding that can dip below zero (fatal once the axis is log10-scaled).
    density!(ax, ua;        boundary=(0.5 * minimum(ua), 2 * maximum(ua)),
             color=(:steelblue, 0.4), label="A: f_2 = 1 (no data)")
    density!(ax, ub_export; boundary=(0.5 * minimum(ub_export), 2 * maximum(ub_export)),
             color=(:tomato, 0.4),    label="B: f_2 sampled")
    axislegend(ax; position=:rt)
    save(figpath("urea_saturation.pdf"), fig)
end
