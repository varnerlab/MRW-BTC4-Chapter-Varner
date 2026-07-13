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

"Draw one parameter set and return the optimal flux vector (or nothing)."
function sample_flux(rng; impute_f2::Bool=false)
    kcat = KCAT0 .* exp.(SIGMA_LN .* randn(rng, 5))
    e0   = lognrand(rng, E0, SIGMA_LN)          # one shared draw per sample
    dG   = DG0 .+ DG_SIGMA .* randn(rng, 5)
    f    = sample_f(rng; impute_f2 = impute_f2)
    solve_flux(urea_cycle_model(; kcat = kcat, e0 = e0, dG = dG, f = f))
end

"Run N draws; return an (kept x nreactions) flux matrix."
function run_ensemble(seed; impute_f2::Bool=false)
    rng = MersenneTwister(seed)
    m0  = urea_cycle_model()
    n   = length(m0.reactions)
    rows = Vector{Vector{Float64}}()
    for _ in 1:N
        v = sample_flux(rng; impute_f2 = impute_f2)
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

# Figure: nominal flux bars with 2.5-97.5 percentile range bars
let fig = Figure()
    ax = Axis(fig[1,1], xticks=(1:nrow(stats), stats.reaction), ylabel="flux",
              xticklabelrotation=pi/4)
    barplot!(ax, 1:nrow(stats), stats.nominal)
    rangebars!(ax, 1:nrow(stats), stats.q025, stats.q975; color=:black, whiskerwidth=8)
    save(figpath("urea_fba.pdf"), fig)
end

b4 = findfirst(==("b4"), m0.reactions)
ua = -F[:, b4]  # urea export magnitude
println("configA_kept=", kept)
println("configA_urea_export mean=", mean(ua), " sd=", std(ua),
        " median=", quantile(ua, 0.5), " ci=[", quantile(ua, 0.025), ",", quantile(ua, 0.975), "]")
println("configA_v5 mean=", mean(F[:, findfirst(==("v5"), m0.reactions)]))
@assert isapprox(quantile(ua, 0.5), 0.0328; rtol=0.15) "Config A urea median should stay near 0.0328"
