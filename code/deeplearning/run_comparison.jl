# run_comparison.jl — S4 vs LSTM head-to-head comparison figure.
# Produces: code/figs/s4_vs_lstm.pdf
#
# Figure layout:
#   Panel 1 (top):    Per-state RMSE normalized by each state's data range
#                     (relative RMSE %), displayed as grouped bars.
#                     Normalization makes all 7 states legible on a single
#                     y-axis despite RMSE spanning ~5 orders of magnitude
#                     (mAb ~35–162 vs X ~0.03).
#   Panel 2 (bottom): Truth vs LSTM vs S4 forecast overlay for mAb (the
#                     headline product state) over the first test window,
#                     in physical units (mg/L).
#
# Reproducibility: uses LOOKBACK=24, HORIZON=12, test_frac=0.2, seed=42 —
# identical to run_lstm.jl and run_s4.jl so the comparison is fair.
#
# Usage: cd code && julia --project=. deeplearning/run_comparison.jl

include(joinpath(@__DIR__, "..", "Include.jl"))

include(joinpath(@__DIR__, "data_prep.jl"))
include(joinpath(@__DIR__, "lstm.jl"))
include(joinpath(@__DIR__, "s4.jl"))

# ── Reproducibility ─────────────────────────────────────────────────────────
Random.seed!(42)

# ── Hyperparameters (must match run_lstm.jl / run_s4.jl exactly) ────────────
const LOOKBACK  = 24
const HORIZON   = 12
const HIDDEN    = 32      # LSTM hidden units
const H_ORDER   = 16      # S4 per-channel LegS basis order
const DT        = 1.0     # S4 bilinear discretization timestep
const EPOCHS    = 300
const LR_LSTM   = 5e-3
const LR_S4     = 1e-2
const BATCHSIZE = 32

const state_cols = ["V", "X", "Glc", "Gln", "Lac", "Amm", "mAb"]
const N_states   = length(state_cols)

# ── Load data ────────────────────────────────────────────────────────────────
df = CSV.read(joinpath(_PATH_TO_DATA, "cho_trajectories.csv"), DataFrame)

# ── Build windows and scaler ─────────────────────────────────────────────────
# Identical call signature to run_lstm.jl and run_s4.jl
Xtrain, Ytrain, Xtest, Ytest, scaler = make_windows(df;
    lookback  = LOOKBACK,
    horizon   = HORIZON,
    test_frac = 0.2,
)

println("Data shapes:")
println("  Xtrain: ", size(Xtrain), "  (states × lookback × train_windows)")
println("  Ytrain: ", size(Ytrain), "  (states × horizon  × train_windows)")
println("  Xtest:  ", size(Xtest),  "  (states × lookback × test_windows)")
println("  Ytest:  ", size(Ytest),  "  (states × horizon  × test_windows)")

# ── Train LSTM ───────────────────────────────────────────────────────────────
Random.seed!(42)
model_lstm = train_lstm(Xtrain, Ytrain;
    epochs    = EPOCHS,
    hidden    = HIDDEN,
    lr        = LR_LSTM,
    batchsize = BATCHSIZE,
    seed      = 42,
)

# ── Train S4 ─────────────────────────────────────────────────────────────────
Random.seed!(42)
model_s4 = train_s4(Xtrain, Ytrain;
    epochs    = EPOCHS,
    h         = H_ORDER,
    Δt        = DT,
    lr        = LR_S4,
    batchsize = BATCHSIZE,
    seed      = 42,
)

# ── Forecast (normalized) ────────────────────────────────────────────────────
Ŷ_lstm_norm = forecast(model_lstm, Xtest)   # N × H × Btest
Ŷ_s4_norm   = forecast(model_s4,   Xtest)   # N × H × Btest

# ── Inverse-transform to physical units ──────────────────────────────────────
function denorm_all(Ŷ_norm, scaler)
    Ŷ_phys = similar(Ŷ_norm)
    for b in axes(Ŷ_norm, 3)
        Ŷ_phys[:, :, b] = Float32.(scaler.denormalize(Float64.(Ŷ_norm[:, :, b])))
    end
    return Ŷ_phys
end

Ŷ_lstm_phys = denorm_all(Ŷ_lstm_norm, scaler)
Ŷ_s4_phys   = denorm_all(Ŷ_s4_norm,   scaler)

Y_phys = similar(Ytest)
for b in axes(Ytest, 3)
    Y_phys[:, :, b] = Float32.(scaler.denormalize(Float64.(Ytest[:, :, b])))
end

# ── Per-state RMSE in physical units ─────────────────────────────────────────
rmse_lstm = Float64[]
rmse_s4   = Float64[]
for s in 1:N_states
    err_lstm = Y_phys[s, :, :] .- Ŷ_lstm_phys[s, :, :]
    err_s4   = Y_phys[s, :, :] .- Ŷ_s4_phys[s, :, :]
    push!(rmse_lstm, sqrt(mean(err_lstm .^ 2)))
    push!(rmse_s4,   sqrt(mean(err_s4   .^ 2)))
end

println("\nPer-state RMSE comparison (physical units):")
println(lpad("State", 8), "  ", lpad("LSTM", 12), "  ", lpad("S4", 12), "  ",
        lpad("S4/LSTM", 10))
for (name, rl, rs) in zip(state_cols, rmse_lstm, rmse_s4)
    println(lpad(name, 8), "  ", lpad(round(rl; digits=4), 12), "  ",
            lpad(round(rs; digits=4), 12), "  ", lpad(round(rs/rl; digits=4), 10))
end

# ── Normalization for bar chart: relative RMSE % ─────────────────────────────
# Divide each state's RMSE by that state's data range (max - min across the
# full trajectory) and multiply by 100 to get relative RMSE %.
# This makes all 7 states legible on a single linear y-axis despite spanning
# ~5 orders of magnitude in absolute RMSE.
raw_phys = Matrix(df[:, [Symbol(s) for s in state_cols]])'  # N × T

data_range = Float64[]
for s in 1:N_states
    rng_s = maximum(raw_phys[s, :]) - minimum(raw_phys[s, :])
    push!(data_range, rng_s > 0 ? rng_s : 1.0)
end

rel_rmse_lstm = (rmse_lstm ./ data_range) .* 100.0   # %
rel_rmse_s4   = (rmse_s4   ./ data_range) .* 100.0   # %

println("\nRelative RMSE % (RMSE / data-range × 100):")
println(lpad("State", 8), "  ", lpad("LSTM%", 10), "  ", lpad("S4%", 10))
for (name, rl, rs) in zip(state_cols, rel_rmse_lstm, rel_rmse_s4)
    println(lpad(name, 8), "  ", lpad(round(rl; digits=3), 10), "  ",
            lpad(round(rs; digits=3), 10))
end

# ── mAb overlay time indices ──────────────────────────────────────────────────
t_all           = df.t
n_train_windows = size(Xtrain, 3)
b               = 1   # first test window
t0_idx          = n_train_windows + b
t_input_idx     = t0_idx:(t0_idx + LOOKBACK - 1)
t_target_idx    = (t0_idx + LOOKBACK):(t0_idx + LOOKBACK + HORIZON - 1)
t_input         = t_all[t_input_idx]
t_target        = t_all[t_target_idx]

mAb_idx = findfirst(==("mAb"), state_cols)
println("\nOverlay drawn for state: ", state_cols[mAb_idx], " (index ", mAb_idx, ")")
println("  first test window, train_windows=", n_train_windows,
        ", t_target = ", t_all[first(t_target_idx)], "–", t_all[last(t_target_idx)], " h")

# ── Build figure ─────────────────────────────────────────────────────────────
fig = Figure(size = (900, 700))

# ── Panel 1 (top): Grouped bars of relative RMSE % ──────────────────────────
ax1 = Axis(fig[1, 1];
    title      = "Per-State Relative RMSE  (RMSE / data-range × 100%)",
    xlabel     = "State",
    ylabel     = "Relative RMSE (%)",
    xticks     = (1:N_states, state_cols),
    titlesize  = 13,
    xlabelsize = 12,
    ylabelsize = 12,
)

bar_width = 0.35
offsets   = [-bar_width / 2, bar_width / 2]

# LSTM bars (left, steel blue)
barplot!(ax1, (1:N_states) .+ offsets[1], rel_rmse_lstm;
    width = bar_width,
    color = :steelblue,
    label = "LSTM",
)
# S4 bars (right, orange)
barplot!(ax1, (1:N_states) .+ offsets[2], rel_rmse_s4;
    width = bar_width,
    color = :darkorange,
    label = "S4",
)

axislegend(ax1; position = :rt, labelsize = 11)

# ── Panel 2 (bottom): mAb forecast overlay ───────────────────────────────────
ax2 = Axis(fig[2, 1];
    title      = "mAb Forecast Overlay — first test window (physical units)",
    xlabel     = "Time (h)",
    ylabel     = "mAb (mg/L)",
    titlesize  = 13,
    xlabelsize = 12,
    ylabelsize = 12,
)

# Full trajectory (grey background)
lines!(ax2, t_all, raw_phys[mAb_idx, :];
    color = (:grey, 0.35), linewidth = 1, label = "full trajectory")

# Input window (blue)
lines!(ax2, t_input, raw_phys[mAb_idx, t_input_idx];
    color = :steelblue, linewidth = 2, linestyle = :solid, label = "input window")

# Ground truth forecast (black dashed)
lines!(ax2, t_target, Y_phys[mAb_idx, :, b];
    color = :black, linewidth = 2, linestyle = :dash, label = "truth")

# LSTM forecast (steel blue, dotted)
lines!(ax2, t_target, Ŷ_lstm_phys[mAb_idx, :, b];
    color = :crimson, linewidth = 2, linestyle = :solid, label = "LSTM")

# S4 forecast (orange)
lines!(ax2, t_target, Ŷ_s4_phys[mAb_idx, :, b];
    color = :darkorange, linewidth = 2, linestyle = :solid, label = "S4")

axislegend(ax2; position = :lt, labelsize = 11)

# ── Figure title ─────────────────────────────────────────────────────────────
Label(fig[0, 1], "S4 vs LSTM — CHO Fed-Batch Forecaster Comparison";
    fontsize = 14, font = :bold)

# ── Save ─────────────────────────────────────────────────────────────────────
save(joinpath(_PATH_TO_FIGS, "s4_vs_lstm.pdf"), fig)
println("\nSaved: figs/s4_vs_lstm.pdf")
println("Overlay state: mAb — confirmed")

println("\nCMP_OK")
