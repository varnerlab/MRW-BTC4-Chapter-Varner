# run_lstm.jl — Train LSTM forecaster on CHO bioreactor data, save figure + metrics.
# Usage: cd code && julia --project=. deeplearning/run_lstm.jl

include(joinpath(@__DIR__, "..", "Include.jl"))

include(joinpath(@__DIR__, "data_prep.jl"))
include(joinpath(@__DIR__, "lstm.jl"))

# ── Reproducibility ─────────────────────────────────────────────────────────
Random.seed!(42)

# ── Hyperparameters ──────────────────────────────────────────────────────────
const LOOKBACK  = 24    # history window (timesteps)
const HORIZON   = 12    # forecast horizon (timesteps)
const HIDDEN    = 32    # LSTM hidden units
const EPOCHS    = 300
const LR        = 5e-3
const BATCHSIZE = 32

# ── Load data ────────────────────────────────────────────────────────────────
df = CSV.read(joinpath(_PATH_TO_DATA, "cho_trajectories.csv"), DataFrame)

# ── Build windows and scaler ─────────────────────────────────────────────────
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

# ── Train ────────────────────────────────────────────────────────────────────
model = train_lstm(Xtrain, Ytrain;
    epochs    = EPOCHS,
    hidden    = HIDDEN,
    lr        = LR,
    batchsize = BATCHSIZE,
    seed      = 42,
)

# ── Forecast ─────────────────────────────────────────────────────────────────
Ŷ_norm = forecast(model, Xtest)    # N × H × Btest  (normalized)

state_cols = ["V", "X", "Glc", "Gln", "Lac", "Amm", "mAb"]
N_states   = length(state_cols)

# Inverse-transform to physical units
# Ytest and Ŷ_norm are (N × H × B); denormalize expects (N × *)
Ŷ_phys = similar(Ŷ_norm)
Y_phys = similar(Ytest)
for b in axes(Ŷ_norm, 3)
    Ŷ_phys[:, :, b] = Float32.(scaler.denormalize(Float64.(Ŷ_norm[:, :, b])))
    Y_phys[:, :, b]  = Float32.(scaler.denormalize(Float64.(Ytest[:, :, b])))
end

# ── Per-state RMSE (physical units, averaged over all test windows & horizon steps) ──
rmse_vals = Float64[]
for s in 1:N_states
    err = Y_phys[s, :, :] .- Ŷ_phys[s, :, :]
    push!(rmse_vals, sqrt(mean(err .^ 2)))
end

println("\nPer-state RMSE (physical units):")
for (name, r) in zip(state_cols, rmse_vals)
    println("  ", lpad(name, 6), " : ", round(r; digits=4))
end

# ── Save metrics CSV ─────────────────────────────────────────────────────────
metrics_df = DataFrame(state = state_cols, rmse = rmse_vals)
CSV.write(joinpath(_PATH_TO_DATA, "lstm_metrics.csv"), metrics_df)
println("\nSaved: data/lstm_metrics.csv")

# ── Figure: truth vs forecast overlay ────────────────────────────────────────
# Use the first test window for illustration
fig = Figure(size = (1200, 900))

# Time axis for the forecast: last LOOKBACK steps + HORIZON steps
t_all = df.t
n_total = nrow(df)
n_test  = size(Xtest, 3)
n_train_windows = size(Xtrain, 3)

# For first test window: find the corresponding raw time indices
b = 1   # first test window
t0_idx = n_train_windows + b   # 1-indexed start of this window in time series
t_input_idx  = t0_idx:(t0_idx + LOOKBACK - 1)
t_target_idx = (t0_idx + LOOKBACK):(t0_idx + LOOKBACK + HORIZON - 1)

t_input  = t_all[t_input_idx]
t_target = t_all[t_target_idx]

# Full trajectory in physical units for reference
raw_phys = Matrix(df[:, [Symbol(s) for s in state_cols]])'  # N × T

for (s, name) in enumerate(state_cols)
    row = (s - 1) ÷ 4 + 1
    col = mod(s - 1, 4) + 1
    ax = Axis(fig[row, col];
        title  = name,
        xlabel = "Time (h)",
        ylabel = name,
    )

    # Full trajectory (grey background)
    lines!(ax, t_all, raw_phys[s, :]; color = (:grey, 0.4), linewidth = 1, label = "trajectory")

    # Input window (blue)
    lines!(ax, t_input, raw_phys[s, t_input_idx]; color = :steelblue, linewidth = 2, label = "input")

    # Ground truth forecast window (black)
    lines!(ax, t_target, Y_phys[s, :, b]; color = :black, linewidth = 2, linestyle = :dash, label = "truth")

    # LSTM forecast (red)
    lines!(ax, t_target, Ŷ_phys[s, :, b]; color = :crimson, linewidth = 2, label = "LSTM")

    axislegend(ax; position = :rt, labelsize = 9)
end

# Empty 8th panel (7 states, 4×2 grid has 8 slots)
ax_empty = Axis(fig[2, 4])
hidedecorations!(ax_empty)
hidespines!(ax_empty)

Label(fig[0, :], "LSTM Forecast vs Truth — CHO Fed-Batch (first test window)";
    fontsize = 14, font = :bold)

save(joinpath(_PATH_TO_FIGS, "lstm_cho.pdf"), fig)
println("Saved: figs/lstm_cho.pdf")

println("\nLSTM_OK")
