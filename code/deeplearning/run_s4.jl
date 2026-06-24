# run_s4.jl — Train S4/HiPPO-LegS forecaster on CHO bioreactor data, save figure + metrics.
# Usage: cd code && julia --project=. deeplearning/run_s4.jl
#
# §5 S4 counterpart to run_lstm.jl.  Uses the IDENTICAL windows/scaler as the
# LSTM (same LOOKBACK / HORIZON) so the two forecasters are directly comparable
# (Task 4.3 head-to-head).

include(joinpath(@__DIR__, "..", "Include.jl"))

include(joinpath(@__DIR__, "data_prep.jl"))
include(joinpath(@__DIR__, "s4.jl"))

# ── Reproducibility ─────────────────────────────────────────────────────────
Random.seed!(42)

# ── Hyperparameters ──────────────────────────────────────────────────────────
# LOOKBACK / HORIZON are BINDING: identical to run_lstm.jl for comparability.
const LOOKBACK  = 24      # history window (timesteps) — matches run_lstm.jl
const HORIZON   = 12      # forecast horizon (timesteps) — matches run_lstm.jl
const H_ORDER   = 16      # per-channel LegS basis order (hidden dim per channel)
const DT        = 1.0     # bilinear discretization timestep
const EPOCHS    = 300
const LR        = 1e-2
const BATCHSIZE = 32

# ── Load data ────────────────────────────────────────────────────────────────
df = CSV.read(joinpath(_PATH_TO_DATA, "cho_trajectories.csv"), DataFrame)

# ── Build windows and scaler (SAME pipeline + params as run_lstm.jl) ─────────
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
model = train_s4(Xtrain, Ytrain;
    epochs    = EPOCHS,
    h         = H_ORDER,
    Δt        = DT,
    lr        = LR,
    batchsize = BATCHSIZE,
    seed      = 42,
)

# ── Forecast ─────────────────────────────────────────────────────────────────
Ŷ_norm = forecast(model, Xtest)    # N × H × Btest  (normalized)

state_cols = ["V", "X", "Glc", "Gln", "Lac", "Amm", "mAb"]
N_states   = length(state_cols)

# Inverse-transform to physical units (same scaler as the LSTM).
Ŷ_phys = similar(Ŷ_norm)
Y_phys = similar(Ytest)
for b in axes(Ŷ_norm, 3)
    Ŷ_phys[:, :, b] = Float32.(scaler.denormalize(Float64.(Ŷ_norm[:, :, b])))
    Y_phys[:, :, b]  = Float32.(scaler.denormalize(Float64.(Ytest[:, :, b])))
end

# ── Per-state RMSE (physical units, over all test windows & horizon steps) ──
rmse_vals = Float64[]
for s in 1:N_states
    err = Y_phys[s, :, :] .- Ŷ_phys[s, :, :]
    push!(rmse_vals, sqrt(mean(err .^ 2)))
end

println("\nPer-state RMSE (physical units):")
for (name, r) in zip(state_cols, rmse_vals)
    println("  ", lpad(name, 6), " : ", round(r; digits=4))
end

@assert all(isfinite, rmse_vals) "Non-finite RMSE detected"

# ── Save metrics CSV ─────────────────────────────────────────────────────────
metrics_df = DataFrame(state = state_cols, rmse = rmse_vals)
CSV.write(joinpath(_PATH_TO_DATA, "s4_metrics.csv"), metrics_df)
println("\nSaved: data/s4_metrics.csv")

# ── Figure: truth vs forecast overlay (first test window) ────────────────────
fig = Figure(size = (1200, 900))

t_all = df.t
n_total = nrow(df)
n_test  = size(Xtest, 3)
n_train_windows = size(Xtrain, 3)

b = 1   # first test window
t0_idx = n_train_windows + b
t_input_idx  = t0_idx:(t0_idx + LOOKBACK - 1)
t_target_idx = (t0_idx + LOOKBACK):(t0_idx + LOOKBACK + HORIZON - 1)

t_input  = t_all[t_input_idx]
t_target = t_all[t_target_idx]

raw_phys = Matrix(df[:, [Symbol(s) for s in state_cols]])'  # N × T

for (s, name) in enumerate(state_cols)
    row = (s - 1) ÷ 4 + 1
    col = mod(s - 1, 4) + 1
    ax = Axis(fig[row, col];
        title  = name,
        xlabel = "Time (h)",
        ylabel = name,
    )

    lines!(ax, t_all, raw_phys[s, :]; color = (:grey, 0.4), linewidth = 1, label = "trajectory")
    lines!(ax, t_input, raw_phys[s, t_input_idx]; color = :steelblue, linewidth = 2, label = "input")
    lines!(ax, t_target, Y_phys[s, :, b]; color = :black, linewidth = 2, linestyle = :dash, label = "truth")
    lines!(ax, t_target, Ŷ_phys[s, :, b]; color = :darkorange, linewidth = 2, label = "S4")

    axislegend(ax; position = :rt, labelsize = 9)
end

ax_empty = Axis(fig[2, 4])
hidedecorations!(ax_empty)
hidespines!(ax_empty)

Label(fig[0, :], "S4 / HiPPO-LegS Forecast vs Truth — CHO Fed-Batch (first test window)";
    fontsize = 14, font = :bold)

save(joinpath(_PATH_TO_FIGS, "s4_cho.pdf"), fig)
println("Saved: figs/s4_cho.pdf")

println("\nS4_OK")
