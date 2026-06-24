# run_hybrid.jl — Mechanistic+ML hybrid demonstration on CHO bioreactor data.
#
# Pedagogical setup
# -----------------
# TRUE data    : code/data/cho_trajectories.csv  (correct §2 model).
# IMPERFECT mechanism : default_cho_params() with mu_max scaled ×0.75 and
#                       Y_X_glc scaled ×0.75 (≈25% under-estimation of growth
#                       rate and yield — plausible identification error).
# THREE forecasts compared on the held-out tail of the trajectory:
#   1. Pure mechanistic  — biased ODE prediction (re-initialised from
#                          the true state at the split point).
#   2. Pure data-driven  — LSTM trained directly on normalised true-data windows.
#   3. Hybrid            — LSTM trained to predict the RESIDUAL
#                          (true − mechanistic) added back to the mechanistic.
#
# Figures saved to : code/figs/hybrid_cho.pdf
# Usage            : cd code && julia --project=. deeplearning/run_hybrid.jl
#
# Expected output  : three RMSE tables + "HYB_OK"

include(joinpath(@__DIR__, "..", "Include.jl"))

include(joinpath(@__DIR__, "data_prep.jl"))
include(joinpath(@__DIR__, "lstm.jl"))
include(joinpath(@__DIR__, "..","kinetics", "cho_model.jl"))

Random.seed!(42)

# ── Hyper-parameters ─────────────────────────────────────────────────────────
const LOOKBACK  = 24
const HORIZON   = 12
const HIDDEN    = 32
const EPOCHS    = 300
const LR        = 5e-3
const BATCHSIZE = 32

const STATE_COLS = ["V", "X", "Glc", "Gln", "Lac", "Amm", "mAb"]
const N_STATES   = length(STATE_COLS)

# ── 1. Load true data ─────────────────────────────────────────────────────────
df_true = CSV.read(joinpath(_PATH_TO_DATA, "cho_trajectories.csv"), DataFrame)

T_total = nrow(df_true)
println("True trajectory length: $T_total timesteps")

# ── 2. Build imperfect mechanistic prediction ─────────────────────────────────
# Perturb mu_max and Y_X_glc by −25% to introduce systematic bias.
p_biased = default_cho_params()
# We can't modify the NamedTuple fields directly, so reconstruct with new values.
p_biased = (
    mu_max    = 0.75 * p_biased.mu_max,   # ×0.75
    K_glc     = p_biased.K_glc,
    K_gln     = p_biased.K_gln,
    K_I_lac   = p_biased.K_I_lac,
    K_I_amm   = p_biased.K_I_amm,
    k_d       = p_biased.k_d,
    KD_lac    = p_biased.KD_lac,
    KD_amm    = p_biased.KD_amm,
    alpha_P   = p_biased.alpha_P,
    beta_P    = p_biased.beta_P,
    Y_X_glc   = 0.75 * p_biased.Y_X_glc,  # ×0.75
    Y_X_gln   = p_biased.Y_X_gln,
    Y_lac_glc = p_biased.Y_lac_glc,
    Y_amm_gln = p_biased.Y_amm_gln,
    S_glc_f   = p_biased.S_glc_f,
    S_gln_f   = p_biased.S_gln_f,
    F_max     = p_biased.F_max,
    Glc_min   = p_biased.Glc_min,
    Glc_max   = p_biased.Glc_max,
    feed_on   = Ref(1.0),                  # feed on at t=0 (Glc0 < Glc_min)
)

println("\nImperfect mechanistic model perturbations:")
println("  mu_max  : default × 0.75  (growth rate under-estimated by 25%)")
println("  Y_X_glc : default × 0.75  (biomass-on-glucose yield under-estimated by 25%)")

# Simulate the imperfect model over the full time span.
t_start = df_true.t[1]
t_end   = df_true.t[end]

# Use the same initial conditions as the true simulation.
u0_true = Float64.([df_true.V[1], df_true.X[1], df_true.Glc[1],
                    df_true.Gln[1], df_true.Lac[1], df_true.Amm[1], df_true.mAb[1]])

prob_biased = ODEProblem(cho_rhs!, u0_true, (t_start, t_end), p_biased)
cb_biased   = _build_feed_callbacks(p_biased)
sol_biased  = solve(prob_biased, Tsit5();
                    callback  = cb_biased,
                    saveat    = 1.0,
                    reltol    = 1e-6,
                    abstol    = 1e-8,
                    isoutofdomain = (u, p, t) -> any(x -> x < -1e-6, u))

# Extract imperfect mechanistic trajectory on same time grid as true data.
t_mech  = sol_biased.t
arr_mec = hcat(sol_biased.u...)'   # T_mech × 7

# Align to the true-data time grid (both saved at saveat=1.0 starting from 0).
# Trim to same length in case solver stopped early.
T_mech = size(arr_mec, 1)
T_use  = min(T_total, T_mech)

arr_true_full = Float64.(Matrix(df_true[1:T_use, Symbol.(STATE_COLS)]))  # T_use × 7
arr_mec_full  = arr_mec[1:T_use, :]                                      # T_use × 7

println("Aligned trajectory length: $T_use timesteps")

# Build a DataFrame for the mechanistic trajectory (same format as df_true).
df_mec = DataFrame(
    t   = df_true.t[1:T_use],
    V   = arr_mec_full[:, 1],
    X   = arr_mec_full[:, 2],
    Glc = arr_mec_full[:, 3],
    Gln = arr_mec_full[:, 4],
    Lac = arr_mec_full[:, 5],
    Amm = arr_mec_full[:, 6],
    mAb = arr_mec_full[:, 7],
)

df_true_use = df_true[1:T_use, :]

# ── 3. Sliding windows: true data for pure-ML, residual for hybrid ────────────
# Residual series: true − mechanistic (physical units, T_use × 7).
arr_resid_full = arr_true_full .- arr_mec_full   # T_use × 7

df_resid = DataFrame(
    t   = df_true.t[1:T_use],
    V   = arr_resid_full[:, 1],
    X   = arr_resid_full[:, 2],
    Glc = arr_resid_full[:, 3],
    Gln = arr_resid_full[:, 4],
    Lac = arr_resid_full[:, 5],
    Amm = arr_resid_full[:, 6],
    mAb = arr_resid_full[:, 7],
)

# make_windows on true data (for pure-ML).
Xtrain_ml, Ytrain_ml, Xtest_ml, Ytest_ml, scaler_ml =
    make_windows(df_true_use; lookback=LOOKBACK, horizon=HORIZON, test_frac=0.2)

# make_windows on residual data (for hybrid residual model).
Xtrain_res, Ytrain_res, Xtest_res, Ytest_res, scaler_res =
    make_windows(df_resid; lookback=LOOKBACK, horizon=HORIZON, test_frac=0.2)

n_train_win = size(Xtrain_ml, 3)
n_test_win  = size(Xtest_ml,  3)
println("\nWindow split: $n_train_win training, $n_test_win test windows")

# ── 4. Train the pure-ML LSTM ─────────────────────────────────────────────────
println("\n=== Training pure-ML LSTM (on true data) ===")
model_ml = train_lstm(Xtrain_ml, Ytrain_ml;
    epochs    = EPOCHS,
    hidden    = HIDDEN,
    lr        = LR,
    batchsize = BATCHSIZE,
    seed      = 42,
)

# ── 5. Train the residual LSTM ────────────────────────────────────────────────
println("\n=== Training residual LSTM (on true − mechanistic residuals) ===")
model_res = train_lstm(Xtrain_res, Ytrain_res;
    epochs    = EPOCHS,
    hidden    = HIDDEN,
    lr        = LR,
    batchsize = BATCHSIZE,
    seed      = 42,
)

# ── 6. Generate forecasts ─────────────────────────────────────────────────────
# Pure-ML forecast (normalized → physical).
Ŷ_ml_norm = forecast(model_ml,  Xtest_ml)   # N × H × Btest (normalized)
Ŷ_ml_phys = similar(Ŷ_ml_norm)
Y_ml_phys  = similar(Ytest_ml)
for b in axes(Ŷ_ml_norm, 3)
    Ŷ_ml_phys[:, :, b] = Float32.(scaler_ml.denormalize(Float64.(Ŷ_ml_norm[:, :, b])))
    Y_ml_phys[:, :, b]  = Float32.(scaler_ml.denormalize(Float64.(Ytest_ml[:, :, b])))
end

# Residual forecast (normalized residual → physical residual).
Ŷ_res_norm = forecast(model_res, Xtest_res)  # N × H × Btest
Ŷ_res_phys = similar(Ŷ_res_norm)
for b in axes(Ŷ_res_norm, 3)
    Ŷ_res_phys[:, :, b] = Float32.(scaler_res.denormalize(Float64.(Ŷ_res_norm[:, :, b])))
end

# Mechanistic prediction over test windows (physical units).
# For each test window b, the mechanistic forecast for the horizon is taken
# directly from arr_mec_full at the corresponding time indices.
# Window b (1-indexed) in test set corresponds to start index n_train_win + b.
Ŷ_mec_phys = zeros(Float32, N_STATES, HORIZON, n_test_win)
for b in 1:n_test_win
    t0 = n_train_win + b   # 1-indexed start of this window in the time series
    # mechanistic horizon: t0+LOOKBACK .. t0+LOOKBACK+HORIZON-1
    h_start = t0 + LOOKBACK
    h_end   = t0 + LOOKBACK + HORIZON - 1
    if h_end <= T_use
        Ŷ_mec_phys[:, :, b] = Float32.(arr_mec_full[h_start:h_end, :]')
    end
end

# Hybrid forecast = mechanistic prediction + learned residual correction.
Ŷ_hyb_phys = Ŷ_mec_phys .+ Ŷ_res_phys

# Ground-truth physical values (aligned with pure-ML Ytest).
# Y_ml_phys already holds denormalized true values.
Y_true_phys = Y_ml_phys   # N × H × Btest

# ── 7. Compute per-state RMSE ─────────────────────────────────────────────────
function per_state_rmse(Ypred, Ytrue)
    # Ypred, Ytrue : N × H × B  Float32 arrays
    rmse = Float64[]
    for s in 1:size(Ytrue, 1)
        err = Float64.(Ytrue[s, :, :]) .- Float64.(Ypred[s, :, :])
        push!(rmse, sqrt(mean(err .^ 2)))
    end
    return rmse
end

rmse_mec = per_state_rmse(Ŷ_mec_phys, Y_true_phys)
rmse_ml  = per_state_rmse(Ŷ_ml_phys,  Y_true_phys)
rmse_hyb = per_state_rmse(Ŷ_hyb_phys, Y_true_phys)

println("\n" * "="^62)
println(lpad("State", 6), "  ",
        lpad("Pure-Mech", 12), "  ",
        lpad("Pure-ML",   12), "  ",
        lpad("Hybrid",    12))
println("-"^62)
units = ["L", "gDW/L", "mM", "mM", "mM", "mM", "mg/L"]
for s in 1:N_STATES
    println(lpad(STATE_COLS[s], 6), " [", rpad(units[s], 5), "]  ",
            lpad(round(rmse_mec[s]; digits=4), 12), "  ",
            lpad(round(rmse_ml[s];  digits=4), 12), "  ",
            lpad(round(rmse_hyb[s]; digits=4), 12))
end
println("="^62)

# Overall (mean across states, in normalised scale for summary).
mean_mec = mean(rmse_mec)
mean_ml  = mean(rmse_ml)
mean_hyb = mean(rmse_hyb)
println("\nMean RMSE across states:")
println("  Pure-Mechanistic : $(round(mean_mec; digits=4))")
println("  Pure-ML (LSTM)   : $(round(mean_ml;  digits=4))")
println("  Hybrid           : $(round(mean_hyb; digits=4))")

@assert rmse_hyb[end] < rmse_mec[end] "Hybrid did not correct the mAb mechanistic bias!"

# ── 8. Figure: truth vs three forecasts (mAb state, first test window) ─────────
MAB_IDX = 7   # index of mAb in state vector

b_plot  = 1   # first test window
t_all   = df_true_use.t

t0_plot   = n_train_win + b_plot
t_in_idx  = t0_plot:(t0_plot + LOOKBACK - 1)
t_ho_idx  = (t0_plot + LOOKBACK):(t0_plot + LOOKBACK + HORIZON - 1)
t_input   = t_all[t_in_idx]
t_horizon = t_all[t_ho_idx]

raw_true = Float64.(Matrix(df_true_use[:, Symbol.(STATE_COLS)])')  # N × T
raw_mec  = arr_mec_full'                                           # N × T

fig = Figure(size = (900, 560))
ax  = Axis(fig[1, 1];
    title  = "CHO Hybrid Model: mAb Forecast",
    xlabel = "Time (h)",
    ylabel = "mAb (mg/L)",
)

# Full true trajectory (light grey)
lines!(ax, t_all, raw_true[MAB_IDX, :];
    color = (:grey, 0.35), linewidth = 1.5, label = "True trajectory")

# Input window (blue)
lines!(ax, t_input, raw_true[MAB_IDX, t_in_idx];
    color = :steelblue, linewidth = 2.5, label = "Input window")

# Ground truth horizon (dashed black)
lines!(ax, t_horizon, Float64.(Y_true_phys[MAB_IDX, :, b_plot]);
    color = :black, linewidth = 2.5, linestyle = :dash, label = "Truth (horizon)")

# Imperfect mechanistic (orange)
lines!(ax, t_horizon, Float64.(Ŷ_mec_phys[MAB_IDX, :, b_plot]);
    color = :darkorange, linewidth = 2.5, label = "Pure mechanistic (biased)")

# Pure ML LSTM (red)
lines!(ax, t_horizon, Float64.(Ŷ_ml_phys[MAB_IDX, :, b_plot]);
    color = :crimson, linewidth = 2.5, label = "Pure ML (LSTM)")

# Hybrid (green)
lines!(ax, t_horizon, Float64.(Ŷ_hyb_phys[MAB_IDX, :, b_plot]);
    color = :darkgreen, linewidth = 2.5, label = "Hybrid (mech + residual ML)")

axislegend(ax; position = :lt, labelsize = 10)

# RMSE annotation box
rmse_text = "RMSE (mAb, mg/L)\n" *
            "  Mech : $(round(rmse_mec[MAB_IDX]; digits=2))\n" *
            "  ML   : $(round(rmse_ml[MAB_IDX];  digits=2))\n" *
            "  Hybrid: $(round(rmse_hyb[MAB_IDX]; digits=2))"
text!(ax, t_horizon[end] - 4, maximum(raw_true[MAB_IDX, :]) * 0.12;
    text = rmse_text, fontsize = 10, align = (:right, :bottom))

# Second panel: RMSE bar chart for all three approaches (mAb only for focus)
ax2 = Axis(fig[1, 2];
    title  = "Per-state RMSE comparison",
    xlabel = "State",
    ylabel = "RMSE (physical units)",
    xticks = (1:N_STATES, STATE_COLS),
    xticklabelrotation = π/6,
)

xs = 1:N_STATES
barplot!(ax2, xs .- 0.25, rmse_mec; width = 0.22, color = :darkorange, label = "Pure-Mech")
barplot!(ax2, xs,          rmse_ml;  width = 0.22, color = :crimson,    label = "Pure-ML")
barplot!(ax2, xs .+ 0.25,  rmse_hyb; width = 0.22, color = :darkgreen,  label = "Hybrid")
axislegend(ax2; position = :lt, labelsize = 9)

fig_path = joinpath(_PATH_TO_FIGS, "hybrid_cho.pdf")
save(fig_path, fig)
println("\nSaved: figs/hybrid_cho.pdf")

println("\nHYB_OK")
