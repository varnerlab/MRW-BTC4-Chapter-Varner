"""
    data_prep.jl — Sliding-window data preparation for deep learning forecasters.

    PRIMARY INTERFACE (stable; reused by Tasks 4.2 S4 and 4.5 hybrid):

        make_windows(df; lookback=24, horizon=12, test_frac=0.2)
            -> (Xtrain, Ytrain, Xtest, Ytest, scaler)

    Array shapes
    ============
    Let N = number of states (7 for CHO: V,X,Glc,Gln,Lac,Amm,mAb)
        L = lookback  (number of history timesteps fed as input)
        H = horizon   (number of future timesteps to predict)
        B = batch size (number of sliding-window samples)

    Xtrain : Float32 array  (N × L × Btrain)   — input windows (normalized)
    Ytrain : Float32 array  (N × H × Btrain)   — target windows (normalized)
    Xtest  : Float32 array  (N × L × Btest)
    Ytest  : Float32 array  (N × H × Btest)

    scaler : NamedTuple with fields
        μ   :: Vector{Float64}   (length N) — per-state mean
        σ   :: Vector{Float64}   (length N) — per-state std
        normalize   :: Function  matrix -> matrix  (apply z-score)
        denormalize :: Function  matrix -> matrix  (invert z-score)

    The scaler `denormalize` function accepts a matrix of shape (N, *) and
    returns values in physical units.  Tasks 4.2/4.5 can reuse the same
    scaler for their own forecast arrays.

    Normalization: z-score per state column computed on the TRAINING windows only.
    A minimum σ floor of 1e-8 guards against constant states.
"""

"""
    make_windows(df; lookback=24, horizon=12, test_frac=0.2)

Slide a window of length (lookback + horizon) over all state columns of `df`.
Normalizes each state with z-score statistics computed on training windows.
Holds out the *tail* of the time series as the test set.

# Arguments
- `df`         : DataFrame with columns [t, V, X, Glc, Gln, Lac, Amm, mAb]
- `lookback`   : number of past timesteps used as model input
- `horizon`    : number of future timesteps to forecast
- `test_frac`  : fraction of total windows reserved for test (tail)

# Returns
`(Xtrain, Ytrain, Xtest, Ytest, scaler)` — see module docstring for shapes.
"""
function make_windows(df; lookback::Int=24, horizon::Int=12, test_frac::Float64=0.2)

    state_cols = [:V, :X, :Glc, :Gln, :Lac, :Amm, :mAb]
    N = length(state_cols)

    # Raw data matrix: N × T  (states as rows, time as columns)
    raw = Float64.(Matrix(df[:, state_cols]))'  # N × T

    T = size(raw, 2)
    window = lookback + horizon

    # Build all windows
    n_windows = T - window + 1
    @assert n_windows > 0 "Trajectory too short for lookback=$lookback + horizon=$horizon"

    # Split train/test by window index (tail = test)
    n_test  = max(1, round(Int, test_frac * n_windows))
    n_train = n_windows - n_test

    @assert n_train > 0 "Not enough windows for training; reduce test_frac or window size."

    # Compute z-score statistics from training windows only
    # Collect all training-window timesteps for stats
    train_vals = hcat([raw[:, i:(i + lookback - 1)] for i in 1:n_train]...)  # N × (n_train*lookback)
    μ = vec(mean(train_vals, dims=2))         # length N
    σ = vec(std(train_vals, dims=2))          # length N
    σ = max.(σ, 1e-8)                         # floor to avoid div-by-zero

    normalize   = Z -> (Z .- μ) ./ σ          # works on N × * matrices
    denormalize = Z -> Z .* σ .+ μ            # inverse of above

    # Normalise full series
    norm_raw = normalize(raw)                  # N × T

    # Allocate arrays:  N × L × B  and  N × H × B
    Xtrain = Array{Float32}(undef, N, lookback, n_train)
    Ytrain = Array{Float32}(undef, N, horizon,  n_train)
    Xtest  = Array{Float32}(undef, N, lookback, n_test)
    Ytest  = Array{Float32}(undef, N, horizon,  n_test)

    for b in 1:n_train
        t0 = b  # window starts at timestep t0 (1-indexed)
        Xtrain[:, :, b] = Float32.(norm_raw[:, t0:(t0 + lookback - 1)])
        Ytrain[:, :, b] = Float32.(norm_raw[:, (t0 + lookback):(t0 + window - 1)])
    end

    for b in 1:n_test
        t0 = n_train + b
        Xtest[:, :, b] = Float32.(norm_raw[:, t0:(t0 + lookback - 1)])
        Ytest[:, :, b] = Float32.(norm_raw[:, (t0 + lookback):(t0 + window - 1)])
    end

    scaler = (
        μ           = μ,
        σ           = σ,
        normalize   = normalize,
        denormalize = denormalize,
    )

    return Xtrain, Ytrain, Xtest, Ytest, scaler
end
