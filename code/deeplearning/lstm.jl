"""
    lstm.jl — Flux LSTM forecaster for multi-step sequence-to-sequence prediction.

    Tested against Flux 0.16.x API.

    Model architecture:
        Chain(LSTM(N_in => hidden), LastStep(), Dense(hidden => N_out*horizon))

    Input  shape: (N_states × lookback × batch)   — matches make_windows Xtrain
    Output shape: (N_states*horizon × batch)       — reshaped to (N_states × horizon × batch)

    PUBLIC INTERFACE
    ================
    train_lstm(Xtrain, Ytrain; epochs, hidden, lr, batchsize, seed) :: Chain
    forecast(model, Xtest)  :: Matrix{Float32}   (N_states × horizon × Btest) reshaped to (N×H×B)
"""

# ── custom layer: select the last time-step of an LSTM 3-D output ──────────
struct LastStep end
(::LastStep)(x::AbstractArray{<:Any,3}) = x[:, end, :]  # (hidden × T × B) -> (hidden × B)
Flux.@layer LastStep

"""
    build_lstm_model(n_states, lookback, horizon; hidden=32) -> Chain

Construct the LSTM model.  The Dense readout produces `n_states * horizon`
outputs which `forecast` reshapes to (n_states × horizon × batch).
"""
function build_lstm_model(n_states::Int, lookback::Int, horizon::Int; hidden::Int=32)
    Chain(
        LSTM(n_states => hidden),
        LastStep(),
        Dense(hidden => n_states * horizon),
    )
end

"""
    train_lstm(Xtrain, Ytrain; epochs=200, hidden=32, lr=1e-3,
               batchsize=32, seed=42) -> Chain

Train an LSTM forecaster on windowed CHO data.

# Arguments
- `Xtrain`    : Float32 array (N × lookback × Btrain)
- `Ytrain`    : Float32 array (N × horizon  × Btrain)
- `epochs`    : number of full passes over training data
- `hidden`    : LSTM hidden size
- `lr`        : ADAM learning rate
- `batchsize` : mini-batch size (number of windows per gradient step)
- `seed`      : RNG seed for reproducibility

# Returns
Trained `Flux.Chain`.
"""
function train_lstm(
    Xtrain::AbstractArray{Float32,3},
    Ytrain::AbstractArray{Float32,3};
    epochs::Int    = 200,
    hidden::Int    = 32,
    lr::Float64    = 1e-3,
    batchsize::Int = 32,
    seed::Int      = 42,
)
    Random.seed!(seed)

    n_states, lookback, n_train = size(Xtrain)
    _,        horizon,  _       = size(Ytrain)

    model = build_lstm_model(n_states, lookback, horizon; hidden=hidden)
    opt   = Adam(lr)
    opt_state = Flux.setup(opt, model)

    # Flatten Ytrain target to (N*H × Btrain) to match Dense output
    Y_flat = reshape(Ytrain, n_states * horizon, n_train)   # (N*H × B)

    # MSE loss: reset LSTM state before each forward pass
    function loss_fn(m, x_batch, y_batch)
        Flux.reset!(m)
        ŷ = m(x_batch)            # (N*H × batchsize)
        return Flux.mse(ŷ, y_batch)
    end

    println("Training LSTM (Flux $(pkgversion(Flux))):  epochs=$epochs, hidden=$hidden, lr=$lr")
    println(lpad("epoch", 8), "  ", lpad("loss", 14))

    for epoch in 1:epochs
        # Shuffle batch indices each epoch
        idx = randperm(n_train)
        epoch_loss = 0.0
        n_batches  = 0

        for start in 1:batchsize:n_train
            batch_idx = idx[start:min(start + batchsize - 1, n_train)]
            x_b = Xtrain[:, :, batch_idx]    # (N × L × b)
            y_b = Y_flat[:, batch_idx]        # (N*H × b)

            l, grads = Flux.withgradient(loss_fn, model, x_b, y_b)
            Flux.update!(opt_state, model, grads[1])

            epoch_loss += l
            n_batches  += 1
        end

        avg_loss = epoch_loss / n_batches

        # Print every 10 epochs
        if epoch == 1 || epoch % 10 == 0 || epoch == epochs
            println(lpad(epoch, 8), "  ", lpad(round(avg_loss; digits=6), 14))
        end
    end

    return model
end

"""
    forecast(model, Xtest) -> Array{Float32,3}   (N_states × horizon × Btest)

Run the trained LSTM on `Xtest` (N × lookback × Btest) and return
predictions shaped (N × horizon × Btest), matching Ytest from `make_windows`.
"""
function forecast(model::Chain, Xtest::AbstractArray{Float32,3})
    n_states, lookback, n_test = size(Xtest)
    Flux.reset!(model)
    ŷ_flat = model(Xtest)                    # (N*horizon × Btest)

    # Infer horizon from model output size
    n_out   = size(ŷ_flat, 1)
    horizon = n_out ÷ n_states
    @assert n_out == n_states * horizon

    return reshape(ŷ_flat, n_states, horizon, n_test)
end
