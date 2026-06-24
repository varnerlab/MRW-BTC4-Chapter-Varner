"""
    s4.jl — Structured State-Space (S4 / HiPPO-LegS) forecaster as a custom Flux layer.

    Tested against Flux 0.16.x API (matches lstm.jl / Task 4.1).

    This is the §5 S4 counterpart to the LSTM forecaster.  It implements a discrete
    linear state-space model with a HiPPO-LegS initialization of the continuous
    (A, B) pair, bilinear (Tustin) discretization to (Ā, B̄), and a recurrence over
    the lookback window followed by a trainable readout.

    HiPPO-LegS construction + bilinear discretization adapted from
        CHEME-5820 Week-14 lab L14d  (src/Compute.jl, src/Types.jl):
        _build_legS_single / build_legS_matrices_mimo / discretize
    and the L14a HiPPO-SSM lecture example.

    Model (MIMO, block-diagonal, one LegS filter per input channel):

        continuous:  dx/dt = A x + B u
        discrete:    x_k   = Ā x_{k-1} + B̄ u_k          (rolled over the lookback)
        readout:     y     = C x_L (+ D u_L)            (final-state readout)

    `A`, `B` (hence `Ā`, `B̄`) are FIXED at their HiPPO-LegS values.  The readout
    `C` (and feedforward `D`) are trained by gradient descent — this is the stable
    S4/SSM training recipe used in L14a, lifted from a closed-form ridge solve to
    a Flux gradient loop so it composes with the rest of the deep-learning §.

    Array shapes (match make_windows / lstm.jl):
        input  : (N × lookback × batch)   — N = d_in = 7 CHO state-channels
        output : (N*horizon × batch)       — reshaped to (N × horizon × batch)

    PUBLIC INTERFACE
    ================
    S4Layer(n_states, lookback, horizon; h, Δt)   :: Flux-compatible layer
    train_s4(Xtrain, Ytrain; epochs, ...)         :: Chain (trained)
    forecast(model, Xtest)                        :: Array{Float32,3} (N × horizon × B)
"""

# ── HiPPO-LegS matrix construction (adapted from L14d src/Compute.jl) ────────
"""
    _build_legS_single(h) -> (A, b)

Single-channel LegS state matrix `A ∈ ℝ^{h×h}` and input vector `b ∈ ℝ^{h}` in
the 1-based sign convention of Gu et al. 2020 with the stable form
`dx/dt = A x + b u`.  Verbatim recipe from L14d.
"""
function _build_legS_single(h::Int)
    h >= 1 || error("h must be positive")
    A = zeros(Float64, h, h)
    b = zeros(Float64, h)
    for i in 1:h
        for k in 1:h
            if i > k
                A[i, k] = -sqrt((2i + 1) * (2k + 1))
            elseif i == k
                A[i, k] = -(i + 1)
            end
        end
        b[i] = sqrt(2i + 1)
    end
    return (A, b)
end

"""
    build_legS_matrices_mimo(h, d_in) -> (A, B)

Block-diagonal MIMO LegS state matrix `A ∈ ℝ^{(h·d_in)×(h·d_in)}` and input
matrix `B ∈ ℝ^{(h·d_in)×d_in}`.  Each input channel drives one independent
`h`-dim LegS filter.  Adapted from L14d.
"""
function build_legS_matrices_mimo(h::Int, d_in::Int)
    h >= 1    || error("h must be positive")
    d_in >= 1 || error("d_in must be positive")
    (A_single, b_single) = _build_legS_single(h)
    H = h * d_in
    A = zeros(Float64, H, H)
    B = zeros(Float64, H, d_in)
    for j in 1:d_in
        rows = ((j - 1) * h + 1):(j * h)
        A[rows, rows] = A_single
        B[rows, j]    = b_single
    end
    return (A, B)
end

"""
    discretize_bilinear(A, B; Δt) -> (Ā, B̄)

Bilinear (Tustin) discretization of the continuous pair `(A, B)` at step `Δt`:

    Ā = (I - Δt/2 A)⁻¹ (I + Δt/2 A)
    B̄ = (I - Δt/2 A)⁻¹ (Δt B)

Adapted from L14d `discretize`.
"""
function discretize_bilinear(A::AbstractMatrix, B::AbstractMatrix; Δt::Float64)
    Δt > 0 || error("Δt must be positive")
    H, W = size(A)
    H == W || error("A must be square, got $(size(A))")
    size(B, 1) == H || error("B has $(size(B, 1)) rows, expected $(H)")
    I_H = Matrix{Float64}(I, H, H)
    M   = I_H - (Δt / 2) .* A
    Ā   = M \ (I_H + (Δt / 2) .* A)
    B̄   = M \ (Δt .* B)
    return (Ā, B̄)
end

# ── S4Layer: Flux-compatible HiPPO-LegS state-space layer ────────────────────
"""
    S4Layer

A discrete HiPPO-LegS state-space layer.  Holds the FIXED discrete operators
`Ā` (H×H) and `B̄` (H×d_in) and the TRAINABLE readout `C` (n_out×H) and
feedforward `D` (n_out×d_in).

Calling the layer on a window `u :: (d_in × L × B)` rolls the recurrence
`x_k = Ā x_{k-1} + B̄ u_k` from `x₀ = 0` over the `L` lookback steps, then reads
out the FINAL state: `y = C x_L + D u_L`, returning `(n_out × B)`.

`n_out = n_states * horizon`, so the layer maps a lookback window directly to a
flattened multi-step forecast (reshaped to (N × horizon × B) by `forecast`).
"""
struct S4Layer{M1,M2,M3,M4}
    Ā::M1      # (H × H)        discrete state matrix      (FIXED, HiPPO)
    B̄::M2      # (H × d_in)     discrete input matrix      (FIXED, HiPPO)
    C::M3      # (n_out × H)    readout                    (TRAINABLE)
    D::M4      # (n_out × d_in) feedforward                (TRAINABLE)
end

# Only C and D are trained; Ā and B̄ stay frozen at their HiPPO-LegS values.
Flux.@layer S4Layer trainable=(C, D)

"""
    S4Layer(n_states, lookback, horizon; h=16, Δt=1.0, seed=42)

Construct an `S4Layer`.  Builds the block-diagonal HiPPO-LegS `(A, B)` for
`d_in = n_states` channels at per-channel order `h`, discretizes bilinearly at
step `Δt`, and initializes the trainable readout `C` (small random) and `D`
(zeros).  Output dimension is `n_states * horizon`.
"""
function S4Layer(n_states::Int, lookback::Int, horizon::Int;
                 h::Int = 16, Δt::Float64 = 1.0, seed::Int = 42)
    d_in  = n_states
    n_out = n_states * horizon
    (A, B) = build_legS_matrices_mimo(h, d_in)
    (Ā, B̄) = discretize_bilinear(A, B; Δt = Δt)
    H = h * d_in

    rng = MersenneTwister(seed)
    C = Float32.(0.1 .* randn(rng, n_out, H) ./ sqrt(H))   # small random readout
    D = zeros(Float32, n_out, d_in)                         # feedforward starts at 0

    return S4Layer(Float32.(Ā), Float32.(B̄), C, D)
end

"""
    (m::S4Layer)(u)

Apply the layer to a batch of lookback windows `u :: (d_in × L × B)`.
Rolls the discrete SSM recurrence over the `L` steps from zero initial state and
reads out the final hidden state.  Returns `(n_out × B)`.

Implemented with batched matrix multiplies so it is `Flux`/`Zygote`
differentiable (gradients flow into `C` and `D`; `Ā`/`B̄` are frozen).
"""
function (m::S4Layer)(u::AbstractArray{<:Real,3})
    d_in, L, B = size(u)
    H = size(m.Ā, 1)
    T = eltype(m.C)
    # roll the recurrence:  x_k = Ā x_{k-1} + B̄ u_k,   x_0 = 0
    x = zeros(T, H, B)
    @inbounds for k in 1:L
        uk = @view u[:, k, :]            # (d_in × B)
        x = m.Ā * x .+ m.B̄ * uk          # (H × B)
    end
    uL = @view u[:, L, :]                # (d_in × B) — last input step
    return m.C * x .+ m.D * uL           # (n_out × B)
end

# Zygote-friendly variant: gather the final input step without a view.
function (m::S4Layer)(u::AbstractMatrix{<:Real})
    # single time step convenience (d_in × B) -> one-step rollout
    return reshape(u, size(u, 1), 1, size(u, 2)) |> m
end

"""
    build_s4_model(n_states, lookback, horizon; h, Δt, seed) -> Chain

Wrap an `S4Layer` in a `Chain`.  The S4Layer already produces the flattened
`n_states*horizon` forecast, so the Chain is a thin wrapper kept for parity with
the LSTM `Chain` interface (Task 4.3 compares `Chain`s).
"""
function build_s4_model(n_states::Int, lookback::Int, horizon::Int;
                        h::Int = 16, Δt::Float64 = 1.0, seed::Int = 42)
    return Chain(S4Layer(n_states, lookback, horizon; h = h, Δt = Δt, seed = seed))
end

"""
    train_s4(Xtrain, Ytrain; epochs=300, h=16, Δt=1.0, lr=1e-2,
             batchsize=32, seed=42) -> Chain

Train the S4/HiPPO-LegS forecaster on windowed CHO data.  The HiPPO `(Ā, B̄)`
state-space operators are FIXED; the readout `C` and feedforward `D` are trained
by Adam gradient descent on the MSE between the final-state readout and the
flattened horizon target.

# Arguments
- `Xtrain`    : Float32 array (N × lookback × Btrain)   — normalized input windows
- `Ytrain`    : Float32 array (N × horizon  × Btrain)   — normalized targets
- `epochs`    : number of full passes over training data
- `h`         : per-channel LegS basis order (hidden dim per channel)
- `Δt`        : bilinear discretization timestep
- `lr`        : Adam learning rate
- `batchsize` : mini-batch size
- `seed`      : RNG seed

# Returns
Trained `Flux.Chain` (single `S4Layer`).
"""
function train_s4(
    Xtrain::AbstractArray{Float32,3},
    Ytrain::AbstractArray{Float32,3};
    epochs::Int    = 300,
    h::Int         = 16,
    Δt::Float64    = 1.0,
    lr::Float64    = 1e-2,
    batchsize::Int = 32,
    seed::Int      = 42,
)
    Random.seed!(seed)

    n_states, lookback, n_train = size(Xtrain)
    _,        horizon,  _       = size(Ytrain)

    model = build_s4_model(n_states, lookback, horizon; h = h, Δt = Δt, seed = seed)
    opt_state = Flux.setup(Adam(lr), model)

    # Flatten Ytrain to (N*H × Btrain) to match the layer output.
    Y_flat = reshape(Ytrain, n_states * horizon, n_train)

    loss_fn(m, x_batch, y_batch) = Flux.mse(m(x_batch), y_batch)

    H_total = h * n_states
    ρ = maximum(abs, eigvals(Float64.(model[1].Ā)))
    println("Training S4/HiPPO-LegS (Flux $(pkgversion(Flux))):  epochs=$epochs, h=$h, " *
            "N=$(H_total), Δt=$Δt, lr=$lr")
    println("  discrete Ā spectral radius (should be < 1): ", round(ρ; digits=6))
    println(lpad("epoch", 8), "  ", lpad("loss", 14))

    losses = Float64[]
    for epoch in 1:epochs
        idx = randperm(n_train)
        epoch_loss = 0.0
        n_batches  = 0
        for start in 1:batchsize:n_train
            batch_idx = idx[start:min(start + batchsize - 1, n_train)]
            x_b = Xtrain[:, :, batch_idx]
            y_b = Y_flat[:, batch_idx]
            l, grads = Flux.withgradient(loss_fn, model, x_b, y_b)
            Flux.update!(opt_state, model, grads[1])
            epoch_loss += l
            n_batches  += 1
        end
        avg_loss = epoch_loss / n_batches
        push!(losses, avg_loss)
        if epoch == 1 || epoch % 10 == 0 || epoch == epochs
            println(lpad(epoch, 8), "  ", lpad(round(avg_loss; digits=6), 14))
        end
    end

    return model
end

"""
    forecast(model, Xtest) -> Array{Float32,3}   (N_states × horizon × Btest)

Run the trained S4 model on `Xtest` (N × lookback × Btest) and return
predictions shaped (N × horizon × Btest), matching `Ytest` from `make_windows`.
"""
function forecast(model::Chain, Xtest::AbstractArray{Float32,3})
    n_states, lookback, n_test = size(Xtest)
    ŷ_flat = model(Xtest)                       # (N*horizon × Btest)
    n_out   = size(ŷ_flat, 1)
    horizon = n_out ÷ n_states
    @assert n_out == n_states * horizon
    return reshape(ŷ_flat, n_states, horizon, n_test)
end
