# cho_model.jl
# Mechanistic fed-batch CHO antibody-production model
#
# Adapted from:
#   CHEME-5820 Week-12 L12d src/{Types,Parameters,Kinetics,MassBalances}.jl
#   and fedbatch-cho-model.md (governing equations and feed policy)
#
# State vector (brief §2 interface):
#   u = [V, X, Glc, Gln, Lac, Amm, mAb]
#   units: V (L), X (gDW/L), Glc/Gln/Lac/Amm (mM), mAb (mg/L)
#
# Feed policy: glucose-triggered hysteresis (per fedbatch-cho-model.md §10)
#   Feed ON  when Glc drops below Glc_min
#   Feed OFF when Glc rises above Glc_max
#
# Rate laws: Monod growth with dual by-product inhibition (lactate, ammonia),
#   by-product-driven death rate, Luedeking-Piret mAb productivity,
#   yield-based substrate uptake. (Xing et al. 2010 parameter basis.)

# ---------------------------------------------------------------------------
# Parameter NamedTuple
# ---------------------------------------------------------------------------

"""
    default_cho_params() -> NamedTuple

Return default kinetic, yield, and feed-policy parameters for the CHO fed-batch
ODE.  All values are adapted from the L12d Parameters.jl defaults (Xing et al.
2010) with a mutable `feed_on` field embedded as a Ref for callback mutation.

Units
- rates:      1/h
- Monod/inhib constants: mM
- yield on biomass:  gDW/mmol
- yield on product:  mg/mmol
- feed rate:  L/h
- feed concs: mM
"""
function default_cho_params()
    return (
        # --- growth kinetics (Xing et al. 2010 Table 2–3) ---
        mu_max   = 0.029,        # max specific growth rate (1/h)
        K_glc    = 0.10,         # Monod constant, glucose (mM)
        K_gln    = 0.05,         # Monod constant, glutamine (mM)
        K_I_lac  = 43.0,         # inhibition constant, lactate (mM)
        K_I_amm  = 6.5,          # inhibition constant, ammonia (mM)

        # --- death rate (Xing et al. 2010 Table 2–3) ---
        k_d      = 0.016,        # maximum specific death rate (1/h)
        KD_lac   = 45.8,         # half-sat for lactate in death term (mM)
        KD_amm   = 6.5,          # half-sat for ammonia in death term (mM)

        # --- Luedeking–Piret mAb productivity ---
        alpha_P  = 100.0,        # growth-associated coefficient (mg/gDW)
        beta_P   = 5.0,          # non-growth-associated coefficient (mg/gDW/h)

        # --- yield coefficients ---
        Y_X_glc  = 0.070,        # biomass yield on glucose  (gDW/mmol)
        Y_X_gln  = 0.210,        # biomass yield on glutamine (gDW/mmol)
        Y_lac_glc = 1.23,        # lactate yield on glucose  (mmol/mmol)
        Y_amm_gln = 0.67,        # ammonia yield on glutamine (mmol/mmol)

        # --- feed concentrations ---
        S_glc_f  = 500.0,        # glucose in feed (mM)
        S_gln_f  = 167.0,        # glutamine in feed (mM)

        # --- feed policy ---
        F_max    = 0.005,        # maximum feed rate (L/h)  ~5 mL/h into a 1 L vessel
        Glc_min  = 10.0,         # glucose threshold: feed turns ON  (mM)
        Glc_max  = 30.0,         # glucose threshold: feed turns OFF (mM)

        # --- mutable feed state (Ref so callbacks can mutate inside NamedTuple) ---
        feed_on  = Ref(0.0),     # 0.0 = off, 1.0 = on
    )
end

# ---------------------------------------------------------------------------
# Rate-law helper functions (adapted from L12d Kinetics.jl, de-typed)
# ---------------------------------------------------------------------------

"""
    _growth_rate(Glc, Gln, Lac, Amm, p) -> Float64

Monod growth with dual by-product inhibition (lactate, ammonia).
"""
function _growth_rate(Glc, Gln, Lac, Amm, p)
    monod_glc = Glc / (p.K_glc + Glc)
    monod_gln = Gln / (p.K_gln + Gln)
    inhib_lac = p.K_I_lac / (p.K_I_lac + Lac)
    inhib_amm = p.K_I_amm / (p.K_I_amm + Amm)
    return p.mu_max * monod_glc * monod_gln * inhib_lac * inhib_amm
end

"""
    _death_rate(Lac, Amm, p) -> Float64

By-product-driven death rate (Xing et al. 2010 form).
"""
function _death_rate(Lac, Amm, p)
    return p.k_d * (Lac / (p.KD_lac + Lac)) * (Amm / (p.KD_amm + Amm))
end

"""
    _product_rate(mu, p) -> Float64

Luedeking–Piret specific mAb productivity (mg/gDW/h).
"""
function _product_rate(mu, p)
    return p.alpha_P * mu + p.beta_P
end

# ---------------------------------------------------------------------------
# ODE right-hand side (brief interface)
# ---------------------------------------------------------------------------

"""
    cho_rhs!(du, u, p, t)

In-place RHS for the 7-state fed-batch CHO system.

State: `u = [V, X, Glc, Gln, Lac, Amm, mAb]`
- V   : reactor volume (L)
- X   : viable biomass concentration (gDW/L)
- Glc : glucose concentration (mM)
- Gln : glutamine concentration (mM)
- Lac : lactate concentration (mM)
- Amm : ammonia concentration (mM)
- mAb : monoclonal antibody concentration (mg/L)

`p` must respond to the fields returned by `default_cho_params()`.
"""
function cho_rhs!(du, u, p, t)

    # unpack state
    V, X, Glc, Gln, Lac, Amm, mAb = u

    # clamp to prevent negative-concentration artefacts
    V   = max(V,   1e-9)
    X   = max(X,   0.0)
    Glc = max(Glc, 0.0)
    Gln = max(Gln, 0.0)
    Lac = max(Lac, 0.0)
    Amm = max(Amm, 0.0)
    mAb = max(mAb, 0.0)

    # feed rate and dilution rate
    F = p.feed_on[] * p.F_max
    D = F / V

    # kinetic rates
    mu    = _growth_rate(Glc, Gln, Lac, Amm, p)
    mu_d  = _death_rate(Lac, Amm, p)
    q_P   = _product_rate(mu, p)
    q_glc = mu / p.Y_X_glc          # specific glucose uptake (mmol/gDW/h)
    q_gln = mu / p.Y_X_gln          # specific glutamine uptake (mmol/gDW/h)
    q_lac = p.Y_lac_glc * q_glc     # specific lactate formation (mmol/gDW/h)
    q_amm = p.Y_amm_gln * q_gln     # specific ammonia formation (mmol/gDW/h)

    # mass balances (fed-batch, per fedbatch-cho-model.md §1)
    du[1] = F                                         # dV/dt  (L/h)
    du[2] = (mu - mu_d - D) * X                       # dX/dt  (gDW/L/h)
    du[3] = D * (p.S_glc_f - Glc) - q_glc * X        # dGlc/dt (mM/h)
    du[4] = D * (p.S_gln_f - Gln) - q_gln * X        # dGln/dt (mM/h)
    du[5] = q_lac * X - D * Lac                       # dLac/dt (mM/h)
    du[6] = q_amm * X - D * Amm                       # dAmm/dt (mM/h)
    du[7] = q_P * X - D * mAb                         # dmAb/dt (mg/L/h)

    return nothing
end

# ---------------------------------------------------------------------------
# Feed callbacks (glucose-triggered hysteresis)
# ---------------------------------------------------------------------------

"""
    _build_feed_callbacks(p) -> CallbackSet

ContinuousCallback pair implementing glucose-triggered square-wave feeding:
- Feed ON  when Glc drops below `p.Glc_min` (downcrossing of u[3] - Glc_min)
- Feed OFF when Glc rises above `p.Glc_max` (upcrossing  of u[3] - Glc_max)

State index 3 = Glc (brief ordering: [V, X, Glc, Gln, Lac, Amm, mAb]).
"""
function _build_feed_callbacks(p)

    # Feed turns ON when glucose drops below Glc_min (downcrossing)
    cond_on(u, t, integrator)  = u[3] - integrator.p.Glc_min
    affect_on!(integrator)     = (integrator.p.feed_on[] = 1.0)
    cb_on  = ContinuousCallback(cond_on, nothing, affect_on!)

    # Feed turns OFF when glucose rises above Glc_max (upcrossing)
    cond_off(u, t, integrator) = u[3] - integrator.p.Glc_max
    affect_off!(integrator)    = (integrator.p.feed_on[] = 0.0)
    cb_off = ContinuousCallback(cond_off, affect_off!, nothing)

    return CallbackSet(cb_on, cb_off)
end

# ---------------------------------------------------------------------------
# Convenience simulation entry-point
# ---------------------------------------------------------------------------

"""
    simulate_cho(; tspan=(0.0, 240.0), saveat=1.0) -> DataFrame

Simulate the CHO fed-batch model with default parameters and return a DataFrame
with columns `t, V, X, Glc, Gln, Lac, Amm, mAb`.

Initial conditions (typical inoculation):
- V₀   = 1.0 L
- X₀   = 0.2 gDW/L  (~4 × 10⁵ cells/mL × 5 × 10⁻¹⁰ gDW/cell)
- Glc₀ = 5.0  mM  (below Glc_min so feed engages at t=0)
- Gln₀ = 4.0  mM
- Lac₀ = 0.0  mM
- Amm₀ = 0.0  mM
- mAb₀ = 0.0  mg/L
"""
function simulate_cho(; tspan=(0.0, 240.0), saveat=1.0)

    p = default_cho_params()
    # Glc₀ < Glc_min ⟹ feed should be ON at t = 0.
    # ContinuousCallback only catches a downcrossing mid-run, so we initialise
    # feed_on = 1 explicitly when starting below the threshold.
    p.feed_on[] = 1.0

    u0 = [
        1.0,   # V   (L)
        0.2,   # X   (gDW/L)
        5.0,   # Glc (mM)  — below Glc_min so feed is ON from t=0
        4.0,   # Gln (mM)
        0.0,   # Lac (mM)
        0.0,   # Amm (mM)
        0.0,   # mAb (mg/L)
    ]

    prob = ODEProblem(cho_rhs!, u0, tspan, p)
    cb   = _build_feed_callbacks(p)
    sol  = solve(prob, Tsit5();
                 callback  = cb,
                 saveat    = saveat,
                 reltol    = 1e-6,
                 abstol    = 1e-8,
                 isoutofdomain = (u, p, t) -> any(x -> x < -1e-6, u))

    # build DataFrame with brief-specified column names
    t_arr   = sol.t
    states  = hcat(sol.u...)'  # n_timepoints × 7

    return DataFrame(
        t   = t_arr,
        V   = states[:, 1],
        X   = states[:, 2],
        Glc = states[:, 3],
        Gln = states[:, 4],
        Lac = states[:, 5],
        Amm = states[:, 6],
        mAb = states[:, 7],
    )
end
