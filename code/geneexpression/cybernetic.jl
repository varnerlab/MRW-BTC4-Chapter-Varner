# cybernetic.jl
# Kompala cybernetic model for diauxic growth on two substrates (glucose + xylose).
# Distilled from: varnerlab/Kompala-Model-LP-Paper
#   code/discrete-cellmass-simulation/src/Kompala.jl  (ODE balances + control laws)
#   code/discrete-cellmass-simulation/conf/Parameters-Fig-6.toml  (parameter set)
#   code/discrete-cellmass-simulation/plotme_Fig_6-Kompala.jl     (initial conditions)
#
# Model state: x = [S1, S2, E1, E2, C]
#   S1 = glucose concentration  (g/L)
#   S2 = xylose concentration   (g/L)
#   E1 = enzyme 1 (glucose pathway) concentration
#   E2 = enzyme 2 (xylose pathway) concentration
#   C  = cell mass concentration (gDW/L)
#
# Cybernetic control law (Kompala matching law):
#   rG_i = μmax_i * (E_i/emax_i) * (S_i / (K_i + S_i))   [substrate-specific growth rate]
#   u_i  = rG_i / Σ rG_j    [enzyme SYNTHESIS allocation — "matching law"]
#   v_i  = rG_i / max(rG_j) [enzyme ACTIVITY modulation]
#   μ    = Σ rG_i * v_i      [net specific growth rate]
#
# Reference: Kompala et al. (1984) Biotechnology & Bioengineering.
# Adapted for DifferentialEquations.jl (Tsit5) interface.

# ── Parameters (Fig-6 two-substrate case: glucose vs xylose) ─────────────────

const _KOMPALA_PARAMS = (
    # Substrate 1: glucose
    μmax1 = 1.08,    # h⁻¹  max specific growth rate
    K1    = 0.01,    # g/L  Monod saturation constant
    Y1    = 0.505,   # gDW/gS  yield
    τ1    = 0.6226053639846744,   # h  enzyme synthesis time constant
    β1    = 0.05,    # h⁻¹  enzyme degradation rate
    emax1 = 1.421375085091899,    # enzyme capacity (normalisation)

    # Substrate 2: xylose
    μmax2 = 0.82,    # h⁻¹
    K2    = 0.20,    # g/L
    Y2    = 0.58,    # gDW/gS
    τ2    = 0.6226053639846744,
    β2    = 0.05,    # h⁻¹
    emax2 = 1.846153846153846,

    # Maintenance / death
    kd    = 0.022,   # h⁻¹  (from utilitybalances death term)
)

# ── Cybernetic ODE right-hand side ───────────────────────────────────────────

function _cybernetic_rhs!(dx, x, p, t)
    S1, S2, E1, E2, C = x

    # Clamp to avoid negative values driving instability
    S1 = max(S1, 0.0)
    S2 = max(S2, 0.0)
    E1 = max(E1, 0.0)
    E2 = max(E2, 0.0)
    C  = max(C,  0.0)

    # ── Substrate-specific growth rates ─────────────────────────────────────
    rG1 = p.μmax1 * (E1 / p.emax1) * (S1 / (p.K1 + S1))
    rG2 = p.μmax2 * (E2 / p.emax2) * (S2 / (p.K2 + S2))

    rG_sum = rG1 + rG2
    rG_max = max(rG1, rG2)

    # ── Cybernetic control variables ─────────────────────────────────────────
    # u_i: fraction of enzyme-synthesis resources allocated to pathway i
    #      (matching law: proportional to relative return rG_i / Σ rG)
    # v_i: relative enzyme activity (proportional to rG_i / max rG)
    if rG_sum < 1e-12
        u1, u2 = 0.5, 0.5   # quiescent: distribute equally (doesn't matter)
    else
        u1 = rG1 / rG_sum
        u2 = rG2 / rG_sum
    end

    if rG_max < 1e-12
        v1, v2 = 0.0, 0.0
    else
        v1 = rG1 / rG_max
        v2 = rG2 / rG_max
    end

    # ── Enzyme synthesis rates (1/τ when substrate present) ─────────────────
    rE1 = (S1 > 1e-6) ? (1.0 / p.τ1) : 0.0
    rE2 = (S2 > 1e-6) ? (1.0 / p.τ2) : 0.0

    # ── Net specific growth rate ─────────────────────────────────────────────
    μ = rG1 * v1 + rG2 * v2

    # ── Balances ─────────────────────────────────────────────────────────────
    # Substrate balances:  dS/dt = -(1/Y) * rG_i * v_i * C
    dx[1] = -(1.0 / p.Y1) * rG1 * v1 * C
    dx[2] = -(1.0 / p.Y2) * rG2 * v2 * C

    # Enzyme balances:  dE/dt = rE_i * u_i - (μ + β_i) * E_i
    dx[3] = rE1 * u1 - (μ + p.β1) * E1
    dx[4] = rE2 * u2 - (μ + p.β2) * E2

    # Cell mass balance:  dC/dt = (μ - kd) * C
    dx[5] = (μ - p.kd) * C

    return nothing
end

# ── Public interface ──────────────────────────────────────────────────────────

"""
    simulate_diauxie(; tspan=(0.0, 20.0)) -> DataFrame

Simulate the Kompala cybernetic model for diauxic growth on glucose (S1) and
xylose (S2).  Returns a DataFrame with columns: t, biomass, S1, S2.

The matching-law cybernetic variables u_i and v_i concentrate enzyme-synthesis
resources on the preferred substrate (glucose, higher μmax / lower K), producing
the classic diauxic phenotype: S1 depletes first, followed by a lag before S2
consumption begins.

Initial conditions follow plotme_Fig_6-Kompala.jl:
  S1=0.5 g/L, S2=2.5 g/L, E1=0.9·emax1, E2=0.18·emax2, C=4e-3 gDW/L
"""
function simulate_diauxie(; tspan=(0.0, 20.0))
    p = _KOMPALA_PARAMS

    # Initial conditions (from plotme_Fig_6-Kompala.jl)
    x0 = [
        0.5,              # S1: glucose (g/L)
        2.5,              # S2: xylose  (g/L)
        0.90 * p.emax1,   # E1: glucose enzyme
        0.18 * p.emax2,   # E2: xylose enzyme
        4e-3,             # C:  cell mass (gDW/L)
    ]

    prob = ODEProblem(_cybernetic_rhs!, x0, tspan, p)
    sol  = solve(prob, Tsit5();
                 abstol=1e-8, reltol=1e-8,
                 saveat=range(tspan[1], tspan[2]; length=2001),
                 isoutofdomain=(x, p, t) -> any(x .< -1e-6))

    t       = sol.t
    biomass = [max(u[5], 0.0) for u in sol.u]
    S1      = [max(u[1], 0.0) for u in sol.u]
    S2      = [max(u[2], 0.0) for u in sol.u]

    return DataFrame(t=t, biomass=biomass, S1=S1, S2=S2)
end
