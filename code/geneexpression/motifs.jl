# code/geneexpression/motifs.jl
#
# Hill input functions and four canonical network-motif ODEs.
# All motif simulators return a DataFrame keyed on `t`.
#
# Motifs covered:
#   1. Negative autoregulation (NAR) vs. simple regulation — response-time speed-up
#   2. Coherent type-1 FFL (AND gate)      — persistence detector / sign-sensitive delay
#   3. Incoherent type-1 FFL               — pulse generator
#   4. Goodwin-style negative-feedback oscillator
#
# Requires: DifferentialEquations (Tsit5 / ODEProblem), DataFrames — loaded by Include.jl.

# ── Hill input functions ──────────────────────────────────────────────────────

"""
    hill_act(x; K=1.0, n=2)

Activating Hill function: x^n / (K^n + x^n).
Returns a value in [0, 1].
"""
hill_act(x; K=1.0, n=2) = x^n / (K^n + x^n)

"""
    hill_rep(x; K=1.0, n=2)

Repressing Hill function: K^n / (K^n + x^n).
Returns a value in [0, 1].
"""
hill_rep(x; K=1.0, n=2) = K^n / (K^n + x^n)

# ── 1. Negative autoregulation vs. simple regulation ─────────────────────────

"""
    simulate_nar(; tspan, α, K, n, Y_st) -> DataFrame(t, simple, nar)

Negative autoregulation (NAR) vs. simple (unregulated) regulation.

The canonical Alon comparison requires **both circuits to share the same steady
state** Y_st so that the speed difference is attributable to the feedback alone.

  Simple:  dY/dt = β_simple - α·Y                  β_simple = α·Y_st
  NAR:     dY/dt = β_nar·hill_rep(Y; K, n) - α·Y   β_nar = α·Y_st / hill_rep(Y_st; K, n)

At steady state both settle to Y_st, but NAR reaches ½·Y_st faster because
the high-production phase at Y≈0 is unattenuated by feedback only in NAR
(the repressor is absent), giving a burst of synthesis that quickly lifts Y
to the half-max threshold before self-repression kicks in.
"""
function simulate_nar(; tspan=(0.0, 10.0), α=1.0, K=0.5, n=2, Y_st=1.0)
    # Derive production rates so both circuits share steady state Y_st
    β_simple = α * Y_st                             # simple: Y_ss = β/α = Y_st
    β_nar    = α * Y_st / hill_rep(Y_st; K=K, n=n) # NAR: β·h(Y_st) = α·Y_st

    f!(du, u, p, t) = begin
        du[1] = β_simple - α * u[1]                          # simple regulation
        du[2] = β_nar * hill_rep(u[2]; K=K, n=n) - α * u[2] # negative autoregulation
    end
    sol = solve(ODEProblem(f!, [0.0, 0.0], tspan), Tsit5();
                saveat=0.02, abstol=1e-9, reltol=1e-9)
    DataFrame(t=sol.t,
              simple=getindex.(sol.u, 1),
              nar=getindex.(sol.u, 2),
              Y_st=fill(Y_st, length(sol.t)))
end

# ── 2. Coherent type-1 FFL (C1-FFL, AND gate) ────────────────────────────────

"""
    simulate_c1ffl(; tspan, βy, βz, αy, αz, Kxy, Kxz, Kyz, n,
                    pulse, step_on) -> DataFrame(t, X, Y, Z)

Coherent type-1 FFL with AND-gate integration at Z:
  X → Y,  X → Z,  Y → Z   (all activating)

Persistence-detector / sign-sensitive delay:
  • A brief input pulse (X ON for < 1/αy) is FILTERED — Y cannot accumulate
    fast enough, so the AND gate blocks Z.
  • A sustained step keeps X ON long enough for Y to rise; once Y > K_yz the
    gate opens and Z turns ON after a delay of ~1/αy.

Two-phase input signal:
  X(t) = 1  for pulse[1] ≤ t ≤ pulse[2]   (transient pulse, filtered)
  X(t) = 1  for t ≥ step_on               (sustained step, passes through)
  X(t) = 0  otherwise
"""
function simulate_c1ffl(; tspan=(0.0, 20.0),
                          βy=2.0, βz=2.0, αy=1.0, αz=1.0,
                          Kxy=0.3, Kxz=0.3, Kyz=0.5,
                          n=2,
                          pulse=(1.0, 1.8),   # short pulse: ~0.8 time units
                          step_on=8.0)
    Xsig(t) = (pulse[1] ≤ t ≤ pulse[2]) || (t ≥ step_on) ? 1.0 : 0.0

    f!(du, u, p, t) = begin
        X = Xsig(t)
        du[1] = βy * hill_act(X; K=Kxy, n=n) - αy * u[1]                          # Y
        du[2] = βz * hill_act(X; K=Kxz, n=n) * hill_act(u[1]; K=Kyz, n=n) - αz * u[2]  # Z (AND)
    end
    sol = solve(ODEProblem(f!, [0.0, 0.0], tspan), Tsit5();
                saveat=0.05, abstol=1e-9, reltol=1e-9)
    DataFrame(t=sol.t,
              X=Xsig.(sol.t),
              Y=getindex.(sol.u, 1),
              Z=getindex.(sol.u, 2))
end

# ── 3. Incoherent type-1 FFL (I1-FFL) ────────────────────────────────────────

"""
    simulate_i1ffl(; tspan, βy, βz, αy, αz, Kxy, Kxz, Kyz, n,
                    step_on) -> DataFrame(t, X, Y, Z)

Incoherent type-1 FFL:
  X → Y,  X → Z,  Y ⊣ Z   (X activates both; Y represses Z)

Pulse-generator: when X turns ON —
  1. X immediately activates Z (direct path, fast rise).
  2. X also activates Y (indirect path, slower rise via 1/αy delay).
  3. Once Y accumulates it represses Z, pulling Z back toward a lower plateau.
  Result: Z exhibits a transient pulse — maximum(Z) > 1.5 · Z[end].

Parameters tuned so Y builds on the time-scale 1/αy = 1 and the pulse is
clearly visible within tspan.
"""
function simulate_i1ffl(; tspan=(0.0, 15.0),
                          βy=2.0, βz=3.0, αy=1.0, αz=1.0,
                          Kxy=0.3, Kxz=0.3, Kyz=0.3,
                          n=2, step_on=1.0)
    Xsig(t) = t ≥ step_on ? 1.0 : 0.0

    f!(du, u, p, t) = begin
        X = Xsig(t)
        du[1] = βy * hill_act(X; K=Kxy, n=n) - αy * u[1]                          # Y
        du[2] = βz * hill_act(X; K=Kxz, n=n) * hill_rep(u[1]; K=Kyz, n=n) - αz * u[2]  # Z
    end
    sol = solve(ODEProblem(f!, [0.0, 0.0], tspan), Tsit5();
                saveat=0.05, abstol=1e-9, reltol=1e-9)
    DataFrame(t=sol.t,
              X=Xsig.(sol.t),
              Y=getindex.(sol.u, 1),
              Z=getindex.(sol.u, 2))
end

# ── 4. Goodwin-style negative-feedback oscillator ─────────────────────────────

"""
    simulate_oscillator(; tspan, a, b, c, d, e, f, K, n) -> DataFrame(t, X, Y, Z)

Three-variable Goodwin (1965) negative-feedback oscillator:
  dX/dt = a · hill_rep(Z; K, n) - b·X
  dY/dt = c·X - d·Y
  dZ/dt = e·Y - f·Z

A high Hill coefficient (n ≥ 9) is required for sustained or strongly-damped
oscillations.  The repressor Z feeds back to suppress X with a sigmoidal
nonlinearity; the three-stage delay (X → Y → Z ⊣ X) enables limit-cycle
behaviour.

Parameters here give sustained oscillations with period ≈ 20–30 time units;
two clear peaks are visible within tspan=(0, 100).
"""
function simulate_oscillator(; tspan=(0.0, 100.0),
                               a=8.0, b=0.5, c=1.0, d=0.2, e=1.0, f=0.2,
                               K=1.0, n=10)
    g!(du, u, p, t) = begin
        du[1] = a * hill_rep(u[3]; K=K, n=n) - b * u[1]  # X
        du[2] = c * u[1] - d * u[2]                        # Y
        du[3] = e * u[2] - f * u[3]                        # Z
    end
    sol = solve(ODEProblem(g!, [0.5, 0.5, 0.5], tspan), Tsit5();
                saveat=0.1, abstol=1e-9, reltol=1e-9)
    DataFrame(t=sol.t,
              X=getindex.(sol.u, 1),
              Y=getindex.(sol.u, 2),
              Z=getindex.(sol.u, 3))
end
