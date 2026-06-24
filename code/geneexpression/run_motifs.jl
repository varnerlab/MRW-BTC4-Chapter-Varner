# code/geneexpression/run_motifs.jl
#
# Simulate the four network-motif models and write a 4-panel figure to
# code/figs/motifs.pdf.
#
# Usage (from code/):
#   julia --project=. geneexpression/run_motifs.jl
#
# Verification checklist printed to stdout:
#   (a) NAR half-rise time < simple half-rise time
#   (b) C1-FFL filters brief X pulse (Z_pulse_max ≈ 0) but passes sustained step
#   (c) I1-FFL: maximum(Z) > 1.5 * Z[end]   (pulse signature)
#   (d) Oscillator: X has ≥ 2 peaks over the window

include(joinpath(@__DIR__, "..", "Include.jl"))
include(joinpath(@__DIR__, "motifs.jl"))

# ── Simulate all four motifs ──────────────────────────────────────────────────

nar = simulate_nar(tspan=(0.0, 10.0), α=1.0, K=0.5, n=2, Y_st=1.0)
c1  = simulate_c1ffl(tspan=(0.0, 20.0), pulse=(1.0, 1.8), step_on=8.0)
i1  = simulate_i1ffl(tspan=(0.0, 15.0))
osc = simulate_oscillator(tspan=(0.0, 100.0))

# ── Signature verification ────────────────────────────────────────────────────

## (a) NAR vs simple: half-rise times toward the shared steady state Y_st
Y_st   = nar.Y_st[1]
half_Y = 0.5 * Y_st

# First time index where each circuit exceeds ½·Y_st
idx_simple = findfirst(y -> y ≥ half_Y, nar.simple)
idx_nar    = findfirst(y -> y ≥ half_Y, nar.nar)

t_half_simple = isnothing(idx_simple) ? Inf : nar.t[idx_simple]
t_half_nar    = isnothing(idx_nar)    ? Inf : nar.t[idx_nar]

println("── (a) NAR vs simple ──────────────────────────────────────")
println("   Shared steady state Y_st  = ", Y_st)
println("   t½ simple = ", round(t_half_simple; digits=3))
println("   t½ NAR    = ", round(t_half_nar;    digits=3))
println("   NAR faster? ", t_half_nar < t_half_simple)
@assert t_half_nar < t_half_simple "NAR half-rise time should be less than simple (got $(t_half_nar) vs $(t_half_simple))"

## (b) C1-FFL persistence detector
# Z during the brief pulse window (t ∈ pulse)
pulse_mask = (c1.t .>= 1.0) .& (c1.t .<= 1.8)
step_mask  = c1.t .>= 10.0   # well after sustained step_on=8.0

Z_during_pulse = c1.Z[pulse_mask]
Z_during_step  = c1.Z[step_mask]

Z_pulse_max = isempty(Z_during_pulse) ? 0.0 : maximum(Z_during_pulse)
Z_step_max  = isempty(Z_during_step)  ? 0.0 : maximum(Z_during_step)

println("\n── (b) C1-FFL (persistence detector) ─────────────────────")
println("   Z_max during brief pulse  = ", round(Z_pulse_max; digits=4), "  (should be ≈ 0)")
println("   Z_max during sustained step = ", round(Z_step_max; digits=4), "  (should be > 0)")
println("   Pulse filtered? ", Z_pulse_max < 0.05)
println("   Step passes?    ", Z_step_max  > 0.1)
@assert Z_pulse_max < 0.05  "C1-FFL should filter brief pulse: Z_pulse_max=$(Z_pulse_max) (expected < 0.05)"
@assert Z_step_max  > 0.1   "C1-FFL should pass sustained step: Z_step_max=$(Z_step_max) (expected > 0.1)"

## (c) I1-FFL pulse generator
Z_max_i1 = maximum(i1.Z)
Z_end_i1 = i1.Z[end]

println("\n── (c) I1-FFL (pulse generator) ───────────────────────────")
println("   max(Z)   = ", round(Z_max_i1; digits=4))
println("   Z[end]   = ", round(Z_end_i1; digits=4))
println("   ratio    = ", round(Z_max_i1 / max(Z_end_i1, 1e-12); digits=3))
println("   Pulse signature? max(Z) > 1.5·Z[end] : ", Z_max_i1 > 1.5 * Z_end_i1)
@assert Z_max_i1 > 1.5 * Z_end_i1 "I1-FFL pulse signature not met: max(Z)=$(Z_max_i1) vs 1.5*Z[end]=$(1.5*Z_end_i1)"

## (d) Oscillator: count peaks in X (local maxima)
X_osc = osc.X
n_peaks = sum(
    X_osc[i] > X_osc[i-1] && X_osc[i] > X_osc[i+1]
    for i in 2:(length(X_osc) - 1)
)

println("\n── (d) Negative-feedback oscillator ───────────────────────")
println("   Number of peaks in X over t ∈ [0, 100] = ", n_peaks)
println("   Oscillation confirmed? (≥ 2 peaks) : ", n_peaks ≥ 2)
@assert n_peaks ≥ 2 "Oscillator should show ≥ 2 peaks in X (got $(n_peaks))"

println("\n── All signature checks passed ─────────────────────────────")

# ── 4-panel CairoMakie figure ─────────────────────────────────────────────────

fig = Figure(size=(920, 660))

# Panel 1: NAR vs simple regulation
let ax = Axis(fig[1, 1],
              title="(A)  Negative autoregulation",
              xlabel="time",
              ylabel="Y")
    lines!(ax, nar.t, nar.simple; label="simple",     linewidth=2)
    lines!(ax, nar.t, nar.nar;    label="NAR",        linewidth=2)
    hlines!(ax, [nar.Y_st[1]]; linestyle=:dash, color=:gray, label="Y*")
    axislegend(ax; position=:rb)
end

# Panel 2: C1-FFL persistence detector
let ax = Axis(fig[1, 2],
              title="(B)  C1-FFL — persistence detector",
              xlabel="time")
    lines!(ax, c1.t, c1.X; label="X (input)",  linewidth=2, linestyle=:dash)
    lines!(ax, c1.t, c1.Y; label="Y",           linewidth=1.5, color=:green)
    lines!(ax, c1.t, c1.Z; label="Z (output)",  linewidth=2)
    axislegend(ax; position=:rt)
end

# Panel 3: I1-FFL pulse generator
let ax = Axis(fig[2, 1],
              title="(C)  I1-FFL — pulse generator",
              xlabel="time")
    lines!(ax, i1.t, i1.X; label="X (input)",  linewidth=2, linestyle=:dash)
    lines!(ax, i1.t, i1.Y; label="Y",           linewidth=1.5, color=:green)
    lines!(ax, i1.t, i1.Z; label="Z (output)",  linewidth=2)
    axislegend(ax; position=:rt)
end

# Panel 4: Negative-feedback oscillator
let ax = Axis(fig[2, 2],
              title="(D)  Negative-feedback oscillator",
              xlabel="time")
    lines!(ax, osc.t, osc.X; label="X", linewidth=2)
    lines!(ax, osc.t, osc.Y; label="Y", linewidth=1.5)
    lines!(ax, osc.t, osc.Z; label="Z", linewidth=1.5)
    axislegend(ax; position=:rt)
end

save(figpath("motifs.pdf"), fig)
println("\nFigure written to: ", figpath("motifs.pdf"))
