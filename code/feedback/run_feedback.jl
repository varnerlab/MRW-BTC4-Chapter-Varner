include(joinpath(@__DIR__, "..", "Include.jl"))
include(joinpath(@__DIR__, "dual_feedback.jl"))

# --------------------------------------------------------------------------- #
# Example 2 driver: integrate the dual-feedback reference, compare it with four
# committed-step bound configurations, and emit the CSV and figures.
# --------------------------------------------------------------------------- #
reference = feedback_reference()
controls = control_factors(reference)

v_reference = reference_metabolic_fluxes(reference)
v_naive = feedback_fba(controls; expression=false, activity=false)
v_expr  = feedback_fba(controls; expression=true,  activity=false)
v_act   = feedback_fba(controls; expression=false, activity=true)
v_both  = feedback_fba(controls; expression=true,  activity=true)

df = DataFrame(reaction=METAB_REACTIONS, reference=v_reference, naive=v_naive,
               expression_only=v_expr, activity_only=v_act,
               both_controls=v_both)
CSV.write(datapath("feedback_fba.csv"), df)

# ---- prose numbers -------------------------------------------------------- #
sp = reference.seq
X3 = reference.Xss[findfirst(==("X3"), reference.species)]
T  = v_reference[findfirst(==("r3"), METAB_REACTIONS)]
capacities = [controls.Vmax0, controls.Vmax0*controls.e_e0,
              controls.Vmax0*controls.θ,
              controls.Vmax0*controls.e_e0*controls.θ]
println("fixed point: X3*=", round(X3,digits=3), "  T*=", round(T,digits=3),
        "  theta*=", round(controls.θ,digits=3),
        "  (e/e0)*=", round(controls.e_e0,digits=3))
println("overshoot factor Vmax0/T* = ", round(controls.Vmax0/T,digits=3))
println("capacities [Vmax0, x(e/e0), x theta, both] = ",
        round.(capacities,digits=3))
println("illustrative elongation times: transcription=",
        round(sp.transcription_time_s,digits=2), " s, translation=",
        round(sp.translation_time_s,digits=2), " s (not used in calibration)")
println("mu fraction of clearance: transcript=", round(sp.mu/(sp.theta_m+sp.mu),digits=3),
        "  protein=", round(sp.mu/(sp.theta_p+sp.mu),digits=3))

# ---- figure: effects of the two controls (wiring is a TikZ schematic) ------ #
# chapter/figures/feedback_network.tex; the two loop colors below match it) --- #
fig = Figure(size=(520, 400))
axL = Axis(fig[1,1], xticks=(1:4, ["V°max","×(e/e°)","×θ","both"]),
           ylabel="committed-step capacity (AU)",
           title="Expression and activity factors multiply")
barplot!(axL, 1:4, capacities; color=[:gray70,:steelblue,:orange,:seagreen])
hlines!(axL, [controls.Vmax0]; color=:gray, linestyle=:dash, label="un-inhibited capacity")
hlines!(axL, [T];        color=:crimson, linestyle=:dot, label="BST reference T*")
axislegend(axL; position=:rt, framevisible=false)

save(figpath("feedback_ledger.pdf"), fig)
println("run_feedback OK: wrote feedback_fba.csv and feedback_ledger.pdf")

# ---- structural robustness: sweep the ASSUMED activity-factor shape (K,n) -- #
# The activity factor's functional form is a modeling choice, unrelated to the
# reference kinetics. The sweep tests whether the dual-control range includes
# the reference value and whether either single-control range does.
sw = theta_shape_sweep(reference, controls)
CSV.write(datapath("feedback_sweep.csv"),
          DataFrame(K = repeat(sw.Ks, outer=length(sw.ns)),
                    n = repeat(sw.ns, inner=length(sw.Ks)),
                    T_both = vec(sw.Tboth), T_act = vec(sw.Tact)))

Tlo, Thi = extrema(sw.Tboth)
println("shape sweep: T_both in [", round(Tlo,digits=2), ", ", round(Thi,digits=2),
        "] ; includes reference ", round(T,digits=2), " = ", Tlo < T < Thi)
println("single-factor reach across window: expr=", round(sw.Texpr,digits=2),
        " (shape-free)  act in [", round(minimum(sw.Tact),digits=2), ", ",
        round(maximum(sw.Tact),digits=2), "] ; both stay above reference = ",
        minimum(sw.Tact) > T && sw.Texpr > T)

# figure: dual-control throughput over the (K,n) grid, color diverging about the
# reference so blue and red show values on either side.
fig2 = Figure(size=(560, 430))
ax2  = Axis(fig2[1,1], xlabel="half-saturation K", ylabel="Hill coefficient n",
            title="An assumed activity factor brackets the BST reference")
Δ  = maximum(abs.(sw.Tboth .- T))
hm = heatmap!(ax2, sw.Ks, sw.ns, sw.Tboth; colormap=:RdBu, colorrange=(T-Δ, T+Δ))
contour!(ax2, sw.Ks, sw.ns, sw.Tboth; levels=[T], color=:black, linewidth=2)
scatter!(ax2, [K_THETA], [N_THETA]; color=:white, strokecolor=:black,
         strokewidth=1.5, markersize=13)
text!(ax2, K_THETA, N_THETA; text="  (K,n)=(5,2)", align=(:left,:center), fontsize=11)
Colorbar(fig2[1,2], hm, label="throughput T (both factors included)")
save(figpath("feedback_robustness.pdf"), fig2)
println("run_feedback OK (robustness): wrote feedback_sweep.csv and feedback_robustness.pdf")
