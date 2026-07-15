include(joinpath(@__DIR__, "..", "Include.jl"))
include(joinpath(@__DIR__, "dual_feedback.jl"))

# --------------------------------------------------------------------------- #
# Example 2 driver: integrate the dual-feedback truth, approximate it through the
# committed-step bound in four escalating cases, and emit the CSV + figure.
# --------------------------------------------------------------------------- #
t  = feedback_truth()
gw = gateway_factors(t)

vtruth  = truth_metabolic_fluxes(t)
v_naive = feedback_fba(gw; expression=false, activity=false)
v_expr  = feedback_fba(gw; expression=true,  activity=false)
v_act   = feedback_fba(gw; expression=false, activity=true)
v_both  = feedback_fba(gw; expression=true,  activity=true)

df = DataFrame(reaction=METAB_REACTIONS, truth=vtruth, naive=v_naive,
               expression_only=v_expr, activity_only=v_act, both_open=v_both)
CSV.write(datapath("feedback_fba.csv"), df)

# ---- prose numbers -------------------------------------------------------- #
sp = t.seq
X3 = t.Xss[findfirst(==("X3"), t.species)]
T  = vtruth[findfirst(==("r3"), METAB_REACTIONS)]
ledger = [gw.Vmax0, gw.Vmax0*gw.e_e0, gw.Vmax0*gw.θ, gw.Vmax0*gw.e_e0*gw.θ]
println("fixed point: X3*=", round(X3,digits=3), "  T*=", round(T,digits=3),
        "  theta*=", round(gw.θ,digits=3), "  (e/e0)*=", round(gw.e_e0,digits=3))
println("overshoot factor Vmax0/T* = ", round(gw.Vmax0/T,digits=3))
println("ledger [Vmax0, x(e/e0), x theta, both] = ", round.(ledger,digits=3))
println("mu fraction of clearance: transcript=", round(sp.mu/(sp.theta_m+sp.mu),digits=3),
        "  protein=", round(sp.mu/(sp.theta_p+sp.mu),digits=3))

# ---- figure: [wiring schematic | capacity ledger] ------------------------- #
fig = Figure(size=(940, 380))

# panel (a): dual-feedback wiring schematic
axS = Axis(fig[1,1], title="Dual end-product control of the committed enzyme")
hidedecorations!(axS); hidespines!(axS); xlims!(axS, 0, 10); ylims!(axS, 0, 6)
nodes = Dict("X1"=>(2.0,1.2), "X2"=>(4.0,1.2), "X3"=>(6.0,1.2),
             "gene"=>(6.0,4.6), "m"=>(4.0,4.6), "E0"=>(2.0,4.6))
for (lbl,(x,y)) in nodes
    scatter!(axS, [x],[y]; markersize=34, color=(:steelblue,0.20), strokecolor=:steelblue, strokewidth=1.5)
    text!(axS, x, y; text=lbl, align=(:center,:center), fontsize=13)
end
# metabolic chain (solid) and expression cascade (solid, rightward gene->m->E0)
arrows!(axS, [0.5,2.3,4.3,6.3], [1.2,1.2,1.2,1.2], [1.1,1.4,1.4,1.1], [0,0,0,0]; color=:black)  # r0 in, X1->X2, X2->X3, r3 out
arrows!(axS, [5.6,3.6],[4.6,4.6],[-1.2,-1.2],[0,0]; color=:seagreen)                             # gene->m, m->E0
arrows!(axS, [2.0],[4.1],[0.0],[-2.2]; color=:seagreen)                                          # E0 catalyzes r0 (down)
text!(axS, 1.0,1.55; text="r0", fontsize=11); text!(axS, 7.1,1.55; text="r3", fontsize=11)
text!(axS, 2.35,3.0; text="catalysis", fontsize=10, color=:seagreen, align=(:left,:center))
# two repression arrows from X3 (dashed, crimson): slow expression + fast activity
lines!(axS, [6.0,6.0],[1.6,4.2]; color=:crimson, linestyle=:dash)                               # X3 -| gene (expression)
lines!(axS, [5.7,2.3],[1.5,4.2]; color=:crimson, linestyle=:dash)                               # X3 -| E0   (activity)
text!(axS, 6.15,3.0; text="(e/e°)\nexpression, slow", color=:crimson, fontsize=10, align=(:left,:center))
text!(axS, 3.7,2.4; text="θ\nactivity, fast", color=:crimson, fontsize=10, align=(:left,:center))

# panel (b): capacity ledger
axL = Axis(fig[1,2], xticks=(1:4, ["V°max","×(e/e°)","×θ","both"]),
           ylabel="committed-step capacity (AU)",
           title="Capacity ledger: the two gateways multiply")
barplot!(axL, 1:4, ledger; color=[:gray70,:steelblue,:orange,:seagreen])
hlines!(axL, [gw.Vmax0]; color=:gray, linestyle=:dash, label="un-inhibited capacity")
hlines!(axL, [T];        color=:crimson, linestyle=:dot, label="BST truth T*")
axislegend(axL; position=:rt, framevisible=false)

save(figpath("feedback_gateway.pdf"), fig)
println("run_feedback OK: wrote feedback_fba.csv and feedback_gateway.pdf")

# ---- structural robustness: sweep the FABRICATED gateway shape (K,n) ------- #
# The activity gateway's functional form is a modeling choice, unrelated to the
# truth's kinetics. Sweeping it shows the recovery brackets the truth rather than
# hitting it, and that only the dual-gateway model reaches the truth at all.
sw = theta_shape_sweep(t, gw)
CSV.write(datapath("feedback_sweep.csv"),
          DataFrame(K = repeat(sw.Ks, outer=length(sw.ns)),
                    n = repeat(sw.ns, inner=length(sw.Ks)),
                    T_both = vec(sw.Tboth), T_act = vec(sw.Tact)))

Tlo, Thi = extrema(sw.Tboth)
println("shape sweep: T_both in [", round(Tlo,digits=2), ", ", round(Thi,digits=2),
        "] ; brackets truth ", round(T,digits=2), " = ", Tlo < T < Thi)
println("single-gateway reach across window: expr=", round(sw.Texpr,digits=2),
        " (shape-free)  act in [", round(minimum(sw.Tact),digits=2), ", ",
        round(maximum(sw.Tact),digits=2), "] ; both stay above truth = ",
        minimum(sw.Tact) > T && sw.Texpr > T)

# figure: both-open throughput over the (K,n) grid, color diverging about the
# truth so blue/red on either side make the bracketing visible at a glance.
fig2 = Figure(size=(560, 430))
ax2  = Axis(fig2[1,1], xlabel="half-saturation K", ylabel="Hill coefficient n",
            title="A fabricated gateway shape still brackets the truth")
Δ  = maximum(abs.(sw.Tboth .- T))
hm = heatmap!(ax2, sw.Ks, sw.ns, sw.Tboth; colormap=:RdBu, colorrange=(T-Δ, T+Δ))
contour!(ax2, sw.Ks, sw.ns, sw.Tboth; levels=[T], color=:black, linewidth=2)
scatter!(ax2, [K_THETA], [N_THETA]; color=:white, strokecolor=:black,
         strokewidth=1.5, markersize=13)
text!(ax2, K_THETA, N_THETA; text="  (K,n)=(5,2)", align=(:left,:center), fontsize=11)
Colorbar(fig2[1,2], hm, label="throughput T (both gateways open)")
save(figpath("feedback_robustness.pdf"), fig2)
println("run_feedback OK (robustness): wrote feedback_sweep.csv and feedback_robustness.pdf")
