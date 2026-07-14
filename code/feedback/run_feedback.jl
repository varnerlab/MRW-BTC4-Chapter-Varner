include(joinpath(@__DIR__, "..", "Include.jl"))
include(joinpath(@__DIR__, "dual_feedback.jl"))

# --------------------------------------------------------------------------- #
# Example 2 driver: integrate the dual-feedback truth, recover it through the
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
text!(axS, 6.15,3.0; text="(e/e0)\nexpression, slow", color=:crimson, fontsize=10, align=(:left,:center))
text!(axS, 3.7,2.4; text="theta\nactivity, fast", color=:crimson, fontsize=10, align=(:left,:center))

# panel (b): capacity ledger
axL = Axis(fig[1,2], xticks=(1:4, ["Vmax0","x(e/e0)","x theta","both"]),
           ylabel="committed-step capacity (AU)",
           title="Capacity ledger: the two gateways multiply")
barplot!(axL, 1:4, ledger; color=[:gray70,:steelblue,:orange,:seagreen])
hlines!(axL, [gw.Vmax0]; color=:gray, linestyle=:dash, label="un-inhibited capacity")
hlines!(axL, [T];        color=:crimson, linestyle=:dot, label="BST truth T*")
axislegend(axL; position=:rt, framevisible=false)

save(figpath("feedback_gateway.pdf"), fig)
println("run_feedback OK: wrote feedback_fba.csv and feedback_gateway.pdf")
