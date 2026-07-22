# Simulated Peer Review — "Mathematical Models in Biotechnology"

**Manuscript type:** Book chapter (invited) for *Comprehensive Biotechnology, 4th Edition* — a review/methods chapter with two original Julia worked examples and two reviewed capstones.

**Note on scope.** Because this is a reference-work chapter rather than a primary research article, "novelty" is judged as the value of the synthesis and pedagogy, not as a new algorithm. The reviewers below nonetheless hold the technical content, the worked examples, and the framing claims to research-grade standards, since the chapter advances a specific thesis ("flux bounds are the gateway"; retain the growth-dilution term `μ`). Recommendations should be read as "readiness for inclusion" rather than journal accept/reject.

---

## Reviewer 1 — Moderate (metabolic engineering / systems biology)

### Summary
The chapter develops flux balance analysis (FBA) as an integrative framework in which the flux bounds, written as a product of thermodynamic, kinetic, expression, and regulatory factors, are the port through which biological information enters a constraint-based prediction. The central pedagogical device (the "gateway" factorization of a bound) is clear, well-motivated, and supported by an honest derivation of the flux constraint and two compact worked examples. Overall this is a strong, well-written chapter that will teach the material effectively; a few claims are over-sold relative to what the examples demonstrate, and these should be tempered.

### Strengths
1. The derivation in "From Mole Balances to the Flux Constraint" is genuinely valuable: deriving `S·v̂ = μx` from an open mole balance and showing `S·v = 0` as a labeled special case (Eqs. 1–5) is more honest than the usual "assume steady state" and correctly flags where the reduction fails (fed-batch, cell-free).
2. The bound factorization in "Bounds as Gateways" (Eq. 7) is an elegant organizing idea. Attaching a distinct data source to each factor (eQuilibrator, BRENDA, BioNumbers, expression profiles) makes the abstract framework concrete and actionable.
3. The two worked examples are well-chosen to be hand-checkable, and the pairing of the urea-cycle network schematic with the optimal-flux bar chart (Fig. 1) is a good expository choice.
4. The writing is dense and disciplined, and the estimation-versus-design framing of the two capstones ("Integration in Practice") is a satisfying way to close the arc.

### Weaknesses
1. **The chapter's signature claim — retain `μ` — is never exercised.** The derivation makes `S·v̂ = μx` (Eq. 4) the intellectual centerpiece, but both worked examples explicitly fall back to `S·v̂ = 0` ("Urea-Cycle Metabolism": "retains the same steady-state constraint `S·v̂ = 0`"; "A Feedback-Inhibited Linear Pathway", Eq. near the LP setup: "the same steady-state constraint `S·v̂ = 0`"). The "Outlook" concedes that estimating against `S·v̂ = μx` is future work. A reader will reasonably ask why the growth-dilution term is emphasized so heavily and then dropped in every demonstration. *Fix:* either add a short worked calculation that actually carries `μx` (even a two-metabolite toy), or soften the framing to "a term worth retaining, illustrated conceptually, and demonstrated in the reviewed cell-free capstone."
2. **The feedback example's "recovers the truth without being told the answer" is over-stated.** The text states `X_3 = T` at steady state (because `r_3` is first-order with unit rate constant), and then sets the open-gateway bound from the *measured* `X_3 ≈ 4.64`. Since `X_3` equals the throughput, the bound is being set from the answer (transformed by `X_3^{-1/2}`). The disclaimer "not copied from the throughput answer itself" is technically true but misleading. *Fix:* reframe as a self-consistency demonstration of *where* regulatory information enters (the bound), and drop the "without ever being told the answer" language.
3. **No code or data availability statement.** The chapter references JuMP/HiGHS and BSTModelKit.jl and regenerable figures, but the manuscript itself gives no repository link or availability statement. For a reference work with reproducible examples this should be explicit. *Fix:* add a short "Code and data availability" note.
4. **The examples are too small to exhibit the problem the framework is motivated by.** "The Linear Program and Its Geometry" motivates gateways by the large null space of genome-scale `S` (n ≫ m), but the urea example has `n − m = 1` and the feedback example `n − m = 1`. In both, conservation plus one binding bound essentially fixes the answer. *Fix:* acknowledge explicitly that the toy examples illustrate the mechanism, not the genome-scale underdetermination, and lean on the capstones for scale.

### Questions for Authors
1. In the urea example, is the reported optimum unique? With 14 exchange reactions and an objective that rewards only urea export, are the cofactor exchange fluxes uniquely determined, or are there alternative optima? (See Requested Analyses.)
2. The reference enzyme abundance `e° = 0.01 mmol/gDW` is applied uniformly to all five enzymes. How sensitive is the "sharp" prediction to this assumption?
3. Why is the threshold `ΔG* = −10 kJ/mol` chosen? The text says it is "fixed by the physiological concentration range," but the number appears without derivation.

### Requested Experiments / Analyses
1. Run flux variability analysis (FVA) on the urea-cycle optimum and report the min/max range for each of the 19 fluxes at the optimal urea-export value. This directly supports (or corrects) the "single, sharply determined answer" claim and is easy with the existing JuMP model.
2. Add a one-paragraph sensitivity sweep on `V°_{max,v2}` (the binding capacity) and on `e°`, showing how the predicted urea flux tracks them.

### Minor Comments
1. Fig. 1 left panel (network schematic) labels are small when compressed to half-width; consider trimming the wide dashed-border margin so the content fills the panel.
2. The abstract and Introduction ¶1 open on very similar sentences ("Predicting the intracellular flux distribution…" / "Biotechnology depends on predicting…"); consider differentiating them slightly.
3. "Integration in Practice" states cell-free metabolite tracking is "on the order of sixty metabolites" — give the actual number from the cited work if available.

### Recommendation
**Minor Revision.** The framework and exposition are strong; the required changes are tempering two over-claims, adding an availability statement, and one confirmatory FVA.

---

## Reviewer 2 — Hard (constraint-based modeling; FBA/FVA/enzyme-constrained methods)

### Summary
The chapter proposes a factorized flux bound (Eq. 7) as a unifying account of how thermodynamic, kinetic, expression, and regulatory information constrain flux. The formulation is clean and the derivation of the growth-aware constraint is a genuine strength. However, the two computational demonstrations are engineered so that the answer is dictated by a single pre-set bound, the claimed "sharp" and "recovered" predictions are presented without the variability and uncertainty analyses the chapter itself argues are necessary, and one control-function identification is mathematically inconsistent with the framework's own definitions.

### Strengths
1. The mole-balance derivation (Eqs. 1–5) is rigorous and the explicit product-rule expansion (Eq. 2 → Eq. 3) is pedagogically excellent; the identification of `S·v = 0` as a special case rather than an axiom is correct and well done.
2. The `V°_{max} = k°_cat · e°` construction (Eq. 9) and the transcript/protein steady-state balances with `(θ + μ)` denominators (Eqs. 14–16) are a correct and compact statement of enzyme-constrained and expression-aware FBA.
3. Fig. 1's pairing of network and flux is a clean way to present a small solved model.

### Weaknesses
1. **`θ(X_3) = X_3^{-1/2}` violates the chapter's own definition `θ_j ∈ [0,1]`.** "Bounds as Gateways" defines the control function as "a dimensionless number in [0,1] reporting the fraction of a reaction's capacity switched on." In the feedback example the control function is set to the power law `X_3^{-1/2}` (Eq. 21 mapping), which is unbounded: for `X_3 < 1` it exceeds 1, so the bound `V°_{max,r0}·θ` would exceed the *uninhibited* capacity, contradicting the interpretation of `θ` as a fraction switched on. The claim that the Boltzmann/partition-function control function (a saturating, Hill-type object bounded in [0,1]) "collapses to the concrete power law `X_3^{-1/2}`" (line ~70) is asserted, not derived, and is in fact inconsistent with the partition-function form, which cannot produce an unbounded power law. *Fix:* either derive an explicit, bounded `θ(X_3)` from the partition function for this effector and show it reproduces the throttling, or state plainly that `X_3^{-1/2}` is borrowed from the S-system rate law for illustration and is not a partition-function control function.
2. **Neither example runs FVA, despite the chapter arguing FVA is essential.** "The Linear Program and Its Geometry" ends by explaining that a single optimal flux vector "cannot report" the ranges that FVA reveals. The urea example then reports exactly such a single vector and calls it "sharply determined." Uniqueness is not shown. In a network with 14 exchange reactions and an objective that fixes only urea export, alternative optima in the cofactor/byproduct exchanges are plausible. *Fix:* report FVA ranges at the optimum; if fluxes are non-unique, say so.
3. **The urea "prediction" reduces to `min` over capacities and does not exercise the framework.** In the closed cycle at steady state the backbone reactions carry equal flux, so maximal urea export equals `min(V°_{max,v1..v4}) = V°_{max,v2} = 0.0328`. The thermodynamic gateway (irreversibility of v3, v4, v5) never binds because the optimum is all-forward, and three of the four gateways are set to unity. The elaborate LP therefore adds nothing to a one-line `min`. *Fix:* choose an example where the branch point (v5) is actually contested at the optimum (e.g., a competing objective or a tighter arginine supply), so the gateways demonstrably change the result.
4. **Uniform `e°` across all enzymes is a strong, unquantified assumption.** `V°_{max}` is proportional to enzyme abundance, which varies over orders of magnitude; setting a single `e° = 0.01 mmol/gDW` for all five enzymes means the "sharp" prediction is really a statement about assumed equal abundance. This is precisely the regime where enzyme-constrained methods (see Weakness 6) use proteomics. *Fix:* add a sensitivity analysis and cite the enzyme-constrained FBA literature.
5. **No uncertainty propagation.** The entire chapter is premised on data scarcity, yet every prediction is a point estimate. `k_cat` values from BRENDA and `ΔG°` from group contribution both carry large, well-documented uncertainties. Propagating even nominal ±1 order-of-magnitude `k_cat` uncertainty through the urea LP would show whether `v2` remains the bottleneck. *Fix:* add a short Monte Carlo or interval analysis on the binding parameters.
6. **Missing positioning against enzyme-constrained and thermodynamic FBA.** The kinetic and thermodynamic gateways are, respectively, enzyme-constrained FBA (e.g., GECKO, Sánchez et al., *Mol Syst Biol* 2017) and thermodynamics-based flux analysis (Henry, Broadbelt & Hatzimanikatis, *Biophys J* 2007). Neither is cited. *Fix:* cite and briefly relate the gateways to these established methods.

### Questions for Authors
1. For the feedback example, what is the explicit partition-function expression that yields `θ(X_3) = X_3^{-1/2}`, and over what range of `X_3` is it bounded in [0,1]? If none exists, on what basis is `X_3^{-1/2}` called a control function?
2. Is the urea optimum unique (all 19 fluxes), or only in the objective value? What does FVA give?
3. Does the thermodynamic gateway change the urea result at all? If every optimal flux is forward, what does assigning `δ_j = 0` to v3/v4/v5 accomplish here?
4. What organism, pH, ionic strength, and temperature were used for the eQuilibrator `ΔG°` and BRENDA `k_cat` lookups? These strongly affect the numbers.

### Requested Experiments / Analyses
1. FVA at the urea optimum, reported as a table or an added error-bar overlay on Fig. 1 (right).
2. A parameter-uncertainty analysis (Monte Carlo over `k_cat`, `e°`, and `ΔG°` within literature ranges) reporting the resulting distribution of the predicted urea flux and the probability that `v2` remains the bottleneck.
3. A revised or additional urea-cycle scenario in which the v5 (NOS) branch is genuinely contested at the optimum, so that opening/closing the thermodynamic and kinetic gateways changes the predicted split.

### Minor Comments
1. Eq. 7 writes the bound symmetrically with `δ_j` on both sides; clarify that `θ_j`, `f_j`, `(e/e°)` are assumed identical for the forward and reverse capacity, which is not obvious for reversible reactions.
2. The `−1000 ≤ v̂ ≤ 1000` default exchange bounds should be stated in the same units and justified as "effectively unconstrained relative to the 0.0328 scale."
3. Fig. 1 right panel y-axis is labeled "flux" with no units on the axis; put units on the axis, not only in the caption.

### Recommendation
**Major Revision.** The formulation is sound but the demonstrations do not yet support the claims; FVA, an uncertainty analysis, resolution of the `θ ∈ [0,1]` inconsistency, and positioning against enzyme-constrained/thermodynamic FBA are required.

---

## Reviewer 3 — Very Hard (author of competing integrated constraint-based methods)

### Summary
The chapter reframes a set of well-established extensions of FBA (enzyme constraints, thermodynamic constraints, metabolism–expression coupling, regulatory constraints) as four "gateways" on a factorized bound. The reframing is pedagogically tidy but is presented as if it were a new integrative insight, while the foundational methods that already implement each gateway go uncited. The two original examples are constructed so that the answer is inserted as a bound and then read back out, so neither validates the framework's predictive value. I recommend the framework be positioned honestly against prior art and that at least one example be made non-trivial.

### Strengths
1. The honest derivation of `S·v̂ = μx` and the explicit statement of the assumptions behind `S·v = 0` (Eqs. 1–5) are the strongest and most original part of the chapter and are worth keeping front and center.
2. The single-sentence unification at the end of "Integration in Practice" (estimation and design as the same operation on the bounds) is a genuinely useful framing.
3. The chapter is unusually careful about what is assumed away, which is rare in FBA pedagogy.

### Weaknesses
1. **The "gateways" are established methods presented without attribution.** Each gateway has a canonical, named predecessor that is not cited: the kinetic gateway (`V°_{max} = k_cat·e°`) is enzyme-constrained FBA / GECKO (Sánchez et al., *Mol Syst Biol* 2017); the expression gateway (transcript/protein balances feeding `(e/e°)`) is metabolism-and-expression (ME) modeling (O'Brien et al., *Mol Syst Biol* 2013; Lloyd et al., COBRAme, *PLoS Comput Biol* 2018); the thermodynamic gateway (`δ_j` from `ΔG`) is thermodynamics-based flux analysis (Henry et al., *Biophys J* 2007); the regulatory gateway (`θ_j`) is regulatory FBA (Covert, Schilling & Palsson, *J Theor Biol* 2001) and PROM (Chandrasekaran & Price, *PNAS* 2010). As written, a reader would think the factorization is new. *Fix:* add a paragraph mapping each gateway to its established method with citations, and state clearly that the contribution is the unifying pedagogy, not the individual constraints.
2. **The feedback example is circular by construction.** The text establishes `X_3 = T` (steady state, unit-rate export), then sets the open-gateway bound to `10·X_3^{-1/2} = 10·T^{-1/2}`, which at the self-consistent `T = 4.64` equals `4.64`; the LP then "returns 4.64." The bound is the answer. Because "the optimal throughput is read directly off the bound" (the authors' own words), the LP performs no inference. This cannot be described as the flux-balance model "recovering the truth"; it is an identity. *Fix:* either predict `X_3` (and hence the bound) from an independent measurement that is not equal to `T`, or present the example honestly as an illustration of *where* regulation enters, not as a validation.
3. **`θ = X_3^{-1/2}` is not a partition-function control function.** The claim that the Boltzmann accounting "collapses to" this power law is false in general: partition-function occupancies are bounded in [0,1] and saturate; a power law that diverges as `X_3 → 0` is not among them. The example simply reuses the S-system's own kinetic order and relabels it `θ`. *Fix:* derive a genuine bounded `θ(X_3)` or retract the identification.
4. **The central "keep μ" thesis is asserted but never demonstrated in an original calculation.** Every original example uses `S·v = 0`. The only setting where the retained term is exercised is the reviewed cell-free capstone, which is prior published work by others. The chapter therefore does not itself demonstrate its headline contribution. *Fix:* add an original calculation (fed-batch or cell-free toy) in which carrying `μx` changes a predicted flux relative to `S·v = 0`.
5. **"Estimation and design are the same operation on the bounds" conflates two different problems.** The estimation capstone tightens continuous bounds from data (a parameter-fitting/LP problem); the design capstone (Wayman) sets bounds to zero via gene knockouts, which is a discrete combinatorial search (typically MILP / bilevel, e.g., OptKnock). Calling these "the same operation" glosses over the combinatorial explosion that makes design hard. *Fix:* qualify the claim, or cite the design-optimization literature (e.g., OptKnock, Burgard et al., *Biotechnol Bioeng* 2003).
6. **No genome-scale demonstration of the authors' own.** The framework is motivated by genome-scale underdetermination, but both original examples have `n − m = 1`. The genome-scale content is entirely in reviewed third-party capstones. *Fix:* acknowledge this limitation prominently, or add a genome-scale gateway demonstration.

### Questions for Authors
1. What does the "gateway" factorization provide that GECKO (kinetic), ME-models (expression), TFA (thermodynamic), and rFBA/PROM (regulatory) do not, other than a common notation?
2. In the feedback example, if `X_3 ≠ T` (e.g., add a nonzero degradation of `X_3` so that the export flux and the pool decouple), does the open gateway still "recover" the throughput, or does the identity break?
3. Since the LP reads the throughput directly off the `r_0` bound, what role does linear programming play in this example at all?
4. Can the authors exhibit any prediction in the chapter that could have come out *wrong* — i.e., where the bound was set from data not equal to the quantity being predicted?

### Requested Experiments / Analyses
1. Break the `X_3 = T` degeneracy in the feedback example (add `X_3` degradation or a second consumer) and repeat, so that opening the regulatory gateway is a genuine prediction rather than an identity. Report whether recovery survives.
2. Add a "related methods" comparison table mapping each gateway to its established method and stating the incremental contribution.
3. Provide one calculation in which retaining `μx` (`S·v̂ = μx`) yields a measurably different flux than `S·v = 0`, to substantiate the chapter's thesis.

### Minor Comments
1. The abstract's "the same operation on the bounds, run in opposite directions" inherits the estimation/design conflation of Weakness 5; soften.
2. "Bounds as Gateways" presents `f_j` (saturation) and `θ_j` (regulation) as independent factors, but substrate saturation and allosteric regulation are often coupled through the same enzyme state; a sentence acknowledging this would prevent over-reading the factorization as physically independent.
3. The uniform statement that constraint-based models "impose conservation of mass at steady state" (Introduction) should note that exchange/biomass pseudo-reactions are not mass-balanced in the usual sense.

### Recommendation
**Major Revision.** As a synthesis the chapter is useful, but it must (i) cite and position against the established gateway-by-gateway methods, (ii) fix or reframe the circular feedback example and the unbounded `θ`, and (iii) either demonstrate the `μ`-retaining thesis in an original calculation or downgrade it from headline claim to motivating remark.

---

## Summary of Actionable Items (consolidated, deduplicated, prioritized)

### Priority 1 — Must fix (claims not supported by content)
1. **Resolve the `θ ∈ [0,1]` vs `θ = X_3^{-1/2}` inconsistency** (R2-W1, R3-W3): either derive a bounded partition-function `θ(X_3)` or explicitly state the power law is borrowed from the S-system rate law and is not a partition-function control function.
2. **Reframe or repair the circular feedback example** (R1-W2, R3-W2/Q2, R2 implicitly): `X_3 = T` makes the "recovery without being told the answer" an identity. Either break the `X_3 = T` degeneracy (add `X_3` degradation / second consumer) so the open gateway makes a genuine prediction, or drop the validation framing and present it as "where regulation enters."
3. **Demonstrate the `μ`-retaining thesis or downgrade it** (R1-W1, R2 summary, R3-W4): add an original calculation in which carrying `μx` changes a flux relative to `S·v = 0`, or reframe "keep `μ`" as a motivating remark illustrated only in the reviewed capstone.
4. **Cite and position against established gateway methods** (R2-W6, R3-W1): GECKO/enzyme-constrained FBA (Sánchez 2017), ME-models (O'Brien 2013; Lloyd/COBRAme 2018), thermodynamics-based flux analysis (Henry 2007), rFBA (Covert 2001) and PROM (Chandrasekaran & Price 2010); add a mapping paragraph or table.

### Priority 2 — Strongly requested analyses
5. **Run FVA at the urea optimum** and report per-flux ranges (R1-A1, R2-W2/A1); confirm or retract "single, sharply determined answer."
6. **Uncertainty / sensitivity analysis** over `k_cat`, `e°`, `ΔG°` for the urea example (R1-A2/Q2, R2-W4/W5/A2); report the distribution of the predicted urea flux and stability of the `v2` bottleneck.
7. **Make at least one example non-trivial** (R2-W3/A3, R3-W6): a urea scenario where the v5 branch is contested at the optimum, so the gateways demonstrably change the result rather than the answer being `min(capacity)`.

### Priority 3 — Presentation and scope
8. **Add a Code and Data availability statement** with the repository and solver versions (R1-W3).
9. **Qualify "estimation and design are the same operation"** (R3-W5/minor, abstract): note the design side is discrete/combinatorial; optionally cite OptKnock (Burgard 2003).
10. **Acknowledge the toy-scale limitation** explicitly (R1-W4, R3-W6): the original examples have `n − m = 1` and do not exhibit the genome-scale underdetermination the framework is motivated by.
11. **Figure and unit fixes** (R1-m1/m3, R2-m3): enlarge the Fig. 1 network panel (trim its dashed-border margin), add flux units to the Fig. 1 right-panel axis, and give the exact cell-free metabolite count.
12. **Report database query conditions** for eQuilibrator/BRENDA (organism, pH, ionic strength, temperature) (R2-Q4).

### Overall disposition
Two reviewers recommend **Major Revision** and one **Minor Revision**. The consensus blocking issues are the unsupported/over-stated claims (circular feedback recovery, unbounded `θ`, undemonstrated `μ` thesis) and the missing positioning against prior integrated-modeling methods. None of the issues is fatal; all are addressable within the existing examples plus one added analysis, and the derivation and framing are strong enough to make the chapter a valuable contribution once the claims are brought into line with what the examples show.
