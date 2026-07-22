# Phase 2: Claims Alignment

**Date:** 2026-07-14
**Source:** `chapter_audit_results.md` + `peer-review.md`, continuing the Phase 1 decomposition.
**Sections touched:** `chapter/sections/introduction.tex`, `linearprogram.tex`, `gateways.tex`,
  `example_urea.tex`, `capstones.tex`
**Code touched:** `code/fba/urea_cycle_uq.jl`
**Depends on:** Phase 1 (merged) — all numbers below assume the corrected units and the
  bounded feedback gateway are already in place.

## Scope

Phase 1 fixed code/math correctness. Phase 2 aligns chapter *claims* with what the corrected
code and cited literature actually support — no further code changes except one (§4, common
random numbers), which is small and independently verified below.

**Already satisfied by Phase 1, dropped from this phase's scope:** "reframe the feedback
example as a genuine (if toy-scale) prediction" — Task 11 of Phase 1 already did this in full
(bounded θ_FBA, softened "recovering...to numerical tolerance" claim).

**Explicitly deferred (unchanged from the Phase 1 decision):** database provenance (BRENDA
query conditions, eQuilibrator organism/pH/ionic-strength/temperature) stays a placeholder —
no new literature sourcing for those specific parameters in this phase either.

---

## §1. Thermodynamic gate: label as a heuristic

**Files:** `chapter/sections/gateways.tex`, `chapter/sections/example_urea.tex`

**Decision:** label the −10 kJ/mol cutoff explicitly as a heuristic directionality rule, not a
result of a shown concentration-range calculation (no new data sourcing, consistent with the
provenance-placeholder decision).

`gateways.tex`, in the thermodynamic-gateway paragraph, the sentence "a reaction whose driving
force cannot plausibly change sign anywhere in that range is declared irreversible" should be
followed by an explicit heuristic label, e.g.: "...is declared irreversible — a heuristic
directionality rule that substitutes for propagating a concentration range through $\Delta G =
\Delta G^{\circ} + RT\ln Q$, not a substitute for it."

`example_urea.tex`, the sentence introducing $\Delta G^{*}_j=-10$ kJ/mol should similarly be
flagged as a heuristic threshold rather than one derived from a shown range calculation.

---

## §2. UQ Algorithm 1: fix the reused-$Z$ pseudocode and show substrate/Km draws

**File:** `chapter/sections/example_urea.tex` (Algorithm 1, `alg:uq`)

**Problem:** the algorithm reuses a single symbol $Z$ for the shared $e^{\circ}$ draw, each
per-reaction $k_{\mathrm{cat},j}$ draw, and each per-reaction $\Delta G_j$ draw, implying
perfectly correlated draws. The actual code (`urea_cycle_uq.jl`) draws:
- `e0`: one shared lognormal draw (`lognrand(rng, E0, SIGMA_LN)`) — correctly shown as shared.
- `kcat`: 5 independent lognormal draws (`KCAT0 .* exp.(SIGMA_LN .* randn(rng, 5))`).
- `dG`: 5 independent normal draws (`DG0 .+ DG_SIGMA .* randn(rng, 5)`).
- `f_j` (where substrate data exists): independent lognormal draws for *both* concentration
  and $K_M$ per reaction (`sample_f`), not shown at all in the current pseudocode as distinct
  quantities.

**Fix:** rewrite Algorithm 1 to use distinct symbols per draw ($Z_e$ shared;
$Z_{k,j}$/$Z_{\Delta G,j}$ independent per reaction; explicit $[S_j]$/$K_{M,j}$ draws shown
where substrate data exists, `else` branch for $f_j=1$). This directly fixes both audit
sub-items ("reuses a single symbol Z" and "does not show distinct draws for substrate
concentration and Km").

---

## §3. UQ prose: soften three overstated empirical claims

**File:** `chapter/sections/example_urea.tex`

Verified against the current (Phase-1-regenerated) `code/data/urea_fba_uq.csv`:
```
v1: mean=156.08965310772047  q50=105.06709989875603
v2: mean=156.08965310772047  q50=105.06709989875603   (identical to v1)
v3: mean=156.04172201526268  q50=105.05381342754373
v4: mean=156.04172201526268  q50=105.05381342754373   (identical to v3, ~0.03% below v1/v2)
v5: mean=0.02396554622889309  std=1.663797720005487    (nonzero — active in rare draws)
```

- **"the nitric oxide synthase branch stays silent across the ensemble"** → false; v5 has a
  small positive mean and nonzero spread. Replace with language acknowledging it is active
  (contributes flux) only in a small fraction of draws, remaining negligible in typical draws.
- **"the sampled bound on v2 alone carries the uncertainty"** → soften "alone" to "dominates",
  since v1/v2 and v3/v4 are two distinct (if very close) pairs, not one identical value across
  the whole backbone.
- **Figure 1 caption, "every nonzero flux equals the single throughput up to sign, so the
  parametric-bootstrap uncertainty...is one quantity common to them all"** → soften to "nearly
  equal" / "essentially one quantity", consistent with the v1/v2 vs. v3/v4 split above.
- **Terminology:** add one clarifying sentence near Algorithm 1's introduction that the
  parameter distributions (lognormal spread ~2×, normal $\Delta G$ spread) are modeling
  assumptions, not distributions fitted to replicate data — addressing the audit's fallback
  option for keeping "parametric bootstrap" as the algorithm's name without renaming it
  everywhere.

---

## §4. UQ code: give Config A/B common random numbers

**File:** `code/fba/urea_cycle_uq.jl`

**Problem:** Config B currently draws from `MersenneTwister(SEED + 1)` — a genuinely different
seed from Config A's `MersenneTwister(SEED)` — so kcat/e0/dG/measured-$f_j$ are independently
resampled between the two configs, not shared. The audit's suggestion ("common random numbers
for the shared parameters, varying only the imputed f2 factor") is not currently true.

**Fix:** change `let rng = MersenneTwister(SEED + 1)` to `let rng = MersenneTwister(SEED)` in
the Config B block. Because Config A and B draw kcat/e0/dG/measured-$f_j$ in the identical
order via the identical code path, and Config A never draws $f_2$ (its `sample_f` call uses
`impute_f2=false`) while Config B does (`impute_f2=true`), this single seed change makes every
draw *up to* the $f_2$ draw bit-identical between A and B — a true paired comparison that
isolates the $f_2$ effect. No other code change needed.

**Verified during planning** (scratch re-run with this one-line change):
```
configB_f2 median=0.4976497320733263        (was 0.5021191499901032)
configB_urea_export median=50.82740723016741 ci=[4.19,400.59]   (was median=51.27, ci=[4.35,397.81])
```
Both round to the same prose values already in the chapter ("51", "0.50") — **no prose numbers
need to change**, only `code/data/urea_uq_sensitivity.csv` and `code/figs/urea_saturation.pdf`
need regenerating.

---

## §5. Write up the FVA result in the urea section

**File:** `chapter/sections/example_urea.tex`

Phase 1 (Task 3) added `fva(m; tol=1e-6)` and a test confirming the nominal optimum is unique
(max range across all 19 reactions < 1e-4). This was never written into the chapter prose. Add
one or two sentences after the nominal-solution paragraph reporting that an independent FVA
confirmed every reaction's minimum and maximum flux coincide at the optimum, substantiating
"single, sharply determined answer" with an actual check rather than an assertion — directly
answering Reviewer 1's/Reviewer 2's requested analysis.

---

## §6. `linearprogram.tex`: two overstated claims

- **"the reaction count $n$ exceeds the metabolite count $m$ by a wide margin, often an order
  of magnitude or more"** (genome-scale reconstructions in general aren't reliably 10×). The
  file already states the *correct*, rank-based claim later in the same section ("the rank of
  $\mathbf S$ is smaller than the reaction count, $\mathrm{rank}(\mathbf S)\ll n$..."). Fix: trim
  the first (overstated, redundant) claim rather than duplicating the rank statement — e.g.
  soften to "the reaction count $n$ typically exceeds the metabolite count $m$" and let the
  later, already-correct rank statement carry the precise quantitative claim.
- **"A reaction whose minimum and maximum coincide is pinned by the network topology and
  bounds regardless of objective choice"** — FVA pins a reaction only within the *specific*
  objective (or near-optimal set) used for that FVA calculation; a different objective can
  produce a different feasible face and a different range for the same reaction. Fix: qualify
  "regardless of objective choice" to something like "at the optimal objective value used for
  this calculation; a different objective can produce a different feasible face and a
  different range."

---

## §7. `introduction.tex`: two overly-exclusive claims about stoichiometry

- **"A genome sequence and the stoichiometry of its annotated reactions are knowable in
  full"** — genome-scale reconstructions carry real uncertainty in reaction membership,
  directionality, compartmentalization, cofactors, GPR rules, exchange definitions, and biomass
  composition. Soften "knowable in full" to acknowledge this.
- **"The bounds, not the stoichiometry, are where biological information enters"** — replace
  with the audit's suggested, more defensible framing: "Flux bounds are a major interface
  through which condition-specific thermodynamic, kinetic, expression, and regulatory
  information enters a constraint-based model" — preserves the gateway framing without making
  it exclusive.

---

## §8. Qualify "estimation and design are the same operation"

**Files:** `chapter/sections/capstones.tex`, `chapter/sections/introduction.tex`

Both files currently make this claim unqualified (introduction.tex's closing sentence
paraphrases capstones.tex's). Once a deletion set is specified, both problems solve an LP with
altered bounds — that part is true. But finding the deletion set is a discrete, often
combinatorial search (unlike continuous bound-fitting for estimation), which the current
wording elides. Add the qualification in both places, and cite OptKnock (Burgard et al. 2003,
new bib entry — see §11) as the canonical example of that combinatorial search, in
`capstones.tex` where gene-deletion search is discussed.

---

## §9. Capstone wording fixes

**File:** `chapter/sections/capstones.tex`

- **"Two published systems"** → the Vilkhovoy2023 cell-free integration is a bioRxiv preprint
  (confirmed in `References.bib`: `journal = {bioRxiv}`, `note = {preprint}`), not a published
  journal article. Replace with "Two reported systems" (matches the audit's suggested fix,
  applies to both capstones without singling one out awkwardly).
- **"on the order of sixty metabolites tracked"** → the audit states the cited work reports the
  exact number, 63; replace "on the order of sixty" with "63."
- **Wayman/OptKnock claim — verified against the actual paper** (Wayman et al. 2019,
  *Metabolic Engineering Communications* 9:e00088): the single Δgnd knockout, tested
  experimentally after being built, *did* raise glycan production nearly 3-fold over wild
  type — more than the Δsdh–Δgnd double-knockout combination (~2.5-fold) — so "the best single
  knockout" framing is factually correct and should be **kept**. What needs fixing is a
  different, more subtle overclaim: the current text says this 3-fold outcome was "an outcome
  the model predicted by re-solving the same linear program under the modified bound." Per the
  paper, the constraint-based model only *selected/ranked candidate knockouts* (by predicted
  growth-coupling); the quantitative 3-fold increase was measured experimentally after the
  strain was built, not itself an LP-predicted number. Fix: rephrase so the model's role is
  "selecting which bound to close" (direction), not forecasting the magnitude, e.g.: "...the
  best single knockout the model identified, chosen for its predicted growth-coupling among
  the candidate deletions, raised glycan production roughly threefold over the unmodified
  strain once built and measured; the linear program re-solved under each candidate's modified
  bound predicted which deletion to make, not the magnitude of the yield increase eventually
  observed."

---

## §10. Park data cross-organism relabeling

**File:** `chapter/sections/example_urea.tex`

The opening sentence, "The urea cycle in HL-60 cells is small enough to make the opposite case
concretely..." implies an HL-60-specific prediction, but the Park et al. saturation data is
explicitly cross-organism (human/mouse/yeast/*E. coli*; see `park_saturation.csv`'s own
`organism` column). The chapter already acknowledges this later in the section but the opening
framing oversells it. Soften the opening claim along the lines the audit suggests: "the urea
cycle in HL-60 cells, parameterized here with cross-organism public data, is small enough to
make the opposite case concretely..." — keeps the worked example's identity (it's still about
the HL-60 urea cycle network) while being honest that the *parameters* are not HL-60-specific.

---

## §11. Related-methods positioning paragraph

**File:** `chapter/sections/gateways.tex`
**New bib entries required in `chapter/References.bib`** (none of these six exist yet —
verified against the current file):

```bibtex
@article{Sanchez2017,
  author  = {S{\'a}nchez, Benjam{\'i}n J. and Zhang, Cheng and Nilsson, Avlant and
             Lahtvee, Petri-Jaan and Kerkhoven, Eduard J. and Nielsen, Jens},
  title   = {Improving the Phenotype Predictions of a Yeast Genome-Scale Metabolic Model by
             Incorporating Enzymatic Constraints},
  journal = {Molecular Systems Biology},
  volume  = {13}, number = {8}, pages = {935}, year = {2017},
  doi     = {10.15252/msb.20167411}
}

@article{OBrien2013,
  author  = {O'Brien, Edward J. and Lerman, Joshua A. and Chang, Roger L. and
             Hyduke, Daniel R. and Palsson, Bernhard {\O}.},
  title   = {Genome-Scale Models of Metabolism and Gene Expression Extend and Refine Growth
             Phenotype Prediction},
  journal = {Molecular Systems Biology},
  volume  = {9}, pages = {693}, year = {2013},
  doi     = {10.1038/msb.2013.52}
}

@article{Henry2007,
  author  = {Henry, Christopher S. and Broadbelt, Linda J. and Hatzimanikatis, Vassily},
  title   = {Thermodynamics-Based Metabolic Flux Analysis},
  journal = {Biophysical Journal},
  volume  = {92}, number = {5}, pages = {1792--1805}, year = {2007},
  doi     = {10.1529/biophysj.106.093138}
}

@article{Covert2001,
  author  = {Covert, Markus W. and Schilling, Christophe H. and Palsson, Bernhard},
  title   = {Regulation of Gene Expression in Flux Balance Models of Metabolism},
  journal = {Journal of Theoretical Biology},
  volume  = {213}, number = {1}, pages = {73--88}, year = {2001},
  doi     = {10.1006/jtbi.2001.2405}
}

@article{Chandrasekaran2010,
  author  = {Chandrasekaran, Sriram and Price, Nathan D.},
  title   = {Probabilistic Integrative Modeling of Genome-Scale Metabolic and Regulatory
             Networks in {Escherichia coli} and {Mycobacterium tuberculosis}},
  journal = {Proceedings of the National Academy of Sciences},
  volume  = {107}, number = {41}, pages = {17845--17850}, year = {2010},
  doi     = {10.1073/pnas.1005139107}
}

@article{Burgard2003,
  author  = {Burgard, Anthony P. and Pharkya, Priti and Maranas, Costas D.},
  title   = {{OptKnock}: A Bilevel Programming Framework for Identifying Gene Knockout
             Strategies for Microbial Strain Optimization},
  journal = {Biotechnology and Bioengineering},
  volume  = {84}, number = {6}, pages = {647--657}, year = {2003},
  doi     = {10.1002/bit.10803}
}
```

**Paragraph placement:** a new paragraph in `gateways.tex`, after the "every gateway just
opened demands data..." paragraph (which introduces the information-free-defaults idea) and
before the sequence-specific expression extension paragraph that currently closes the section.
Style: continuous prose (matches the chapter's established convention — no tables), one
paragraph mapping each gateway to its established predecessor:
- thermodynamic gateway ($\delta_j$) ↔ thermodynamics-based flux analysis (Henry2007)
- kinetic gateway ($V^{\circ}_{\max}$, $f_j$) ↔ enzyme-constrained FBA / GECKO (Sanchez2017)
- expression gateway ($e/e^{\circ}$) ↔ metabolism-and-expression (ME) models (OBrien2013)
- regulatory gateway ($\theta_j$) ↔ regulatory FBA (Covert2001) and PROM (Chandrasekaran2010)
State plainly, per both the audit and all three peer reviewers, that the chapter's
contribution is the unifying pedagogical notation and workflow, not the invention of the
individual constraints.

---

## §12. Toy-scale ($n-m=1$) limitation acknowledgment

**File:** `chapter/sections/capstones.tex`

The urea example has $m=18$, $n=19$ ($n-m=1$); the feedback example's metabolic submatrix has
$m=3$, $n=4$ ($n-m=1$). Both are far from the genome-scale underdetermination
(`linearprogram.tex`'s $\mathrm{rank}(\mathbf S)\ll n$) that motivates the whole gateway
framework. `capstones.tex`'s opening paragraph already contrasts "networks small enough that a
hand calculation could check the result" against the two at-scale capstones — add one clause
making the $n-m=1$ point explicit (e.g., "...on networks small enough that a hand calculation
could check the result against a known truth — one degree of freedom in each case, $n-m=1$, far
short of the genome-scale null space §\ref{sec:lp} describes...") rather than opening a new
paragraph.

---

## Sequencing

1. §11 (new bib entries) must land before §11's positioning paragraph is written (citations
   must resolve) — otherwise independent of everything else.
2. §4 (UQ common random numbers, code) must land before §2/§3/§5 (UQ prose fixes) since §3
   quotes CSV-derived numbers that could in principle shift (verified above that they don't,
   but the regenerated CSV/figure should exist before the prose that describes it is finalized).
3. Everything else (§1, §6, §7, §8, §9, §10, §12) is independent prose editing across five
   files and can proceed in any order.
4. Full rebuild + both Julia test suites at the end.
