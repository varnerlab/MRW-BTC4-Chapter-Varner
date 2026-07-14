# Phase 1: Correctness Fixes at the Model Boundary

**Date:** 2026-07-14
**Source:** `chapter_audit_results.md` (technical/narrative audit) + `peer-review.md`
  (three simulated reviewers), reconciled during brainstorming.
**Sections touched:** `chapter/sections/gateways.tex`, `chapter/sections/example_feedback.tex`
**Code touched:** `code/fba/urea_cycle.jl`, `code/fba/urea_cycle_uq.jl`, `code/feedback/dual_feedback.jl`,
  `code/feedback/run_feedback.jl`, tests in both directories

## Why this phase is first

The audit's findings span code correctness, prose claims, and narrative structure. Prose in
Phase 2/3 quotes numbers this phase produces (urea flux values, feedback ledger values, FVA
results), so the code and math have to be right before anything downstream is rewritten.
This is sub-project 1 of a larger decomposition; Phases 2-4 are out of scope here.

## Explicitly out of scope for this phase

**The mole-balance derivation (`chapter/sections/derivation.tex`) is not being changed.**
The audit (item #1) claimed the $d\bar V/dt=0$ assumption is inconsistent with ordinary batch
culture. This was checked against the canonical source
([`CHEME-5430-Advanced-Derivation-FluxBalanceAnalysis`](https://github.com/varnerlab/Lecture-5430-FluxBalanceAnalysis/blob/main/CHEME-5430-Advanced-Derivation-FluxBalanceAnalysis-Spring-2026.ipynb))
and confirmed to be a misreading: $\bar V$ (the specific/per-biomass volume) is the fixed
quantity by construction, and $V=B\bar V$ moves as $B$ grows — exactly what the chapter and the
notebook both derive. The chapter's derivation already mirrors the notebook's algebra step for
step (same $V=B\bar V$ split, same product-rule expansion, same three assumptions), differing
only in the deliberate final step where the chapter retains $\mu x$ instead of dropping it —
which is the chapter's stated thesis. No action item here.

Database provenance (audit item, BRENDA/eQuilibrator query conditions) is deferred as a
placeholder TODO per an earlier scope decision — not part of this phase.

---

## 1. Urea-cycle seconds-to-hours unit fix

**File:** `code/fba/urea_cycle.jl`

**Problem:** `KCAT0` is documented and used as `s^-1` (BRENDA turnover numbers), `E0` is
`mmol/gDW`, and `Vmax = kcat .* e0` (line 148) is computed with no conversion. The result is in
`mmol gDW^-1 s^-1`, but every downstream consumer — the worked-example prose, the bound
comments (`mmol gDW^-1 h^-1`), the UQ code, the figures — treats it as `mmol gDW^-1 h^-1`. Every
reported urea flux is low by 3600x. Worked example: argininosuccinate lyase,
$3.28\ \text{s}^{-1}\times0.01\ \text{mmol gDW}^{-1} = 0.0328\ \text{mmol gDW}^{-1}\text{s}^{-1}
= 118.08\ \text{mmol gDW}^{-1}\text{h}^{-1}$.

**Fix:** `Vmax = kcat .* e0 .* 3600.0` (or introduce a named `SECONDS_PER_HOUR = 3600.0`
constant and multiply by it), with a one-line comment noting the conversion explicitly rather
than leaving it implicit. Add a unit test asserting the conversion factor directly (hardening
item, see §4) rather than only checking hardcoded downstream bound values.

**Ripple effects to regenerate/update:**
- `code/fba/test_urea_cycle.jl`: `expected_lb`/`expected_ub` literals (currently `[-0.1,
  -0.0328, 0.0, 0.0, 0.0]` / `[0.1, 0.0328, 1.9, 4.1, 0.1]`) and the `v[b4] ≈ 0.0328` assertion
  all scale by 3600x.
- `code/fba/urea_cycle_uq.jl`: regenerate `urea_fba_uq.csv`, `urea_uq_sensitivity.csv`, the
  `configA_urea_export ≈ 0.0328` assertion (rtol 0.15), figures `urea_fba.pdf`,
  `urea_saturation.pdf`.
- `chapter/sections/example_urea.tex` (and any figure captions): every quoted flux value,
  including the "sharp answer" `0.0328` and any UQ interval numbers, must be updated to the
  corrected (3600x) values. A grep for `0.0328`, `kcat`, and other quoted flux literals across
  `chapter/sections/` is needed at implementation time to find every occurrence — not
  enumerated exhaustively here.
- The qualitative result (v2/argininosuccinate lyase remains the bottleneck) is unchanged; only
  magnitudes and their uncertainty intervals move.

---

## 2. Align the general protein balance with the feedback implementation

**File:** `chapter/sections/gateways.tex`, Eq.~\eqref{eq:protein} and Eq.~\eqref{eq:expression-ss}

**Problem:** The general balance is
$$\dot p_j = r_{L,j}\,w_j - (\theta_{p,j}+\mu)\,p_j,$$
with no dependence on transcript abundance $m_j$ — translation is driven only by a dimensionless
control variable $w_j$. But `code/feedback/dual_feedback.jl`'s `Dual-Feedback.toml` declares
`rTL::{m}` with kinetic order 1 (translation rate $\propto m$, asserted explicitly:
`model.G[trow("m"), idx["rTL"]] == 1.0`), i.e. the actual implementation makes translation
first-order in the transcript, not merely gated by a dimensionless switch.

**Fix:** Add the missing $m_j$ dependence:
$$\dot p_j = r_{L,j}\,w_j\,m_j - (\theta_{p,j}+\mu)\,p_j,$$
and update the steady-state form in Eq.~\eqref{eq:expression-ss} accordingly:
$$p^{\star}_j = \frac{r_{L,j}w_j\,m^{\star}_j}{\theta_{p,j}+\mu}.$$
This is a more standard ME-model-style translation term (rate proportional to ribosome capacity,
a control variable, and available transcript) and matches what's actually implemented. Contained
edit to `gateways.tex` only — `example_feedback.tex` cites Eq.~\eqref{eq:protein} but doesn't
restate its RHS, so no contradiction is introduced there. A grep for other `\eqref{eq:protein}`
or `\eqref{eq:expression-ss}` citations across `chapter/sections/` should be checked at
implementation time to confirm nothing else assumes the old form.

---

## 3. Feedback example: bounded, independently-motivated $\theta(X_3)$

**Files:** `code/feedback/dual_feedback.jl`, `code/feedback/run_feedback.jl`,
`code/feedback/test_dual_feedback.jl`, `chapter/sections/example_feedback.tex`

**Problem (audit item #5 / peer-review Weakness 1 + 2):** Two distinct issues, both rooted in
the same line. `gateway_factors()` computes
`θ = rate_with / rate_noallo`, which evaluates to exactly `X3^{-a}` — algebraically identical to
the S-system truth model's own kinetic term for the allosteric loop, evaluated at the exact
converged truth trajectory `Xss`. This is (a) unbounded ($X_3^{-a}>1$ for $X_3<1$, violating the
chapter's own definition of $\theta\in[0,1]$ as "the fraction of capacity switched on"), and (b)
circular: since the FBA bound `ub[r0] = Vmax0 * e_e0 * θ` is built from the *exact same
functional form* evaluated at the *exact same point* that generates the true rate, `ub[r0]`
equals the true `r0` rate by construction, not by inference — the LP reports back an identity,
not a prediction.

**Fix:** Replace `θ = X3^{-a}` with a genuinely bounded two-state Hill/partition-function
occupancy, independent of the S-system's power-law kinetics:
$$\theta(X_3) = \frac{K^n}{K^n + X_3^{\,n}}.$$
No leak floor: unlike the transcriptional loop (which retains a basal leak $\lambda$, justified
in the existing prose — "the gene is never silenced completely"), a two-state allosteric
on/off switch has no such floor, so $\theta\to0$ as $X_3\to\infty$ is the correct asymptote.
This is the same Hill-type occupancy form already used for $f_s$ in the partition-function
machinery (Eq.~\eqref{eq:partition}), so it stays inside the chapter's own established
vocabulary rather than introducing a new device.

**Proposed starting parameters (explicitly provisional, to iterate on):** $K=5$, $n=2$ — a mild
cooperativity, chosen independently of the truth model's $a=0.4$ exponent, not fit to reproduce
it. The truth model's own S-system kinetics (`feedback_truth()`, `a=0.4`) are **unchanged** —
only the FBA-layer gateway function used to interpret the measured $X_3$ changes.

**Resulting numbers at these parameters** (Vmax0=10, e/e0=0.464 and $T^\star\approx2.66$
unchanged, since the truth model is untouched):

| case | old θ (X3^-a) | old ub | new θ (Hill, K=5,n=2) | new ub |
|---|---|---|---|---|
| naive (both off) | 1 | 10 | 1 | 10 |
| expression-only | — | 4.64 | — | 4.64 (unchanged) |
| activity-only | 0.574 | 5.74 | 0.610 | 6.10 |
| both open | 0.574 | **2.66** (exact) | 0.610 | **≈2.83** (~6% over truth) |

The escalation story (naive overshoots most; single gateways bracket; only the product
approaches truth) survives. "Both open" no longer equals truth to numerical tolerance — it
approximates it, which is the point. Prose in `example_feedback.tex` needs to drop "recovering
the BST truth to numerical tolerance" language in favor of something like "approximating the
true throughput to within about 6%," consistent with both the audit's and peer review's request
to soften "recovered"/"truth" language.

**Code changes:**
- `gateway_factors()`: replace the `_powerlaw`-based `rate_with`/`rate_noallo`/`θ` computation
  with a direct evaluation of the Hill formula at `truth.Xss[X3]`. This also fully removes this
  function's dependency on `BSTModelKit._powerlaw` (a separate hardening item, §4) — the only
  remaining use of `_powerlaw` in this file is in `reaction_fluxes()`, which legitimately needs
  the truth model's own rate law to report integrated fluxes.
- `run_feedback.jl`: ledger numbers and figure update automatically once `gateway_factors`
  changes; no structural change expected.
- `test_dual_feedback.jl` Task 3: numeric targets change (`gw.θ ≈ 0.574` → `≈0.610` at the
  proposed parameters) and the "both-open must recover truth" assertion loosens from `< 1e-2`
  absolute to a stated relative tolerance (e.g. 10%), since exact recovery is no longer the
  claim being tested. The bracketing assertions (`T < Ee < N`, `T < Aa < N`, naive overshoot)
  should still hold and stay as tight checks.
- Side benefit: the Hill form is well-behaved at $X_3=0$ ($\theta\to1$), removing the current
  code's implicit requirement that $X_3$ stay strictly positive to avoid a blowup under
  $X_3^{-a}$.

---

## 4. Code hardening pass

Full scope, per earlier decision. None of these change chapter claims; all reduce risk of a
silent numerical error going undetected before submission.

- **Input validation:** vector-length checks, finite-value checks, nonnegative saturation
  factors, `lb <= ub` in `urea_cycle_model()` and `feedback_fba()`.
- **Explicit solve-failure handling:** `solve_fba()` currently builds a `DataFrame` from
  `solve_flux()`'s result even when that's `nothing` on an infeasible/unbounded solve. Raise or
  return a clearly-tagged failure instead of silently constructing a DataFrame of `nothing`.
- **UQ termination-status recording:** `urea_cycle_uq.jl`'s `run_ensemble()` currently just
  `continue`s past a `nothing` (failed) draw and reports `kept`. Record termination status per
  draw so the "probability `v2` remains the bottleneck" style claims have a documented failure
  rate, not just a silently-dropped count.
- **Load `park_saturation.csv` properly:** `urea_cycle_uq.jl`'s `SAT_NOMINAL` hardcodes the
  already-reduced (concentration, Km) pairs "reduced from `code/data/park_saturation.csv` by
  the aggregation rule" in a comment, rather than actually loading and reducing the CSV at run
  time as the original UQ design spec (`2026-07-13-urea-cycle-bootstrap-uq-design.md`)
  specifies. Load the CSV and apply the stated aggregation rule (least-saturated substrate;
  prefer *Homo sapiens*, else geometric mean across organisms) in code.
- **Stop relying on `BSTModelKit._powerlaw`:** mostly resolved by the §3 redesign
  (`gateway_factors` no longer needs it). The remaining use in `reaction_fluxes()` should be
  wrapped in a locally-named function with a comment pinning the `BSTModelKit` version in
  `code/Manifest.toml`, so a future package change fails at the wrapper rather than silently
  inside a private call.
- **New tests:** a dimensional-conversion test for the §1 unit fix (assert `Vmax` uses the
  3600x factor directly, not just that downstream bounds match a hardcoded literal); an FVA
  test/routine for the urea model asserting the nominal optimum's uniqueness (the audit's
  "Verification performed" section already confirmed this manually — codify it as a repeatable
  test so it stays checked if the network changes).

---

## Sequencing within this phase

1. §1 (urea unit fix) and §2 (gateways.tex equation) are independent and can be done in either
   order.
2. §3 (feedback redesign) is independent of §1/§2.
3. §4 (hardening) touches the same files as §1 and §3, so do it after those land to avoid
   rebasing validation/error-handling code around numbers that are about to change.
4. Full rebuild (`cd chapter && make`) and full test run (`test_urea_cycle.jl`,
   `test_dual_feedback.jl`) at the end of the phase, before handing off to Phase 2.
