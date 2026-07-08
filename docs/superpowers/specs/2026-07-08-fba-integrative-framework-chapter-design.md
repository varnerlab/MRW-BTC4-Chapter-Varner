# Design: FBA as an Integrative Framework — "Bounds as Gateways"

**Date:** 2026-07-08
**Author:** Jeffrey D. Varner (with Claude)
**Deliverable:** *Comprehensive Biotechnology, 4th Edition* (Elsevier MRW) chapter, solo-author, ~35 pp, **due 2026-07-24**. Assigned slot title: *Mathematical Models in Biotechnology* (kept as the official title; the FBA-integrative framing below is the chapter's argument/through-line). Branch: `chapter-draft`.

**Supersedes:** the six-pillar "mechanistic → data-driven arc" design at
`docs/superpowers/specs/2026-06-23-mathematical-models-biotechnology-chapter-design.md`.
That draft is complete on `chapter-draft` and remains fully recoverable in git history; this pivot rewrites the chapter around flux balance analysis (FBA).

---

## 1. Thesis and intellectual signature

This is **not** a survey of FBA. It is a rigorous, self-contained, opinionated account with two load-bearing ideas:

1. **Bounds are the integrative gateway.** The FBA flux bound is not a throwaway box constraint — it is a *model* that carries the biology. In the general form
   $$-\delta_j\big[\,V^{\circ}_{\max,j}\,(e/e^{\circ})\,\theta_j(\cdot)\,f_j(\cdot)\,\big]\ \le\ \hat v_j\ \le\ V^{\circ}_{\max,j}\,(e/e^{\circ})\,\theta_j(\cdot)\,f_j(\cdot)$$
   every factor is a gateway for a distinct class of information: `δ_j` ← thermodynamics; `V°_max = k°_cat·e°` and `f_j` ← kinetics / metabolite state; `(e/e°)` ← gene expression; `θ_j` ← regulation. Flux estimation and the integration of metabolic/kinetic/regulatory data are *the same problem* — deciding where and how to set the bounds.

2. **Keep the growth-dilution term `μ` everywhere.** The standard literature drops it twice, unexamined. We do not:
   - Metabolite balance → `S·v̂ = μx` (honest); `S·v̂ = 0` is a *labeled special case* with stated validity **and stated failure** (fed-batch, cell-free).
   - Expression balance → the `(θ+μ)` denominators in the steady-state `m*`, `p*`.
   The single consistent move — retain `μ` — is the chapter's signature and the reason it reads as rigorous rather than folk-standard.

**Tone / references policy:** minimal citations, cite only what we build on, no survey-padding, no genuflecting to incremental literature. Present *an alternative, correct, complete account*.

**Source materials (all Varner-lab):**
- CHEME-5430 FBA lecture + **advanced derivation** notebook (the keep-`μ` derivation).
- CHEME-5450 lecture notes: L5c/L6a flux-bounds model, L6a gene-expression sub-model, L6c **Boltzmann/partition-function** control function.
- Vilkhovoy 2018 (*ACS Synth Biol* 7:1844) — sequence-specific constraint-based CFPS.
- Allen & Palsson 2003 (*J Theor Biol* 220:1) — sequence-based metabolic demands for protein synthesis.
- Adhikari 2020 (*Front Bioeng Biotechnol* 8:539081, PMID 33324619) — **the** effective biophysical (= partition-function/Boltzmann) cell-free TX/TL model. **Adhikari and "Boltzmann" are the same object, not two.**
- Vilkhovoy 2023 (bioRxiv 2023.02.10.528035) — integrated *dynamic* constraint-based cell-free model with real time-resolved intracellular data (capstone).
- Wayman 2019 (*Metab Eng Commun* 9:e00088) — model-guided metabolic engineering / designer glycans (capstone).

---

## 2. Chapter spine

1. **Introduction.** Flux is the phenotype; estimation is the problem; bounds are the gateway. Promise a rigorous, self-contained account.
2. **From mole balances to the constraint (keep `μ`).** Open species mole balance → concentration balance → specific units `V = B·V̄` → `S·v̂ = μx`; `S·v̂ = 0` as special case with validity/failure stated; exchange reactions for openness.
3. **The linear program and its geometry.** The FBA LP; objective choice (biomass vs. product); null-space structure of `S` via SVD (right null space = feasible flux modes; left null space = conserved moieties); FVA for what the constraints leave undetermined.
4. **Bounds as gateways** (the heart). The product bound, walked factor by factor: thermodynamic `δ` → kinetic capacity + saturation `V°_max, f` → expression `(e/e°)` sub-model (with kept `μ`) → regulation `θ` via the **partition-function control function** (Adhikari/Boltzmann). Simplified bounds model as the tractable baseline; gateways opened one at a time.
5. **Worked example 1 — urea-cycle metabolism** (BUILT, steady-state). Baseline LP with simplified bounds (`δ`, `V°_max`); then *open a gateway* and show flux redistribution. Exercises the thermodynamic + kinetic gateways from **databases** (BRENDA `k_cat`, eQuilibrator `ΔG`) — honest, no intracellular measurement needed.
6. **Worked example 2 — regulated sequence-specific cell-free expression** (BUILT, steady-state). The expression + regulation gateways made concrete: sequence-specific TX/TL with precursors supplied by exchange (no metabolic resource limitation, à la Allen–Palsson), transcription regulated by a partition-function control function (Adhikari/Boltzmann). Exercises the `(e/e°)` and `θ` gateways.
7. **Capstones (reviewed, figure + prose, not rebuilt).** Vilkhovoy 2023 — §5 and §6 coupled back together, dynamically, with real intracellular data (cell-free = where intracellular *is* the reactor, so the kept term is honestly measurable). Wayman 2019 — bounds as engineering design levers.
8. **Outlook.** Brief, rigorous, non-crap.
- **Appendix.** Long derivations: full open-system mole balance; SVD null-space structure; steady-state expression `m*`, `p*`.

**Structural logic:** §5 is the metabolism half, §6 is the expression half — both steady-state built examples exercising different gateways; §7 shows them united in the real dynamic system.

---

## 3. Key mathematics the chapter must contain

**Derivation (keep `μ`), §2 / appendix:**
$$\sum_{s}d_s C_{i,s}\dot V_s + \sum_j \sigma_{ij}\hat v_j V = \frac{d}{dt}(C_i V), \quad V=B\bar V$$
Under steady intracellular pools, constant culture volume, no physical transport →
$$C_i\mu = \sum_j \sigma_{ij}\hat v_j \iff \mathbf{S}\hat{\mathbf v}=\mu\mathbf{x};\qquad \text{drop } \mu\mathbf{x}\Rightarrow \mathbf{S}\hat{\mathbf v}=\mathbf 0 \ \text{(Palsson, special case)}.$$
State explicitly: valid when growth is slow vs. metabolism and pools are dilute; **fails** for fed-batch (`dV̄/dt≠0`) and cell-free (`μ=0` but `dx/dt≠0`, transport ≠ 0).

**LP, §3:** `max cᵀv̂  s.t.  S·v̂ = μx (or 0),  L ≤ v̂ ≤ U`. SVD `S = UΣVᵀ`; `N(S)` from zero-singular-value columns of `V`; conserved moieties from left null space.

**Bounds gateway, §4:** the product bound above, plus:
- `δ_j ∈ {0,1}` from `sign(ΔG° − ΔG*)` or `K_eq` cutoff.
- `V°_max,j = k°_cat,j·e°`; `f_j(·)` substrate saturation (Michaelis–Menten form — re-homed from old kinetics section).
- `(e/e°)` from expression sub-model.

**Expression sub-model, §4 / appendix (keep `μ`):**
$$\dot m_j = r_{X,j}u_j - (\theta_{m,j}+\mu)m_j + \lambda_j,\qquad \dot p_j = r_{L,j}w_j - (\theta_{p,j}+\mu)p_j$$
$$m^{\star}_j = \frac{r_{X,j}u_j+\lambda_j}{\theta_{m,j}+\mu},\qquad p^{\star}_j = \frac{r_{L,j}w_j}{\theta_{p,j}+\mu};\qquad (e/e^\circ)=p^\star_j/e^\circ \ \text{via GPR}.$$

**Partition-function control function `θ`/`u`, §4 (Adhikari/Boltzmann):**
$$p_s = \frac{f_s e^{-\beta\epsilon_s}}{Z},\quad Z=\sum_{s}f_s e^{-\beta\epsilon_s};\qquad \bar u = \underbrace{\sum_{s\in\mathcal A}p_s}_{u}+\underbrace{\sum_{s\in\mathcal B}p_s}_{u^\dagger},\quad \lambda \equiv r_X u^\dagger.$$
`f_s` Hill-type in inducer/TF concentration; `e^{-βε_s}` estimated from data.

**Sequence-specific TX/TL, §6 (Allen–Palsson / Vilkhovoy 2018):**
- Transcription: `Σ_{x∈{A,U,G,C}} n_{x,G}·NTP_x → mRNA_G + PPi(·)` — coefficients from the gene's nucleotide composition.
- Translation: `Σ_a n_{a,P}·(charged aa_a) + ~2·L_P·GTP + charging energy → protein_P + tRNA + GDP + Pi` — coefficients from the protein's residue composition and length.
- Kinetic limits `r_{X,G}`, `r_{L,G}` (elongation-rate based) set the expression `V°_max`.
- Precursors (NTP, aa, GTP/ATP energy) supplied by **exchange reactions** — no metabolic network modeled.

---

## 4. Built examples (detailed)

### Example 1 — Urea cycle (metabolic, steady-state) — REUSE existing code
- Network already in repo (`code/fba/`, 5 enzymatic + 14 exchange reactions, HL-60). `S ∈ ℝ^{18×19}`.
- Objective: maximize urea export. Solver: JuMP + HiGHS (keep) or GLPK (5430 uses GLPK) — pick one, keep consistent.
- **Baseline:** simplified bounds (`δ` from ΔG, `V°_max` from `k_cat`). Reproduce the argininosuccinate-lyase bottleneck.
- **Open a gateway:** impose a saturation `f_j` or activity `θ_j` factor on one reaction and show the optimal flux redistribute — demonstrating the bound carries biology.
- Figure: flux bar chart (adapt existing `urea_fba.pdf`).

### Example 2 — Regulated sequence-specific cell-free expression (steady-state) — NEW code
- **Product:** deGFP (one of the two proteins Vilkhovoy 2018 had data for). Compute TX/TL stoichiometry from the actual deGFP nucleotide + amino-acid composition → genuinely sequence-specific and reproducible.
- **Network:** exchanges (NTP, aa, GTP/ATP) → `TX_deGFP` → mRNA → `TL_deGFP` → deGFP → export. Precursors unbounded (no resource limitation).
- **Objective:** maximize deGFP synthesis (translation) flux — protein productivity.
- **Regulation (the gateway):** transcription upper bound `= r_{X}·ū(I)`, with `ū(I)` from a small (2–3 microstate) partition-function promoter model whose `f_s` are Hill functions of an inducer `I`. Sweep `I` → `ū(I)` → transcription bound → deGFP productivity. Shows: *sequence sets the demand (NTP/aa/energy per deGFP), regulation sets the achievable rate through the bound.*
- **`μ` note:** cell-free ⇒ `μ=0`, so `(θ+μ)→θ`; the general kept-`μ` derivation still stands, cell-free is its `μ→0` limit. Consistency, not contradiction.
- Solver stack identical to Example 1 (one `code/Project.toml`, one `Include.jl`).
- Figure(s): productivity vs. inducer (gateway opening); sequence-specific demand breakdown.

### Reviewed capstones (§7) — NOT built
- **Vilkhovoy 2023:** dynamic integrated cell-free model; §5⊕§6 coupled with real 63-metabolite + mRNA/protein/enzyme-activity time-series; the honest home of the kept term. One reproduced/adapted figure + prose.
- **Wayman 2019:** designer-glycan model-guided engineering; bounds as design levers. One figure + prose.

---

## 5. Disposition of existing content

| Asset | Action |
|---|---|
| `sections/fba.tex`, `code/fba/` (urea cycle) | **KEEP / adapt** — becomes Example 1; rewrite prose around keep-`μ` derivation + gateways. |
| `sections/kinetics.tex` (Michaelis–Menten) | **RE-HOME** → the `f_j` saturation gateway in §4. Delete as standalone pillar. |
| `sections/geneexpression.tex` (motifs) | **CUT** prose + `code/geneexpression/`. Gene-expression *sub-model* re-homed to `(e/e°)` gateway; network motifs removed. |
| `sections/deeplearning.tex`, `code/deeplearning/` (LSTM/S4) | **CUT** entirely (prose + code). |
| `sections/hybrid.tex` | **CUT** entirely. |
| `sections/introduction.tex`, `appendix.tex` | **REWRITE** to the new spine. |
| `code/kinetics/` (CHO ODE) | **PARK** — only needed if a dFBA example is later added; not used by Examples 1–2. Keep in tree for now, unwired. |
| Old six-pillar spec | Leave in place; superseded, recoverable via git. |
| `code/Include.jl`, `Project.toml`, `Manifest.toml` | **KEEP** — rewire to the two built examples; still commit `Manifest.toml`. |

---

## 6. Conventions (unchanged from prior locked decisions)
- LaTeX: plain `article` class; build via `chapter/Makefile` (pdflatex → bibtex → pdflatex ×2).
- Full derivations inline (text or `appendix.tex`).
- Julia in `code/`, one `Include.jl` + `Project.toml`, commit `Manifest.toml`. Figures generated into `code/figs/`, copied into `chapter/figures/`.
- One reproducible run per example: `cd code && julia --project=. <example>/run_*.jl`.

---

## 7. Open / parked decisions (require author input at spec review)
1. **dFBA example — PARKED.** Whether to build any dynamic example at all. Current plan builds none; all dynamics live in the reviewed Vilkhovoy 2023 capstone. If added later, it needs intracellular `x(t)` (synthetic, twin-experiment) and reuses `code/kinetics/`.
2. **Example 2 specifics — proposed above, confirm at review:** deGFP as product; 2–3 microstate inducible promoter as the circuit; two lumped sequence-specific TX/TL reactions (not elongation-resolved). Alternative circuit available: the local gluconate-sensor structured promoter.
3. **Solver:** standardize on HiGHS (existing chapter) vs. GLPK (5430 notebooks). Recommend HiGHS.
4. **Official title:** keep assigned "Mathematical Models in Biotechnology" vs. request an FBA-specific title from the editors.

---

## 8. Risks / de-risking cut order (if time tight before 2026-07-24)
1. Make the Wayman capstone review a single paragraph + figure.
2. Simplify Example 2's promoter to a 2-microstate ON/OFF (still partition-function, fewer parameters).
3. Reduce the gateway-opening in Example 1 to a described (not separately coded) variant.
