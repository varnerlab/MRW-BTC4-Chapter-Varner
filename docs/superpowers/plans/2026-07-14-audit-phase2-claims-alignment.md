# Phase 2 Claims Alignment Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Align chapter prose with what the (Phase-1-corrected) code and the cited literature
actually support: label the thermodynamic gate as a heuristic, fix the UQ pseudocode/prose and
give Config A/B common random numbers, write up the FVA result, soften several overstated
claims, add a related-methods positioning paragraph with verified citations, and fix capstone
wording (including a literature-verified correction to the Wayman knockout claim).

**Architecture:** Almost entirely LaTeX prose edits across five `chapter/sections/*.tex` files
plus one bib file. One small Julia code change (a single-line RNG-seed fix) with a
verified-in-advance regeneration step.

**Tech Stack:** Julia 1.12 (JuMP + HiGHS, only for the one UQ script re-run), plain
`article`-class LaTeX built via `chapter/Makefile`.

## Global Constraints

- Spec: `docs/superpowers/specs/2026-07-14-audit-phase2-claims-alignment-design.md`.
- Depends on Phase 1 (merged to `main` and `chapter-draft`) — all numbers quoted below assume
  the Phase 1 unit fix and bounded feedback gateway are already in place.
- Julia invocation is always `julia --project=code <script path>` from the repo root; do not
  `cd` into `code/` first.
- Figures must be copied from `code/figs/` to `chapter/figures/` before `make` picks them up.
- The one code-level fix in this plan (Task 4) was verified in advance in a scratch copy during
  planning: it changes `code/data/urea_uq_sensitivity.csv` and `code/figs/urea_saturation.pdf`
  but does **not** change any number already quoted in the chapter prose (both round to the
  same values). `code/data/urea_fba_uq.csv` and `code/figs/urea_fba.pdf` are unaffected (Config
  A's seed is untouched) and should come out byte-identical to what's already committed.

---

## Task 1: Add six new bibliography entries

**Files:**
- Modify: `chapter/References.bib`

**Interfaces:**
- Produces: six new BibTeX keys — `Sanchez2017`, `OBrien2013`, `Henry2007`, `Covert2001`,
  `Chandrasekaran2010`, `Burgard2003` — consumed by Task 2 (positioning paragraph) and Task 10
  (estimation/design qualification, cites `Burgard2003`).

- [ ] **Step 1: Append the six entries to the end of the file**

`chapter/References.bib` currently ends (line 154) with the closing `}` of the `Park2016`
entry. Append, after that closing brace:

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

- [ ] **Step 2: Verify the bib file still parses correctly**

Run: `cd chapter && bibtex Chapter 2>&1 | tail -20 && cd ..`
Expected: no `"..." is undefined` or syntax errors from the six new entries specifically (an
unrelated pre-existing "There were N warning(s)" summary line, if present, is not a new
regression from this edit — only investigate if it names one of the six new keys).

- [ ] **Step 3: Commit**

```bash
git add chapter/References.bib
git commit -m "Add related-methods citations (GECKO, ME-models, TFA, rFBA, PROM, OptKnock)"
```

---

## Task 2: Add the related-methods positioning paragraph to `gateways.tex`

**Files:**
- Modify: `chapter/sections/gateways.tex`

**Interfaces:**
- Consumes: the six bib keys from Task 1.
- Produces: no interface consumed by later tasks.

- [ ] **Step 1: Insert the new paragraph**

In `chapter/sections/gateways.tex`, insert a new paragraph between the "Every gateway just
opened demands data..." paragraph (ending "...which gateway must be reopened to relax each of
them.") and the "The expression and regulatory gateways lumped transcription and translation
into aggregate capacities..." paragraph that currently closes the section. That is, insert
immediately after:

```latex
Equation~\eqref{eq:simplified-bound} retains only the thermodynamic gate
and the capacity scale, and it is exactly the reversible or irreversible
$\pm V_{\max}$ box used throughout genome-scale flux balance analysis.
It is the tractable baseline, the honest statement of a bound when
nothing beyond a capacity estimate is known, and
Equation~\eqref{eq:general-bound} records precisely which assumptions
that baseline makes and which gateway must be reopened to relax each of
them.
```

and immediately before:

```latex
The expression and regulatory gateways lumped transcription and
translation into aggregate capacities $r_{X,j}$ and $r_{L,j}$, but these
```

insert this new paragraph:

```latex
Each gateway just described has a named, established predecessor in the
constraint-based modeling literature, and this chapter's contribution is
the unifying notation and workflow that connects them, not the invention
of any individual constraint. The thermodynamic gateway $\delta_j$ is the
central device of thermodynamics-based flux analysis, which propagates
standard free energies and metabolite activity ranges through the whole
network to rule out infeasible flux and concentration profiles together
\cite{Henry2007}. The kinetic gateway, $V^{\circ}_{\max,j}$ built from a
turnover number and an enzyme abundance together with the saturation
factor $f_j$, is the organizing idea behind enzyme-constrained flux
balance analysis, most fully realized in the GECKO framework, which
folds a proteome-wide capacity constraint into a genome-scale
reconstruction \cite{Sanchez2017}. The expression gateway $(e/e^{\circ})$,
transcript and protein balances entered as a modeled quantity rather
than a fixed parameter, is the metabolism-and-expression (ME) model
construction, which couples the transcription and translation machinery
directly to the metabolic network at genome scale
\cite{OBrien2013}. The regulatory gateway $\theta_j$ has two established
lines of descent: regulatory flux balance analysis, which gates
reactions on and off with Boolean rules derived from a transcriptional
regulatory network \cite{Covert2001}, and PROM, which replaces the
Boolean rule with a continuous, data-derived probability of an
interaction being active \cite{Chandrasekaran2010}, the same move from a
binary switch to a graded control fraction that the partition-function
construction of Equations~\eqref{eq:partition}
and~\eqref{eq:control-split} makes explicit. What the gateway
factorization of Equation~\eqref{eq:general-bound} adds is not a fifth
constraint class but a common accounting of where each of these four
already-established devices sits inside one bound, so that a single
notation and a single worked derivation carry across all four rather
than requiring a separate formalism, and a separate paper, for each.
```

- [ ] **Step 2: Rebuild the chapter and confirm all six new citations resolve**

Run:
```bash
cp code/figs/*.pdf chapter/figures/ 2>/dev/null
cd chapter && make && cd ..
```

Run: `grep -c "Sanchez2017\|OBrien2013\|Henry2007\|Covert2001\|Chandrasekaran2010" chapter/Chapter.log`
Expected: this grep is checking the LOG, not the source — instead run:
`grep -i "undefined" chapter/Chapter.log | grep -E "Sanchez2017|OBrien2013|Henry2007|Covert2001|Chandrasekaran2010"`
Expected: no output (empty — none of these five citations used in this paragraph are
undefined). Note `Burgard2003` is not cited in this paragraph — it is used in Task 10.

- [ ] **Step 3: Commit**

```bash
git add chapter/sections/gateways.tex
git commit -m "Add related-methods positioning paragraph to Bounds as Gateways"
```

---

## Task 3: Label the thermodynamic gate as a heuristic

**Files:**
- Modify: `chapter/sections/gateways.tex`
- Modify: `chapter/sections/example_urea.tex`

**Interfaces:** none (prose-only, no cross-task dependency).

- [ ] **Step 1: Update `gateways.tex`'s general description of the threshold**

Replace:

```latex
The assignment follows the sign
of $\Delta G^{\circ}_j-\Delta G^{*}_j$, the standard free energy of
reaction measured against a threshold $\Delta G^{*}_j$ fixed by the
physiological concentration range, or equivalently a cutoff on the
equilibrium constant $K_{\mathrm{eq},j}$: a reaction whose driving force
cannot plausibly change sign anywhere in that range is declared
irreversible, and its reverse flux is forbidden. Standard free energies
for this test are drawn from group- and component-contribution
estimates, for which eQuilibrator provides genome-scale coverage
\cite{Beber2022}.
```

with:

```latex
The assignment follows the sign
of $\Delta G^{\circ}_j-\Delta G^{*}_j$, the standard free energy of
reaction measured against a threshold $\Delta G^{*}_j$ meant to
represent the physiological concentration range, or equivalently a
cutoff on the equilibrium constant $K_{\mathrm{eq},j}$: a reaction whose
driving force cannot plausibly change sign anywhere in that range is
declared irreversible, and its reverse flux is forbidden. This is a
heuristic directionality rule, a stand-in for propagating a stated
concentration range through $\Delta G=\Delta G^{\circ}+RT\ln Q$ and
checking directly whether the sign can flip, not a substitute for that
calculation where the range is actually known. Standard free energies
for this test are drawn from group- and component-contribution
estimates, for which eQuilibrator provides genome-scale coverage
\cite{Beber2022}.
```

- [ ] **Step 2: Update `example_urea.tex`'s specific cutoff value**

Replace:

```latex
The threshold
$\Delta G^{*}_j$ of Equation~\eqref{eq:thermo-gate} was fixed here at
$-10$ kJ/mol: a reaction whose standard free energy falls below that
cutoff is declared irreversible and its reverse flux forbidden
($\delta_j=0$), while a reaction whose free energy sits above it is
left reversible ($\delta_j=1$).
```

with:

```latex
The threshold
$\Delta G^{*}_j$ of Equation~\eqref{eq:thermo-gate} was fixed here at
$-10$ kJ/mol, a heuristic cutoff rather than one derived from a shown
physiological concentration range: a reaction whose standard free energy
falls below that
cutoff is declared irreversible and its reverse flux forbidden
($\delta_j=0$), while a reaction whose free energy sits above it is
left reversible ($\delta_j=1$).
```

- [ ] **Step 3: Rebuild and commit**

Run: `cp code/figs/*.pdf chapter/figures/ 2>/dev/null; cd chapter && make && cd ..`
Expected: `chapter/Chapter.pdf` produced, no new undefined-reference warnings.

```bash
git add chapter/sections/gateways.tex chapter/sections/example_urea.tex
git commit -m "Label the -10 kJ/mol thermodynamic cutoff explicitly as a heuristic"
```

---

## Task 4: Give Config A/B common random numbers in the UQ script

**Files:**
- Modify: `code/fba/urea_cycle_uq.jl`
- Generated: `code/data/urea_uq_sensitivity.csv`, `code/figs/urea_saturation.pdf` (regenerated)

**Interfaces:**
- Produces: `code/data/urea_uq_sensitivity.csv` and `code/figs/urea_saturation.pdf` regenerated
  with Config A/B now sharing kcat/e0/dG/measured-$f_j$ draws. `code/data/urea_fba_uq.csv` and
  `code/figs/urea_fba.pdf` are NOT expected to change (Config A's seed is untouched).

- [ ] **Step 1: Change Config B's seed**

Replace:

```julia
let rng = MersenneTwister(SEED + 1)
```

with:

```julia
let rng = MersenneTwister(SEED)  # common random numbers with Config A: vary only f2
```

- [ ] **Step 2: Run the script and confirm output**

Run: `julia --project=code code/fba/urea_cycle_uq.jl`
Expected (values verified during planning in a scratch copy with this exact change):
```
configA_kept=10000
configA_failed=0
configA_urea_export mean=156.04172201526268 sd=157.22177306463496 median=105.05381342754373 ci=[17.060725222471447,613.5765361248797]
configA_v5 mean=0.02396554622889309
configCap_kept=10000
configCap_failed=0
configCap_urea_export mean=167.33961130464067 sd=169.8734191808569 median=112.91950708196956 ci=[17.596821995102992,672.00858977199]
configB_f2 median=0.4976497320733263
configB_urea_export median=50.82740723016741 ci=[4.187265966550702,400.59159073853664]
```
(`configA_*` and `configCap_*` values are unchanged from before this task, since Config A's
seed was not touched — only `configB_*` values shift slightly.)

- [ ] **Step 3: Confirm exactly the expected files changed**

Run: `git status --short`
Expected: only `code/data/urea_uq_sensitivity.csv` and `code/figs/urea_saturation.pdf` show as
modified. If `code/data/urea_fba_uq.csv` or `code/figs/urea_fba.pdf` also show as modified,
stop and check whether Config A's output actually changed (it should not have) before
proceeding — that would indicate the seed change had an unintended effect.

- [ ] **Step 4: Copy the regenerated figure into the chapter**

Run: `cp code/figs/urea_saturation.pdf chapter/figures/urea_saturation.pdf`

- [ ] **Step 5: Commit**

```bash
git add code/fba/urea_cycle_uq.jl code/data/urea_uq_sensitivity.csv code/figs/urea_saturation.pdf chapter/figures/urea_saturation.pdf
git commit -m "Give Config A/B common random numbers in the UQ ensemble"
```

---

## Task 5: Fix Algorithm 1's pseudocode and the "parametric bootstrap" terminology note

**Files:**
- Modify: `chapter/sections/example_urea.tex`

**Interfaces:** none.

- [ ] **Step 1: Rewrite Algorithm 1**

Replace the entire `algorithmic` body:

```latex
  \begin{algorithmic}[1]
    \Require nominal $k^{\circ}_{\mathrm{cat}}$, $e^{\circ}$, $\Delta G^{\circ}$;
      saturation inputs $[S_j], K_{M,j}$; spreads $\sigma_{\ln}, \sigma_{\Delta G}$;
      threshold $\Delta G^{*}$; draws $N$; seed
    \State seed the random generator
    \For{$i = 1$ to $N$}
      \State $e \gets e^{\circ}\exp(\sigma_{\ln} Z)$
        \Comment{one shared abundance draw, $Z\sim\mathcal N(0,1)$}
      \For{each enzymatic reaction $j$}
        \State $k_{\mathrm{cat},j} \gets k^{\circ}_{\mathrm{cat},j}\exp(\sigma_{\ln} Z)$,
          \quad $\Delta G_j \gets \Delta G^{\circ}_j + \sigma_{\Delta G} Z$
        \State $\delta_j \gets \mathbb{1}[\Delta G_j > \Delta G^{*}]$
        \State $f_j \gets [S_j]/(K_{M,j}+[S_j])$ with $[S_j], K_{M,j}$ sampled,
          or $f_j \gets 1$ if no substrate data
        \State $V_{\max,j} \gets k_{\mathrm{cat},j}\,e$;
          \quad $u_j \gets V_{\max,j} f_j$,
          \quad $\ell_j \gets -\delta_j V_{\max,j} f_j$
      \EndFor
      \State $\hat{\mathbf v}^{(i)} \gets
        \arg\max\{\mathbf c^{\top}\mathbf v : \mathbf S\mathbf v = \mathbf 0,\
        \boldsymbol\ell \le \mathbf v \le \mathbf u\}$
    \EndFor
    \State \Return per-reaction mean, standard deviation, and $2.5/50/97.5$
      percentiles of $\{\hat{\mathbf v}^{(i)}\}_{i=1}^{N}$
  \end{algorithmic}
```

with:

```latex
  \begin{algorithmic}[1]
    \Require nominal $k^{\circ}_{\mathrm{cat}}$, $e^{\circ}$, $\Delta G^{\circ}$;
      nominal saturation inputs $\overline{[S_j]}, \overline{K_{M,j}}$ where measured;
      spreads $\sigma_{\ln}, \sigma_{\Delta G}$;
      threshold $\Delta G^{*}$; draws $N$; seed
    \State seed the random generator
    \For{$i = 1$ to $N$}
      \State $e \gets e^{\circ}\exp(\sigma_{\ln} Z_e)$
        \Comment{one shared abundance draw, $Z_e\sim\mathcal N(0,1)$}
      \For{each enzymatic reaction $j$}
        \State $k_{\mathrm{cat},j} \gets k^{\circ}_{\mathrm{cat},j}\exp(\sigma_{\ln} Z_{k,j})$,
          \quad $\Delta G_j \gets \Delta G^{\circ}_j + \sigma_{\Delta G} Z_{\Delta G,j}$
          \Comment{independent $Z_{k,j}, Z_{\Delta G,j}\sim\mathcal N(0,1)$ per reaction}
        \State $\delta_j \gets \mathbb{1}[\Delta G_j > \Delta G^{*}]$
        \If{substrate data measured for reaction $j$}
          \State $[S_j] \gets \overline{[S_j]}\exp(\sigma_{\ln} Z_{S,j})$,
            \quad $K_{M,j} \gets \overline{K_{M,j}}\exp(\sigma_{\ln} Z_{K,j})$
            \Comment{independent per reaction}
          \State $f_j \gets [S_j]/(K_{M,j}+[S_j])$
        \Else
          \State $f_j \gets 1$
        \EndIf
        \State $V_{\max,j} \gets k_{\mathrm{cat},j}\,e$;
          \quad $u_j \gets V_{\max,j} f_j$,
          \quad $\ell_j \gets -\delta_j V_{\max,j} f_j$
      \EndFor
      \State $\hat{\mathbf v}^{(i)} \gets
        \arg\max\{\mathbf c^{\top}\mathbf v : \mathbf S\mathbf v = \mathbf 0,\
        \boldsymbol\ell \le \mathbf v \le \mathbf u\}$
    \EndFor
    \State \Return per-reaction mean, standard deviation, and $2.5/50/97.5$
      percentiles of $\{\hat{\mathbf v}^{(i)}\}_{i=1}^{N}$
  \end{algorithmic}
```

- [ ] **Step 2: Add the terminology clarification where Algorithm 1 is first introduced**

Replace:

```latex
and re-solving the linear program over ten thousand
parameter draws by the parametric bootstrap of Algorithm~\ref{alg:uq}, spreads
urea export over a right-skewed band with median $105$ and mean $156$
mmol\,gDW$^{-1}$\,h$^{-1}$ and a central ninety-five percent interval of
$17$ to $614$.
```

with:

```latex
and re-solving the linear program over ten thousand
parameter draws by the procedure of Algorithm~\ref{alg:uq}, a Monte Carlo
propagation of assumed parameter distributions rather than a bootstrap
resampled from replicate measurements, spreads
urea export over a right-skewed band with median $105$ and mean $156$
mmol\,gDW$^{-1}$\,h$^{-1}$ and a central ninety-five percent interval of
$17$ to $614$.
```

- [ ] **Step 3: Rebuild and commit**

Run: `cp code/figs/*.pdf chapter/figures/ 2>/dev/null; cd chapter && make && cd ..`
Expected: clean build.

```bash
git add chapter/sections/example_urea.tex
git commit -m "Fix Algorithm 1 pseudocode (distinct draws) and clarify UQ terminology"
```

---

## Task 6: Soften the NOS-branch, v2-uncertainty, and figure-caption overclaims

**Files:**
- Modify: `chapter/sections/example_urea.tex`

**Interfaces:** none.

- [ ] **Step 1: Soften the "v2 alone / NOS stays silent" claim**

Replace:

```latex
so the sampled bound on $v_2$ alone carries the uncertainty and the
nitric oxide synthase branch stays silent across the ensemble.
```

with:

```latex
so the sampled bound on $v_2$ dominates the uncertainty in most draws,
and the nitric oxide synthase branch, though it carries a small positive
mean and a nonzero spread in the rare draws where it turns on, remains
negligible overall.
```

- [ ] **Step 2: Soften the Figure 1 caption's "identical" claim**

Replace:

```latex
Because the cycle is linear,
    every nonzero flux equals the single throughput up to sign, so the
    parametric-bootstrap uncertainty of Algorithm~\ref{alg:uq} is one quantity
    common to them all; it is drawn as the $2.5$--$97.5$ percentile whisker on
    the urea-export flux $b_4$, the reported objective, rather than repeated
    identically on every bar.
```

with:

```latex
Because the cycle is linear,
    every nonzero backbone flux nearly equals the single throughput up to
    sign ($v_1$/$v_2$ and $v_3$/$v_4$ differ from each other only in the
    rare draws where the nitric oxide synthase branch activates), so the
    parametric-bootstrap uncertainty of Algorithm~\ref{alg:uq} is
    essentially one quantity common to them all; it is drawn as the
    $2.5$--$97.5$ percentile whisker on
    the urea-export flux $b_4$, the reported objective, rather than repeated
    identically on every bar.
```

- [ ] **Step 3: Soften the closing paragraph's "silent competing branch"**

Replace:

```latex
a single hard
bottleneck at $v_2$ and a silent competing branch at $v_5$, both
explainable entirely from the reversibility and capacity assigned
above.
```

with:

```latex
a single hard
bottleneck at $v_2$ and a competing branch at $v_5$ that remains
negligible in typical draws, both
explainable entirely from the reversibility and capacity assigned
above.
```

- [ ] **Step 4: Rebuild and commit**

Run: `cp code/figs/*.pdf chapter/figures/ 2>/dev/null; cd chapter && make && cd ..`
Expected: clean build.

```bash
git add chapter/sections/example_urea.tex
git commit -m "Soften overstated NOS-branch and per-flux-uncertainty claims"
```

---

## Task 7: Write up the FVA uniqueness result

**Files:**
- Modify: `chapter/sections/example_urea.tex`

**Interfaces:**
- Consumes: `fva()` and its test from Phase 1 (already committed, already confirms the nominal
  optimum's uniqueness to < 1e-4).

- [ ] **Step 1: Extend the nominal-solution paragraph**

Replace:

```latex
and every exchange flux carrying carbon or
nitrogen into or out of the cycle, carbamoyl phosphate and aspartate
on the uptake side, urea and fumarate on the secretion side, scaled to
match it, giving the full optimal flux distribution
(right panel, Fig.~\ref{fig:urea}).
```

with:

```latex
and every exchange flux carrying carbon or
nitrogen into or out of the cycle, carbamoyl phosphate and aspartate
on the uptake side, urea and fumarate on the secretion side, scaled to
match it, giving the full optimal flux distribution
(right panel, Fig.~\ref{fig:urea}). An independent flux-variability
analysis confirmed this optimum is unique rather than one vertex among
several tied solutions: minimizing and then maximizing each of the
nineteen fluxes in turn, with the objective held at its optimal value,
gave an identical minimum and maximum for every reaction to numerical
tolerance, so ``a single, sharply determined answer'' is a checked
property of this model rather than an assumption.
```

- [ ] **Step 2: Rebuild and commit**

Run: `cp code/figs/*.pdf chapter/figures/ 2>/dev/null; cd chapter && make && cd ..`
Expected: clean build.

```bash
git add chapter/sections/example_urea.tex
git commit -m "Write up the FVA uniqueness result in the urea example"
```

---

## Task 8: Relabel the urea example's Park data as cross-organism

**Files:**
- Modify: `chapter/sections/example_urea.tex`

**Interfaces:** none.

- [ ] **Step 1: Soften the opening sentence**

Replace:

```latex
The general bound of Equation~\eqref{eq:general-bound} carries four
gateways, and a genome-scale reconstruction may in principle open all
four at once. The urea cycle in HL-60 cells is small enough to make the
opposite case concretely: two gateways, thermodynamic and kinetic,
filled entirely from public databases, already pin down a sharp
prediction without a single measurement taken on the cells themselves.
```

with:

```latex
The general bound of Equation~\eqref{eq:general-bound} carries four
gateways, and a genome-scale reconstruction may in principle open all
four at once. The urea cycle in HL-60 cells, parameterized here with
cross-organism public data rather than HL-60-specific measurements, is
small enough to make the
opposite case concretely: two gateways, thermodynamic and kinetic,
filled entirely from public databases, already pin down a sharp
prediction without a single measurement taken on the cells themselves.
```

- [ ] **Step 2: Rebuild and commit**

Run: `cp code/figs/*.pdf chapter/figures/ 2>/dev/null; cd chapter && make && cd ..`
Expected: clean build.

```bash
git add chapter/sections/example_urea.tex
git commit -m "Relabel Park saturation data as cross-organism, not HL-60-specific"
```

---

## Task 9: Fix the overstated genome-scale-dimension and FVA-interpretation claims

**Files:**
- Modify: `chapter/sections/linearprogram.tex`

**Interfaces:** none.

- [ ] **Step 1: Soften the "order of magnitude" claim**

Replace:

```latex
A genome-scale reconstruction lists every reaction
annotated to a genome, and the reaction count $n$ exceeds the metabolite
count $m$ by a wide margin, often an order of magnitude or more, because
a single metabolite typically participates in several reactions while a
single reaction typically touches only a handful of metabolites.
```

with:

```latex
A genome-scale reconstruction lists every reaction
annotated to a genome, and the reaction count $n$ typically exceeds the
metabolite count $m$, because
a single metabolite typically participates in several reactions while a
single reaction typically touches only a handful of metabolites.
```

- [ ] **Step 2: Fix the FVA-interpretation sentence**

Replace:

```latex
A reaction
whose minimum and maximum coincide is pinned by the network topology and
bounds regardless of objective choice, while a reaction with a wide range
is left essentially free, and the pattern of which fluxes fall into which
category is exactly the information that a single optimal flux vector,
taken alone, cannot report.
```

with:

```latex
A reaction
whose minimum and maximum coincide is pinned by the network topology and
bounds at the optimal objective value used for that calculation, while a
reaction with a wide range
is left essentially free; a different objective can produce a different
feasible face and a different range for the same reaction, and the
pattern of which fluxes fall into which
category is exactly the information that a single optimal flux vector,
taken alone, cannot report.
```

- [ ] **Step 3: Rebuild and commit**

Run: `cp code/figs/*.pdf chapter/figures/ 2>/dev/null; cd chapter && make && cd ..`
Expected: clean build.

```bash
git add chapter/sections/linearprogram.tex
git commit -m "Soften overstated genome-scale dimension and FVA-interpretation claims"
```

---

## Task 10: Fix the "stoichiometry knowable in full" claims and qualify estimation/design

**Files:**
- Modify: `chapter/sections/introduction.tex`

**Interfaces:**
- Consumes: `Burgard2003` bib key from Task 1.

- [ ] **Step 1: Soften the "knowable in full" claim**

Replace:

```latex
A genome sequence and the
stoichiometry of its annotated reactions are knowable in full, but the
internal fluxes, the kinetic parameters that would set them, and the
metabolite concentrations that drive them are, at genome scale, largely
unmeasured.
```

with:

```latex
A genome sequence and the
stoichiometry of its annotated reactions are knowable in broad outline,
though real reconstructions carry their own uncertainty in reaction
membership, directionality, compartmentalization, and biomass
composition; the internal fluxes, the kinetic parameters that would set
them, and the metabolite concentrations that drive them are, at genome
scale, far less constrained still.
```

- [ ] **Step 2: Replace the "not the stoichiometry" claim**

Replace:

```latex
The bounds, not the stoichiometry, are where biological information
enters, and the question this chapter takes up is how it enters:
```

with:

```latex
Flux bounds are a major interface through which condition-specific
thermodynamic, kinetic, expression, and regulatory information enters a
constraint-based model, and the question this chapter takes up is how
it enters:
```

- [ ] **Step 3: Fix "Two published systems" and qualify the estimation/design claim**

Replace:

```latex
Two published systems then show
the same gateways opened together at a scale no hand calculation could
reach, one estimating fluxes in a dynamic cell-free reaction
\cite{Vilkhovoy2023} and one designing a glycan-producing strain by
model-guided gene deletion \cite{Wayman2019}. Read across these cases,
estimating a flux distribution from data and designing one by
intervention are the same operation on the bounds, run in opposite
directions.
```

with:

```latex
Two reported systems then show
the same gateways opened together at a scale no hand calculation could
reach, one estimating fluxes in a dynamic cell-free reaction
\cite{Vilkhovoy2023} and one designing a glycan-producing strain by
model-guided gene deletion \cite{Wayman2019}. Read across these cases,
estimating a flux distribution from data and designing one by
intervention evaluate the bounds the same way once an intervention is
chosen; design differs in also requiring a combinatorial search, over
candidate deletions or expression changes, for which bounds to move
\cite{Burgard2003}.
```

- [ ] **Step 4: Rebuild and confirm the new citation resolves**

Run: `cp code/figs/*.pdf chapter/figures/ 2>/dev/null; cd chapter && make && cd ..`

Run: `grep -i "undefined" chapter/Chapter.log | grep "Burgard2003"`
Expected: no output.

- [ ] **Step 5: Commit**

```bash
git add chapter/sections/introduction.tex
git commit -m "Soften stoichiometry claims and qualify the estimation/design equivalence"
```

---

## Task 11: Fix capstone wording (published/reported, metabolite count, Wayman claim, toy-scale note)

**Files:**
- Modify: `chapter/sections/capstones.tex`

**Interfaces:**
- Consumes: `Burgard2003` bib key from Task 1 (for the same qualification pattern as Task 10,
  applied here to `capstones.tex`'s own closing paragraph).

- [ ] **Step 1: Add the toy-scale ($n-m=1$) note to the opening paragraph**

Replace:

```latex
The urea-cycle and feedback examples opened the gateways of
Equation~\eqref{eq:general-bound} one at a time, on networks small
enough that a hand calculation could check the result against a
known truth. Two published systems show the same gateways opened
together, at a scale no hand calculation could reach, and driven by
measurement rather than by a modeler's guess at what a bound should
be.
```

with:

```latex
The urea-cycle and feedback examples opened the gateways of
Equation~\eqref{eq:general-bound} one at a time, on networks small
enough that a hand calculation could check the result against a
known truth: one degree of freedom in each case, $n-m=1$, far short of
the large right null space that motivates the gateway framework at
genome scale (Section~\ref{sec:lp}). Two reported systems show the same gateways opened
together, at a scale no hand calculation could reach, and driven by
measurement rather than by a modeler's guess at what a bound should
be.
```

- [ ] **Step 2: Fix the exact metabolite count**

Replace:

```latex
This integration is checked against denser
data: time-resolved measurement of the
intracellular state over the course of a batch reaction, on the
order of sixty metabolites tracked through the reaction, together
with time courses of mRNA, protein abundance, and enzyme activity.
```

with:

```latex
This integration is checked against denser
data: time-resolved measurement of the
intracellular state over the course of a batch reaction, 63 metabolites tracked through the reaction, together
with time courses of mRNA, protein abundance, and enzyme activity.
```

- [ ] **Step 3: Fix the Wayman knockout-prediction claim**

Replace:

```latex
The deletions this campaign screened were sufficient on their own: the
best single knockout identified by the model raised glycan production
roughly threefold over the unmodified strain once built, an outcome the
model predicted by re-solving the same linear program under the
modified bound rather than by any change to the stoichiometry $\mathbf
S$. The bounds are, in this use, not a report of what the cell is doing
but a lever for what the cell is made to do.
```

with:

```latex
The deletions this campaign screened were sufficient on their own: the
best single knockout the model identified, chosen for its predicted
growth-coupling among the candidate deletions, raised glycan production
roughly threefold over the unmodified strain once built and measured;
the linear program re-solved under each candidate's modified bound
predicted which deletion to make, not the magnitude of the yield
increase eventually observed. The bounds are, in this use, not a report
of what the cell is doing but a lever for what the cell is made to do.
```

- [ ] **Step 4: Qualify the closing "same operation" claim and cite OptKnock**

Replace:

```latex
Read together, the two capstones are the same operation on
Equation~\eqref{eq:general-bound} run in opposite directions. In the
cell-free system, measurement narrows the gateways and the flux
distribution is estimated; in the glycan strain, the engineer narrows
a gateway and the flux distribution is designed. Nothing about the
linear program distinguishes the two uses,
since both amount to fixing the bounds $\boldsymbol\ell$ and
$\mathbf u$ and solving the same constrained optimization; only the
direction in which information moves through the gateway, from data
inward for estimation, from intention outward for design, changes.
Estimation and design are, at bottom, the same operation on the
bounds, run in two directions.
```

with:

```latex
Read together, the two capstones run the same operation on
Equation~\eqref{eq:general-bound} in opposite directions, once the
intervention itself has been chosen. In the
cell-free system, measurement narrows the gateways and the flux
distribution is estimated; in the glycan strain, the engineer narrows
a gateway and the flux distribution is designed. Once a set of bound
changes is specified, both problems solve the identical linear program
under the modified box, since fixing $\boldsymbol\ell$ and $\mathbf u$
and re-solving the same constrained optimization is all either use
requires; only the direction in which information moves through the
gateway, from data inward for estimation, from intention outward for
design, changes. Choosing which bounds to change is not itself part of
that shared operation, though: for design, searching a combinatorial
space of candidate deletions or expression changes for one that
couples the desired product to growth is a discrete optimization
problem in its own right, solved by frameworks such as bilevel
programming \cite{Burgard2003} rather than by the continuous linear
program each candidate is then checked against. Estimation and design
are, at bottom, the same evaluation of the bounds, run in two
directions; only design additionally requires searching for which
bounds to move.
```

- [ ] **Step 5: Rebuild and confirm the citation resolves**

Run: `cp code/figs/*.pdf chapter/figures/ 2>/dev/null; cd chapter && make && cd ..`

Run: `grep -i "undefined" chapter/Chapter.log | grep -E "Burgard2003|sec:lp"`
Expected: no output.

- [ ] **Step 6: Commit**

```bash
git add chapter/sections/capstones.tex
git commit -m "Fix capstone wording: preprint status, metabolite count, Wayman claim, toy-scale note"
```

---

## Task 12: Final full rebuild and verification

**Files:** none modified; this task only runs and verifies.

**Interfaces:**
- Consumes: every prior task in this plan.

- [ ] **Step 1: Full chapter rebuild from a clean state**

Run:
```bash
cp code/figs/*.pdf chapter/figures/
cd chapter && make clean && make && cd ..
```
Expected: `chapter/Chapter.pdf` produced, no fatal `pdflatex` errors.

Run: `grep -i "undefined" chapter/Chapter.log`
Expected: no output referencing any label or citation touched by this plan (`Sanchez2017`,
`OBrien2013`, `Henry2007`, `Covert2001`, `Chandrasekaran2010`, `Burgard2003`, `sec:lp`).
Pre-existing unrelated undefined references, if any, are out of scope.

- [ ] **Step 2: Confirm no stale claims survive**

Run:
```bash
grep -rn "an order of magnitude or more\|regardless of objective choice\|knowable in full\|not the stoichiometry, are where\|stays silent across the ensemble\|silent competing branch\|Two published systems\|on the order of sixty" chapter/sections/
```
Expected: no output — every one of these phrases was rewritten by this plan's tasks.

- [ ] **Step 3: Confirm the new citations and content are present**

Run: `grep -c "Sanchez2017\|OBrien2013\|Henry2007\|Covert2001\|Chandrasekaran2010\|Burgard2003" chapter/sections/gateways.tex chapter/sections/introduction.tex chapter/sections/capstones.tex`
Expected: nonzero counts distributed across the three files (5 in `gateways.tex`, `Burgard2003`
in both `introduction.tex` and `capstones.tex`).

Run: `grep -c "63 metabolites" chapter/sections/capstones.tex`
Expected: `1`.

- [ ] **Step 4: Run the one Julia script this plan touched, one more time, to confirm determinism**

Run: `julia --project=code code/fba/urea_cycle_uq.jl`
Expected: identical printed output to Task 4 Step 2 (deterministic given the fixed seeds).

Run: `git status --short`
Expected: no changes (the script's outputs should already match what Task 4 committed). If
anything shows as modified, investigate before proceeding — it would mean the script is not
fully deterministic across runs, which was not the case for any Phase 1 script.

- [ ] **Step 5: Commit anything Step 1's figure copy produced (if any)**

Run: `git status --short`

If `chapter/Chapter.pdf` or any `chapter/figures/*.pdf` shows as modified from the Step 1
rebuild, commit it:

```bash
git add chapter/Chapter.pdf chapter/figures
git commit -m "Rebuild chapter PDF with all Phase 2 corrections"
```

If nothing changed, skip this step.
