# Bootstrap Uncertainty Quantification for the Urea-Cycle Example

**Date:** 2026-07-13
**Section:** `chapter/sections/example_urea.tex` (Sec.~\ref{sec:urea})
**Code:** `code/fba/`

## Goal

Extend the urea-cycle worked example with an uncertainty quantification (UQ)
layer: treat the parameters that set the flux bounds as uncertain, propagate
that uncertainty through the linear program by parametric Monte Carlo, and
report a mean and spread for every estimated flux under the objective of
maximizing urea export. Three gateways of Eq.~\eqref{eq:general-bound} carry
parameters here: the thermodynamic switch `delta_j` (via `dG`), the capacity
scale `Vmax_j = kcat_j * e0`, and the saturation factor `f_j` (via substrate
concentration and `Km`). The expression and regulatory gateways stay at unity.

Produce (a) a runnable Julia UQ script, (b) an updated figure with uncertainty
bars plus a second sensitivity figure, (c) a pseudo-code algorithm block for the
text, and (d) two prose paragraphs integrating the result.

## Method: parametric Monte Carlo ("parametric bootstrap")

There is no replicate dataset to resample, so the UQ is a *parametric*
bootstrap: each bound-setting parameter is assigned a probability distribution
centered (median) on its nominal value, `N = 10_000` parameter sets are drawn,
the bounds are rebuilt from each draw, the same linear program of
Eq.~\eqref{eq:fba-lp} is re-solved, and the resulting flux vectors form an
ensemble. A fixed RNG seed makes the run reproducible.

### Parameter distributions (moderate, factor ~2 spread)

| Parameter | Distribution | Spread |
|---|---|---|
| `kcat[j]`, j=1..5 (enzyme turnover) | Lognormal, median = nominal | `sigma_ln = 0.69` (geometric SD ~2), independent per enzyme |
| `e0` (reference enzyme abundance) | Lognormal, median = 0.01 | `sigma_ln = 0.69`, **one shared draw per sample**, applied to all five capacities |
| `dG[j]`, j=1..5 (Gibbs free energy) | Normal, mean = nominal | `sigma = 2.0` kJ/mol, independent per reaction |
| `conc[j]`, `Km[j]` (saturation inputs) | Lognormal, median = Park value | `sigma_ln = 0.69` each; see saturation gateway below |

Nominal kinetic/thermodynamic values (unchanged from the existing model):
`kcat = [10.0, 3.28, 190.0, 410.0, 10.0]` (1/s), `e0 = 0.01` (mmol/gDW),
`dG = [-4.3, -5.5, -51.0, -30.3, -1220.0]` (kJ/mol).

The lognormal median equals the nominal value (`mu = ln(nominal)`), so draws are
multiplicatively unbiased about the reported figure; the ensemble mean sits
slightly above nominal because of lognormal skew, reported honestly.

## Saturation gateway `f_j` (new)

The chapter defines `f_j` as the single-substrate Michaelis-Menten fraction
(`gateways.tex:99`), `f_j = [S_j] / (K_M,j + [S_j]) in [0,1]`, so each reaction
needs exactly **one** characteristic substrate concentration and **one** `Km`.

### Data source

Park et al. 2016 (*Nat. Chem. Biol.*, PMC4912430), Supplementary Data Set
(file `data/41589_2016_BFnchembio2077_MOESM585_ESM.xlsx`, "Comparison of
absolute concentrations to Km"), which pairs, per EC number, an intracellular
metabolite concentration with a BRENDA `Km` and the source organism. Extracted
rows for the five urea-cycle EC numbers (concentration and `Km` in M):

| Reaction | EC | Substrate row(s) present | conc / Km (M) | implied f | organism |
|---|---|---|---|---|---|
| v1 ASS | 6.3.4.5 | ATP; aspartate | 4.673e-3 / 3.923e-4; 1.492e-2 / 1.543e-4 | 0.923; 0.990 | Mus musculus |
| v2 ASL | 4.3.2.1 | **only products** (arginine, fumarate) | — | **none** | Homo sapiens |
| v3 arginase | 3.5.3.1 | arginine (human); arginine (yeast) | 2.555e-4 / 1.546e-3; 2.182e-2 / 1.570e-2 | 0.142; 0.582 | H. sapiens; S. cerevisiae |
| v4 OTC | 2.1.3.3 | ornithine (yeast); ornithine (E. coli) | 4.489e-3 / 1.600e-3; 1.010e-5 / 8.500e-4 | 0.737; 0.012 | S. cerevisiae; E. coli |
| v5 NOS | 1.14.13.39 | arginine | 2.555e-4 / 3.497e-6 | 0.986 | Mus musculus |

Exhaustive name search of all three supplement files confirms
**argininosuccinate, citrulline, and carbamoyl phosphate are absent** (the
"carbamoyl-aspartate" hit is a pyrimidine-pathway metabolite). So the `v2`
forward substrate (argininosuccinate) is genuinely unmeasured; the table's `v2`
rows are its reverse-direction products.

### Single-substrate rule and nominal `f_j`

Two aggregation rules, applied in order: (i) across the distinct substrates of
one reaction (e.g. `v1` has ATP and aspartate), take the **least-saturated**
substrate (smallest `f`), since that is what limits the capacity ceiling; (ii)
across multiple organism rows for the *same* substrate, prefer the *Homo
sapiens* row (this is a human/HL-60 context), and if none is human take the
geometric mean of the available rows' `conc` and `Km`. Nominal `f_j`:
`f1 = 0.923` (ATP, Mus musculus only), `f3 = 0.142` (human arginine row),
`f4 = 0.154` (ornithine, geometric mean of yeast and E. coli; no human row),
`f5 = 0.986` (arginine, Mus musculus only), `f2 = 1.0` (no substrate data;
default).

### The `f_2` gap decides everything (main result)

Multiplying each nominal `f_j` into its capacity gives effective ceilings
`0.10*0.923 = 0.092` (v1), `1.9*0.142 = 0.27` (v3), `4.1*0.154 = 0.63` (v4),
`0.10*0.986 = 0.099` (v5) -- **all above the v2 bottleneck at `0.0328`**. So
opening the saturation gateway with every substrate concentration Park provides
leaves urea export **unchanged at 0.0328**: the only reaction that binds is the
one whose substrate nobody measured. This is the quantitative form of the
"you need data on the binding constraint" lesson.

### Two Monte Carlo configurations

- **Config A (main):** sample `f1, f3, f4, f5` as capacity multipliers -- for
  each draw, `conc[j] ~ Lognormal(median = nominal conc, sigma_ln = 0.69)` and
  `Km[j] ~ Lognormal(median = nominal Km, sigma_ln = 0.69)`, then
  `f_j = conc/(Km+conc)`. Keep `f2 = 1`. Expected: the four data-informed
  factors are inert (v2 remains the bottleneck in essentially every draw), and
  the urea-export distribution equals that of `kcat[2]*e0`.
- **Config B (f_2 sensitivity):** additionally impute `f2`. Draw argininosuccinate
  concentration `[ASA] ~ Lognormal(median = 5e-6 M, sigma_ln = 1.0)` (spans
  ~1-100 uM; normal intracellular is low-uM, argininosuccinic aciduria is
  diagnosed at plasma 5-110 uM per HMDB0000052 / GeneReviews) and
  `Km_ASL ~ Lognormal(median = 5e-6 M, sigma_ln = 0.7)` (a few uM to a few tens
  of uM, BRENDA EC 4.3.2.1; exact central value read from BRENDA at
  implementation). `f2 = [ASA]/(Km_ASL+[ASA])` then multiplies `Vmax[2]`, pulling
  the bottleneck below 0.0328 and adding an argininosuccinate-dominated band.
  This is an explicitly illustrative sweep, labeled as such; it commits no
  argininosuccinate point value to the model.

## Per-draw bound construction

For each enzymatic reaction `j in 1..5`:
- `delta[j] = (dG[j] > dG_threshold) ? 1 : 0`, `dG_threshold = -10.0` kJ/mol
- `Vmax[j] = kcat[j] * e0`
- `cap[j] = Vmax[j] * f[j]` (f from the config; `f2 = 1` in A, sampled in B)
- `lb[j] = -delta[j] * cap[j]`, `ub[j] = cap[j]`

Exchange reactions `b1..b14` keep `lb = -1000`, `ub = +1000`. Objective
unchanged: `c[b4] = -1` (maximize urea export). The LP is always feasible
(`v = 0`) and bounded (capacities cap the backbone); any infeasible/unbounded
draw is counted and skipped.

## Expected results to confirm

1. **Thermodynamic noise is inert.** At the optimum the backbone runs forward,
   the lower bounds never bind, and `dG` flips of `v1`/`v2` to irreversible
   (~1% and ~0.2% of draws) do not move the optimum.
2. **The four data-informed `f_j` are inert (Config A).** Urea export tracks
   `min` of the cycle capacities, and v2 stays the minimum in essentially every
   draw, so urea export equals `kcat[2]*e0`: median ~0.0328, right-skewed,
   CV ~120-130% (two lognormal factors multiplying).
3. **`v5` stays robustly at zero** (diverting arginine can only lower urea).
4. **`f_2` is the sole mover (Config B).** Imputing argininosuccinate drags the
   urea median below 0.0328 and widens the band; the width is governed by the
   argininosuccinate/`Km` uncertainty, i.e. by the one missing measurement.

These are expectations to confirm from the ensemble, not asserted results.

## Statistics and outputs

Per reaction: nominal flux, ensemble mean, ensemble SD, CV, and 2.5/50/97.5
percentiles. Terminology for the text: the ensemble SD *is* the bootstrap
standard error of the flux (SD of the bootstrap distribution); the Monte-Carlo
error of the mean, `SD/sqrt(N)`, is separate and negligible at `N = 10_000`.

Outputs (in `code/data/`):
- `park_saturation.csv` -- extracted `[S]/Km` rows for the five ECs, with a
  provenance header (source file, sheet, columns, citation). ~9 rows.
- `urea_fba_uq.csv` -- Config A per-reaction stats
  (`reaction, nominal, mean, std, cv, q025, q50, q975`).
- `urea_uq_sensitivity.csv` -- Config B urea-export (`b4`) summary across a grid
  or ensemble of `f2` draws, plus the marginal `f2` distribution.

## Components

### 1. `code/data/park_saturation.csv` (new, committed)

Small derived CSV of the five-EC substrate rows above, extracted from the
supplement with a documented one-off step (Python + openpyxl; the raw supplement
stays out of git). Header comment records provenance and the Park et al.
citation. Columns: `reaction, ec, substrate, conc_M, km_M, organism`.

### 2. `code/fba/urea_cycle.jl` (refactor, behavior-preserving)

Signature `urea_cycle_model(; kcat=KCAT0, e0=E0, dG=DG0, dG_threshold=-10.0,
f=ones(5))`, default constants = nominal values, `f` the five saturation
multipliers (default all ones). Bounds derive from `(kcat, e0, dG,
dG_threshold, f)`. The zero-argument call `urea_cycle_model()` must reproduce
the current model exactly: `lb = [-0.1, -0.0328, 0, 0, 0, -1000...]`,
`ub = [0.1, 0.0328, 1.9, 4.1, 0.1, 1000...]`, optimum still `0.0328`,
`max|Sv| < 1e-6`. `S`, `reactions`, `metabolites`, `c` unchanged.

Also move `solve_fba(m)` out of `run_fba.jl` and into `urea_cycle.jl` (own
docstring) so both `run_fba.jl` and the UQ script reuse one solver. `run_fba.jl`
keeps only its top-level driver code; behavior unchanged.

### 3. `code/fba/urea_cycle_uq.jl` (new)

- `include` `Include.jl` and `urea_cycle.jl` (not `run_fba.jl`, so no nominal
  file I/O). Load `park_saturation.csv` and reduce to per-reaction nominal
  `(conc, Km)` (geometric mean across organism rows; least-saturated substrate).
- Set seed, `N`, nominal params, spreads.
- **Config A loop:** sample `kcat`, `e0`, `dG`, and `f1,f3,f4,f5` (from the
  Park nominal `conc`/`Km` with lognormal jitter), `f2 = 1`; build model, solve,
  store the flux row. Reduce to stats; write `urea_fba_uq.csv`.
- **Config B loop:** as Config A but also sample `f2` from imputed
  argininosuccinate `[ASA]`/`Km_ASL`; record urea export (`b4`) and `f2`; write
  `urea_uq_sensitivity.csv`.
- Regenerate `code/figs/urea_fba.pdf`: the flux bar chart with asymmetric error
  bars spanning the 2.5-97.5 percentile band (bar = nominal flux), reusing the
  existing axis styling (Config A).
- New `code/figs/urea_saturation.pdf`: the Config-A vs Config-B urea-export
  distributions (overlaid densities/histograms), showing the leftward shift and
  widening when `f2` is imputed.
- Print urea-export summaries for both configs to stdout for the prose numbers.

### 4. Figure integration

Copy `code/figs/urea_fba.pdf` and `code/figs/urea_saturation.pdf` into
`chapter/figures/`. `urea_fba.pdf` keeps its filename (existing
`\includegraphics` reference and placement unchanged); `urea_saturation.pdf` is
a new figure.

### 5. `chapter/Chapter.tex` (preamble)

Add `\usepackage{algorithm}` and `\usepackage{algpseudocode}`. Confirm `make`
still compiles.

### 6. `chapter/sections/example_urea.tex` (prose + pseudo-code + figure)

- Add one `algorithm` float, `\label{alg:uq}`, with the Monte-Carlo pseudo-code:
  inputs (nominal params, spreads, N, seed), the sample-rebuild-solve loop
  including the `f_j` construction, and the mean/SD/percentile reduction.
- Add **two long paragraphs** after the current result paragraph. First
  (Config A): the parametric-bootstrap setup, the factor-2 assumption, urea
  export as mean +/- SD with 95% interval (~unchanged 0.0328), and the finding
  that the four Park-informed `f_j` are inert while the thermodynamic noise is
  inert too, so the whole band flows through `kcat[2]*e0` at the v2 bottleneck.
  Second (Config B): the argininosuccinate gap, the imputed sweep, the resulting
  reduced-and-widened urea band, and the binding-constraint lesson.
- Add the new figure `urea_saturation.pdf` with `\label{fig:urea_uq}` and a
  caption describing the two configurations.
- Update the `\caption` of Fig.~\ref{fig:urea} to describe the percentile bars.
- No results table (chapter is figure-driven; full stats live in the CSVs).

## Data provenance and git

- Add `data/` (the raw Park supplement: 1 PDF + 3 xlsx, third-party) to
  `.gitignore`; it is not committed.
- Commit only the small derived `code/data/park_saturation.csv` with a
  provenance/citation header.

## Style / conventions to honor

- No em dashes; no subsection headings in the example; no "this section..."
  self-references; long complete paragraphs (merge, don't fragment).
- Numeric results in prose tie to a figure (the error bars / the sensitivity
  figure) or are structural setup numbers (N, sigma, threshold, the Park
  `conc`/`Km` inputs), which are exempt.
- Reproduce the nominal solve exactly before layering UQ on top.
- `f_2 = 1` is an honest "no data" default in Config A; Config B is labeled an
  illustrative sensitivity sweep, not a committed prediction.

## Out of scope

- Varying the expression ratio `e/e0` (beyond the shared `e0` draw) or the
  regulatory factor `theta_j` (those gateways stay at unity).
- Any measurement on HL-60 cells; all inputs are database-derived and assumed,
  with a cross-species caveat (Park is E. coli / yeast / mouse iBMK, not human).
- Changes to the feedback example or other sections.

## Verification

- `urea_cycle_model()` zero-arg output byte-identical to current `lb`/`ub`;
  nominal solve returns `0.0328` and passes `max|Sv| < 1e-6`.
- With nominal `f` from Park applied (f2=1), the solve still returns `0.0328`
  (confirms the four factors are inert and v2 binds).
- UQ script runs to completion, writes both CSVs with all reactions, and
  regenerates both figure PDFs.
- `make` in `chapter/` compiles without errors with the new packages and prose.
- Reported prose numbers match the CSV / stdout summaries.
