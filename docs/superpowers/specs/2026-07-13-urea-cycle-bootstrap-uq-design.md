# Bootstrap Uncertainty Quantification for the Urea-Cycle Example

**Date:** 2026-07-13
**Section:** `chapter/sections/example_urea.tex` (Sec.~\ref{sec:urea})
**Code:** `code/fba/`

## Goal

Extend the urea-cycle worked example with an uncertainty quantification (UQ)
layer: treat the parameters that set the flux bounds as uncertain, propagate
that uncertainty through the linear program by parametric Monte Carlo, and
report a mean and spread for every estimated flux under the objective of
maximizing urea export. Produce (a) a runnable Julia UQ script, (b) an updated
figure with uncertainty bars, (c) a pseudo-code algorithm block for the text,
and (d) a prose paragraph integrating the result.

## Method: parametric Monte Carlo ("parametric bootstrap")

There is no replicate dataset to resample here, so the UQ is a *parametric*
bootstrap: each bound-setting parameter is assigned a probability distribution
centered (median) on its nominal chapter value, `N = 10_000` parameter sets are
drawn, the bounds are rebuilt from each draw, the same linear program of
Eq.~\eqref{eq:fba-lp} is re-solved, and the resulting flux vectors form an
ensemble. A fixed RNG seed makes the run reproducible.

### Parameter distributions (moderate, factor ~2 spread)

| Parameter | Distribution | Spread |
|---|---|---|
| `kcat[j]`, j=1..5 (enzyme turnover) | Lognormal, median = nominal | `sigma_ln = 0.69` (geometric SD ~2), independent per enzyme |
| `e0` (reference enzyme abundance) | Lognormal, median = 0.01 | `sigma_ln = 0.69`, **one shared draw per sample**, applied to all five capacities |
| `dG[j]`, j=1..5 (Gibbs free energy) | Normal, mean = nominal | `sigma = 2.0` kJ/mol, independent per reaction |

Nominal values (unchanged from the existing model):
`kcat = [10.0, 3.28, 190.0, 410.0, 10.0]` (1/s), `e0 = 0.01` (mmol/gDW),
`dG = [-4.3, -5.5, -51.0, -30.3, -1220.0]` (kJ/mol).

The lognormal median equals the nominal value (`mu = ln(nominal)`), so the
draws are multiplicatively unbiased about the reported figure; the ensemble
mean sits slightly above nominal because of lognormal skew, which is reported
honestly rather than hidden.

### Per-draw bound construction (identical to the nominal recipe)

For each enzymatic reaction `j in 1..5`:
- `delta[j] = (dG[j] > dG_threshold) ? 1 : 0` with `dG_threshold = -10.0` kJ/mol
- `Vmax[j] = kcat[j] * e0`
- `lb[j] = -delta[j] * Vmax[j]`, `ub[j] = Vmax[j]`

Exchange reactions `b1..b14` keep `lb = -1000`, `ub = +1000`. Objective is
unchanged: `c[b4] = -1` (maximize urea export). The LP is always feasible
(`v = 0` is feasible) and bounded (capacities cap the backbone), so no draw is
expected to fail; any infeasible/unbounded draw is counted and skipped.

## What the result is expected to show

1. **The thermodynamic gateway's noise is inert for this objective.** At the
   max-urea optimum the backbone runs forward, so the lower bounds never bind;
   `dG` noise flipping `v1`/`v2` to irreversible (probability ~1% and ~0.2% per
   draw at `sigma = 2`) does not move the optimum. All flux uncertainty flows
   through the capacity `Vmax = kcat * e0`.
2. **Urea export tracks the minimum cycle capacity, dominated by `v2`.** Because
   `v1..v4` carry equal flux around the cycle, urea export equals the minimum of
   the sampled cycle capacities; `v2` (nominal 0.0328, ~3x below its neighbors)
   is the bottleneck in the large majority of draws. The urea-export
   distribution is therefore close to a minimum-of-lognormals, right-skewed,
   with CV ~120-130% (two lognormal factors, `kcat[2]` and `e0`, multiplying).
3. **The NOS branch `v5` stays robustly at zero** across draws, since diverting
   arginine away from urea can only lower the objective.

These are expectations to confirm from the ensemble, not asserted results.

## Statistics reported

Per reaction: nominal flux, ensemble mean, ensemble standard deviation (SD),
coefficient of variation, and the 2.5 / 50 / 97.5 percentiles. Terminology note
for the text: the ensemble SD *is* the bootstrap standard error of the flux
estimate (the SD of the bootstrap distribution of that flux); the Monte-Carlo
error of the mean, `SD / sqrt(N)`, is a separate quantity that is negligible at
`N = 10_000` and is not what the error bars represent.

Output: `code/data/urea_fba_uq.csv` with columns
`reaction, nominal, mean, std, cv, q025, q50, q975`.

## Components

### 1. `code/fba/urea_cycle.jl` (refactor, behavior-preserving)

Change the signature to
`urea_cycle_model(; kcat=KCAT0, e0=E0, dG=DG0, dG_threshold=-10.0)` where the
default constants are the nominal values above. The bounds `lb`/`ub` are
derived from `(kcat, e0, dG, dG_threshold)` instead of being hardcoded. The
zero-argument call `urea_cycle_model()` must reproduce the current model
exactly: `lb = [-0.1, -0.0328, 0, 0, 0, -1000...]`,
`ub = [0.1, 0.0328, 1.9, 4.1, 0.1, 1000...]`. Verified by re-solving and
checking the optimum is still `0.0328` and `Sv = 0`. `S`, `reactions`,
`metabolites`, `c` are unchanged.

Also move `solve_fba(m)` out of `run_fba.jl` and into `urea_cycle.jl` (with its
own docstring) so both `run_fba.jl` and the new UQ script reuse one solver and
neither duplicates it. `run_fba.jl` keeps only its top-level driver code
(build model, `solve_fba`, write CSV, plot); its behavior is unchanged.

### 2. `code/fba/urea_cycle_uq.jl` (new)

- `include` `Include.jl` and `urea_cycle.jl`; reuse the `solve_fba(m)` now
  living in `urea_cycle.jl` (does not `include` `run_fba.jl`, so no nominal-run
  file I/O is triggered).
- Set seed, `N`, nominal params, `sigma_ln`, `dG_sigma`.
- Loop `N` times: sample params, `m = urea_cycle_model(; kcat, e0, dG)`, solve,
  store the flux row. Track feasibility.
- Reduce the ensemble to per-reaction mean/SD/CV/percentiles; write the CSV.
- Regenerate `code/figs/urea_fba.pdf`: same bar chart of nominal flux per
  reaction, now with asymmetric error bars spanning the 2.5-97.5 percentile band
  (bar height = nominal flux). Reuse the existing axis styling.
- Print the urea-export summary (`b4`) to stdout for the prose numbers.

### 3. Figure integration

Copy `code/figs/urea_fba.pdf` to `chapter/figures/urea_fba.pdf` (same filename;
the `\includegraphics{urea_fba.pdf}` reference and figure placement are
unchanged).

### 4. `chapter/Chapter.tex` (preamble)

Add `\usepackage{algorithm}` and `\usepackage{algpseudocode}` (standard TeX
Live). Confirm the chapter still compiles with `make`.

### 5. `chapter/sections/example_urea.tex` (prose + pseudo-code)

- Add one `algorithm` float, `\label{alg:uq}`, with the Monte-Carlo pseudo-code:
  inputs (nominal params, spreads, N, seed), the sample-rebuild-solve loop, and
  the mean/SD/percentile reduction. Referenced as Algorithm~\ref{alg:uq}.
- Add **one long paragraph** after the current result paragraph (the one ending
  at Fig.~\ref{fig:urea}): introduce the parametric-bootstrap setup and the
  factor-2 assumption, report urea export as mean +/- SD with its 95% interval,
  and give the interpretation (v2 bottleneck carries the uncertainty, the
  thermodynamic gateway's noise is inert because the reversibility bounds do not
  bind, v5 stays at zero). Long-form prose, no em dashes, no subsection heading,
  no section self-reference, consistent with house style.
- Update the `\caption` of Fig.~\ref{fig:urea} to describe the percentile error
  bars.
- No results table (the chapter is figure-driven; full stats live in the CSV).

## Style / conventions to honor

- No em dashes; no subsection headings in the example; no "this section..."
  self-references; long complete paragraphs (merge, don't fragment).
- Numeric results in prose must tie to the figure (the error bars) or be
  structural setup numbers (N, sigma, threshold), which are exempt.
- Reproduce the nominal solve exactly before layering UQ on top.

## Out of scope

- Varying the saturation factor `f_j`, expression `e/e0` beyond the reference
  draw, or regulation `theta_j` (the other two gateways stay at unity).
- Any measurement on HL-60 cells; all uncertainty is database-derived and
  assumed.
- Changes to the feedback example or other sections.

## Verification

- `urea_cycle_model()` zero-arg output byte-identical to current `lb`/`ub`;
  nominal solve returns `0.0328` and passes `max|Sv| < 1e-6`.
- UQ script runs to completion, writes the CSV with all 19 reactions, and
  regenerates the figure PDF.
- `make` in `chapter/` compiles without errors with the new packages and prose.
- Reported prose numbers match the CSV / stdout summary.
