# Current Chapter Audit

**Chapter:** *Mathematical Models in Biotechnology*

**Audit date:** 2026-07-21

**Branch:** chapter-draft

## Current status

The chapter-wide language and paragraph-structure revision is complete. The
prose is simple, concise, and direct, and related sentences are grouped into
sustained paragraphs rather than short blocks. The balance derivation,
modeling terminology, examples, captions, and outlook use consistent language.
The chapter compiles with no broken references, overfull boxes, or visible
layout defects.

The language work does not close every technical and reproducibility issue.
The remaining items are listed below.

## Resolved

### Balance and notation

- Reactor volume is written as $\bar V$, not $V$.
- The biomass basis is $\mathcal B=B\bar V$, where $B$ is biomass
  concentration in gDW/L, $\bar V$ is reactor volume in L, and $\mathcal B$
  is total biomass in gDW.
- Intracellular concentration $C_i$ has units mmol/gDW, and intracellular
  flux $\hat v_j$ has units mmol/gDW/h.
- A distinct stream concentration $c_{i,s}$ is used for the hypothetical
  physical-flow term. Active transport remains part of
  $\mathbf S\hat{\mathbf v}$.
- The derivation expands $\mathcal B=B\bar V$ before imposing a fixed-volume
  assumption. The $d\bar V/dt$ term therefore remains visible for fed-batch
  operation.
- The conventional FBA approximation
  $\mathbf S\hat{\mathbf v}=\mathbf0$ is distinguished from the
  growth-dilution form $\mathbf S\hat{\mathbf v}=\mu\mathbf x$.

### Model scope and terminology

- The introduction now distinguishes unstructured extracellular models,
  structured kinetic models, reaction-level kinetic models, and
  constraint-based models.
- Shuler and coworkers' structured cell models are included.
- Bounds are described as a major route for condition-specific information,
  not as the only place where biology enters a model.
- The phrase “retained-cell culture” and other unclear reactor terminology
  have been removed.

### Linear program

- The LP section now states directly that the balance equations admit many
  flux distributions and that bounds and an objective narrow this set.
- The growth-aware formulation is qualified: it remains linear when
  $\mathbf x$ is measured and $\mu$ is fixed or enters linearly. If both are
  unknown, $\mu\mathbf x$ is bilinear.
- The chapter notes that explicit dilution demands and a biomass
  pseudo-reaction must be defined consistently to avoid double counting.
- Flux variability is interpreted on the selected optimal face.

### Urea-cycle example

- The example is identified as a cross-organism public-data parameterization,
  not an HL-60-specific measurement set.
- The seconds-to-hours factor of 3600 is explicit in both the calculation and
  Algorithm 1.
- Uncertainty is consistently described as Monte Carlo propagation, not a
  bootstrap.
- The nominal optimum, uncertainty summaries, and flux-variability result are
  stated directly.
- Metaphors such as “content poured into” and “left over to spend” have been
  removed.

### Feedback example

- The kinetic simulation is described as a BST reference, not independent
  experimental validation or “truth.”
- The FBA calculation is described as a cross-model comparison that receives
  steady-state $X_3$ and enzyme abundance from the kinetic simulation.
- The BST power law is written with the dimensionless ratio
  $X_3/X_3^\circ$.
- Sequence lengths and elongation rates are described as illustrative
  expression timescales, not as values that determine the calibrated steady
  state.
- Metaphors such as “load-bearing,” “capacity ledger,” “landing at,” and
  “recovery” have been removed from the chapter.

### Chapter-wide language

- The abstract, introduction, derivation, LP section, integrative-bounds
  section, both examples, capstones, captions, and outlook were revised.
- Long sentences were split where needed, but related ideas remain grouped in
  full paragraphs.
- Short setup, explanation, and interpretation blocks were merged across all
  sections. Most prose paragraphs now contain roughly 100--300 words, depending
  on their equations and figures.
- Promotional or absolute wording such as “fully known,” “the same operation
  in opposite directions,” and “independent recovery” has been removed.
- The sentence “Linear programming addresses the underdetermination of the
  constraint set” has been replaced with a direct explanation.
- Every figure reference is integrated into the opening sentence of the
  paragraph that discusses it; standalone constructions such as “Figure 1
  shows these steps” have been removed. The chapter currently has no tables.

### Code reproducibility and terminology

- The FBA solver now exposes the optimization termination status while
  preserving the existing flux-only interface.
- The urea Monte Carlo calculation writes a complete draw log with draw number,
  solver status, active bounds, and urea export. Failed draws are retained in
  the log rather than silently discarded.
- The capacity-only ensemble and the saturation comparison also record
  draw-level status and active bounds.
- The imputed $f_2$ calculation uses a separate random-number stream, so the
  other sampled inputs remain matched between Configurations A and B.
- Obsolete “parametric bootstrap,” identical-band, and “gateway” comments were
  removed from the UQ implementation.
- Feedback code and tests now use reference and control terminology instead of
  truth and gateway terminology.
- Sequence elongation times are named and reported as illustrative metadata,
  with explicit documentation that they do not calibrate the expression
  coefficients.

### Availability and source verification

- A dedicated code-and-data availability statement now links to the chapter
  repository and identifies the pinned Julia project and manifest files.
- The repository README now describes the current urea-cycle and feedback
  examples, their generated outputs, tests, and the chapter build process.
- The Wayman et al. paper confirms that the single $\Delta gnd$ deletion,
  drawn from the model-identified strain families, increased measured glycan
  production nearly threefold relative to wild type. The chapter now separates
  the model-identified perturbation from the magnitude measured experimentally.

## Remaining technical and reproducibility items

### Parameter provenance

The chapter identifies eQuilibrator, BRENDA, BioNumbers, and Park et al. as
sources, but it does not provide a complete record for each parameter. A
reproducible parameter table should include the database record or query,
organism, substrate, assay direction, temperature, pH, original units,
conversion, selection rule, and access date.

The thermodynamic calculations also need their transformed-state conditions,
including pH, ionic strength, temperature, magnesium assumptions, and
uncertainty. The reported nitric oxide synthase value of
$-1220$ kJ/mol warrants a separate stoichiometric and thermodynamic check.

## Verification

- latexmk -pdf -interaction=nonstopmode -halt-on-error Chapter.tex: passed.
- git diff --check: passed.
- PDF text scan found none of the previously flagged phrases.
- All 15 rendered pages were inspected. No clipped text, overlapping elements,
  broken figures, unreadable captions, or margin overflow was found.
- The urea and feedback Julia test suites passed after the code changes.
- The full 10,000-draw UQ calculations completed with zero failed solves in
  each ensemble, and every successful draw has at least one labeled active
  bound.

## Conclusion

The language and paragraph-structure audits are closed. The chapter now uses
direct, professional prose in sustained paragraphs throughout. The audit
remains open only for parameter provenance and the separate thermodynamic
check identified above.
