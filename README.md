# Mathematical Models in Biotechnology — Chapter Repository

This repository contains the manuscript and reproducible code for the chapter
"Mathematical Models in Biotechnology" by Jeffrey D. Varner, to appear in
*Comprehensive Biotechnology*, 4th edition (Elsevier/MRW-BTC4). The chapter
surveys the spectrum of mathematical modeling in biotechnology — from mechanistic
ODE systems grounded in conservation laws and reaction kinetics, through
flux-balance analysis and gene-regulatory circuit motifs, to pure data-driven
deep sequence architectures (LSTM, S4/HiPPO) and hybrid mechanistic+ML models —
using a CHO fed-batch monoclonal-antibody bioreactor as a unified worked example
throughout all six sections.

---

## Repository layout

```
chapter/                  LaTeX manuscript
  Chapter.tex             Root document (title, abstract, \input{sections/*})
  sections/               One .tex file per section
    introduction.tex
    kinetics.tex
    fba.tex
    geneexpression.tex
    deeplearning.tex
    hybrid.tex
    appendix.tex
  figures/                PDF figures (copied from code/figs/ before building)
  References.bib          BibTeX database
  Makefile                pdflatex + bibtex build recipe

code/                     Julia project
  Project.toml            Package manifest (pinned environment)
  Manifest.toml           Resolved dependency tree
  Include.jl              Shared entry point (activates project, loads src/)
  src/
    Runtime.jl            Shared utilities (data loading, plotting helpers)
  kinetics/
    cho_model.jl          CHO fed-batch ODE model
    run_cho.jl            -> figs/cho_kinetics.pdf
  fba/
    urea_cycle.jl         Urea-cycle stoichiometry (JuMP/HiGHS)
    run_fba.jl            -> figs/urea_fba.pdf
  geneexpression/
    motifs.jl             NAR, C1-FFL, I1-FFL, negative-feedback ODE circuits
    run_motifs.jl         -> figs/motifs.pdf
    cybernetic.jl         Cybernetic diauxic-growth model
    run_cybernetic.jl     -> figs/cybernetic_diauxie.pdf
  deeplearning/
    data_prep.jl          Sliding-window data preparation for CHO trajectories
    lstm.jl               LSTM implementation (Flux)
    s4.jl                 S4/HiPPO-LegS implementation (Flux)
    run_lstm.jl           -> figs/lstm_cho.pdf
    run_s4.jl             -> figs/s4_cho.pdf
    run_comparison.jl     -> figs/s4_vs_lstm.pdf
    run_hybrid.jl         -> figs/hybrid_cho.pdf
  data/                   CSV data files (CHO trajectories, pre-computed metrics)
  figs/                   Generated PDF figures
```

---

## Reproducing the figures

Install the Julia environment once (Julia 1.12+ required):

```bash
julia --project=code -e 'using Pkg; Pkg.instantiate()'
```

Then run each script from the repo root. Scripts are listed below grouped by
chapter section; each exits 0 and writes a single PDF to `code/figs/`.

### Section 2 — Kinetics

```bash
julia --project=code code/kinetics/run_cho.jl        # -> code/figs/cho_kinetics.pdf
```

### Section 3 — Flux-Balance Analysis

```bash
julia --project=code code/fba/run_fba.jl             # -> code/figs/urea_fba.pdf
```

### Section 4 — Gene Expression

```bash
julia --project=code code/geneexpression/run_motifs.jl     # -> code/figs/motifs.pdf
julia --project=code code/geneexpression/run_cybernetic.jl # -> code/figs/cybernetic_diauxie.pdf
```

### Sections 5 & 6 — Deep Learning and Hybrid Models

The four scripts below train neural networks (fixed seeds; ~2–5 min each on a
laptop CPU):

```bash
julia --project=code code/deeplearning/run_lstm.jl       # -> code/figs/lstm_cho.pdf
julia --project=code code/deeplearning/run_s4.jl         # -> code/figs/s4_cho.pdf
julia --project=code code/deeplearning/run_comparison.jl # -> code/figs/s4_vs_lstm.pdf
julia --project=code code/deeplearning/run_hybrid.jl     # -> code/figs/hybrid_cho.pdf
```

---

## Building the chapter PDF

Requires a LaTeX distribution (TeX Live 2023+ recommended):

```bash
cp code/figs/*.pdf chapter/figures/
cd chapter && make
```

`make` runs `pdflatex` + `bibtex` + two additional `pdflatex` passes and
produces `chapter/Chapter.pdf`.

---

## Requirements

| Requirement | Version |
|---|---|
| Julia | 1.12+ |
| LaTeX | TeX Live (any recent distribution) |

Julia package versions are pinned in `code/Project.toml` and
`code/Manifest.toml`. Running `Pkg.instantiate()` (see above) installs exactly
those versions.

---

## Design specification

The full chapter design specification (section outlines, figure descriptions,
notation conventions, and implementation notes) is at:

```
docs/superpowers/specs/2026-06-23-mathematical-models-biotechnology-chapter-design.md
```
