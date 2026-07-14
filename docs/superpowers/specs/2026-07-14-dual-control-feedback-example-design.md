# Dual-Control Feedback: Expression and Activity on One Committed Enzyme

**Date:** 2026-07-14
**Section:** `chapter/sections/example_feedback.tex` (Sec.~\ref{sec:feedback})
**Code:** `code/feedback/`

## Goal

Rebuild the linear-feedback worked example so it exercises **two** gateways of
Eq.~\eqref{eq:general-bound} on a single committed enzyme, driven by the same end
product: the **regulatory** factor `theta_j` (fast allosteric end-product
inhibition, already present) *and* the **expression** ratio `(e/e0)` (slow
transcriptional end-product repression, new). The expression feedback is cast as
an Alon **negative-feedback loop**; its transcription and translation capacities
are computed from sequence-derived quantities with realistic-but-fake values; and
the committed-step bound becomes the product `Vmax0 * (e/e0) * theta`, showing the
two gateways multiply. This is the trp-operon dual-control architecture: one
effector, two loops, two timescales.

The current example opens only `theta` (the `(e/e0)` gateway is left at its
default of unity). This rewrite makes both load-bearing and keeps the example's
signature move: generate a self-consistent "truth" and recover it *through the
bound*, never by copying fluxes.

Produce (a) an extended BSTModelKit S-system as the truth, (b) a rewritten Julia
script running four escalating FBA cases, (c) one two-panel figure, and (d) the
rewritten prose (~5 long paragraphs).

## The system

Same linear chain, now with the committed enzyme `E0` and its transcript `m`
modeled explicitly:

```
        (transcription, repressed by X3)      (allosteric, inhibited by X3)
   {} --> m --> E0  ......................> catalyzes ...............>
                                                              r0
   {} --r0--> X1 --r1--> X2 --r2--> X3 --r3--> {}
```

`X3` closes two loops on `r0`:
- **Activity (fast):** `X3` allosterically inhibits `E0` -> `theta`.
- **Expression (slow):** `X3`, through a repressor, transcriptionally represses
  `gene(E0)` -> `(e/e0)`. Alon negative-feedback loop.

`r1`, `r2` are fast first-order steps that never limit; `r3` is a first-order
export that pins `X3 = T/k3` at steady state (`T` = throughput). Because the chain
is linear, every metabolic flux equals the throughput `T`.

## Truth model: one extended BST S-system

Extend the BSTModelKit model to five dynamic species: `X1, X2, X3, m, E0`, GMA
power-law throughout. The expression states are **nondimensional ratios**
normalized to their reference (de-repressed) level so both sit at 1 when
repression is off: `mhat = m/m_ref`, `ehat = E0/e0`, and `(e/e0) = ehat`. This
keeps the metabolic chain in the same arbitrary units as the existing example
(reference capacity `Vmax0 = 10`) while the expression cascade contributes a pure
dimensionless ratio.

| reaction | connection | kinetic order | rate |
|---|---|---|---|
| `rTX`  | `{} --> m`      | `X3^{-b}`, coeff `(theta_m+mu)(1-phi)` | repressed transcription |
| `rTXb` | `{} --> m`      | constant, coeff `(theta_m+mu) phi`     | basal leak `lambda` |
| `rMdeg`| `m --> {}`      | `m^{1}`, coeff `(theta_m+mu)`          | mRNA turnover + dilution |
| `rTL`  | `{} --> E0`     | `m^{1}`, coeff `(theta_p+mu)`          | translation (m catalytic) |
| `rEdeg`| `E0 --> {}`     | `E0^{1}`, coeff `(theta_p+mu)`         | enzyme turnover + **dilution** |
| `r0`   | `{} --> X1`     | `E0^{1} X3^{-a}`, coeff `Vmax0`        | committed step: expression x activity |
| `r1`   | `X1 --> X2`     | `X1^{1}`, coeff `100`                  | fast |
| `r2`   | `X2 --> X3`     | `X2^{1}`, coeff `100`                  | fast |
| `r3`   | `X3 --> {}`     | `X3^{1}`, coeff `k3`                   | export |

Coefficients on `rTX`/`rMdeg` are tied so the reference (de-repressed, `X3^{-b} ->
1`) steady state is `mhat = 1`; likewise `rTL`/`rEdeg` give reference `ehat = 1`.
At steady state:

```
mhat* = X3^{-b} (1 - phi) + phi          (control fraction, with basal leak phi)
ehat* = mhat*                            (translation SS; (e/e0) = ehat*)
v_r0  = Vmax0 * ehat* * X3^{-a}          (= Vmax0 * (e/e0) * theta)
```

Power-law repression `X3^{-b}` is the **concrete operating-point instance** of the
Hill/partition control function of Eqs.~\eqref{eq:partition}-\eqref{eq:control-split}
-- the exact move the chapter already makes for `theta` ("collapses to the concrete
power law"). Prose states this explicitly so the two gateways stay consistent with
`gateways.tex`.

**mu convention.** Metabolite balances keep `S v = 0` (growth dilution of
metabolites negligible at these fluxes, as in the existing example and the urea
example). The `+mu` term lives in the **expression denominators** `(theta_m+mu)`,
`(theta_p+mu)`, where it is the chapter's visible signature and where -- for a
stable enzyme -- dilution is the dominant clearance term (see below). State this
modeling choice in one sentence.

### Chosen parameters and the fixed point

| symbol | value | meaning |
|---|---|---|
| `a` | `0.4` | allosteric feedback strength (activity) |
| `b` | `0.6` | transcriptional repression strength (expression) |
| `phi` | `0.05` | basal transcription (leak) fraction |
| `Vmax0` | `10.0` | reference committed-step capacity (`= kcat0 * e0`, AU) |
| `k3` | `0.666` | export rate constant, chosen to place `X3* = 4` |
| `alpha(r1), alpha(r2)` | `100` | fast internal steps |

Fixed point (solve `T = Vmax0 * ehat*(X3) * X3^{-a}`, `X3 = T/k3`):

```
X3*      = 4.00
T*       = 2.66
theta*   = X3*^{-a}                 = 4^{-0.4} ~ 0.574   (activity)
(e/e0)*  = X3*^{-b}(1-phi)+phi      ~ 0.464              (expression)
product  = theta* * (e/e0)*         ~ 0.266
naive cap (both = 1)                = 10.0
overshoot factor = Vmax0 / T*       = X3*^{a+b} (phi=0)  ~ 3.8x
```

The overshoot factor is `~X3*^{a+b}`, set only by the operating point and the two
exponents -- independent of the absolute capacity scale. These are the analytic
targets; the code confirms them from the BST integration (a fixed-point solve is
not needed at runtime -- integration to steady state finds `X3*` directly, and the
analytic fixed point is derived inline in the prose and checked against the run).

## Sequence-specific capacities (realistic-but-fake)

`E0` is a fake ~330-aa committed-step enzyme; its gene is 990 nt. Capacities and
timescales are computed from length + elongation, with **no precursors entered in
`S`** (the full precursor-in-`S` construction is the cell-free capstone's job; add
a one-line forward-pointer). These numbers set the reference abundance `e0` (hence
`Vmax0`) and the residence times that carry `mu`; the dimensionless `(e/e0)`
response is carried by the exponent `b`.

| quantity | value | source/derivation |
|---|---|---|
| `L_prot` | `330` aa | fake enzyme length |
| `L_gene` | `990` nt | `= 3 * L_prot` |
| `e_X` (txn elongation) | `60` nt/s | E. coli range 40-70 |
| `e_L` (tln elongation) | `16` aa/s | E. coli range 12-20 |
| `tau_X = L_gene/e_X` | `16.5` s | transcript elongation delay |
| `tau_L = L_prot/e_L` | `20.6` s | protein elongation delay |
| `r_X` | `~6` transcripts/min | initiation-limited loading x copy number |
| `r_L` | `~10` proteins/(transcript min) | ribosome loading |
| mRNA half-life | `2.5` min | -> `theta_m = ln2/2.5 = 0.277` /min |
| protein half-life | `~35` min | stable enzyme -> `theta_p ~ 0.020` /min |
| doubling time | `60` min | -> `mu = ln2/60 = 0.0116` /min |
| `kcat0` | `30` /s | committed-enzyme turnover |
| `e0` | `= r_L m_ref/(theta_p+mu)` | reference enzyme abundance |

**The mu signature, made concrete.** `theta_m + mu = 0.277 + 0.012 = 0.289` /min
(mu ~ 4% of transcript clearance) but `theta_p + mu = 0.020 + 0.012 = 0.032` /min
(**mu ~ 36% of enzyme clearance**). For the stable enzyme, growth dilution is a
dominant clearance term -- the strongest single illustration of why the chapter
carries `mu` in every intracellular balance. Emphasize this in the prose.

`Vmax0 = kcat0 * e0` is the reference capacity; it is fixed at `10` in the
arbitrary units of the metabolic chain (the sequence numbers calibrate the scale,
they do not need a molecule-to-AU conversion). The dimensionless `(e/e0)` is what
enters the bound.

## Control scheme (the made-up illustrative regulation)

- **Repressor senses X3.** A repressor TF, corepressed by the end product `X3`,
  represses `gene(E0)`. Steady-state transcription control fraction
  `u(X3) = X3^{-b}(1-phi) + phi` (the power-law instance of the repressive Hill /
  partition control), floored at the basal leak `phi` so the gene is never fully
  off -- `lambda = r_X * phi` is exactly the leak source of Eq.~\eqref{eq:mrna} and
  the basal occupancy `u_dagger` of Eq.~\eqref{eq:control-split}.
- **Allosteric activity.** `X3` inhibits `E0` activity, `theta(X3) = X3^{-a}`.
- Strengths `a = 0.4`, `b = 0.6` chosen so **each** gateway is individually
  load-bearing (`theta* ~ 0.57`, `(e/e0)* ~ 0.46`) and neither alone reaches the
  truth: expression alone throttles the cap to `4.64`, activity alone to `5.74`,
  together to `2.66`. Expression is the stronger throttle (`b > a`), foregrounding
  the new gateway, while the two remain visibly distinct.

## FBA recovery + four escalating cases

Load-bearing identity (proved inline): `Vmax0 * (e/e0)(X3*) * theta(X3*)` computed
from the model at the measured `X3*` **equals** the BST committed-step flux
`v_r0 = kcat0 * E0* * X3*^{-a}` -- nothing copied from `truth.vss`. Both factors
are obtained the same way `theta` already is (rate-with / rate-without at the
measured state), fully parallel to the existing `capacity_and_control`:

- `theta`   = `v_r0(a active) / v_r0(a zeroed)` at `X3*`, reference enzyme -> `X3*^{-a}`.
- `(e/e0)`  = `ehat*(b active) / ehat*(b zeroed)` at `X3*` -> `mhat*` -> the
  expression ratio. `e0` = reference enzyme (`b` zeroed); `Vmax0 = kcat0 * e0`.

Four cases on the `r0` bound (every other bound identical between runs):

| case | `ub[r0]` | value | reads |
|---|---|---|---|
| naive (both closed) | `Vmax0` | `10.0` | feedback-blind overshoot |
| expression only | `Vmax0 * (e/e0)` | `4.64` | enzyme throttled, activity blind |
| activity only | `Vmax0 * theta` | `5.74` | activity throttled, abundance blind |
| both open | `Vmax0 * (e/e0) * theta` | `2.66` | recovers `T*` |

Objective: maximize export `v_r3`. Downstream steps keep a generous constant
capacity (`100`, well above `Vmax0`) so `r0` is the sole throughput determinant,
exactly as in the S-system.

## Figure

One two-panel figure, `code/figs/feedback_gateway.pdf` (keeps the filename so the
existing `\includegraphics` reference resolves), regenerated by the script:

- **(a) Dual-feedback wiring schematic** (Makie: labeled nodes + arrows). The
  chain `X1->X2->X3`, the gene/transcript/enzyme `gene(E0)->m->E0` feeding `r0`,
  and two inhibition arrows from `X3`: a dashed **slow** arrow to transcription
  (`(e/e0)`) and a dashed **fast** arrow to `E0` activity (`theta`). Annotate the
  two timescales.
- **(b) Capacity ledger** (bar chart): descending bars `Vmax0 = 10` ->
  `x(e/e0) = 4.64` -> `xtheta = 5.74` -> `both = 2.66`, with the un-inhibited
  capacity as a dashed line at `10` and the BST truth `T* = 2.66` marked so the
  both-open bar visibly lands on it. This is "the bound is a ledger" made visual:
  each gateway whittles the committed-step capacity, and only both together reach
  the truth.

If the Makie schematic proves fiddly, fall back to a compact TikZ diagram in the
`.tex` for panel (a) and make `feedback_gateway.pdf` the ledger panel alone; note
this fallback in the implementation plan.

## Components

### 1. `code/feedback/Dual-Feedback.toml` (new; replaces `Linear-Feedback.toml`)

Five dynamic species `X1,X2,X3,m,E0`; the nine connection + kinetics records of
the truth-model table. Delete `Linear-Feedback.toml`. Metadata description updated
to "linear pathway under dual end-product control (expression + activity)".

### 2. `code/feedback/run_feedback.jl` (substantial rewrite)

- `feedback_truth()`: build the extended model; set `a, b, phi, Vmax0, k3`, the
  fast/export constants, and the sequence-derived `theta_m, theta_p, mu` (with the
  reference-normalizing coefficients); integrate to steady state; return
  `(species, Xss, reactions, vss, model)` plus the expression states. Keep the
  by-name lookups + ordering asserts already in the file so a package reordering
  fails loudly.
- `capacity_and_control(...)`: generalize to return `Vmax0`, `theta`, and `(e/e0)`
  from the measured state via rate-with/rate-without on the `a` and `b` orders
  (never reading `vss`).
- `feedback_fba(truth; expression::Bool, activity::Bool)`: one bound builder,
  four case calls. `ub[r0] = Vmax0 * (expression ? e_ratio : 1) * (activity ? theta : 1)`.
- Sequence-specific block: compute `tau_X, tau_L, r_X, r_L, theta_m, theta_p, mu,
  e0, Vmax0` from the length/rate/half-life inputs; print them (prose numbers).
- TDD gates (asserts):
  1. BST steady state satisfies the linear-chain mass balance (`r0=r1=r2=r3`
     within `rtol`), and `S_FEEDBACK == model.S` for the metabolic submatrix.
  2. `X3* > 1`, `theta* < 0.9`, `(e/e0)* < 0.9` (both gateways genuinely sub-unity,
     guarding against inert-feedback regressions).
  3. both-open FBA recovers the BST truth: `max|v_open - vss| < 1e-2`.
  4. naive overshoots: `v_naive[r3] - vss[r3] > 0.5`.
  5. the two partial cases are distinct and each strictly between truth and naive
     (`truth < expr_only < naive`, `truth < activity_only < naive`,
     `expr_only != activity_only`).
- Write `code/data/feedback_fba.csv`:
  `reaction, truth, naive, expression_only, activity_only, both_open`.
- Regenerate `code/figs/feedback_gateway.pdf` (two panels above).
- Print the fixed point (`X3*, T*, theta*, (e/e0)*`, overshoot) and the
  sequence-derived timescales for the prose.

### 3. Figure integration

Copy `code/figs/feedback_gateway.pdf` into `chapter/figures/` (filename unchanged;
placement and `\includegraphics` reference unchanged).

### 4. `chapter/sections/example_feedback.tex` (rewrite, ~5 long paragraphs)

Honor the long-paragraph preference (merge, do not fragment). Arc:

1. **Dual control motivation.** The urea example opened thermodynamic + kinetic;
   the first feedback pass opened only the regulatory `theta`. A committed
   biosynthetic step is controlled at two levels by its own end product -- fast
   allosteric inhibition of enzyme activity and slow transcriptional repression of
   the enzyme gene -- regulated as the trp operon is. This example opens both
   gateways at once on one enzyme (light trp anchor, abstract symbols kept).
2. **Extended BST truth.** The cascade `gene(E0)->m->E0` plus the metabolic chain
   as one S-system; `X3` represses transcription (slow) and enzyme activity
   (fast); integrate to the coupled steady state `(T*, X3*, E0*, m*)`. Power-law
   repression as the concrete instance of the partition-function control; the
   `S v = 0` / `+mu`-in-expression convention stated.
3. **Sequence-specific capacities + the mu signature.** Compute `r_X, r_L` from
   gene/protein length and elongation rates; the `(theta+mu)` steady states give
   `m*, p*, (e/e0)`; the stable enzyme's clearance is dominated by dilution
   (`mu ~ 36%` of `theta_p+mu`) -- the visible payoff of carrying `mu` everywhere.
4. **Coupled steady state + control scheme.** The repressor/leak scheme; the
   self-consistency in `X3` solved analytically and matched to the integration
   ("two ways, checked"); `theta* ~ 0.57`, `(e/e0)* ~ 0.45`.
5. **FBA escalation + the ledger payoff.** The bound `Vmax0 (e/e0) theta` equals
   the BST committed-step flux; the four cases (naive overshoot `~3.8x` -> each
   partial -> both recover `2.66`); nothing about either loop appears in `S`; both
   gateways are load-bearing and they multiply. Reference Fig.~\ref{fig:feedback}.

Update the figure caption to describe the two panels (wiring; capacity ledger).

## Style / conventions to honor

- No em dashes; no subsection headings in the example; no self-references ("this
  section..."); long complete paragraphs.
- Numeric results in prose tie to the figure (the ledger bars) or are structural
  setup numbers (exponents, lengths, elongation rates, half-lives, `mu`, `phi`),
  which are exempt.
- Keep the "never copy the truth" discipline: `Vmax0`, `(e/e0)`, `theta` all
  computed from the model at the measured state.
- Keep the by-name BSTModelKit lookups + ordering asserts.

## Out of scope

- Precursor accounting (NTP/aa columns in `S`) -- reserved for the cell-free
  capstone; only a one-line forward-pointer here.
- Translation-level regulation `w_j` (kept at 1; only transcription is repressed).
- Any change to the urea example or other sections.
- Dynamic (time-course) analysis; the example stays steady-state, with the
  timescales reported only to motivate the `mu` signature.

## Verification

- `run_feedback.jl` runs to completion; all TDD asserts pass; `feedback_fba.csv`
  and `feedback_gateway.pdf` regenerate.
- BST integration reproduces the analytic fixed point (`X3* ~ 4.0`, `T* ~ 2.66`,
  `theta* ~ 0.574`, `(e/e0)* ~ 0.464`) within tolerance.
- Both-open FBA recovers `T*`; naive overshoots to `~10`; the two partials are
  distinct and bracketed.
- `make` in `chapter/` compiles cleanly; the figure reference and label resolve;
  reported prose numbers match the script's stdout / CSV.
