# Technical and Narrative Audit: *Mathematical Models in Biotechnology*

## Overall assessment

The chapter has a strong overall arc: balance, linear-program geometry, biological bounds, worked examples, and larger applications. The central idea that flux bounds provide an interface for biological information is pedagogically useful.

The draft needs technical revision before sentence-level polishing. Two issues materially affect the mathematical or numerical claims:

1. The mole-balance derivation mixes reactor-volume and biomass-specific bases.
2. The urea-cycle capacity calculation omits the conversion from seconds to hours.

Several other claims describe stronger biological provenance, prediction, or uncertainty results than the code currently supports.

## Highest-priority technical issues

### 1. The growth-dilution derivation mixes incompatible bases

In `chapter/sections/derivation.tex`, `V` is introduced as the reactor volume. The reaction term `v_j V` therefore treats `v_j` as a volumetric reaction rate. The worked examples later use biomass-specific fluxes in mmol gDW^-1 h^-1.

There is also a mathematical problem with the step leading to

```text
S v = mu x.
```

The chapter defines the specific culture volume as

```text
Vbar = V / B.
```

For an ordinary fixed-volume batch culture, biomass grows while reactor volume stays fixed, so

```text
(1 / Vbar) dVbar/dt = -mu.
```

In Equation 3, this term cancels the displayed `mu C_i` term. Thus Equation 4 does not follow for a fixed-volume culture under the definitions currently given.

A cleaner biomass-specific derivation is

```text
d(B x)/dt = B S v,
```

and therefore

```text
dx/dt + mu x = S v,
```

where `x` is metabolite amount per unit biomass and `v` is biomass-specific flux. Under balanced growth, `dx/dt = 0`, giving

```text
S v = mu x.
```

The text should choose one basis and retain it throughout. If intracellular concentration rather than amount per biomass is desired, the relevant intracellular or total-cell volume must be distinguished clearly from reactor volume.

The fed-batch discussion should then be reconsidered. Changes in reactor working volume and dilution of intracellular pools are different effects and should not enter through an ambiguously defined `Vbar`.

### 2. The urea-cycle capacities omit a factor of 3600

The chapter reports `kcat` in s^-1, `e0` in mmol gDW^-1, and `Vmax` in mmol gDW^-1 h^-1. The implementation in `code/fba/urea_cycle.jl` computes

```julia
Vmax = kcat .* e0
```

without converting seconds to hours.

For argininosuccinate lyase,

```text
3.28 s^-1 * 0.01 mmol gDW^-1
    = 0.0328 mmol gDW^-1 s^-1
    = 118.08 mmol gDW^-1 h^-1.
```

Therefore every reported urea flux is low by a factor of 3600 unless the `kcat` inputs were actually meant to be h^-1. The qualitative result that `v2` is the nominal bottleneck remains unchanged, but all dimensional flux values and their uncertainty intervals change.

This issue occurs in the prose, code comments, generated CSV files, and figure captions. The fix should be made at the model boundary, with an explicit conversion and a unit test.

### 3. The thermodynamic gate is too categorical for the calculation shown

The urea example classifies reaction direction from standard Gibbs free energy using a fixed -10 kJ/mol cutoff. Actual reaction direction depends on

```text
Delta G = Delta G degree + R T ln(Q),
```

not on standard free energy alone.

The chapter says the cutoff represents physiological concentration ranges, but neither the range nor the calculation is shown. The text should do one of the following:

- propagate stated physiological activity or concentration ranges and determine whether `Delta G` can change sign; or
- label the -10 kJ/mol rule explicitly as a heuristic directionality assignment.

The eQuilibrator conditions should also be reported: transformed or untransformed standard state, pH, ionic strength, temperature, magnesium assumptions, and uncertainty. The value reported for nitric oxide synthase is especially large and deserves a stoichiometric and thermodynamic audit.

### 4. The feedback example does not use its claimed sequence-derived capacities

The prose says gene length, protein length, and polymerase and ribosome elongation rates determine the transcription and translation capacities `r_X` and `r_L`.

The code calculates `tau_X` and `tau_L`, but these quantities are never used in the kinetic model. Instead, the transcription, translation, and degradation coefficients are calibrated directly from the clearance rates so that the unrepressed reference state has normalized abundance one.

The chapter should not say that the expression capacities were computed from sequence unless the implementation actually uses the elongation times to set those capacities. The current calculation uses biologically plausible sequence metadata to describe timescales, but those metadata do not determine the reported steady state.

There is a related equation mismatch. In `chapter/sections/gateways.tex`, the general protein balance is

```text
dp_j/dt = r_L,j w_j - (theta_p,j + mu) p_j,
```

which contains no dependence on transcript abundance. The feedback implementation correctly makes translation proportional to `m`. The general equation should include `m_j`, or `w_j` must be redefined so it is not merely a dimensionless control variable.

### 5. The feedback example demonstrates consistency, not independent recovery

The kinetic reference model defines the committed-step rate as

```text
v_r0 = Vmax0 * (E0/e0) * X3^(-a).
```

The FBA model then reads `E0` and `X3` from that same simulated steady state and places their product on the binding upper bound. Because the linear chain maximizes throughput until that upper bound binds, the matching result is built into the construction.

This is still a useful example: it shows exactly where simulated expression and regulatory information enters an FBA model. It is not an independent prediction or validation. Recommended wording:

- replace “truth” with “reference simulation”;
- replace “recovered” with “reproduced after transferring the simulated state into the bound”; and
- state directly that the exercise checks bookkeeping consistency.

The power law should also be normalized:

```text
theta = (X3 / X3_ref)^(-a).
```

As written, raising a dimensional concentration to a fractional power is invalid. In addition, `X3^(-a)` exceeds one when `X3 < 1`, which conflicts with the earlier definition of `theta` as a fraction in `[0,1]`. A bounded inhibitory control law would align better with the partition-function discussion.

The claim that a partition-function occupancy “collapses” to this power law is too strong without a local log-linear derivation, reference state, and stated range of validity.

### 6. The uncertainty prose and pseudocode do not exactly match the computation

Several details should be corrected:

- Algorithm 1 reuses a single symbol `Z`, implying perfectly correlated draws. The code uses independent normal draws for individual `kcat` and `Delta G` values, plus one shared enzyme-abundance draw.
- The algorithm does not show distinct draws for substrate concentration and `Km`, although the code samples both.
- The chapter says the nitric oxide synthase branch stays silent across the ensemble. In `code/data/urea_fba_uq.csv`, `v5` has a small positive mean and nonzero standard deviation, so it is active in rare samples.
- The figure caption says every nonzero flux has identical bootstrap uncertainty. The CSV shows small differences between `v1/v2` and `v3/v4`, again because the branch becomes active in rare samples.
- The statement that `v2` alone carries the uncertainty should be softened to “dominates the uncertainty in most draws.”

Because the parameter distributions are assumed rather than estimated from replicate data, “Monte Carlo uncertainty propagation” or “probabilistic sensitivity analysis” is clearer than “parametric bootstrap.” If “parametric bootstrap” is retained, the text should explain that the distributions are modeling assumptions, not fitted sampling distributions.

The Config A and Config B comparison would also be cleaner with common random numbers for the shared parameters, varying only the imputed `f2` factor.

## Additional technical corrections

### Stoichiometry and bounds both contain biological information

The statements that reaction-network stoichiometry is “fully known” and that bounds, rather than stoichiometry, are where biological information enters are too broad.

Genome-scale reconstructions contain uncertainty in reaction membership, directionality, compartments, cofactors, GPR rules, exchange definitions, and biomass composition. The stoichiometric matrix and biomass reaction therefore contain substantial biological assumptions.

A more defensible central claim is:

> Flux bounds are a major interface through which condition-specific thermodynamic, kinetic, expression, and regulatory information enters a constraint-based model.

This preserves the useful gateway framing without making it exclusive.

### The growth-aware form needs qualifications in the LP

The LP writes the balance as `S v = 0 (or mu x)`. The growth-aware form remains linear only if `x` is fixed and `mu` is either fixed or introduced linearly as a scalar decision variable. If both metabolite pools and growth are unknown, their product is nonlinear.

The chapter should state what is treated as measured, fixed, or optimized.

It should also explain the relation between `mu x` dilution demands and the conventional biomass pseudo-reaction so that growth-associated precursor drains are not counted twice.

### The genome-scale dimension claim is overstated

`chapter/sections/linearprogram.tex` says reaction count often exceeds metabolite count by an order of magnitude. Genome-scale models generally have more reactions than metabolites, but not reliably ten times more. Replace this with a rank-based statement:

```text
rank(S) < n, often leaving a large right null space.
```

### The FVA interpretation is incorrect in one sentence

The chapter says that a zero FVA range means a reaction is pinned “regardless of objective choice.” FVA pins a reaction only within the optimal or near-optimal objective set used for that FVA calculation. A different objective can produce a different feasible face and different range.

The nominal urea result itself is sound on this point: an independent FVA calculation at optimal urea export gave identical minimum and maximum values for all 19 fluxes. Thus the nominal optimum is unique to numerical tolerance, and the phrase “single, sharply determined answer” is supported for that particular model.

### Database provenance is incomplete

The code comments describe `kcat = 10` for `v1` and `v5` as a default, while the chapter says turnover numbers for all five enzymes were taken from BRENDA. Defaults, database records, and assumed values must be distinguished.

For each database-derived parameter, record at least:

- database accession or query;
- enzyme and organism;
- substrate and assay direction;
- assay temperature and pH where available;
- units and conversion;
- selected record and aggregation rule; and
- access date.

The same concern applies to the common enzyme abundance `e0 = 0.01 mmol/gDW`. It is an assumed reference scale, not an HL-60 measurement. Since a shared `e0` scales all capacities together, it sets the absolute flux scale directly.

### The Park data are cross-organism rather than HL-60-specific

The saturation data combine human, mouse, yeast, and *E. coli* values. Park et al. measured absolute metabolite concentrations and fluxes in *E. coli*, yeast, and a mammalian cell line, but these values are not measurements from HL-60 cells.

The chapter acknowledges the lack of HL-60 measurements late in the example, but the opening claim that this is a prediction “in HL-60 cells” remains stronger than the inputs justify. “A urea-cycle network parameterized with cross-organism public data” would be more accurate.

### Capstone wording should be tightened

The chapter calls the two capstones “published systems,” but the cited 2023 integrated cell-free study is listed in the bibliography as a bioRxiv preprint. Use “reported studies” or identify it explicitly as a preprint.

The cell-free paper reports 63 time-resolved metabolites, so the chapter can replace “on the order of sixty” with the exact number.

The Wayman paper supports roughly threefold improvement from predicted knockout designs. The stronger claim that the best single knockout alone produced that improvement should be verified against the paper or softened to “predicted knockouts and knockout combinations improved production roughly threefold.”

Finally, flux estimation and strain design are not literally the same computational operation. Once a deletion set is specified, both problems solve an LP with altered bounds. Finding the deletion set is a discrete, often combinatorial optimization problem. The conclusion should retain this qualification.

### Position the gateway framework against established methods

The gateway framing unifies several established model classes. The chapter would benefit from a short positioning paragraph connecting:

- thermodynamic bounds to thermodynamics-based flux analysis;
- `kcat` and enzyme abundance to enzyme-constrained FBA and GECKO-like models;
- expression balances to metabolism-and-expression models;
- regulatory factors to regulatory FBA and probabilistic regulatory constraints; and
- knockout design to strain-design and bilevel optimization methods.

The original contribution can then be described as a common pedagogical notation and workflow rather than as the invention of the individual constraints.

## Code-quality observations

The code is compact and readable, and the core nominal and feedback tests pass. The following improvements would make it safer and more reproducible:

- Add explicit physical-unit handling or conversion constants around `kcat` and flux capacities.
- Add input validation for vector lengths, finite values, nonnegative saturation factors, and lower bounds not exceeding upper bounds.
- Make `solve_fba` handle an unsuccessful solve explicitly instead of constructing a `DataFrame` from `nothing`.
- Record optimization termination status in ensemble output rather than silently skipping failed samples.
- Load or generate the saturation aggregation from `park_saturation.csv`; the current UQ implementation hardcodes the reduced values despite the design specification saying the CSV should be loaded and reduced.
- Record which reaction is binding in every UQ draw. This would directly support statements about the probability that `v2` remains the bottleneck.
- Add a test that checks dimensional conversion, not only the current hardcoded bounds.
- Add an FVA routine or test so the uniqueness claim remains checked if the network changes.
- Avoid relying on private package internals such as `BSTModelKit._powerlaw`, or pin the package version and wrap the dependency in a local function.
- Mark the feedback parameters and arbitrary units clearly as illustrative rather than biologically calibrated.

## Professional-language audit

### Assessment

The chapter's structure is sound, but the prose is not consistently suitable for a technical reference work. It frequently uses dramatic, conversational, or evaluative language where a direct statement would be clearer. The sentence “The escalation is the whole point” is one example of a broader pattern.

The current source contains approximately:

- 67 uses of “gateway”;
- 18 uses of “honest,” including equation-label references;
- 12 uses of “exactly”;
- 9 uses of “entirely”; and
- 9 uses of “truth.”

These repetitions create an argumentative tone. The text often tells the reader that a result is important instead of stating the result and its implication.

### Editorial standard to apply

Use the following rules throughout the revision:

1. State the model, result, and limitation directly.
2. Do not characterize an equation or assumption as “honest,” “natural,” or “obvious.”
3. Reserve “truth” for independently observed ground truth. Use “reference simulation” for simulated data.
4. Introduce the “gateway” metaphor once. Thereafter use the mathematical names: reversibility factor, capacity, expression factor, regulatory factor, lower bound, and upper bound.
5. Remove sentences that announce importance, such as “the whole point,” “the result itself,” or “this repays stating precisely.”
6. Avoid conversational metaphors: “poured into,” “load-bearing,” “miss the mark,” “marked down,” “trickle,” “natural home,” “bites hardest,” and “lever.”
7. Avoid unsupported intensifiers: “exactly,” “entirely,” “manifestly,” “genuinely,” “sharp,” and “hard.” Retain them only when they have a defined mathematical meaning.
8. Avoid anthropomorphism. Models do not “see,” bounds do not “carry weight,” and cells do not “spend” flux.
9. Prefer one technical claim per sentence. Several paragraphs currently use sentences longer than 50 words with multiple qualifications and conclusions.
10. End paragraphs with a result or limitation, not a slogan.

### Terms to remove or restrict

#### “Honest”

Remove all evaluative uses of “honest.” The term implies that the standard formulation or competing formulations are dishonest. It also does not identify the mathematical distinction.

Use:

- “balanced-growth balance” for `S v = mu x`;
- “steady-state approximation” for `S v = 0`;
- “factorized bound” for the general bound; and
- “baseline bound” for the simplified bound.

The LaTeX label `eq:honest-constraint` should also be renamed, for example to `eq:balanced-growth-constraint`, so the loaded term does not persist in cross-references.

#### “Truth”

Replace “truth,” “kinetic truth,” and “true throughput” in the feedback example with:

- “reference kinetic model”;
- “reference simulation”; or
- “reference throughput.”

The BST model is a constructed simulation with selected parameters. Calling it “truth” overstates its status and makes the comparison sound adversarial.

#### “Gateway”

Keep “gateway” in the section title and once in the introductory definition. Replace most later occurrences with the factor's name. For example:

- “opening the thermodynamic gateway” -> “setting the reversibility factor”;
- “the kinetic gateway” -> “the capacity and saturation terms”;
- “opening both gateways” -> “including both factors”; and
- “a poorly constrained gateway” -> “a poorly constrained bound.”

The current open/closed language is also logically inconsistent. A factor set to one is described as information-free, physically fully open, and not yet opened by data. Direct mathematical language avoids this ambiguity.

### Required sentence-level revisions

#### Abstract: `chapter/Chapter.tex:17-44`

The abstract is one dense paragraph with repeated metaphor and promotional phrasing. It should state scope, method, examples, and conclusion without “textbook form discards,” “opening,” “truth,” or “run in opposite directions.”

Recommended replacement:

> Intracellular fluxes, kinetic parameters, and metabolite concentrations are sparsely measured at genome scale. Flux balance analysis estimates feasible flux distributions from network stoichiometry, flux bounds, and an optimization objective. This chapter derives the flux constraint from a mole balance and identifies the assumptions that lead to the steady-state form. It then expresses each reaction bound in terms of reversibility, enzyme capacity, substrate saturation, enzyme abundance, and regulation. A urea-cycle model illustrates the use of thermodynamic and kinetic data, and a feedback-inhibited pathway illustrates the effects of expression and allosteric regulation on a reaction bound. Two larger studies show how related constraints have been used in cell-free protein synthesis and strain design. Together, these examples show how condition-specific data and proposed interventions alter flux predictions through reaction bounds.

This version is shorter, factual, and does not claim that the examples establish an independent “truth.”

#### Introduction: `chapter/sections/introduction.tex:4-20`

Current:

> each is, at bottom, a question about the intracellular flux distribution

Replace with:

> each requires an estimate of the intracellular flux distribution.

Current:

> the rates at which the reactions of metabolism actually run

Replace with:

> the rates of metabolic reactions under a specified condition.

Current:

> The problem is therefore to predict a metabolic phenotype from what can be known rather than from what a direct calculation would require.

Replace with:

> The modeling problem is therefore to estimate phenotype from incomplete measurements.

#### Introduction: `chapter/sections/introduction.tex:22-46`

Current:

> they resolve the dynamics in full but demand a rate constant, an affinity, or a kinetic order for every step

Replace with:

> they represent dynamics explicitly but require kinetic parameters for each modeled reaction.

Current:

> Constraint-based models make the opposite trade.

Replace with:

> Constraint-based models require fewer kinetic parameters and provide less dynamic detail.

Current:

> Flux balance analysis is the constraint-based workhorse

Replace with:

> Flux balance analysis is a widely used constraint-based method.

Current:

> it fixes a whole subspace of them

Replace with:

> it defines a feasible set that generally contains multiple flux distributions.

#### Introduction: `chapter/sections/introduction.tex:48-73`

This paragraph repeats the abstract and should be cut by at least half. It should not use “pinned,” “truth,” “opened,” or “no hand calculation could reach.”

Recommended replacement:

> This chapter first derives the flux constraint from a mole balance and then formulates FBA as a linear program. It represents each reaction bound using factors for reversibility, enzyme capacity, substrate saturation, expression, and regulation. Two small examples show how these factors affect a predicted flux distribution. Two larger studies illustrate related applications in cell-free protein synthesis and strain design.

#### Derivation: `chapter/sections/derivation.tex:4-9`

Current:

> FBA rests on an algebraic constraint rather than a differential equation, and that constraint is a mole balance in disguise. Writing the balance out, instead of asserting metabolic steady state by fiat, makes explicit what is assumed away...

Replace with:

> The FBA constraint follows from a mole balance. Deriving it identifies the assumptions required to obtain the steady-state form.

“In disguise” and “by fiat” are rhetorical and unnecessary.

#### Derivation: `chapter/sections/derivation.tex:26-28`

Current:

> Equation 1 is exact: it presumes nothing about growth rate, feeding strategy, or steady state...

Replace with:

> Equation 1 is the unreduced balance. No steady-state or constant-volume assumption has yet been applied.

This states what has and has not been assumed without the absolute claim “presumes nothing.”

#### Derivation: `chapter/sections/derivation.tex:70-73`

Current:

> Equation 4 is the honest flux constraint: net metabolic production ... is balanced not against zero, but against the rate at which that species is diluted by its own growing culture.

Replace with:

> Equation 4 is the balanced-growth form. Net production balances growth dilution for each intracellular species.

#### Derivation: `chapter/sections/derivation.tex:75-101`

Remove “discarding ... outright,” “cannot see,” “manifestly,” “fails from the other direction,” “corner case,” and “formal nuisance.”

Recommended replacement for the argumentative portion:

> The conventional FBA form sets the right-hand side to zero. This approximation is appropriate when growth dilution is small relative to metabolic turnover and intracellular pool sizes. It does not describe accumulation in a dynamic cell-free reaction, where `mu = 0` but `dx/dt` is nonzero. Fed-batch operation also requires explicit treatment of changing volume and external flows.

#### Derivation: `chapter/sections/derivation.tex:103-122`

Current:

> A closed algebraic statement can thereby represent a genuinely open system: openness becomes entirely a matter of which reactions are admitted into S and how their bounds are set...

Replace with:

> Exchange reactions represent uptake and secretion within the stoichiometric model. Their inclusion and bounds define the modeled system boundary.

#### Linear program: `chapter/sections/linearprogram.tex:4-17`

Current:

> is a statement of conservation and nothing more: it does not by itself pick out any particular flux distribution

Replace with:

> imposes conservation but generally does not determine a unique flux distribution.

Current:

> a box of flux bounds

Use “lower and upper flux bounds” after the first geometric explanation. Repeated use of “box” is informal and unnecessary.

#### Linear program: `chapter/sections/linearprogram.tex:38-56`

Current:

> The linear-programming machinery is necessary, rather than a matter of convenience...

Replace with:

> Linear programming selects an optimum from the feasible flux set.

Current:

> a whole subspace, not a point

Replace with:

> a nontrivial subspace.

Current:

> does not go away merely because an optimum has been found

Replace with:

> An optimum may therefore be nonunique.

#### Bounds section: `chapter/sections/gateways.tex:4-30`

The opening paragraph should define the factorization without claiming that bounds are “an honest ledger” or describing factors as repeatedly opened and closed.

Recommended replacement:

> Conservation and flux bounds define the feasible set, and the objective selects an optimum from that set. Reaction bounds can incorporate several forms of condition-specific information. Equation 7 factors an upper bound into a reference capacity and dimensionless terms for reversibility, enzyme abundance, substrate saturation, and regulation. Setting a factor from measurements or a mechanistic model changes the feasible flux range for that reaction.

#### Bounds section: repeated metaphor

Use direct subsection transitions, even if subsection headings are not added:

- “The thermodynamic gateway is the switch...” -> “The reversibility factor `delta_j` determines whether reverse flux is allowed.”
- “When the reaction is allowed to proceed, its capacity is set...” -> “The reference capacity is `Vmax,j = kcat,j e0`.”
- “empty gateway” / “fully open one” -> “ranges from zero to one.”
- “opening this gateway makes the bound...” -> “Because `f_j` depends on substrate concentration, the bound becomes state dependent.”
- “natural port for measured data” -> “Measured transcript or protein abundance can be used to estimate this factor.”
- “visible signature of the growth-aware bookkeeping” -> “The growth rate appears in both steady-state abundance expressions.”
- “Every gateway just opened demands data...” -> “These factors require thermodynamic, kinetic, expression, or regulatory data.”
- “honest statement of a bound” -> “baseline bound used when additional data are unavailable.”

#### Urea example: `chapter/sections/example_urea.tex:4-16`

Current:

> two gateways ... filled entirely from public databases, already pin down a sharp prediction

Replace with:

> the model uses public thermodynamic and kinetic estimates to define the reaction bounds.

Current:

> the content poured into S, l, u, and c

Replace with:

> the stoichiometric matrix, bounds, and objective specified for the application.

#### Urea example: result language

Replace the following throughout:

- “pinned irreversible” -> “assigned a zero lower bound”;
- “left shut” -> “set to one”;
- “drives urea secretion” -> “maximizes urea secretion”;
- “single, sharply determined answer” -> “unique optimum”;
- “entire enzymatic backbone saturated” -> “reactions `v1` through `v4` carried the limiting flux”;
- “none of it was left over to spend” -> “the objective assigned no flux to the competing branch”;
- “the branch stays silent” -> “the branch flux was zero in most samples”;
- “the saturation gateway barely moves the prediction is itself the result” -> “Including the available saturation data changed the median only slightly”; and
- “single hard bottleneck” -> “binding capacity constraint.”

These replacements state the optimization result without personifying flux or announcing significance.

#### Feedback example: `chapter/sections/example_feedback.tex:4-34`

The opening is too long and uses “throttles,” “currencies,” “honest capacity,” and “the whole of the construction.”

Recommended replacement:

> The feedback example adds two controls to the committed reaction `r0`. The end product `X3` inhibits the existing enzyme and represses transcription of its gene. The first effect changes enzyme activity; the second changes enzyme abundance. The model represents these effects with separate regulatory and expression factors in the upper bound on `r0`.

#### Feedback example: `chapter/sections/example_feedback.tex:71-99`

Current:

> The growth rate enters those balances in a way that repays stating precisely, because it is easy to misread.

Replace with:

> The same growth rate `mu` appears in both clearance terms.

Current:

> The consequence is a clean signature.

Delete. State the two percentages directly.

Current:

> though the number diluting them is one and the same, and it is exactly this asymmetry that justifies carrying `mu`...

Replace with:

> Growth dilution contributes about 4% of total mRNA clearance and 37% of total protein clearance. It should therefore be retained explicitly in both balances.

#### Feedback example: `chapter/sections/example_feedback.tex:115-127`

Replace “throttles transcription downward,” “trickle of enzyme,” “load-bearing,” and “miss the mark by a visible margin.”

Recommended replacement:

> Product-dependent repression lowers transcription toward the basal expression rate. The selected feedback exponents make both the expression and activity factors contribute to the reduction in `r0` capacity. Models that include only one factor therefore overestimate the reference throughput.

#### Feedback example: `chapter/sections/example_feedback.tex:161-176`

This is the passage containing “The escalation is the whole point.” Replace the entire result paragraph with:

> The four FBA calculations differed only in the upper bound on `r0`. With both factors set to one, the predicted throughput was 10. Including only the expression factor reduced it to 4.64, and including only the activity factor reduced it to 6.10. Including both factors gave 2.83, compared with 2.66 from the reference kinetic simulation, a relative difference of about 6%. The stoichiometric matrix was unchanged across the four calculations; the differences resulted from the expression and activity terms in the `r0` bound.

This paragraph contains the complete result. The sentences about escalation, a ledger, markdowns, and “the throughput the cell actually ran” add tone but no information.

#### Feedback figure caption: `chapter/sections/example_feedback.tex:181-193`

Replace “capacity ledger,” “gateways together approach,” “truth,” and “landing at” with:

> Right: predicted throughput under four specifications of the `r0` upper bound. The uninhibited capacity is 10; the expression-only and activity-only values are 4.64 and 6.10, respectively; and the combined value is 2.83. The reference kinetic model gives 2.66.

#### Capstones: `chapter/sections/capstones.tex:4-14`

Replace “known truth,” “at a scale no hand calculation could reach,” and “modeler's guess” with:

> The preceding examples use small networks with one degree of freedom. The following studies apply related constraints to larger models using experimental data.

#### Capstones: cell-free discussion

Replace:

- “the natural home” -> “a useful application”;
- “manifestly not at steady state” -> “is dynamic”;
- “fails from the direction opposite to slow growth” -> “does not satisfy the pseudo-steady-state assumption”;
- “one and the same volume” -> “share a common reaction volume”; and
- “for once, under direct experimental observation” -> “can be estimated directly from time-resolved measurements.”

Recommended condensed paragraph:

> Cell-free protein synthesis provides a useful setting for dynamic balances. The growth rate is zero, but transcript, protein, and metabolite concentrations change during the batch reaction. The pseudo-steady-state approximation therefore does not apply. Because the measured species share the reaction volume, their accumulation can be estimated directly from time-resolved concentration data.

#### Capstones: design discussion

Current:

> The bounds are, in this use, not a report of what the cell is doing but a lever for what the cell is made to do.

Replace with:

> In this application, bound changes represent proposed interventions rather than measured conditions.

The final paragraph repeats the same estimation/design comparison several times. Replace it with:

> Flux estimation and strain design both evaluate an FBA model after its bounds have been specified. Estimation derives bounds from measurements; design evaluates proposed bound changes. Design also requires a separate discrete search over candidate interventions, often formulated as a bilevel optimization problem.

#### Outlook: `chapter/sections/outlook.tex:4-17`

Current:

> a control function built by hand

Replace with:

> a specified control function.

Current:

> learned directly from the data they are meant to explain

Replace with:

> estimated from time-resolved data.

#### Outlook: `chapter/sections/outlook.tex:19-31`

Current:

> What changes with scale is the weight the bounds are asked to carry.

Replace with:

> Larger null spaces increase the number of flux distributions consistent with stoichiometry alone.

Current:

> does the most damage to a prediction. Closing gateways with data matters most exactly where the null-space argument bites hardest.

Replace with:

> Poorly constrained bounds then permit wider flux ranges and increase prediction uncertainty.

#### Outlook: `chapter/sections/outlook.tex:33-46`

Replace the entire paragraph with:

> Time-resolved cell-free data support models that retain accumulation terms rather than imposing metabolic steady state. The same formulation can be used for growing cultures when intracellular concentrations and growth rates are measured with sufficient temporal resolution. In those settings, fluxes can be estimated from the dynamic balance instead of the steady-state approximation.

### Phrases that should not appear in the revised chapter

Unless part of a direct quotation, remove these expressions:

- “The escalation is the whole point.”
- “That ... is itself the result.”
- “by fiat”
- “honest constraint” / “honest ledger” / “honest capacity”
- “kinetic truth” / “true throughput”
- “content poured into”
- “left shut” / “reopened”
- “none left over to spend”
- “stays silent”
- “load-bearing”
- “miss the mark”
- “read as a ledger”
- “marked down” / “markdowns”
- “trickle of enzyme”
- “landing on”
- “natural home”
- “formal nuisance”
- “manifestly”
- “for once”
- “does the most damage”
- “bites hardest”
- “at bottom”
- “run in opposite directions”

### Recommended structural tightening

The chapter should retain its current technical sequence but reduce repetition:

1. Introduction: define the problem and scope.
2. Balances: derive the dynamic and steady-state forms.
3. FBA: define bounds, objective, nonuniqueness, and FVA.
4. Factorized bounds: define the four factors and relate them to established methods.
5. Urea example: model definition, unique optimum, and concise uncertainty result.
6. Feedback example: factor comparison against a reference simulation.
7. Applications and outlook: cell-free dynamics, strain design, and limitations.

The final paragraph of the introduction should preview this sequence once. Later sections should begin with their immediate technical purpose rather than restating the full gateway framework.

The chapter title is broader than its contents. A more descriptive title is:

> **Flux Balance Analysis as an Integrative Framework for Biotechnology**

## Recommended revision order

1. Fix the balance basis and redefine all units.
2. Correct the seconds-to-hours conversion and regenerate results and figures.
3. Audit the thermodynamic and kinetic provenance.
4. Align the general expression equations with the feedback implementation.
5. Reframe the feedback example as a consistency demonstration.
6. Correct the UQ pseudocode and empirical claims.
7. Add a short related-methods paragraph and qualify the central bounds claim.
8. Tighten repeated gateway language, the abstract, and the capstone conclusion.
9. Add a code-and-data availability statement with environment and solver versions.

## Verification performed

- `code/fba/test_urea_cycle.jl`: passed.
- `code/feedback/test_dual_feedback.jl`: passed.
- LaTeX `make` in `chapter/`: completed successfully.
- Independent FVA at the nominal optimal urea-export value: every reaction had equal minimum and maximum flux to numerical tolerance, confirming a unique nominal optimum.
- Generated UQ CSVs were compared against the numerical statements in the chapter.
- The existing user edit to `chapter/sections/example_feedback.tex` was not changed.

No source changes were made as part of the audit.
