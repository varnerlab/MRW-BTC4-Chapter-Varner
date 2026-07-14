"""
    urea_cycle_model() -> NamedTuple

Return a pure-data FBA model for the urea cycle in HL-60 cells.

Stoichiometry, flux bounds, and objective are transcribed from:
  CHEME-5450-Example-Solution-UreaCycle-S2026.ipynb (varnerlab/Lecture-5430-FluxBalanceAnalysis)

Network (VFF format, from data/Network.net):
  v1  EC 6.3.4.5  ATP + Citrulline + Aspartate → AMP + Diphosphate + N-(L-Arginino)succinate  [reversible]
  v2  EC 4.3.2.1  N-(L-Arginino)succinate → Fumarate + Arginine                               [reversible]
  v3  EC 3.5.3.1  Arginine + H2O → Ornithine + Urea                                            [irreversible]
  v4  EC 2.1.3.3  Carbamoyl_phosphate + Ornithine → Orthophosphate + Citrulline                [irreversible]
  v5  EC 1.14.13.39  2 Arginine + 4 O2 + 3 NADPH + 3 H → 2 NO + 2 Citrulline + 3 NADP + 4 H2O [reversible, but ΔG ≪ −10 kJ/mol → treated irreversible]
  b1–b14  exchange reactions (±1000 default)

Flux bounds:
  Reversibility δⱼ from eQuilibrator (threshold −10 kJ/mol):
    v1: ΔG = −4.3  kJ/mol  → δ = 1 (reversible)
    v2: ΔG = −5.5  kJ/mol  → δ = 1 (reversible)
    v3: ΔG = −51.0 kJ/mol  → δ = 0 (irreversible)
    v4: ΔG = −30.3 kJ/mol  → δ = 0 (irreversible)
    v5: ΔG = −1220 kJ/mol  → δ = 0 (irreversible)
  Vmax = kcat × eₒ × 3600  (eₒ = 0.01 mmol/gDW; kcat from BRENDA in s⁻¹, ×3600 converts to h⁻¹):
    v1: kcat = 10.0  → Vmax = 360
    v2: kcat =  3.28 → Vmax = 118.08
    v3: kcat = 190.0 → Vmax = 6840
    v4: kcat = 410.0 → Vmax = 14760
    v5: kcat = 10.0  → Vmax = 360
  lb = −δVmax, ub = Vmax for enzymatic reactions; ±1000 for exchange reactions.

Objective: maximise urea export through b4.  Under the secretion-positive
  convention (each exchange written M_i → [], positive flux = secretion), b4
  reads Urea → [], so urea export is positive flux.  The objective coefficient
  is set to c[b4] = +1.0 so that maximising c·v drives v[b4] as positive as
  possible, i.e. maximises the rate of urea export.

Fields:
  S            stoichiometric matrix  (18 metabolites × 19 reactions)
  reactions    reaction names
  metabolites  metabolite names (rows of S, alphabetically sorted)
  lb           lower flux bounds
  ub           upper flux bounds
  c            objective coefficients (c[b4] = +1.0; maximises urea export)
"""
const KCAT0 = [10.0, 3.28, 190.0, 410.0, 10.0]   # 1/s, from BRENDA (v1..v5)
const E0    = 0.01                                # mmol/gDW, reference enzyme abundance
const DG0   = [-4.3, -5.5, -51.0, -30.3, -1220.0] # kJ/mol, from eQuilibrator (v1..v5)

function urea_cycle_model(; kcat=KCAT0, e0=E0, dG=DG0, dG_threshold=-10.0, f=ones(5))::NamedTuple

    # ------------------------------------------------------------------ #
    # Metabolites (alphabetically sorted, 18 total)
    # ------------------------------------------------------------------ #
    metabolites = [
        "M_AMP_c",                      #  1
        "M_ATP_c",                      #  2
        "M_Carbamoyl_phosphate_c",      #  3
        "M_Diphosphate_c",              #  4
        "M_Fumarate_c",                 #  5
        "M_H2O_c",                      #  6
        "M_H_c",                        #  7
        "M_L-Arginine_c",               #  8
        "M_L-Aspartate_c",              #  9
        "M_L-Citrulline_c",             # 10
        "M_L-Ornithine_c",              # 11
        "M_N-(L-Arginino)succinate_c",  # 12
        "M_NADPH_c",                    # 13
        "M_NADP_c",                     # 14
        "M_Nitric_oxide_c",             # 15
        "M_Orthophosphate_c",           # 16
        "M_Oxygen_c",                   # 17
        "M_Urea_c",                     # 18
    ]

    # ------------------------------------------------------------------ #
    # Reactions (19 total: 5 enzymatic + 14 exchange)
    # ------------------------------------------------------------------ #
    reactions = [
        "v1", "v2", "v3", "v4", "v5",
        "b1", "b2", "b3", "b4", "b5",
        "b6", "b7", "b8", "b9", "b10",
        "b11", "b12", "b13", "b14",
    ]

    # ------------------------------------------------------------------ #
    # Stoichiometric matrix S (18 × 19)
    # Rows = metabolites (order above), Columns = reactions (order above)
    # Built row-by-row; zeros omitted for clarity.
    # ------------------------------------------------------------------ #
    #       v1    v2    v3    v4    v5    b1    b2    b3    b4    b5    b6    b7    b8    b9   b10   b11   b12   b13   b14
    S = Float64[
     #  1  AMP_c
        1     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0  ;
     #  2  ATP_c
       -1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0  ;
     #  3  Carbamoyl_phosphate_c
        0     0     0    -1     0     1     0     0     0     0     0     0     0     0     0     0     0     0     0  ;
     #  4  Diphosphate_c
        1     0     0     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0     0  ;
     #  5  Fumarate_c
        0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0  ;
     #  6  H2O_c
        0     0    -1     0     4     0     0     0     0     0     0     0     0     0     0     0     0     0     1  ;
     #  7  H_c
        0     0     0     0    -3     0     0     0     0     0     0     0     0     0     0     1     0     0     0  ;
     #  8  L-Arginine_c
        0     1    -1     0    -2     0     0     0     0     0     0     0     0     0     0     0     0     0     0  ;
     #  9  L-Aspartate_c
       -1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0     0     0  ;
     # 10  L-Citrulline_c
       -1     0     0     1     2     0     0     0     0     0     0     0     0     0     0     0     0     0     0  ;
     # 11  L-Ornithine_c
        0     0     1    -1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0  ;
     # 12  N-(L-Arginino)succinate_c
        1    -1     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0  ;
     # 13  NADPH_c
        0     0     0     0    -3     0     0     0     0     0     0     0     0     0     1     0     0     0     0  ;
     # 14  NADP_c
        0     0     0     0     3     0     0     0     0     0     0     0     0     0     0     0     0     1     0  ;
     # 15  Nitric_oxide_c
        0     0     0     0     2     0     0     0     0     0     0     0     0     0     0     0     1     0     0  ;
     # 16  Orthophosphate_c
        0     0     0     1     0     0     0     0     0     0     0     0     1     0     0     0     0     0     0  ;
     # 17  Oxygen_c
        0     0     0     0    -4     0     0     0     0     0     0     0     0     1     0     0     0     0     0  ;
     # 18  Urea_c
        0     0     1     0     0     0     0     0     1     0     0     0     0     0     0     0     0     0     0  ;
    ]

    # ------------------------------------------------------------------ #
    # Exchange convention: secretion-positive (standard COBRA / Orth 2010).
    # The literal above lists each exchange column b1..b14 in the {} -> M_i
    # direction (a +1 in the metabolite row). Negate those 14 columns so each
    # exchange reads M_i -> {}: a positive exchange flux is then SECRETION and a
    # negative flux UPTAKE, matching the outward boundary arrows and c[b4]=+1.
    # ------------------------------------------------------------------ #
    S[:, 6:end] = -S[:, 6:end]

    # ------------------------------------------------------------------ #
    # Input validation
    # ------------------------------------------------------------------ #
    length(kcat) == 5 || throw(ArgumentError("kcat must have length 5, got $(length(kcat))"))
    length(dG)   == 5 || throw(ArgumentError("dG must have length 5, got $(length(dG))"))
    length(f)    == 5 || throw(ArgumentError("f must have length 5, got $(length(f))"))
    all(isfinite, kcat)   || throw(ArgumentError("kcat must be finite: $kcat"))
    all(isfinite, dG)     || throw(ArgumentError("dG must be finite: $dG"))
    all(isfinite, f)      || throw(ArgumentError("f must be finite: $f"))
    isfinite(e0)          || throw(ArgumentError("e0 must be finite: $e0"))
    all(f .>= 0)          || throw(ArgumentError("saturation factors f must be nonnegative: $f"))

    # ------------------------------------------------------------------ #
    # Flux bounds, built from parameters (Eq. general-bound with e/e0 = theta = 1)
    #   delta_j = 1 if dG_j > threshold (reversible), else 0 (irreversible)
    #   Vmax_j  = kcat_j * e0 * 3600 ;  cap_j = Vmax_j * f_j
    #   lb_j = -delta_j * cap_j ,  ub_j = cap_j   (enzymatic)
    #   exchange reactions: lb = -1000, ub = 1000
    # ------------------------------------------------------------------ #
    delta = Float64[g > dG_threshold ? 1.0 : 0.0 for g in dG]
    Vmax  = kcat .* e0 .* 3600.0   # kcat is s^-1 (BRENDA); convert to h^-1 to match e0 (mmol/gDW)
    cap   = Vmax .* f
    n_ex  = 14
    lb = vcat(-delta .* cap, fill(-1000.0, n_ex))
    ub = vcat(cap,           fill( 1000.0, n_ex))

    # ------------------------------------------------------------------ #
    # Objective: maximise urea *export* through b4 (column 9).
    # Under the secretion-positive convention above, exchange b4 reads
    # M_Urea_c -> {}, so a positive flux through b4 is urea secretion.
    # Maximising export is therefore maximising v[b4] directly: set c[b4] = +1.
    # ------------------------------------------------------------------ #
    c = zeros(Float64, length(reactions))
    c[findfirst(==("b4"), reactions)] = 1.0

    # sanity checks
    @assert size(S) == (length(metabolites), length(reactions))
    @assert length(lb) == length(reactions)
    @assert length(ub) == length(reactions)
    @assert length(c)  == length(reactions)
    all(lb .<= ub) || throw(ArgumentError("lower bounds must not exceed upper bounds: lb=$lb, ub=$ub"))

    return (
        S           = S,
        reactions   = reactions,
        metabolites = metabolites,
        lb          = lb,
        ub          = ub,
        c           = c,
    )
end

"""
    solve_flux(m) -> Union{Vector{Float64},Nothing}

Solve the FBA linear program for model `m` (Max c'v s.t. Sv=0, lb<=v<=ub).
Return the optimal flux vector, or `nothing` if the solve is not optimal/feasible.
"""
function solve_flux(m)
    model = Model(HiGHS.Optimizer); set_silent(model)
    n = length(m.reactions)
    @variable(model, m.lb[i] <= v[i=1:n] <= m.ub[i])
    @constraint(model, m.S * v .== 0)
    @objective(model, Max, sum(m.c[i] * v[i] for i in 1:n))
    optimize!(model)
    is_solved_and_feasible(model) ? value.(v) : nothing
end

"""
    solve_fba(m) -> DataFrame

Solve `m` and return a DataFrame with columns `reaction, flux`.
"""
function solve_fba(m)
    v = solve_flux(m)
    v === nothing && error("solve_fba: linear program did not solve to optimality (infeasible or unbounded)")
    DataFrame(reaction = m.reactions, flux = v)
end

"""
    fva(m; tol=1e-6) -> DataFrame

Flux variability analysis: for each reaction, minimize and maximize its flux subject to
`Sv=0`, `lb<=v<=ub`, and the objective held within `tol` of its optimal value. Returns a
DataFrame with columns `reaction, vmin, vmax`.
"""
function fva(m; tol=1e-6)
    vopt = solve_flux(m)
    vopt === nothing && error("fva: nominal model did not solve to optimality")
    zopt = sum(m.c[i] * vopt[i] for i in eachindex(m.c))
    n = length(m.reactions)
    vmin = zeros(n); vmax = zeros(n)
    for i in 1:n
        for sense in (:min, :max)
            model = Model(HiGHS.Optimizer); set_silent(model)
            @variable(model, m.lb[k] <= v[k=1:n] <= m.ub[k])
            @constraint(model, m.S * v .== 0)
            @constraint(model, sum(m.c[k] * v[k] for k in 1:n) >= zopt - tol)
            if sense == :min
                @objective(model, Min, v[i])
            else
                @objective(model, Max, v[i])
            end
            optimize!(model)
            val = is_solved_and_feasible(model) ? value(v[i]) : NaN
            sense == :min ? (vmin[i] = val) : (vmax[i] = val)
        end
    end
    DataFrame(reaction = m.reactions, vmin = vmin, vmax = vmax)
end
