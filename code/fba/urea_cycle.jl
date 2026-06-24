"""
    urea_cycle_model() -> NamedTuple

Return a pure-data FBA model for the urea cycle in HL-60 cells.

Stoichiometry, flux bounds, and objective are transcribed from:
  CHEME-5450-Example-Solution-UreaCycle-S2026.ipynb (varnerlab/Lecture-5430-FluxBalanceAnalysis)

Network (VFF format, from data/Network.net):
  v1  EC 6.3.4.5  ATP + Citrulline + Aspartate → AMP + Diphosphate + N-(L-Arginino)succinate  [irreversible]
  v2  EC 4.3.2.1  N-(L-Arginino)succinate → Fumarate + Arginine                               [irreversible]
  v3  EC 3.5.3.1  Arginine + H2O → Ornithine + Urea                                            [irreversible]
  v4  EC 2.1.3.3  Carbamoyl_phosphate + Ornithine → Orthophosphate + Citrulline                [irreversible]
  v5  EC 1.15.13.39  2 Arginine + 4 O2 + 3 NADPH + 3 H → 2 NO + 2 Citrulline + 3 NADP + 4 H2O [reversible, but ΔG ≪ −10 kJ/mol → treated irreversible]
  b1–b14  exchange reactions (±1000 default)

Flux bounds:
  Reversibility δⱼ from eQuilibrator (threshold −10 kJ/mol):
    v1: ΔG = −4.3  kJ/mol  → δ = 1 (reversible)
    v2: ΔG = −5.5  kJ/mol  → δ = 1 (reversible)
    v3: ΔG = −51.0 kJ/mol  → δ = 0 (irreversible)
    v4: ΔG = −30.3 kJ/mol  → δ = 0 (irreversible)
    v5: ΔG = −1220 kJ/mol  → δ = 0 (irreversible)
  Vmax = kcat × eₒ  (eₒ = 0.01 mmol/gDW; kcat from BRENDA, default 10 1/s):
    v1: kcat = 10.0  → Vmax = 0.100
    v2: kcat =  3.28 → Vmax = 0.0328
    v3: kcat = 190.0 → Vmax = 1.9
    v4: kcat = 410.0 → Vmax = 4.1
    v5: kcat = 10.0  → Vmax = 0.100
  lb = −δVmax, ub = Vmax for enzymatic reactions; ±1000 for exchange reactions.

Objective: maximise flux through b4 (urea export), c[b4] = 1.0.

Fields:
  S            stoichiometric matrix  (18 metabolites × 19 reactions)
  reactions    reaction names
  metabolites  metabolite names (rows of S, alphabetically sorted)
  lb           lower flux bounds
  ub           upper flux bounds
  c            objective coefficients (maximise b4)
"""
function urea_cycle_model()::NamedTuple

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
    # Flux bounds
    # Enzymatic reactions: lb = -δ*Vmax, ub = Vmax
    # Exchange reactions:  lb = -1000, ub = 1000
    # ------------------------------------------------------------------ #
    lb = Float64[
        -0.1,     # v1  reversible (δ=1), Vmax=0.1
        -0.0328,  # v2  reversible (δ=1), Vmax=0.0328
         0.0,     # v3  irreversible (δ=0), Vmax=1.9
         0.0,     # v4  irreversible (δ=0), Vmax=4.1
         0.0,     # v5  irreversible (δ=0), Vmax=0.1
        -1000.0,  # b1
        -1000.0,  # b2
        -1000.0,  # b3
        -1000.0,  # b4
        -1000.0,  # b5
        -1000.0,  # b6
        -1000.0,  # b7
        -1000.0,  # b8
        -1000.0,  # b9
        -1000.0,  # b10
        -1000.0,  # b11
        -1000.0,  # b12
        -1000.0,  # b13
        -1000.0,  # b14
    ]

    ub = Float64[
         0.1,     # v1
         0.0328,  # v2
         1.9,     # v3
         4.1,     # v4
         0.1,     # v5
         1000.0,  # b1
         1000.0,  # b2
         1000.0,  # b3
         1000.0,  # b4
         1000.0,  # b5
         1000.0,  # b6
         1000.0,  # b7
         1000.0,  # b8
         1000.0,  # b9
         1000.0,  # b10
         1000.0,  # b11
         1000.0,  # b12
         1000.0,  # b13
         1000.0,  # b14
    ]

    # ------------------------------------------------------------------ #
    # Objective: maximise urea *export* through b4 (column 9).
    # Exchange reaction b4 is defined as [] → M_Urea_c (positive flux = import).
    # Urea export corresponds to negative flux through b4.
    # With Max objective in JuMP, we set c[b4] = -1.0 so that maximising
    # sum(c*v) drives v[b4] as negative as possible (maximum export).
    # This matches the notebook convention (objective[b4] = -1, JuMP Min).
    # ------------------------------------------------------------------ #
    c = zeros(Float64, length(reactions))
    c[findfirst(==("b4"), reactions)] = -1.0

    # sanity checks
    @assert size(S) == (length(metabolites), length(reactions))
    @assert length(lb) == length(reactions)
    @assert length(ub) == length(reactions)
    @assert length(c)  == length(reactions)

    return (
        S           = S,
        reactions   = reactions,
        metabolites = metabolites,
        lb          = lb,
        ub          = ub,
        c           = c,
    )
end
