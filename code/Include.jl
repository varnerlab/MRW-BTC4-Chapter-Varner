const _ROOT = @__DIR__
const _PATH_TO_SRC  = joinpath(_ROOT, "src")
const _PATH_TO_DATA = joinpath(_ROOT, "data")
const _PATH_TO_FIGS = joinpath(_ROOT, "figs")

using Pkg
Pkg.activate(_ROOT)

using CSV, DataFrames, JSON, Statistics, LinearAlgebra, Random
using JuMP, HiGHS, BSTModelKit, CairoMakie

include(joinpath(_PATH_TO_SRC, "Runtime.jl"))
