using Distributed

addprocs(20)

@everywhere cd(pwd())
@everywhere using Distributions, Optim, StatsFuns, Random, JLD2, Printf, SharedArrays
@everywhere using LinearAlgebra, Parameters, Interpolations, DelimitedFiles,BenchmarkTools
@everywhere Itp = Interpolations
@everywhere import Base.show

@everywhere include.((
"types.jl",
"utilities.jl",
"parameters.jl",
"fcn_income.jl",
"fcn_utility.jl",
"vfi.jl"))
include.((
"montecarlo.jl",
"main.jl"))


@everywhere P = Params()
@everywhere IG = InitialGuess()
@time main(IG,P)

#@save "results.jld2" mc vfi
