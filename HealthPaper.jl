module HealthPaper

using Distributions, Optim, StatsFuns, Random, JLD2, Printf, SharedArrays
using LinearAlgebra, Parameters, Interpolations, DelimitedFiles,BenchmarkTools
const Itp = Interpolations
import Base.show

include("types.jl")
include("utilities.jl")
inlcude("parameters.jl")
include("fcn_income.jl")
include("fcn_utility.jl")
include("vfi.jl")
include("montecarlo.jl")

end
