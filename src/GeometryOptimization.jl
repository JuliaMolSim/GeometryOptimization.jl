module GeometryOptimization

using LinearAlgebra
using StaticArrays
using Optimization, Optim, LineSearches
using OptimizationOptimJL
using AtomsBase
using AtomsCalculators
using Unitful
using UnitfulAtomic

include("atomsbase_interface.jl")
include("optimization.jl")

end
