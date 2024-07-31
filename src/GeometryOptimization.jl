module GeometryOptimization

using LinearAlgebra
using StaticArrays
using Optimization
using AtomsBase
using AtomsCalculators
using Unitful
using UnitfulAtomic

# Make sure Optim is always available
using OptimizationOptimJL
using LineSearches

AC = AtomsCalculators

include("atomsbase_interface.jl")
include("optimization.jl")

end
