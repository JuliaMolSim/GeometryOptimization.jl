module GeometryOptimization

using LinearAlgebra
using StaticArrays
using Optimization
using OptimizationOptimJL
using AtomsBase
using AtomsCalculators
using Unitful
using UnitfulAtomic

include("atomsbase_interface.jl")
include("strain.jl")
include("optimization.jl")

end
