module GeometryOptimization

using StaticArrays
using ComponentArrays
using Optimization
using OptimizationOptimJL
using AtomsBase
using AtomsCalculators
using Unitful
using UnitfulAtomic

include("atomsbase_interface.jl")
include("optimization.jl")

end
