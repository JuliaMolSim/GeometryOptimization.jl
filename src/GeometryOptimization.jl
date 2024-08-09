module GeometryOptimization

using DocStringExtensions
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

@template METHODS =
"""
$(TYPEDSIGNATURES)

$(DOCSTRING)
"""

include("clamping_updating_positions.jl")
include("optimization.jl")
include("printing.jl")

export minimize_energy!
export GeoOptDefaultPrint

end
