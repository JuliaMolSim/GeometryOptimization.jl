module GeometryOptimization

using AtomsBase
using AtomsCalculators
using DocStringExtensions
using LinearAlgebra
using Optimization
using StaticArrays
using Unitful
using UnitfulAtomic

# Make sure Optim is always available
using OptimizationOptimJL
using LineSearches

# Useful shortcuts
using AtomsCalculators: Energy, Forces, Virial, calculate
AC = AtomsCalculators

@template METHODS =
"""
$(TYPEDSIGNATURES)

$(DOCSTRING)
"""

include("dof_management.jl")
include("optimization.jl")
include("callbacks.jl")

export minimize_energy!
export GeoOptDefaultCallback

end
