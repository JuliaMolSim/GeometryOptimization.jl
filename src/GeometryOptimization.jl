module GeometryOptimization

using AtomsBase
using AtomsCalculators
using DocStringExtensions
using LinearAlgebra
using LineSearches
using Optim
using StaticArrays
using Unitful
using UnitfulAtomic

# Useful shortcuts
using AtomsCalculators: Energy, Forces, Virial
AC = AtomsCalculators

@template METHODS =
"""
$(TYPEDSIGNATURES)

$(DOCSTRING)
"""

include("dof_management.jl")
include("optimization.jl")
include("optim.jl")
include("callbacks.jl")

export minimize_energy!
export GeoOptDefaultCallback
export OptimLBFGS, OptimCG, OptimSD

end
