module GeometryOptimization

using AtomsBase
using AtomsCalculators
using DocStringExtensions
using LinearAlgebra
using LineSearches
using StaticArrays
using Unitful
using UnitfulAtomic

using AtomsCalculators: Energy, Forces, Virial
const AC = AtomsCalculators

@template METHODS =
"""
$(TYPEDSIGNATURES)

$(DOCSTRING)
"""

export minimize_energy!
export GeoOptDefaultCallback
export Autoselect, OptimLBFGS, OptimCG, OptimSD

include("dof_management.jl")
include("minimize_energy.jl")
include("optim.jl")
include("callbacks.jl")

end
