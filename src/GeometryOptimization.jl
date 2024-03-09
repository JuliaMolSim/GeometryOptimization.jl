module GeometryOptimization

using LinearAlgebra
using StaticArrays
using Optimization
using OptimizationOptimJL
using AtomsBase
using AtomsCalculators
using Unitful
using UnitfulAtomic
using ComponentArrays

export update_positions, update_not_clamped_positions, clamp_atoms
include("atomsbase_interface.jl")
export compute_voigt_strain, voigt_to_full
include("strain.jl")
export minimize_energy!
include("optimization.jl")

end
