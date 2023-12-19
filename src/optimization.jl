#=
For the moment, optimization is delegated to Optim.jl (hardcoded dependency). 
In order for Optim.jl to use gradients that are returned directly upon function 
evaluation, one has to use the fix described at: 
https://github.com/JuliaNLSolvers/Optim.jl/blob/master/docs/src/user/tipsandtricks.md

This is the reason for the cumbersome implementation of optimization with gradients.

Note that by default all particles in the system are assumed optimizable.

IMPORTANT: Note that we always work in cartesian coordinates.

=#

export optimize_geometry


""" 
    By default we work in cartesian coordinaes.
    Note that internally, when optimizing the cartesian positions, atomic units 
    are used.

"""
function Optimization.OptimizationFunction(system::AbstractSystem, calculator; kwargs...)
    mask = get_optimizable_mask(system) # mask is assumed not to change during optim.

    f = function(x::AbstractVector{<:Real}, p)
        new_system = update_optimizable_positions(system, x * u"bohr")
        energy = AtomsCalculators.potential_energy(new_system, calculator; kwargs...)
        austrip(energy)
    end

    g! = function(G::AbstractVector{<:Real}, x::AbstractVector{<:Real}, p)
        new_system = update_optimizable_positions(system, x * u"bohr")
        energy = AtomsCalculators.potential_energy(new_system, calculator; kwargs...)

        forces = AtomsCalculators.forces(new_system, calculator; kwargs...)
        # Translate the forces vectors on each particle to a single gradient for the optimization parameter.
        forces_concat = collect(Iterators.flatten(forces[mask]))

        # NOTE: minus sign since forces are opposite to gradient.
        G .= - austrip.(forces_concat)
    end
    OptimizationFunction(f; grad=g!)
end

function optimize_geometry(system::AbstractSystem, calculator; solver=Optim.LBFGS(), kwargs...)
    # Use current system parameters as starting positions.
    x0 = Vector(austrip.(get_optimizable_positions(system))) # Optim modifies x0 in-place, so need a mutable type.
    f_opt = OptimizationFunction(system, calculator)
    problem = OptimizationProblem(f_opt, x0, nothing) # Last argument needed in Optimization.jl.
    solve(problem, solver; kwargs...)
end
