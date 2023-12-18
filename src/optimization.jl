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
function construct_optimization_function(system::AbstractSystem, calculator; kwargs...)
    f = function(x::AbstractVector{<:Real})
            x = 1u"bohr" .* x # Work in atomic units.
            new_system = update_optimizable_coordinates(system, x)
            austrip(AtomsCalculators.potential_energy(new_system, calculator; kwargs...))
    end
    return f
end

function construct_optimization_function_w_gradients(system::AbstractSystem, calculator; kwargs...)
    mask = get_optimizable_mask(system) # mask is assumed not to change during optim.

    f = function(x::AbstractVector{<:Real}, p)
        x = 1u"bohr" .* x # Work in atomic units.
        new_system = update_optimizable_coordinates(system, x)
        energy = AtomsCalculators.potential_energy(new_system, calculator; kwargs...)
        return austrip(energy)
    end

    g! = function(G::AbstractVector{<:Real}, x::AbstractVector{<:Real}, p)
        x = 1u"bohr" .* x # Work in atomic units.
        new_system = update_optimizable_coordinates(system, x)
        energy = AtomsCalculators.potential_energy(new_system, calculator; kwargs...)

        forces = AtomsCalculators.forces(new_system, calculator; kwargs...)
        # Translate the forces vectors on each particle to a single gradient for the optimization parameter.
        forces_concat = collect(Iterators.flatten(forces[mask]))

        # NOTE: minus sign since forces are opposite to gradient.
        G .= - austrip.(forces_concat)
    end
    return (f, g!)
end

function optimize_geometry(system::AbstractSystem, calculator;
        no_gradients=false, solver=Optim.LBFGS(), kwargs...)

    # Use current system parameters as starting positions.
    x0 = Vector(austrip.(get_optimizable_coordinates(system))) # Optim modifies x0 in-place, so need a mutable type.

    if no_gradients
        f = construct_optimization_function(system, calculator)
        f_opt = OptimizationFunction(f)
    else
        (f, g!) = construct_optimization_function_w_gradients(system, calculator)
        f_opt = OptimizationFunction(f; grad=g!)
    end
    problem = OptimizationProblem(f_opt, x0, nothing) # Last argument needed in Optimization.jl.
    solve(problem, solver; kwargs...)
end
