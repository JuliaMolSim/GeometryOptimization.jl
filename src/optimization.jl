#=
# Note that by default all particles in the system are assumed optimizable.
# IMPORTANT: Note that we always work in cartesian coordinates.
=#
using DFTK
export minimize_energy!


function update_state(original_system, new_system, state)
    if bounding_box(original_system) != bounding_box(new_system)
        return DFTK.DFTKState()
    else
        return state
    end
end
    
""" 
By default we work in cartesian coordinates.
Note that internally, when optimizing the cartesian positions, atomic units 
are used.
"""
function Optimization.OptimizationFunction(system, calculator; pressure=0.0, kwargs...)
    mask = not_clamped_mask(system)  # mask is assumed not to change during optim.

    # TODO: Note that this function will dispatch appropriately when called with 
    # a Component vector.
    f = function(x, p)
        new_system = update_not_clamped_positions(system, x * u"bohr")
        state = update_state(system, new_system, calculator.state)
        energy = AtomsCalculators.potential_energy(new_system, calculator; state, kwargs...)
        austrip(energy) + pressure * calculator.state.scfres.basis.model.unit_cell_volume
    end

    g! = function(G::AbstractVector{<:Real}, x::AbstractVector{<:Real}, p)
        new_system = update_not_clamped_positions(system, x * u"bohr")
        energy = AtomsCalculators.potential_energy(new_system, calculator; kwargs...)

        forces = AtomsCalculators.forces(new_system, calculator; kwargs...)
        # Translate the forces vectors on each particle to a single gradient for the optimization parameter.
        forces_concat = collect(Iterators.flatten(forces[mask]))

        # NOTE: minus sign since forces are opposite to gradient.
        G .= - austrip.(forces_concat)
    end
    OptimizationFunction(f; grad=g!)
end

function minimize_energy!(system, calculator; pressure=0.0, solver=Optim.LBFGS(), kwargs...)
    # Use current system parameters as starting positions.
    x0 = austrip.(not_clamped_positions(system)) # Optim modifies x0 in-place, so need a mutable type.
    f_opt = OptimizationFunction(system, calculator; pressure)
    problem = OptimizationProblem(f_opt, x0, nothing)  # Last argument needed in Optimization.jl.
    solve(problem, solver; kwargs...)
end
