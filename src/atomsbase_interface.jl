#
# Interface between AtomsBase and GeomOpt that provides 
# utility functions for manipulating systems.
#
# This interface helps defining which particles in the system are considered 
# optimizable. Note that by default all particles are assumed to be optimizable. 
# The user can call clamp_atoms to fix atoms whose position should not be optimized.
#
# IMPORTANT: Note that we always work in cartesian coordinates.
#
export update_positions, update_optimizable_coordinates, clamp_atoms


@doc raw"""
    update_postions(system::AbstractSystem, positions::Vector{<:AbstractVector{<:Real}}) where {L <: Unitful.Length}

Creates a new system based on ``system`` but with atoms positions updated to the ones specified in `positions`, 
using cartesian coordinates.
"""
function update_positions(system::AbstractSystem, positions::AbstractVector{<:AbstractVector{L}}) where {L <: Unitful.Length}
    particles = [
                 Atom(atom, position=SVector{3,L}(position)) 
                 for (atom, position) in zip(system.particles, positions)]
    AbstractSystem(system, particles=particles)
end

@doc raw"""
    update_optimizable_coordinates(system::AbstractSystem, positions::Vector{<:AbstractVector{<:L}}) where {L <: Unitful.Length}

Creates a new system based on ``system`` with the optimizable coordinates are
updated to the ones provided (in the order in which they appear in system.particles), cartesian coordinates version..
"""
function update_optimizable_coordinates(system::AbstractSystem, positions::AbstractVector{<:L}) where {L <: Unitful.Length}
    mask = get_optimizable_mask(system)
    new_positions = Vector.(position(system)) # make mutable.
    new_positions[mask] = [eachcol(reshape(positions, 3, :))...]
    update_positions(system, new_positions)
end

@doc raw"""
    set_optimizable_mask(system::AbstractSystem, mask::AbstractVector{<:Bool})

Sets the mask defining which coordinates of the system can be optimized. 
The mask is a vector of booleans, specifying which of the atoms can be optimized.

By default (when no mask is specified), all particles are assumed optimizable.

"""
function set_optimizable_mask(system::AbstractSystem, mask::AbstractVector{<:Bool})
    particles = [Atom(atom; optimizable=m) for (atom, m) in zip(system.particles, mask)]
    AbstractSystem(system, particles=particles)
end

@doc raw"""
    get_optimizable_mask(system::AbstractSystem) -> AbstractVector{<:Bool}

Returns the optimizable mask of the system (see documentation for `set_optimizable_mask`.
"""
function get_optimizable_mask(system::AbstractSystem)
    # If flag not set, the atom is considered to be optimizable.
    [haskey(a, :optimizable) ? a[:optimizable] : true for a in system.particles]
end

function get_optimizable_coordinates(system::AbstractSystem)
    mask = get_optimizable_mask(system)
    return collect(Iterators.flatten(system[mask, :position]))
end

@doc raw"""
    Clamp given atoms if the system. Clamped atoms are fixed and their positions 
    will not be optimized. The atoms to be clamped should be given as a list of 
    indies corresponding to their positions in system.particles.

    """
function clamp_atoms(system::AbstractSystem, clamped_indexes::Union{AbstractVector{<:Integer},Nothing})
    mask = trues(length(system.particles))
    mask[clamped_indexes] .= false
    clamped_system = set_optimizable_mask(system, mask)
    return clamped_system
end
