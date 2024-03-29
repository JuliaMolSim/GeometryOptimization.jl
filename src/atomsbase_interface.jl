#
# Interface between AtomsBase.jl and GeometryOptimization.jl that provides 
# utility functions for manipulating systems.
#
# IMPORTANT: Note that we always work in cartesian coordinates.
#
export update_positions, update_not_clamped_positions, clamp_atoms


@doc raw"""
Creates a new system based on ``system`` but with atoms positions updated 
to the ones provided.

"""
function update_positions(system, positions::AbstractVector{<:AbstractVector{<:Unitful.Length}})
    particles = [Atom(atom; position) for (atom, position) in zip(system, positions)]
    AbstractSystem(system; particles)
end

@doc raw"""
Creates a new system based on ``system`` where the non clamped positions are
updated to the ones provided (in the order in which they appear in the system).
"""
function update_not_clamped_positions(system, positions::AbstractVector{<:Unitful.Length})
    mask = not_clamped_mask(system)
    new_positions = deepcopy(position(system))
    new_positions[mask] = reinterpret(reshape, SVector{3, eltype(positions)},
                                      reshape(positions, 3, :))
    update_positions(system, new_positions)
end

@doc raw"""
Returns a mask for selecting the not clamped atoms in the system.

"""
function not_clamped_mask(system)
    # If flag not set, the atom is considered not clamped.
    [haskey(a, :clamped) ? !a[:clamped] : true for a in system]
end

function not_clamped_positions(system)
    mask = not_clamped_mask(system)
    Iterators.flatten(system[mask, :position])
end

@doc raw"""
    Clamp given atoms in the system. Clamped atoms are fixed and their positions 
    will not be optimized. The atoms to be clamped should be given as a list of 
    indices corresponding to their positions in the system.

    """
function clamp_atoms(system, clamped_indexes::Union{AbstractVector{<:Integer},Nothing})
    clamped = falses(length(system))
    clamped[clamped_indexes] .= true
    particles = [Atom(atom; clamped=m) for (atom, m) in zip(system, clamped)]
    AbstractSystem(system, particles=particles)
end
