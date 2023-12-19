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
export update_positions, update_optimizable_positions, clamp_atoms


@doc raw"""
Creates a new system based on ``system`` but with atoms positions updated 
to the ones specified in `positions`, using cartesian coordinates.
"""
function update_positions(system, positions::AbstractVector{<:AbstractVector{<:Unitful.Length}})
    particles = [Atom(atom; position) for (atom, position) in zip(system, positions)]
    AbstractSystem(system; particles)
end

@doc raw"""
Creates a new system based on ``system`` with the optimizable coordinates are
updated to the ones provided (in the order in which they appear in the system),
cartesian coordinates version..
"""
function update_optimizable_positions(system, positions::AbstractVector{<:Unitful.Length})
    mask = get_optimizable_mask(system)
    new_positions = deepcopy(position(system))
    new_positions[mask] = reinterpret(reshape, SVector{3, eltype(positions)},
                                      reshape(positions, 3, :))
    update_positions(system, new_positions)
end

@doc raw"""
Sets the mask defining which coordinates of the system can be optimized. 
The mask is a vector of booleans, specifying which of the atoms can be optimized.

By default (when no mask is specified), all particles are assumed optimizable.

"""
function set_optimizable_mask(system, mask::AbstractVector{<:Bool})
    particles = [Atom(atom; optimizable=m) for (atom, m) in zip(system, mask)]
    AbstractSystem(system, particles=particles)
end

@doc raw"""
    get_optimizable_mask(system) -> AbstractVector{<:Bool}

Returns the optimizable mask of the system (see documentation for `set_optimizable_mask`.
"""
function get_optimizable_mask(system)
    # If flag not set, the atom is considered to be optimizable.
    [haskey(a, :optimizable) ? a[:optimizable] : true for a in system]
end

function get_optimizable_positions(system)
    mask = get_optimizable_mask(system)
    collect(Iterators.flatten(system[mask, :position]))
end

@doc raw"""
    Clamp given atoms if the system. Clamped atoms are fixed and their positions 
    will not be optimized. The atoms to be clamped should be given as a list of 
    indies corresponding to their positions in the system.

    """
function clamp_atoms(system, clamped_indexes::Union{AbstractVector{<:Integer},Nothing})
    mask = trues(length(system))
    mask[clamped_indexes] .= false
    set_optimizable_mask(system, mask)  # Clamp by setting the mask.
end
