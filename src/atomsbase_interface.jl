#
# Interface between AtomsBase.jl and GeometryOptimization.jl that provides 
# utility functions for manipulating systems.
#
# IMPORTANT: Note that we always work in cartesian coordinates.
# BUG: Unitful does not play well with arrays in Julia 1.10. This is why 
# collect statements are needed to restor array type (otherwise we get Vector{Any}).
#


@doc raw"""
Creates a new system based on ``system`` but with atoms positions updated 
to the ones provided. Can also update lattice vectors if `bounding_box` is provided.

"""
function update_positions(system, positions::AbstractVector{<:AbstractVector{<:Unitful.Length}},
        bounding_box=bounding_box(system))
    particles = [Atom(atom; position) for (atom, position) in zip(system, positions)]
    AbstractSystem(system; particles, bounding_box)
end

# Method for natural units.
function update_positions(system, positions::AbstractVector{<:AbstractVector{<:Real}},
        bounding_box=bounding_box(system))
    particles = [Atom(atom; position=position * u"bohr") for (atom, position) in zip(system, positions)]
    AbstractSystem(system; particles, bounding_box)
end

@doc raw"""
Creates a new system based on ``system`` but with atoms positions 
updated to the ones provided and lattice vectors deformed according 
to the provided strain.
New generalized coordinates should be provided in a ComponentArray with 
component `atoms` and `strain`.

"""
function update_positions(system, positions::ComponentVector)
    deformation_tensor = I + voigt_to_full(positions.strain)
    # TODO: Do we want to apply the strain to the atoms too?
    strained_positions = [collect((x for x in deformation_tensor * position))
			  for position in positions.atoms]
    bbox = matrix_to_bbox(deformation_tensor * bbox_to_matrix(bounding_box(system)))
    update_positions(system, strained_positions, bbox)
end

function set_masked_positions(system, positions_flat)
    mask = not_clamped_mask(system)
    new_positions = [austrip.(x) for x in deepcopy(position(system))]
    positions_flat = austrip.(positions_flat)
    new_positions[mask] = reinterpret(reshape, SVector{3, eltype(positions_flat)},
                                      reshape(positions_flat, 3, :))
    return new_positions
end

@doc raw"""
Creates a new system based on ``system`` where the non clamped positions are
updated to the ones provided (in the order in which they appear in the system).
"""
function update_not_clamped_positions(system, positions::AbstractVector)
    new_positions = set_masked_positions(system, positions)
    update_positions(system, new_positions)
end

@doc raw"""
Creates a new system based on ``system`` where the non clamped positions and 
lattice vectors are updated to the ones provided.
Note that the `atoms`component of `positions` should be a vector of 
coordinates and that the `strain` component should be a 6-vector.

"""
function update_not_clamped_positions(system, positions::ComponentVector)
    new_positions = set_masked_positions(system, positions.atoms)
    update_positions(system, ComponentVector(; atoms=new_positions, positions.strain))
end

@doc raw"""
=======
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
