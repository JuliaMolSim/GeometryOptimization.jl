#
# Interface between AtomsBase.jl and GeometryOptimization.jl that provides 
# utility functions for manipulating systems.
#

"""
Convert the input system to a kind of system we can work with, i.e. one where
we have the `particles` keyword argument for setting new particles and positions.
"""
function convert_to_updatable(system)
    if length(system) > 0 && !(system[1] isa Atom)
        particles = [Atom(; pairs(atom)...) for atom in system]
        AbstractSystem(system; particles)
    else
        system
    end
end

"""
Creates a new system based on ``system`` but with atoms positions updated 
to the ones provided.
"""
function update_positions(system, positions::AbstractVector{<:AbstractVector{<:Unitful.Length}})
    particles = [Atom(atom; position) for (atom, position) in zip(system, positions)]
    AbstractSystem(system; particles)
end

"""
Creates a new system based on ``system`` where the non clamped positions are
updated to the ones provided (in the order in which they appear in the system).
"""
function update_not_clamped_positions(system, positions::AbstractVector{<:Unitful.Length})
    mask = [!get(atom, :clamped, false) for atom in system]
    new_positions = deepcopy(position(system, :))
    new_positions[mask] = reinterpret(reshape, SVector{3, eltype(positions)},
                                      reshape(positions, 3, :))
    update_positions(system, new_positions)
end

"""
Clamp given atoms in the system. Clamped atoms are fixed and their positions
will not be optimized. The atoms to be clamped should be given as a list of
indices corresponding to their positions in the system.

!!! warn "Experimental"
    This is this a very experimental interface and will likely change
    in the future.
"""
function clamp_atoms(system, clamped_indexes::Union{AbstractVector{<:Integer},Nothing})
    system = convert_to_updatable(system)
    clamped = falses(length(system))
    clamped[clamped_indexes] .= true
    particles = [Atom(atom; clamped=msk) for (atom, msk) in zip(system, clamped)]
    AbstractSystem(system; particles)
end
