#= Test Geometry Optimization on an aluminium supercell.
=#
@testitem "AtomsBase interface: check positions+cell updating" setup=[TestCases] begin
    using AtomsBase
    using Unitful
    using UnitfulAtomic
    using GeometryOptimization
    using Random
    using ComponentArrays

    silicon_supercell = TestCases.silicon_supercell

    Random.seed!(1234)
    orig_positions = position(silicon_supercell)
    orig_lattice = bounding_box(silicon_supercell)
    σ = 0.1u"angstrom"
    new_positions = orig_positions + [σ * rand(Float64, size(v)) for v in orig_positions]
    new_lattice = orig_lattice + [σ * rand(Float64, size(v)) for v in orig_lattice]
    new_general_pos = ComponentVector(atoms = new_positions, bounding_box = new_lattice)
    new_system = update_positions(silicon_supercell, new_general_pos)
    
    @test position(new_system) == new_positions
    @test bounding_box(new_system) == new_lattice
end

@testitem "AtomsBase interface: check positions+cell updating (with mask)" setup=[TestCases] begin
    using AtomsBase
    using Unitful
    using UnitfulAtomic
    using GeometryOptimization
    using Random
    using ComponentArrays

    silicon_supercell = TestCases.silicon_supercell

    Random.seed!(1234)
    orig_positions = position(silicon_supercell)
    orig_lattice = bounding_box(silicon_supercell)
    σ = 0.1u"angstrom"
    new_positions = orig_positions + [σ * rand(Float64, size(v)) for v in orig_positions]
    new_not_clamped_positions = reduce(vcat, new_positions)
    new_lattice = orig_lattice + [σ * rand(Float64, size(v)) for v in orig_lattice]
    new_general_pos = ComponentVector(atoms = new_not_clamped_positions,
                                      bounding_box = reduce(vcat, new_lattice))
    new_system = update_not_clamped_positions(silicon_supercell, new_general_pos)
    
    @test position(new_system) == new_positions
    @test bounding_box(new_system) == new_lattice
end
