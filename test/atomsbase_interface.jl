#= Test Geometry Optimization on an aluminium supercell.
=#
@testitem "AtomsBase interface: check positions+cell updating" setup=[TestCases] begin
    using AtomsBase
    using Unitful
    using UnitfulAtomic
    using GeometryOptimization
    using Random

    silicon_supercell = TestCases.silicon_supercell

    Random.seed!(1234)
    orig_positions = position(silicon_supercell)
    orig_lattice = bounding_box(silicon_supercell)
    σ = 0.1u"angstrom"
    new_positions = orig_positions + [σ * rand(Float64, size(v)) for v in orig_positions]
    new_lattice = orig_lattice + [σ * rand(Float64, size(v)) for v in orig_lattice]
    new_system = update_positions(silicon_supercell, new_positions, new_lattice)
    
    @test position(new_system) == new_positions
    @test bounding_box(new_system) == new_lattice
end

@testitem "AtomsBase interface: cell deformation under strain" setup=[TestCases] begin
    using AtomsBase
    using Unitful
    using UnitfulAtomic
    using GeometryOptimization
    using ComponentArrays

    silicon = TestCases.silicon

    # Test elongation along first component.
    strain = Vector([1., 0, 0, 0, 0, 0.])
    new_general_pos = ComponentVector(;atoms = position(silicon), strain)
    new_system = update_positions(silicon, new_general_pos)

    @test [v[1] for v in bounding_box(new_system)] ≈ 2. * [v[1] for v in bounding_box(silicon)] rtol=1e-3
end

@testitem "Stain interface: atoms positions deformation under strain" setup=[TestCases] begin
    using AtomsBase
    using Unitful
    using UnitfulAtomic
    using GeometryOptimization
    using ComponentArrays

    silicon = TestCases.silicon

    # Test elongation along first component.
    strain = Vector([1., 0, 0, 0, 0, 0.])
    new_general_pos = ComponentVector(atoms = collect.(position(silicon)); strain)
    new_system = update_positions(silicon, new_general_pos)
    @test [v[1] for v in position(new_system)] ≈ 2. * [v[1] for v in position(silicon)] rtol=1e-3
end

@testitem "AtomsBase interface: check positions+cell updating (with mask)" setup=[TestCases] begin
    using AtomsBase
    using Unitful
    using UnitfulAtomic
    using GeometryOptimization
    using Random
    using LinearAlgebra
    using ComponentArrays

    silicon_supercell = TestCases.silicon_supercell

    Random.seed!(1234)
    orig_positions = position(silicon_supercell)
    σ = 0.1u"angstrom"
    new_positions = orig_positions + [σ * rand(Float64, size(v)) for v in orig_positions]
    new_not_clamped_positions = reduce(vcat, new_positions)
    
    # Test elongation along first component.
    strain = Vector([1., 0, 0, 0, 0, 0.])
    new_general_pos = ComponentVector(;atoms = new_not_clamped_positions, strain)
    new_system = update_not_clamped_positions(silicon_supercell, new_general_pos)
    
    new_positions_strained = [(I + voigt_to_full(strain)) * x for x in position(silicon_supercell)]
    @test [v[1] for v in position(new_system)] ≈ 2. * [v[1] for v in new_positions] rtol=1e-3
    @test [v[1] for v in bounding_box(new_system)] ≈ 2. * [v[1] for v in bounding_box(silicon_supercell)] rtol=1e-3
end

@testitem "Stain interface: strain from deformation" setup=[TestCases] begin
    bbox = [[1., 0, 0],
                   [0., 2., 0],
                   [0., 0, 3]]
    bbox_pert = [[1., 0, 0],
                   [0., 4., 0],
                   [0., 0, 3]]
    strain = compute_voigt_strain(bbox, bbox_pert)
    @test strain ≈ [0., 1., 0, 0, 0, 0] rtol=1e-6
    
    # Test with y-z rotation.
    bbox_pert = [[6,0,0],
                 [0,2,1],
                 [0,1.5,3]]
    strain = compute_voigt_strain(bbox, bbox_pert)
    @test strain ≈ [5., 0, 0, 1., 0, 0] rtol=1e-6
end
