#= Test Geometry Optimization on an aluminium supercell.
=#
@testitem "Stain interface: cell deformation under strain" setup=[TestCases] begin
    using AtomsBase
    using Unitful
    using UnitfulAtomic
    using GeometryOptimization
    using Random
    using ComponentArrays

    silicon = TestCases.silicon

    # Test elongation along first vector.
    strain = Vector([1., 0, 0, 0, 0, 0.])
    new_silicon = apply_voigt_strain(silicon, strain)
    @test [v[1] for v in bounding_box(new_silicon)] ≈ 2. * [v[1] for v in bounding_box(silicon)] rtol=1e-3
end

@testitem "Stain interface: atoms positions deformation under strain" setup=[TestCases] begin
    using AtomsBase
    using Unitful
    using UnitfulAtomic
    using GeometryOptimization
    using Random
    using ComponentArrays

    silicon = TestCases.silicon

    # Test elongation along first vector.
    strain = Vector([1., 0, 0, 0, 0, 0.])
    new_silicon = apply_voigt_strain(silicon, strain)
    @test [v[1] for v in position(new_silicon)] ≈ 2. * [v[1] for v in position(silicon)] rtol=1e-3
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
