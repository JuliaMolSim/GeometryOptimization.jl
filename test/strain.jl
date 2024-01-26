#= Test Geometry Optimization on an aluminium supercell.
=#
@testitem "Stain interface: test cell deformation under strain" setup=[TestCases] begin
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

@testitem "Stain interface: test atoms positions deformation under strain" setup=[TestCases] begin
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
