#= Test Geometry Optimization on an aluminium supercell.
=#
@testitem "Test AtomsCalculators energy and forces interface" tags=[:dont_test_mpi, :dont_test_windows] setup=[TestCalculators] begin
    using Unitful
    using UnitfulAtomic
    using OptimizationOptimJL
    using AtomsBase
    using GeometryOptimization
    
    
    bounding_box = 10.0u"angstrom" .* [[1, 0, 0.], [0., 1, 0], [0., 0, 1]]
    atoms = [:H => [0, 0, 0.0]u"bohr", :H => [0, 0, 1.9]u"bohr"]
    system = periodic_system(atoms, bounding_box)
    system = clamp_atoms(system, [1])
    
    calculator = TestCalculators.DummyCalculator()
    
    solver = OptimizationOptimJL.LBFGS()
    optim_options = (f_tol=1e-6, iterations=4, show_trace=false)
    
    results = optimize_geometry(system, calculator; solver=solver, optim_options...)
    @test isapprox(results.u[1], 0; atol=1e-5)
end
