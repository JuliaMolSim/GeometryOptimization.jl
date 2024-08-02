@testitem "Test GeometryOptimization with AtomsCalculators interface" #=
        =# setup=[TestCalculators] begin
    # This is a fake test. It does not really test anything,
    # just that the interface code compiles and works properly.

    using GeometryOptimization
    using AtomsBuilder
    GO = GeometryOptimization

    system = bulk(:Al)
    system = clamp_atoms(system, [1])

    calculator = TestCalculators.DummyCalc()
    for solver in (GO.Autoselect(), GO.OptimCG(), GO.OptimLBFGS(), GO.OptimSD())
        results = minimize_energy!(system, calculator, solver)
        @test norm(position(results.system[1]) - position(system[1])) < 1e-5u"bohr"
    end
end
