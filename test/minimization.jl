@testitem "Test silicon StillingerWeber fixed cell minimisation" begin
    using AtomsBase
    using AtomsBuilder
    using AtomsCalculators
    using EmpiricalPotentials
    using GeometryOptimization
    using LinearAlgebra
    using Unitful
    using UnitfulAtomic
    AC = AtomsCalculators
    GO = GeometryOptimization
    calculator = StillingerWeber()

    # Get a minimised reference
    silicon_final = minimize_energy!(bulk(:Si; cubic=true), calculator;
                                     tol_forces=1e-10, maxiters=200).system
    silicon_init = rattle!(bulk(:Si; cubic=true), 0.2u"Å")
    energy_final = AC.potential_energy(silicon_final, calculator)
    energy_init  = AC.potential_energy(silicon_init,  calculator)
    @assert energy_final ≤ energy_init

    for solver in (GO.Autoselect(), GO.OptimCG(), GO.OptimLBFGS(), GO.OptimSD())
        results = minimize_energy!(silicon_init, calculator, solver;
                                   tol_forces=1e-8, maxiters=500, verbosity=1)
        @test results.energy ≤ energy_init
        @test austrip(abs(results.energy - energy_final)) < 1e-12

        delta = position(results.system, 1) - position(silicon_final, 1)
        for i in 2:length(silicon_final)
            @test norm(  position(results.system, i) - delta
                       - position(silicon_final, i)) < 1e-5u"bohr"
        end
    end
end


@testitem "Test silicon StillingerWeber variable cell minimisation" begin
    using AtomsBase
    using AtomsBuilder
    using AtomsCalculators
    using EmpiricalPotentials
    using GeometryOptimization
    using LinearAlgebra
    using Unitful
    using UnitfulAtomic
    AC = AtomsCalculators
    GO = GeometryOptimization
    calculator = StillingerWeber()

    # Get a minimised reference
    silicon_rattle = rattle!(bulk(:Si, cubic=true) * (2, 2, 1), 0.1u"Å")
    silicon_final  = minimize_energy!(silicon_rattle, calculator;
                                      variablecell=true, maxiters=200,
                                      tol_virial=1e-8u"hartree", tol_forces=1e-10).system

    # Rattle again a little
    silicon_rattle = rattle!(AbstractSystem(silicon_final), 0.1u"Å")
    F = I + 1e-3randn(3, 3)
    new_cell_vectors = tuple([F * v for v in cell_vectors(silicon_rattle)]...)
    silicon_init = AbstractSystem(silicon_rattle; cell_vectors=new_cell_vectors)

    energy_final = AC.potential_energy(silicon_final, calculator)
    energy_init  = AC.potential_energy(silicon_init,  calculator)
    @assert energy_final ≤ energy_init

    println()
    println()
    println()

    for solver in (GO.Autoselect(), GO.OptimCG(), GO.OptimLBFGS())
        results = minimize_energy!(silicon_init, calculator, solver;
                                   variablecell=true, maxiters=500, verbosity=1,
                                   tol_virial=1e-7u"hartree", tol_forces=1e-8)
        @test maximum(x -> austrip(maximum(abs, x)),
                      AC.virial(results.system, calculator)) < 1e-6
        @test maximum(x -> austrip(maximum(abs, x)),
                      AC.forces(results.system, calculator)) < 1e-7

        @test results.energy ≤ energy_init
        @test austrip(abs(results.energy - energy_final)) < 1e-12
    end
end
