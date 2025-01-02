@testitem "Test aluminium EMT fixed cell minimisation" begin
    using ASEconvert
    using AtomsBuilder
    using AtomsCalculators
    using GeometryOptimization
    using LinearAlgebra
    using PythonCall
    using Unitful
    using UnitfulAtomic
    AC = AtomsCalculators
    GO = GeometryOptimization

    ase_emt = pyimport("ase.calculators.emt")
    calculator = ASEcalculator(ase_emt.EMT())

    # Get a minimised reference
    aluminium_final = minimize_energy!(bulk(:Al; cubic=true), calculator;
                                       tol_forces=1e-10).system
    aluminium_init  = rattle!(bulk(:Al; cubic=true), 0.2u"Å")
    energy_final = AC.potential_energy(aluminium_final, calculator)
    energy_init  = AC.potential_energy(aluminium_init,  calculator)
    @assert energy_final ≤ energy_init

    for solver in (GO.Autoselect(), GO.OptimCG(), GO.OptimLBFGS(), GO.OptimSD())
        results = minimize_energy!(aluminium_init, calculator, solver;
                                   tol_forces=1e-8, maxiters=500)
        @test results.energy ≤ energy_init
        @test austrip(abs(results.energy - energy_final)) < 1e-12

        delta = position(results.system, 1) - position(aluminium_final, 1)
        for i in 2:length(aluminium_final)
            @test norm(  position(results.system, i) - delta
                       - position(aluminium_final, i)) < 1e-6u"bohr"
        end
    end
end


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
    silicon_init  = rattle!(bulk(:Si; cubic=true), 0.2u"Å")
    energy_final = AC.potential_energy(silicon_final, calculator)
    energy_init  = AC.potential_energy(silicon_init,  calculator)
    @assert energy_final ≤ energy_init

    for solver in (GO.Autoselect(), GO.OptimCG(), GO.OptimLBFGS())
        # TODO OptimSD is *really slow* here
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
