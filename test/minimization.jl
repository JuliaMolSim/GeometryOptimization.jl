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


@testitem "Test aluminium LJ fixed cell minimisation" begin
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

    rcut  = 6.0u"Å"
    zAl   = atomic_number(ChemicalSpecies(:Al))
    emins = Dict( (zAl, zAl) => -1.0u"meV", )
    rmins = Dict( (zAl, zAl) => 3.1u"Å",    )
    calculator = LennardJones(emins, rmins, rcut)

    # Get a minimised reference
    aluminium_final = minimize_energy!(bulk(:Al; cubic=true), calculator;
                                       tol_forces=1e-10, maxiters=200).system
    aluminium_init  = rattle!(bulk(:Al; cubic=true), 0.2u"Å")
    energy_final = AC.potential_energy(aluminium_final, calculator)
    energy_init  = AC.potential_energy(aluminium_init,  calculator)
    @assert energy_final ≤ energy_init

    for solver in (GO.Autoselect(), GO.OptimCG(), GO.OptimLBFGS())
        # TODO OptimSD is *really slow* here
        results = minimize_energy!(aluminium_init, calculator, solver;
                                   tol_forces=1e-8, maxiters=500, verbosity=1)
        @test results.energy ≤ energy_init
        @test austrip(abs(results.energy - energy_final)) < 1e-12

        delta = position(results.system, 1) - position(aluminium_final, 1)
        for i in 2:length(aluminium_final)
            @test norm(  position(results.system, i) - delta
                       - position(aluminium_final, i)) < 1e-5u"bohr"
        end
    end
end
