@testitem "Test NLopt solver via Optimization.jl interface" begin
    using AtomsBase
    using AtomsBuilder
    using AtomsCalculators
    using EmpiricalPotentials
    using GeometryOptimization
    using OptimizationNLopt
    using Unitful
    using UnitfulAtomic
    AC = AtomsCalculators
    calculator = StillingerWeber()
    solver = NLopt.LD_TNEWTON()

    # Get a minimised reference
    reference = minimize_energy!(bulk(:Si; cubic=true), calculator;
                                 tol_forces=1e-10, maxiters=200)

    # Equilibrate a perturbed Si crystal using a solver from Optimization.jl
    silicon_init = rattle!(bulk(:Si; cubic=true), 0.2u"Å")
    energy_init  = AC.potential_energy(silicon_init,  calculator)
    results = minimize_energy!(silicon_init, calculator, solver;
                               tol_energy=1e-10, tol_forces=1e-4u"eV/Å", maxeval=100)
    @test results.energy ≤ energy_init
    @test austrip(abs(results.energy - reference.energy)) < 1e-10
end
