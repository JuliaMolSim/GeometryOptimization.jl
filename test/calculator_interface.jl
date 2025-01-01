@testitem "Test GeometryOptimization with AtomsCalculators interface" begin
    # This is a fake test. It does not really test anything,
    # just that the interface code compiles and works properly.

    using AtomsBuilder
    using AtomsCalculators
    using GeometryOptimization
    using LinearAlgebra
    using Unitful
    using UnitfulAtomic
    AC = AtomsCalculators
    GO = GeometryOptimization

    struct DummyCalc end
    AC.energy_unit(::DummyCalc) = u"hartree"
    AC.length_unit(::DummyCalc) = u"bohr"
    AC.@generate_interface function AC.potential_energy(system, calc::DummyCalc; kwargs...)
        AC.zero_energy(system, calc)
    end
    AC.@generate_interface function AC.forces(system, calc::DummyCalc; kwargs...)
        AC.zero_forces(system, calc)
    end

    system = GO.clamp_atoms(bulk(:Al; cubic=true), [1])
    for solver in (GO.Autoselect(), GO.OptimCG(), GO.OptimLBFGS(), GO.OptimSD())
        results = minimize_energy!(system, DummyCalc(), solver; verbosity=1)
        @test norm(position(results.system[1]) - position(system[1])) < 1e-5u"bohr"
    end
end
