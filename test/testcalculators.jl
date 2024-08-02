# Create a dummy AtomsCalculator to test the geometry optimization interfcae.
@testsetup module TestCalculators
    using AtomsCalculators
    using Unitful
    using UnitfulAtomic
    AC = AtomsCalculators

    struct DummyCalc end
    AC.energy_unit(::DummyCalc) = u"hartree"
    AC.length_unit(::DummyCalc) = u"bohr"

    AC.@generate_interface function AC.potential_energy(system, calc::DummyCalc; kwargs...)
        AC.zero_energy(calc)
    end

    AC.@generate_interface function AC.forces(system, calc::DummyCalc; kwargs...)
        AC.zero_forces(system, calc)
    end
end
