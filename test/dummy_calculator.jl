# Create a dummy AtomsCalculator to test the geometry optimization interfcae.
@testsetup module TestCalculators
    using AtomsCalculators
    using Unitful
    using UnitfulAtomic
    
    struct DummyCalculator
        state 
    end
    DummyCalculator() = DummyCalculator(nothing)
    
    AtomsCalculators.@generate_interface function AtomsCalculators.potential_energy(
            system, calculator::DummyCalculator; kwargs...)
        0.0u"eV"
    end
        
    AtomsCalculators.@generate_interface function AtomsCalculators.forces(
            system, calculator::DummyCalculator; kwargs...)
        AtomsCalculators.zero_forces(system, calculator)
    end
end
