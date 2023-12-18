# Create a dummy AtomsCalculator to test the geometry optimization interfcae.
@testsetup module TestCalculators
    using AtomsCalculators
    using Unitful
    using UnitfulAtomic
    
    struct DummyCalculator
    end
    
    AtomsCalculators.@generate_interface function AtomsCalculators.potential_energy(system, calculator::DummyCalculator; kwargs...)
        return 0.0u"eV"
    end
        
    AtomsCalculators.@generate_interface function AtomsCalculators.forces(system, calculator::DummyCalculator; kwargs...)
        return AtomsCalculators.zero_forces(system, calculator)
    end
end
