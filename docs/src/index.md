```@meta
CurrentModule = GeometryOptimization
```

# GeometryOptimization

A geometry optimization package for
[AtomsBase structures](https://github.com/JuliaMolSim/AtomsBase.jl)
and [AtomsCalculator calculators](https://github.com/JuliaMolSim/AtomsCalculators.jl).
The source code can be found [on github](https://github.com/JuliaMolSim/GeometryOptimization.jl).

## Motivating example

We consider the optimisation of the bondlength of a hydrogen
molecule using a simple Lennard Jones potential:

```julia
## TODO: Should run as @example once EmpiricalPotentials is compatible
using AtomsBase
using EmpiricalPotentials
using GeometryOptimization
using Unitful
using UnitfulAtomic

# Setup system and calculator
system = isolated_system([:H => [0, 0, 0.0]u"bohr",
                          :H => [0, 0, 1.9]u"bohr"])
calc = LennardJones(-1.17u"hartree", 0.743u"angstrom", 1, 1, 0.6u"nm")

# Run the geometry optimisation (using verbosity=1 to print the progress)
results = minimize_energy!(system, calc; verbosity=1)

# Inspect the results
optimised_system = results.system
optimised_bondlength = norm(position(optsystem[1]) - position(optsystem[2]))
```

The idea is that
any [AtomsBase](https://github.com/JuliaMolSim/AtomsBase.jl)-compatible
structure can be employed as a `system`
any [AtomsCalculators](https://github.com/JuliaMolSim/AtomsCalculators.jl)-compatible
calculator as a `calc`. See the list of examples to get an overview of possible
calculators to employ.

Note that [`minimize_energy!`](@ref) supports further arguments to fine-tune the
convergence tolerance or customise the solver selection.
See the API documentation or the examples to explore this.
