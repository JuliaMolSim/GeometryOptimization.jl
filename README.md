# GeometryOptimization

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaMolSim.github.io/GeometryOptimization.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaMolSim.github.io/GeometryOptimization.jl/dev/)
[![Build Status](https://github.com/JuliaMolSim/GeometryOptimization.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaMolSim/GeometryOptimization.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaMolSim/GeometryOptimization.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaMolSim/GeometryOptimization.jl)

A package for optimising the structural parameters of an atomistic system,
i.e. the step usually referred to as
**geometry optimisation** or **structural relaxation**
in electronic structure theory and atomistic modelling.

The package is generic in the datastructures used to represent the geometry,
the calculator used to evaluate energies and forces as well as the solver algorithm.
Generally all
[AtomsBase structures](https://github.com/JuliaMolSim/AtomsBase.jl)
and [AtomsCalculator calculators](https://github.com/JuliaMolSim/AtomsCalculators.jl)
should work out of the box.

See the [documentation](https://JuliaMolSim.github.io/GeometryOptimization.jl/stable/)
for examples and further details.

## Motivating example

We consider the optimisation of the bondlength of a hydrogen
molecule using a simple Lennard Jones potential:

```julia
using AtomsBase
using EmpiricalPotentials
using GeometryOptimization
using Unitful
using UnitfulAtomic

# Setup system and calculator
system = isolated_system([:H => [0, 0, 0.0]u"bohr",
                          :H => [0, 0, 1.9]u"bohr"])
calc = LennardJones(-1.17u"hartree", 0.743u"angstrom", 1, 1, 0.6u"nm")

# Run the geometry optimisation
results = minimize_energy!(system, calc)

# Inspect the results
optimised  = results.system
bondlength = norm(position(optimised[1]) - position(optimised[2]))
```
