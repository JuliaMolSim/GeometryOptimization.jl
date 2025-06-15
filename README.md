# GeometryOptimization

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaMolSim.github.io/GeometryOptimization.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaMolSim.github.io/GeometryOptimization.jl/dev/)
[![Build Status](https://github.com/JuliaMolSim/GeometryOptimization.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaMolSim/GeometryOptimization.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaMolSim/GeometryOptimization.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaMolSim/GeometryOptimization.jl)

A package for optimising the structural parameters of an atomistic system,
i.e. the step usually referred to as
**geometry optimisation** or **structural relaxation**
in electronic structure theory and atomistic modelling.
Both relaxing **atomic positions** as well as the **unit cell** is supported.

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
using LinearAlgebra
using Unitful
using UnitfulAtomic

# Setup system and calculator
cell_vectors = ([10.0, 0.0, 0.0]u"Å", [0.0, 10.0, 0.0]u"Å", [0.0, 0.0, 10.0]u"Å")
system = periodic_system([:H => [0, 0, 0.0]u"bohr",
                          :H => [0, 0, 1.9]u"bohr"],
                         cell_vectors)
zH = 1
emins = Dict((zH, zH) => -1.17u"hartree", )
rmins = Dict((zH, zH) =>  0.743u"Å",      )
calc = LennardJones(emins, rmins, 5.0u"Å")

# Run the geometry optimisation
results = minimize_energy!(system, calc)

# Inspect the results
optsystem = results.system
optimised_bondlength = norm(position(optsystem[1]) - position(optsystem[2]))
```
