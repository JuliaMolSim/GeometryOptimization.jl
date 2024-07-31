# Geometry optimisation of a hydrogen molecule
using AtomsBase
using EmpiricalPotentials
using GeometryOptimization
using Unitful
using UnitfulAtomic

system = isolated_system([:H => [0, 0, 0.0]u"bohr", :H => [0, 0, 1.9]u"bohr"])
lj = LennardJones(-1.17u"hartree", 0.743u"angstrom", 1, 1, 0.6u"nm")
results = minimize_enery!(system, lj)
optsystem = results.system
println("Bond length: $(norm(position(optsystem[1]) - position(optsystem[2])))")
