#= Test Geometry Optimization on an aluminium supercell.
=#
using LinearAlgebra
using EmpiricalPotentials
using Unitful
using UnitfulAtomic
using OptimizationOptimJL
using AtomsBase

using GeometryOptimization


bounding_box = 10.0u"angstrom" .* [[1, 0, 0.], [0., 1, 0], [0., 0, 1]]
atoms = [:H => [0, 0, 0.0]u"bohr", :H => [0, 0, 1.9]u"bohr"]
system = periodic_system(atoms, bounding_box)

lj = LennardJones(-1.17u"hartree", 0.743u"angstrom", 1, 1, 0.6u"nm")

solver = OptimizationOptimJL.LBFGS()
optim_options = (f_tol=1e-6, iterations=100, show_trace=false)

results = minimize_energy!(system, lj; solver=solver, optim_options...)
println("Bond length: $(norm(results.minimizer[1:3] - results.minimizer[4:end])).")
