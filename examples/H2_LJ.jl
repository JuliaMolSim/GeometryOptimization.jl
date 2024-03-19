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
jk
optim_options = (; solver, f_tol=1e-10, g_tol=1e-5, iterations=30,
                 show_trace=true, store_trace = true, allow_f_increases=true)

results = minimize_energy!(system, calculator; solver, procedure="relax", optim_options...)
println(results)
@printf "Bond length: %3f bohrs.\n" norm(results.minimizer[1:end])
