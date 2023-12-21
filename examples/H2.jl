#= Test Geometry Optimization on an aluminium supercell.
=#
using Printf
using LinearAlgebra
using DFTK
using Unitful
using UnitfulAtomic
using OptimizationOptimJL

using GeometryOptimization


a = 10.                  # Big box around the atoms.
lattice = a * I(3)
H = ElementPsp(:H; psp=load_psp("hgh/lda/h-q1"));
atoms = [H, H];
positions = [[0, 0, 0], [0, 0, .16]]
system = periodic_system(lattice, atoms, positions)

# Set everything to optimizable.
system = clamp_atoms(system, [1])

# Create a simple calculator for the model.
model_kwargs = (; functionals = [:lda_x, :lda_c_pw])
basis_kwargs = (; kgrid = [1, 1, 1], Ecut = 10.0)
scf_kwargs = (; tol = 1e-7)
calculator = DFTKCalculator(system; model_kwargs, basis_kwargs, scf_kwargs)

solver = OptimizationOptimJL.LBFGS()
optim_options = (f_tol=1e-32, iterations=20, show_trace=true)

results = optimize_geometry(system, calculator; solver=solver, optim_options...)
println(results)
@printf "Bond length: %3f bohrs.\n" norm(results.minimizer[1:end])
