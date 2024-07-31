# Geometry optimisation of a hydrogen molecule
using AtomsBase
using EmpiricalPotentials
using GeometryOptimization
using Unitful
using UnitfulAtomic

bounding_box = 10.0u"angstrom" .* [[1, 0, 0.], [0., 1, 0], [0., 0, 1]]
atoms = [:H => [0, 0, 0.0]u"bohr", :H => [0, 0, 1.9]u"bohr"]
system = periodic_system(atoms, bounding_box)

# Set everything to optimizable.
# system = clamp_atoms(system, [1])

# Create a simple DFT calculator
model_kwargs = (; functionals = [:lda_x, :lda_c_pw])
basis_kwargs = (; kgrid = [1, 1, 1], Ecut = 10.0)
scf_kwargs   = (; tol = 1e-7)
lda = DFTKCalculator(system; model_kwargs, basis_kwargs, scf_kwargs)

results = minimize_enery!(system, lda)
println(results)
optsystem = results.system
println("Bond length: $(norm(position(optsystem[1]) - position(optsystem[2])))")
