#= Test Geometry Optimization on an aluminium supercell.
=#
using Printf
using LinearAlgebra
using DFTK
using ASEconvert
using LazyArtifacts
using AtomsCalculators
using Unitful
using UnitfulAtomic
using Random
using OptimizationOptimJL

using GeometryOptimization

# Get PseudoDojo pseudopotential.
psp_upf  = load_psp(artifact"pd_nc_sr_lda_standard_0.4.1_upf/Al.upf");

function build_al_supercell(rep=1)
    a = 7.65339 # true lattice constant.
    lattice = a * Matrix(I, 3, 3)
    Al = ElementPsp(:Al; psp=psp_upf)
    atoms     = [Al, Al, Al, Al]
    positions = [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
    unit_cell = periodic_system(lattice, atoms, positions)

    # Make supercell in ASE:
    # We convert our lattice to the conventions used in ASE, make the supercell
    # and then convert back ...
    supercell_ase = convert_ase(unit_cell) * pytuple((rep, 1, 1))
    supercell     = pyconvert(AbstractSystem, supercell_ase)

    # Unfortunately right now the conversion to ASE drops the pseudopotential information,
    # so we need to reattach it:
    supercell = attach_psp(supercell; Al=artifact"pd_nc_sr_lda_standard_0.4.1_upf/Al.upf")
    return supercell
end;

al_supercell = build_al_supercell(1)

# Create a simple calculator for the model.
model_kwargs = (; functionals = [:lda_x, :lda_c_pw], temperature = 1e-4)
basis_kwargs = (; kgrid = [8, 8, 8], Ecut = 30.0)
scf_kwargs = (; tol = 1e-6)
calculator = DFTKCalculator(al_supercell; model_kwargs, basis_kwargs, scf_kwargs, verbose_scf=true)

energy_true = AtomsCalculators.potential_energy(al_supercell, calculator)

# Starting position is a random perturbation of the equilibrium one.
Random.seed!(1234)
x0 = vcat(position(al_supercell)...)
σ = 0.5u"angstrom"; x0_pert = x0 + σ * rand(Float64, size(x0))
al_supercell = update_optimizable_coordinates(al_supercell, x0_pert)
energy_pert = AtomsCalculators.potential_energy(al_supercell, calculator)

@printf "Initial guess distance (norm) from true parameters %.3e bohrs.\n" austrip(norm(x0 - x0_pert))
@printf "Initial regret %.3e.\n" energy_pert - energy_true

optim_options = (f_tol=1e-6, iterations=6, show_trace=true)

results = optimize_geometry(al_supercell, calculator; optim_options...)
println(results)
