using LinearAlgebra
using DFTK
using ASEconvert
using LazyArtifacts
using AtomsBase
using AtomsCalculators
using Unitful
using UnitfulAtomic
using Random
using OptimizationOptimJL
using ComponentArrays

using GeometryOptimization


function build_al_supercell(rep=1)
    pseudodojo_psp = artifact"pd_nc_sr_lda_standard_0.4.1_upf/Al.upf"
    a = 7.65339 # true lattice constant.
    lattice = a * Matrix(I, 3, 3)
    Al = ElementPsp(:Al; psp=load_psp(pseudodojo_psp))
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
    supercell = attach_psp(supercell; Al=pseudodojo_psp)
    return supercell
end;

al_supercell = build_al_supercell(1)

# Create a simple calculator for the model.
model_kwargs = (; functionals = [:lda_x, :lda_c_pw], temperature = 1e-2)
basis_kwargs = (; kgrid = [2, 2, 2], Ecut = 10.0)
scf_kwargs = (; tol = 1e-3)
calculator = DFTKCalculator(; model_kwargs, basis_kwargs, scf_kwargs, verbose=true)

optim_options = (f_tol=1e-6, iterations=6, show_trace=true)
results = minimize_energy!(al_supercell, calculator; procedure="vc_relax", optim_options...)
