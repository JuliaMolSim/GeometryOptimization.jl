using DFTK: occupation
#= Test Geometry Optimization on an aluminium supercell.
=#
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
model_kwargs = (; functionals = [:lda_x, :lda_c_pw], temperature = 1e-4)
basis_kwargs = (; kgrid = [6, 6, 6], Ecut = 30.0)
scf_kwargs = (; tol = 1e-6)
calculator = DFTKCalculator(; model_kwargs, basis_kwargs, scf_kwargs, verbose=true)

energy_true = AtomsCalculators.potential_energy(al_supercell, calculator)

# Starting position is a random perturbation of the equilibrium one.
Random.seed!(1234)
x0 = vcat(position(al_supercell)...)
σ = 0.5u"angstrom"; x0_pert = x0 + σ * rand(Float64, size(x0))
al_supercell = update_not_clamped_positions(al_supercell, x0_pert)
energy_pert = AtomsCalculators.potential_energy(al_supercell, calculator)

println("Initial guess distance (norm) from true parameters $(norm(x0 - x0_pert)).")
println("Initial regret $(energy_pert - energy_true).")

optim_options = (f_tol=1e-6, iterations=6, show_trace=true)

results = minimize_energy!(al_supercell, calculator; optim_options...)
println(results)

""" Returns the index of the highest occupied band. 
    `atol` specifies the (fractional) occupations below which 
    a band is considered unoccupied.
"""
function valence_band_index(occupations; atol=1e-36)
    filter = x -> isapprox(x, 0.; atol)
    maximum(maximum.(findall.(!filter, occupations)))
end

function band_gaps(scfres)
    vi = valence_band_index(occupations; atol)
    # If the system is metallic, by convenction band gap is zero.
    if DFTK.is_metal(scfres.eigenvalues, scfres.εF; tol=1e-4)
        return (; εMax_valence=0.0u"hartree", εMin_conduction=0.0u"hartree",
                direct_bandgap=0.0u"hartree", valence_band_index=vi)
    else
        εMax_valence = maximum([εk[vi] for εk in scfres.eigenvalues]) * u"hartree"
        εMin_conduction = minimum([εk[vi + 1] for εk in scfres.eigenvalues]) * u"hartree"
        direct_bandgap = minimum([εk[vi + 1] - εk[vi] for εk in scfres.eigenvalues]) * u"hartree"
        return (; εMax_valence, εMin_conduction, direct_bandgap, valence_band_index=vi)
    end
end
