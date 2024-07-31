# Geometry Optimization of an aluminium supercell

using AtomsBase
using AtomsBuilder
using AtomsCalculators
using DFTK
using GeometryOptimization
using LazyArtifacts
using Unitful
AC = AtomsCalculators
GO = GeometryOptimization
setup_threading(n_blas=1)

# Create a simple LDA-based calculator
model_kwargs = (; functionals = [:lda_x, :lda_c_pw], temperature = 1e-3)
basis_kwargs = (; kgrid = [6, 6, 6], Ecut = 30.0)
basis_kwargs = (; kgrid = [4, 4, 4], Ecut = 10.0)
scf_kwargs   = (; tol=1e-7)
calc = DFTKCalculator(; model_kwargs, basis_kwargs, scf_kwargs, verbose=true)

# Make an aluminium system and compute its energy
system = bulk(:Al, cubic=true)
system = attach_psp(system;
                    Al=artifact"pd_nc_sr_lda_standard_0.4.1_upf/Al.upf")
energy_reference = AC.potential_energy(system, calc)

# Rattle the system
rattled = attach_psp(rattle!(AbstractSystem(system), 0.3u"Å");
                    Al=artifact"pd_nc_sr_lda_standard_0.4.1_upf/Al.upf")

energy_pert = AC.potential_energy(rattled, calc)
println("Initial guess distance    ",
        maximum(norm, position(rattled) - position(system)))
println("Initial energy difference $(energy_pert - energy_reference).")

END
results = minimize_energy!(rattled, calc, GO.AutoLBFGS(); show_trace=true)
println(results)
