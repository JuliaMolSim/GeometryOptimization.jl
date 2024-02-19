# Verify correctness of lattice optimization on silicon.
using DFTK
using AtomsBase
using AtomsCalculators
using ComponentArrays
using LazyArtifacts
using Random
using Unitful
using UnitfulAtomic
using OptimizationOptimJL
using GeometryOptimization
    
    
lattice = [0.0  5.131570667152971 5.131570667152971;
           5.131570667152971 0.0 5.131570667152971;
           5.131570667152971 5.131570667152971  0.0]
hgh_lda_family = artifact"hgh_lda_hgh"
psp_hgh = joinpath(hgh_lda_family, "si-q4.hgh")

positions = [ones(3)/8, -ones(3)/8]
atoms = fill(ElementPsp(:Si; psp=load_psp(psp_hgh)), 2)
system = periodic_system(lattice, atoms, positions)

# Create a simple calculator for the model.
model_kwargs = (; functionals = [:lda_x, :lda_c_pw], temperature = 1e-5)
basis_kwargs = (; kgrid = [5, 5, 5], Ecut = 30.0)
scf_kwargs = (; tol = 1e-6)
calculator = DFTKCalculator(; model_kwargs, basis_kwargs, scf_kwargs, verbose=true)

# Perturb unit cell.
Random.seed!(1234)
σ = 0.2u"bohr"
bounding_box_pert = [v + σ * rand(Float64, size(v)) for v in bounding_box(system)]
system_pert = update_positions(system, position(system); bounding_box=bounding_box_pert)


using LineSearches
linesearch =  BackTracking(c_1= 1e-4, ρ_hi= 0.8, ρ_lo= 0.1, order=2, maxstep=Inf)
solver = OptimizationOptimJL.LBFGS(; linesearch)
optim_options = (; solver, f_tol=1e-10, g_tol=1e-5, iterations=30,
                 show_trace=true, store_trace = true, allow_f_increases=true)

results = minimize_energy!(system_pert, calculator; procedure="vc_relax", optim_options...)
