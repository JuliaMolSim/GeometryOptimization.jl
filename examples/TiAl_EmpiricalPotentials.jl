using AtomsBase
using AtomsCalculators
using EmpiricalPotentials
using ExtXYZ
using Unitful
using UnitfulAtomic
using OptimizationOptimJL

using GeometryOptimization

fname = joinpath(pkgdir(EmpiricalPotentials), "data/", "TiAl-1024.xyz")
data = ExtXYZ.load(fname) |> FastSystem

lj = LennardJones(-1.0u"meV", 3.1u"Å",  13, 13, 6.0u"Å")

# Convert to AbstractSystem, so we have the `particles` attribute.
particles = map(data) do atom
    Atom(; pairs(atom)...)
end
system = AbstractSystem(data; particles)

solver = OptimizationOptimJL.LBFGS()
optim_options = (f_tol=1e-8, g_tol=1e-8, iterations=10, show_trace=true)

results = minimize_energy!(system, lj; solver, optim_options...)
println(results)
