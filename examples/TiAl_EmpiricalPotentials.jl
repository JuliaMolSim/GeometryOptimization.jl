using AtomsBase
using AtomsIO
using EmpiricalPotentials
using GeometryOptimization
using Unitful
using UnitfulAtomic
GO = GeometryOptimization

system = load_system("/home/mfh/.julia/packages/EmpiricalPotentials/APr5T/data/TiAl-1024.xyz")

# Convert to AbstractSystem, so we have the `particles` attribute.
particles = map(system) do atom
    Atom(; pairs(atom)...)
end
system = AbstractSystem(data; particles)

lj = LennardJones(-1.0u"meV", 3.1u"Å", 13, 13, 6.0u"Å")
results = minimize_energy!(system, lj; GO.AutoLBFGS(), maxiters=10, show_trace=true)
println(results)
