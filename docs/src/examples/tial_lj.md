# TiAl structure optimisation

TODO Write some text motivating this example

Setup system:
```@example tial
using AtomsBase
using AtomsIO
using EmpiricalPotentials
using GeometryOptimization
GO = GeometryOptimization

system = load_system(joinpath(pkgdir(EmpiricalPotentials), "data/TiAl-1024.xyz")

# Convert to AbstractSystem, so we have the `particles` attribute.
particles = map(system) do atom
    Atom(; pairs(atom)...)
end
system = AbstractSystem(data; particles)
```

Setup calculator:
```@example tial
using Unitful
using UnitfulAtomic
calc = LennardJones(-1.0u"meV", 3.1u"Å", 13, 13, 6.0u"Å")
```

Minimise energy:
```@example tial
results = minimize_energy!(system, calc, GO.OptimCG(); maxiters=10, show_trace=true)
results.energy
```

Final structure:
```@example tial
results.system
```
