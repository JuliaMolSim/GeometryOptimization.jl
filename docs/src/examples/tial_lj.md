# TiAl structure optimisation

TODO Write some text motivating this example

Setup system:
```@example tial
using AtomsIO
using EmpiricalPotentials

system = load_system(joinpath(pkgdir(EmpiricalPotentials), "data/TiAl-1024.xyz"))
nothing
```

Setup calculator:
```@example tial
using Unitful
using UnitfulAtomic

calc = LennardJones(-1.0u"meV", 3.1u"Å", 13, 13, 6.0u"Å")
nothing
```

Minimise energy:
```@example tial
using GeometryOptimization
GO = GeometryOptimization

results = minimize_energy!(system, calc, GO.OptimCG(); maxiters=10, verbosity=1)
results.energy
```

Final structure:
```@example tial
results.system
```
