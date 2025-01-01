# TiAl structure optimisation

TODO Write some text motivating this example

Setup system:
```@example tial
using AtomsIO
using EmpiricalPotentials

system = load_system(joinpath(pkgdir(EmpiricalPotentials), "data/TiAl-1024.xyz"))
nothing  # hide
```

Setup calculator:
```@example tial
using Unitful
using UnitfulAtomic

# Note: These are completely made up parameters,
#       please do not use in production
rcut = 5.0u"Å"
zAl = atomic_number(ChemicalSpecies(:Al))
zTi = atomic_number(ChemicalSpecies(:Ti))
emins = Dict( (zAl, zAl) => -1.0u"eV",
              (zAl, zTi) => -1.234u"eV",
              (zTi, zTi) => -0.345u"eV" )
rmins = Dict( (zAl, zAl) => 2.7u"Å",
              (zAl, zTi) => 3.2u"Å",
              (zTi, zTi) => 3.0u"Å" )
calc = LennardJones(emins, rmins, rcut)
nothing  # hide
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
