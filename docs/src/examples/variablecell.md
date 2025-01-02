# Silicon variable cell relaxation

In this example we optimise both atomic positions and lattice geometry
of a rattled silicon structure using a Stillinger-Weber potential.

We build an initial structure by first rattling the positions
and than the lattice geometry:

```@example silicon
using AtomsBase
using AtomsBuilder
using Unitful
using UnitfulAtomic

silicon_posrattle = rattle!(bulk(:Si, cubic=true) * (2, 2, 2), 0.1u"Å")
F = I + 1e-3randn(3, 3)
new_cell_vectors = tuple([F * v for v in cell_vectors(silicon_posrattle)]...)
silicon_rattle = AbstractSystem(silicon_posrattle; cell_vectors=new_cell_vectors)
```

and optimise with `variablecell=true`

```@example silicon
using EmpiricalPotentials
using GeometryOptimization

sw = StillingerWeber()
silicon = minimize_energy!(silicon_posrattle, sw; variablecell=true).system

```
