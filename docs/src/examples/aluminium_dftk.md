# Aluminium supercell using density-functional theory

In this example we will optimise the structure of a rattled
aluminium system using density-functional theory.

First we build a rattled aluminium system:

```@example dftk-aluminium
using AtomsBuilder
using Unitful

system = rattle!(bulk(:Al; cubic=true), 0.2u"Å")
```

Next we create a calculator employing the
[density-functional toolkit](https://dftk.org/)
to compute energies and forces at using the LDA density functional.
```@example dftk-aluminium
using DFTK

model_kwargs = (; functionals=[:lda_x, :lda_c_pw], temperature=1e-3)
basis_kwargs = (; kgrid=(3, 3, 3), Ecut=10.0)
scf_kwargs   = (; mixing=KerkerMixing())
calc = DFTKCalculator(; model_kwargs, basis_kwargs, scf_kwargs)
nothing
```

We attach pseudopotentials to the aluminium system,
i.e. we tell DFTK, that each aluminium atom should be modelled using
a pseudopotential rather than the full Coulomb potential.

```@example dftk-aluminium
system = attach_psp(system; Al="hgh/lda/al-q3")
nothing
```

!!! info "Crude computational parameters"
    Note, that these numerical parameters are chosen rather crudely in order
    to give a fast runtime on CI systems. For production calculations one would
    require larger computational parameters.

We perform the structure optimisation using the LBFGS solver
from Optim with solver parameters adapted for our geometry optimisation setting.
This is selected by passing the [GeometryOptimization.OptimLBFGS](@ref)
solver as the third argument. The `verbosity=2` flag makes sure we get
output from both the geometry optimisation as well as the inner SCF solver.

```@example dftk-aluminium
using GeometryOptimization
GO = GeometryOptimization

results = minimize_energy!(system, calc, GO.OptimLBFGS();
                           tol_force=1e-4u"eV/Å", verbosity=2)
nothing
```

!!! tip "Automatically adapted calculator parameters"
    Some calculators (such as DFTK) are able to adapt to the keyword arguments
    and parameters passed to `minimize_energy!`. In this case the SCF tolerance
    is automatically adapted according to the convergence parameters
    (here `tol_force`) passed to `minimize_energy!`.

The final energy is
```@example dftk-aluminium
results.energy
```

We can view the final structure
```@example dftk-aluminium
results.system
```

Some statistics about the optimisation
```@example dftk-aluminium
results.stats
```
or the details about the selected algorithm:
```@example dftk-aluminium
results.alg
```

The final state of the calculator object is also accessible
via `results.state` and could be employed for postprocessing
using the framework of the calculator. E.g. in the case
of `DFTK`, the `results.state` is what `DFTK` calls an `scfres`
and could just be used to plot a density of states or plot
bands or compute response properties.
