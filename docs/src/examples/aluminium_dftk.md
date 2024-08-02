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
scf_kwargs   = (; tol=1e-5, mixing=KerkerMixing())
calc = DFTKCalculator(; model_kwargs, basis_kwargs, scf_kwargs, verbose=true)
```

We attach pseudopotentials to the aluminium system,
i.e. we tell DFTK, that each aluminium atom should be modelled using
a pseudopotential rather than the full Coulomb potential.

```@example dftk-aluminium
system = attach_psp(system; Al="hgh/lda/al-q3");
```

!!! info "Crude computational parameters"
    Note, that these numerical parameters are chosen rather crudely in order
    to give a fast runtime on CI systems. For production calculations one would
    require larger computational parameters.

We perform the structure optimisation using the LBFGS solver
from Optim with solver parameters adapted for our geometry optimisation setting.
This is selected by passing the [GeometryOptimization.OptimLBFGS](@ref)
solver as the third argument.

```@example dftk-aluminium
using GeometryOptimization
GO = GeometryOptimization

results = minimize_energy!(system, calc, GO.OptimLBFGS();
                           tol_force=1e-4u"eV/Å",
                           show_trace=true);
```

The final structure is
```@example dftk-aluminium
system
```
