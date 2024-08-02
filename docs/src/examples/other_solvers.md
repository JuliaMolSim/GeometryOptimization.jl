# Using a different Optimization.jl compatible solver

In this example we perform the simplistic optimisation
the bond length of a Hydrogen molecule using a trust region
quasi-Newton method from
[NLopt](https://github.com/JuliaOpt/NLopt.jl).

We create a calculator employing the
[density-functional toolkit](https://dftk.org/)
to compute energies and forces at using the LDA density functional.

```@example other-solvers
using DFTK

model_kwargs = (; functionals=[:lda_x, :lda_c_pw])
basis_kwargs = (; kgrid=(1, 1, 1), Ecut=20.0)
scf_kwargs   = (; tol=1e-6)
calc = DFTKCalculator(; model_kwargs, basis_kwargs, scf_kwargs, verbose=true)
nothing
```

and we build the hydrogen molecular system,
where we attach pseudopotential information for DFTK:

```@example other-solvers
using AtomsBuilder
using Unitful
using UnitfulAtomic

bounding_box = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]u"Å"
system = periodic_system([:H => [0, 0, 1.]u"bohr",
                          :H => [0, 0, 3.]u"bohr"],
                         bounding_box)
system = attach_psp(system; H="hgh/lda/h-q1")
nothing
```

We now run [`GeometryOptimization.minimize_energy!`](@ref), but notably
pass the `NLopt.LD_TNEWTON` solver from `NLopt` as the
third argument to employ this solver. Extra keyword argument to `NLopt`
can be added, e.g. here the `maxevel=100`, which limits the solver to
100 function evaluations:
```@example other-solvers
using GeometryOptimization
using OptimizationNLopt
solver = NLopt.LD_TNEWTON()

results = minimize_energy!(system, calc, solver;
                           tol_force=1e-4u"eV/Å", maxeval=100)
nothing
```

The final hydrogen bond length is:

```@example other-solvers
using LinearAlgebra
norm(position(results.system[1]) - position(results.system[2]))
```
