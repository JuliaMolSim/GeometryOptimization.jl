# Using a Optimization.jl compatible solver

In this example we perform the simplistic optimisation
the bond length of a Hydrogen molecule using a trust region
quasi-Newton method from
[NLopt](https://github.com/JuliaOpt/NLopt.jl).

We create a calculator employing the
[density-functional toolkit](https://dftk.org/)
to compute energies and forces at using the LDA density functional.

```@example other-solvers
using DFTK
using PseudoPotentialData

pseudopotentials = PseudoFamily("dojo.nc.sr.lda.v0_4_1.oncvpsp3.standard.upf")
model_kwargs = (; functionals=LDA(), pseudopotentials)
basis_kwargs = (; kgrid=(1, 1, 1), Ecut=20.0)
calc = DFTKCalculator(; model_kwargs, basis_kwargs)
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
                           tol_forces=1e-4u"eV/Å", verbosity=1,
                           maxeval=100)
nothing
```

The final hydrogen bond length is:

```@example other-solvers
using LinearAlgebra
norm(position(results.system[1]) - position(results.system[2]))
```

!!! info ""
    While in principle all *first-order* solvers supported
    by Optimization.jl can be employed, only few have been tested with this
    package so far. Given the complexity of the Optimization.jl ecosystem
    we expect that minor changes of our integration may be needed to make
    specific solvers work well. We welcome any PRs.
