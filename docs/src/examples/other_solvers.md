# Using a Optimization.jl compatible solver

In this example we optimize a bulk silicon structure
using a trust region quasi-Newton method from
[NLopt](https://github.com/JuliaOpt/NLopt.jl).

We create a Stillinger-Weber calculator

```@example other-solvers
using EmpiricalPotentials
sw = StillingerWeber()
nothing
```

and we build a slightly rattled silicon structure

```@example other-solvers
using AtomsBuilder
using Unitful
system = rattle!(bulk(:Si, cubic=true) * (2, 2, 2), 0.1u"Å")
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

results = minimize_energy!(system, sw, solver;
                           tol_forces=1e-4u"eV/Å", verbosity=1,
                           maxeval=100)
nothing
```

!!! info ""
    While in principle all *first-order* solvers supported
    by Optimization.jl can be employed, only few have been tested with this
    package so far. Given the complexity of the Optimization.jl ecosystem
    we expect that minor changes of our integration may be needed to make
    specific solvers work well. We welcome any PRs.
