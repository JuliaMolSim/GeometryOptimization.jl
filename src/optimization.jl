#
# Note that by default all particles in the system are assumed optimizable.
# We always work in cartesian coordinates.
#

export minimize_energy!

using AtomsCalculators: Energy, Forces, calculate

mutable struct GeometryOptimizationState
    calc_state  # Reference to the current calculator state
    energy      # Most recent energy value
    forces      # Most recent force value
    virial      # Most recent virial value
end
function GeometryOptimizationState(system, calculator)
    GeometryOptimizationState(AC.get_state(calculator),
                              AC.zero_energy(system, calculator),
                              AC.zero_forces(system, calculator),
                              AC.zero_virial(system, calculator))
end

"""
    Optimization.OptimizationProblem(system, calculator geoopt_state; kwargs...)

Turn `system`, `calculator` and `geoopt_state::GeometryOptimizationState`
into a SciML-compatible `OptimizationProblem`. Note that the `system` is not updated
automatically and that internally atomic units are used.
"""
function Optimization.OptimizationProblem(system, calculator, geoopt_state; kwargs...)
    mask = not_clamped_mask(system)  # mask is assumed not to change during optimisation

    # Use current system parameters as starting positions.
    # Some optimisers modify x0 in-place, so need a mutable type.
    x0 = austrip.(not_clamped_positions(system))

    f = function(x::AbstractVector{<:Real}, ps)
        new_system = update_not_clamped_positions(system, x * u"bohr")
        res = calculate(Energy(), new_system, calculator, ps, geoopt_state.calc_state)
        geoopt_state.calc_state = res.state
        geoopt_state.energy = res.energy
        austrip(res.energy), geoopt_state
    end
    g! = function(G::AbstractVector{<:Real}, x::AbstractVector{<:Real}, ps)
        new_system = update_not_clamped_positions(system, x * u"bohr")
        res = calculate(Forces(), new_system, calculator, ps, geoopt_state.calc_state)
        geoopt_state.calc_state = res.state
        geoopt_state.forces .= res.forces

        # Translate the forces vectors on each particle to a single gradient for the optimization parameter.
        forces_concat = Iterators.flatten(res.forces[mask])
        G .= - austrip.(forces_concat)  # Minus sign since forces are opposite to gradient.
        G
    end
    f_opt = OptimizationFunction(f; grad=g!)

    # TODO Automatically put constraints on positions for periodic systems
    #      to avoid optimising unneccessarily over R^n, but only over the unit cell.

    OptimizationProblem(f_opt, x0, AC.get_parameters(calculator); kwargs...)
end

"""
Minimise the energy of a system using the specified calculator. For now only optimises
atomic positions. Returns a named tuple including the optimised system as first entry.
Under the hood this constructs an `Optimization.OptimizationProblem` and uses
Optimization.jl to solve it using the passed `solver`.

Typical arguments passed as `solver` are
[`GeometryOptimization.Autoselect()`](@ref) (the default),
[`GeometryOptimization.OptimLBFGS()`](@ref),
[`GeometryOptimization.OptimCG()`](@ref),
[`GeometryOptimization.OptimSD()`](@ref).
These automatically choose some heuristics for setting up the solvers,
which we found to work well in practice.
Beyond that any other `solver`
[compatible with Optimization.jl](https://docs.sciml.ai/Optimization/stable/#Overview-of-the-Optimizers)
can also be employed here.

## Keyword arguments:
- `maxiters`: Maximal number of iterations
- `tol_energy`: Tolerance in the energy to stop the minimisation (all `tol_*` need to be satisfied)
- `tol_force`:  Tolerance in the force  to stop the minimisation (all `tol_*` need to be satisfied)
- `tol_virial`: Tolerance in the virial to stop the minimisation (all `tol_*` need to be satisfied)
- `maxstep`: Maximal step size (in AU) to be taken in a single optimisation step
  (not supported for all `solver`s)
- `callback`: A custom callback, which obtains the pair `(optimization_state, geoopt_state)` and is
  expected to return `false` (continue iterating) or `true` (halt iterations).
- `kwargs`: All other keyword arguments are passed to the call to `solve`. Note, that
  if special `kwargs` should be passed to the `Optimization.OptimizationProblem` the user
  needs to setup the problem manually (e.g. `OptimizationProblem(system, calculator)`)
"""
function minimize_energy!(system, calculator, solver;
                          maxiters=100,
                          tol_energy=Inf*u"eV",
                          tol_force=1e-4u"eV/Å",  # VASP default
                          tol_virial=1e-6u"eV",   # TODO How reasonable ?
                          callback=(x,y) -> false,
                          kwargs...)
    system = convert_to_updatable(system)

    geoopt_state = GeometryOptimizationState(system, calculator)
    problem = OptimizationProblem(system, calculator, geoopt_state; sense=Optimization.MinSense)
    converged = false

    Eold = AC.zero_energy(system, calculator)
    function inner_callback(optim_state, ::Any, geoopt_state)
        halt = callback(optim_state, geoopt_state)
        halt && return true

        energy_converged = abs(geoopt_state.energy - Eold) < tol_energy
        force_converged  = maximum(norm, geoopt_state.forces) < tol_force
        virial_converged = maximum(abs, geoopt_state.virial) < tol_virial

        Eold = geoopt_state.energy
        converged = energy_converged && force_converged && virial_converged
        return converged
    end

    optimres = solve(problem, solver; maxiters, callback=inner_callback, kwargs...)
    (; system=update_not_clamped_positions(system, optimres.u * u"bohr"), converged,
       energy=optimres.objective, state=geoopt_state.calc_state,
       optimres.stats, optimres.alg, optimres)
end

"""Use a heuristic to automatically select the minimisation algorithm
(Currently mostly [`OptimCG`](@ref))"""
struct Autoselect end
function minimize_energy!(system, calculator, ::Autoselect=Autoselect(); kwargs...)
    minimize_energy!(system, calculator, OptimCG(); kwargs...)
end

"""Use Optim's LBFGS implementation with some good defaults."""
struct OptimLBFGS end
function minimize_energy!(system, calculator, ::OptimLBFGS; 
                          maxstep=0.8u"bohr", kwargs...)
    maxstep = austrip(maxstep)
    solver  = Optim.LBFGS(; alphaguess=LineSearches.InitialHagerZhang(),
                            linesearch=LineSearches.BackTracking(; order=2, maxstep))
    minimize_energy!(system, calculator, solver; kwargs...)
end

"""Use Optim's ConjugateGradient implementation with some good defaults."""
struct OptimCG end
function minimize_energy!(system, calculator, ::OptimCG;
                          maxstep=0.8u"bohr", kwargs...)
    maxstep = austrip(maxstep)
    solver  = Optim.ConjugateGradient(;
        linesearch=LineSearches.BackTracking(; order=2, maxstep))
    minimize_energy!(system, calculator, solver; kwargs...)
end

"""Use Optim's GradientDescent (Steepest Descent) implementation with some good defaults."""
struct OptimSD end
function minimize_energy!(system, calculator, ::OptimSD;
                          maxstep=0.8u"bohr", kwargs...)
    maxstep = austrip(maxstep)
    solver  = Optim.GradientDescent(;
        linesearch=LineSearches.BackTracking(; order=2, maxstep))
    minimize_energy!(system, calculator, solver; kwargs...)
end
