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
Turn `system` and `calculator` into a SciML-compatible `OptimizationProblem`.
Optionally pass a reference to the state of the calculator in order to be able
to extract the updated state. Note that the `system` is not updated and that
internally atomic units are used.
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

    OptimizationProblem(f_opt, x0, AC.get_parameters(calculator); kwargs...)
end

"""
Minimise the energy of a system using the specified calculator. For now only optimises
atomic positions. Returns a named tuple including the optimised system as first entry.
Under the hood this constructs an `Optimization.OptimizationProblem` and uses
Optimization.jl to solve it using the passed `solver`. `solver` can be any solver
compatible with Optimization.jl or [`Auto()`](@ref), [`AutoLBFGS()`](@ref), [`AutoCG()`](@ref),
[`AutoSD()`](@ref)


## Keyword arguments:
- `maxiters`: Maximal number of iterations
- `tol_energy`: Tolerance in the energy to stop the minimisation (all `tol_*` need to be satisfied)
- `tol_force`:  Tolerance in the force  to stop the minimisation (all `tol_*` need to be satisfied)
- `tol_virial`: Tolerance in the virial to stop the minimisation (all `tol_*` need to be satisfied)
- `maxstep`: Maximal step size (in AU) to be taken in a single optimisation step
  (not considered in all `solver`s)
- `callback`: A custom callback, which obtains the pair `(optimization_state, geoopt_state)` and is
  expected to return `false` (continue iterating) or `true` (halt iterations).
- `kwargs`: All other keyword arguments are passed to the call to `solve`. Note, that
  if special `kwargs` should be passed to the `Optimization.OptimizationProblem` the user
  needs to setup the problem manually (e.g. `OptimizationProblem(system, calculator)`)
"""
function minimize_energy!(system, calculator, solver;
                          maxiters=100,
                          # TODO Check good defaults here
                          tol_energy=Inf*u"eV",
                          tol_force=1e-5u"eV/Å",
                          tol_virial=1e-6u"eV",
                          callback=(x,y) -> false,
                          kwargs...)
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
        @show uconvert(u"eV/Å", maximum(norm, geoopt_state.forces))

        Eold = geoopt_state.energy
        converged =  energy_converged && force_converged && virial_converged
        return converged
    end

    optimres = solve(problem, solver; maxiters, callback=inner_callback, kwargs...)
    (; system=update_not_clamped_positions(system, optimres.u * u"bohr"), converged,
       energy=optimres.objective, state=geoopt_state.calc_state,
       optimres.stats, optimres)
end

"""Use a heuristic to automatically select the minimisation algorithm"""
struct Auto end
function minimize_energy!(system, calculator, ::Auto=Auto(); kwargs...)
    minimize_energy!(system, calculator, AutoCG(); kwargs...)
end

"""Use Optim's LBFGS implementation"""
struct AutoLBFGS end
function minimize_energy!(system, calculator, ::AutoLBFGS; maxstep=Inf, kwargs...)
    solver = Optim.LBFGS(; alphaguess=LineSearches.InitialHagerZhang(),
                           linesearch=LineSearches.BackTracking(; order=2, maxstep))
    minimize_energy!(system, calculator, solver; kwargs...)
end

"""Use Optim's ConjugateGradient implementation"""
struct AutoCG end
function minimize_energy!(system, calculator, ::AutoCG; maxstep=Inf, kwargs...)
    solver = Optim.ConjugateGradient(;linesearch=LineSearches.BackTracking(;order=2, maxstep))
    minimize_energy!(system, calculator, solver; kwargs...)
end

"""Use Optim's GradientDescent (Steepest Descent) implementation"""
struct AutoSD end
function minimize_energy!(system, calculator, ::AutoSD; maxstep=Inf, kwargs...)
    solver = Optim.GradientDescent(;linesearch=LineSearches.BackTracking(;order=2, maxstep))
    minimize_energy!(system, calculator, solver; kwargs...)
end
