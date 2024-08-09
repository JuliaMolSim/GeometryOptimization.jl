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
    if isempty(system)
        throw(ArgumentError("Cannot optimise a system without atoms."))
    end

    # Get the mask of all atoms which are not clamped (i.e. which are optimised)
    mask = [!get(atom, :clamped, false) for atom in system]
    if !any(mask)
        throw(ArgumentError("Cannot optimise systems where all atoms are clamped."))
    end

    # Use current system parameters as starting positions.
    # Some optimisers modify x0 in-place, so need a mutable type.
    LT = eltype(system[1, :position])
    nonclamped_positions = reinterpret(reshape, LT, system[mask, :position])
    x0 = austrip.(vec(nonclamped_positions))

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
function minimize_energy!(system, calculator, solver=Autoselect(); kwargs...)
    _minimize_energy!(system, calculator, solver; kwargs...)
end

struct GeoOptTabularPrint
end
function (cb::GeoOptTabularPrint)(optim_state, geoopt_state)
    show_force = true  # TODO Set to false if have no forces
    show_time  = true

    if optim_state.iter == 0
        label_force = show_force ? ("   max(Force)", "   ----------") : ("", "")
        label_time  = show_time ? ("   Δtime", "   ------") : ("", "")
        println("n     Energy            log10(ΔE)", label_force[1], label_time[1])
        println("---   ---------------   ---------", label_force[2], label_time[2])
    end
    return false

    fstr = 

    tstr = " "^9
    if show_time
        deltatime_ns = 0
        tstr = @sprintf "   % 6s" TimerOutputs.prettytime(deltatime_ns)
    end


    # iteration
    # change in energy (logscale)
    # current force
    # current virial
    # maximal atomic displacement
    # trace of lattice deformation matrix
    # timing
    # other flags (e.g. 



    return false
end

# Function that does all the work. The idea is that calculator implementations
# can provide more specific methods for minimize_energy! calling the function below.
# This allows to adjust (based on the calculator type) the default parameters,
# do some additional calculator-specific setup (e.g. callbacks) and so on. Then
# by calling this function the actual minimisation is started off.
function _minimize_energy!(system, calculator, solver;
                           maxiters=100,
                           tol_energy=Inf*u"eV",
                           tol_force=1e-4u"eV/Å",  # VASP default
                           tol_virial=1e-6u"eV",   # TODO How reasonable ?
                           maxstep=0.8u"bohr",
                           callback=(x,y) -> false,
                           kwargs...)
    solver = setup_solver(system, calculator, solver; maxstep)
    system = convert_to_updatable(system)

    geoopt_state = GeometryOptimizationState(system, calculator)
    problem = OptimizationProblem(system, calculator, geoopt_state; sense=Optimization.MinSense)
    converged = false

    Eold = AC.zero_energy(system, calculator)
    function inner_callback(optim_state, ::Any, geoopt_state)
        halt = callback(optim_state, geoopt_state)
        halt && return true

        energy_converged = austrip(abs(geoopt_state.energy - Eold))    < austrip(tol_energy)
        force_converged  = austrip(maximum(norm, geoopt_state.forces)) < austrip(tol_force)
        virial_converged = austrip(maximum(abs,  geoopt_state.virial)) < austrip(tol_virial)

        Eold = geoopt_state.energy
        converged = energy_converged && force_converged && virial_converged
        return converged
    end

    optimres = solve(problem, solver; maxiters, callback=inner_callback, kwargs...)
    (; system=update_not_clamped_positions(system, optimres.u * u"bohr"), converged,
       energy=optimres.objective, geoopt_state.forces, geoopt_state.virial,
       state=geoopt_state.calc_state, optimres.stats, optimres.alg, optimres)
end

# Default setup_solver function just passes things through
setup_solver(system, calculator, solver::Any; kwargs...) = solver

"""Use a heuristic to automatically select the minimisation algorithm
(Currently [`OptimCG`](@ref), but this may change silently)"""
struct Autoselect end
function setup_solver(system, calculator, ::Autoselect; kwargs...)
    setup_solver(system, calculator, OptimCG(); kwargs...)
end

"""Use Optim's LBFGS implementation with some good defaults."""
struct OptimLBFGS end
function setup_solver(system, calculator, ::OptimLBFGS; maxstep, kwargs...)
    # TODO Maybe directly specialise the method on ::Optim.LBFGS and don't
    #      provide the GeometryOptimisation.OptimLBFGS() marker struct at all,
    #      similar for CG and SD ?
    linesearch = LineSearches.BackTracking(; order=2, maxstep=austrip(maxstep))
    Optim.LBFGS(; linesearch, alphaguess=LineSearches.InitialHagerZhang())
end

"""Use Optim's ConjugateGradient implementation with some good defaults."""
struct OptimCG end
function setup_solver(system, calculator, ::OptimCG; maxstep, kwargs...)
    linesearch = LineSearches.BackTracking(; order=2, maxstep=austrip(maxstep))
    Optim.ConjugateGradient(; linesearch)
end

"""Use Optim's GradientDescent (Steepest Descent) implementation with some good defaults."""
struct OptimSD end
function setup_solver(system, calculator, ::OptimSD; maxstep, kwargs...)
    linesearch = LineSearches.BackTracking(; order=2, maxstep=austrip(maxstep))
    Optim.GradientDescent(; linesearch)
end
