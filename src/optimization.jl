#
# Note that by default all particles in the system are assumed optimizable.
# We always work in cartesian coordinates.
#

mutable struct GeometryOptimizationState
    calculator
    calc_state          # Reference to the most recent calculator state
    energy              # Current energy value
    forces              # Current force value
    virial              # Current virial value
    energy_change       # Change in energy since previous step
    start_time::UInt64  # Time when the calculation started from time_ns()
end
function GeometryOptimizationState(system, calculator; start_time=time_ns())
    GeometryOptimizationState(
        calculator,
        AC.get_state(calculator),
        AC.zero_energy(system, calculator),
        AC.zero_forces(system, calculator),
        Array(AC.zero_virial(system, calculator)),
        AC.zero_energy(system, calculator),
        start_time,
    )
end

"""
    Optimization.OptimizationProblem(system, calculator, dofmgr, geoopt_state; kwargs...)

Turn `system`, `calculator` and `geoopt_state::GeometryOptimizationState`
into a SciML-compatible `OptimizationProblem`. Note that the `system` is not updated
automatically and that internally atomic units are used.
"""
function Optimization.OptimizationProblem(system, calculator, dofmgr, geoopt_state; kwargs...)
    if isempty(system)
        throw(ArgumentError("Cannot optimise a system without atoms."))
    end

    f = function(x::AbstractVector{<:Real}, ps)
        res = energy_dofs(system, calculator, dofmgr, x, ps, geoopt_state.calc_state)
        geoopt_state.calc_state = res.state
        res.energy_unitless, geoopt_state
    end
    g! = function(G::AbstractVector{<:Real}, x::AbstractVector{<:Real}, ps)
        res = gradient_dofs(system, calculator, dofmgr, x, ps, geoopt_state.calc_state)
        geoopt_state.calc_state = res.state
        haskey(res, :forces) && (geoopt_state.forces .= res.forces)
        haskey(res, :virial) && (geoopt_state.virial .= res.virial)
        copy!(G, res.grad)
    end
    f_opt = OptimizationFunction(f; grad=g!)

    # TODO Automatically put constraints on positions for periodic systems
    #      to avoid optimising unneccessarily over R^n, but only over the unit cell.

    # Note: Some optimisers modify Dofs x0 in-place, so x0 needs to be mutable type.
    x0 = get_dofs(system, dofmgr)
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
- `variablecell`: Determines whether the cell is fixed or allowed to change during optimization
- `maxiters`: Maximal number of iterations
- `maxtime`:  Maximal allowed runtime (in seconds)
- `tol_energy`: Tolerance in the energy to stop the minimisation (all `tol_*` need to be satisfied)
- `tol_forces`:  Tolerance in the force  to stop the minimisation (all `tol_*` need to be satisfied)
- `tol_virial`: Tolerance in the virial to stop the minimisation (all `tol_*` need to be satisfied)
- `maxstep`: Maximal step size (in AU or length units) to be taken in a single optimisation step
  (not supported for all `solver`s)
- `verbosity`: Printing level. The idea is that `0` is silent, `1` displays the optimisation
  progress and `≥ 2` starts displaying things from the calculator as well (e.g SCF iterations).
- `callback`: A custom callback, which obtains the pair `(optimization_state, geoopt_state)`
  and is expected to return `false` (continue iterating) or `true` (halt iterations). Note
  that specifying this overwrites the default printing callback. The calculation thus becomes
  silent unless a [`GeoOptDefaultCallback`](@ref) is included in the callback.
- `kwargs`: All other keyword arguments are passed to the call to `solve`. Note, that
  if special `kwargs` should be passed to the `Optimization.OptimizationProblem` the user
  needs to setup the problem manually (e.g. `OptimizationProblem(system, calculator)`)
"""
function minimize_energy!(system, calculator, solver=Autoselect(); kwargs...)
    _minimize_energy!(system, calculator, solver; kwargs...)
end

# Function that does all the work. The idea is that calculator implementations
# can provide more specific methods for minimize_energy! calling the function below.
# This allows to adjust (based on the calculator type) the default parameters,
# do some additional calculator-specific setup (e.g. callbacks) and so on. Then
# by calling this function the actual minimisation is started off.
function _minimize_energy!(system, calculator, solver;
                           variablecell=false,
                           maxiters::Integer=100,
                           maxtime::Integer=60*60*24*365,  # 1 year
                           tol_energy=Inf*u"eV",
                           tol_forces=1e-4u"eV/Å",  # VASP default
                           tol_virial=1e-6u"eV",    # TODO How reasonable ?
                           maxstep=0.8u"bohr",
                           verbosity::Integer=0,
                           callback=GeoOptDefaultCallback(verbosity;
                                                          show_virial=variablecell),
                           kwargs...)
    solver = setup_solver(system, calculator, solver; maxstep)
    system = convert_to_updatable(system)

    # TODO Think carefully whether this is the best interface and integration
    clamp = [iatom for (iatom, atom) in enumerate(system) if get(atom, :clamp, false)]
    dofmgr = DofManager(system; variablecell, clamp)
    geoopt_state = GeometryOptimizationState(system, calculator)
    problem = OptimizationProblem(system, calculator, dofmgr, geoopt_state;
                                  sense=Optimization.MinSense)
    converged = false

    Eold = 0  # Atomic units
    function inner_callback(optim_state, ::Any, geoopt_state)
        geoopt_state.energy        = optim_state.objective * u"hartree"
        geoopt_state.energy_change = (optim_state.objective - Eold) * u"hartree"

        halt = callback(optim_state, geoopt_state)
        halt && return true

        energy_converged = abs(optim_state.objective - Eold) < austrip(tol_energy)
        if optim_state.iter < 1
            force_converged = false
        else
            force_converged  = austrip(maximum(norm, geoopt_state.forces)) < austrip(tol_forces)
        end

        if variablecell  # Also check virial to determine convergence
            virial_converged = maximum(abs, geoopt_state.virial) < tol_virial
        else
            virial_converged = true
        end

        Eold = optim_state.objective
        converged = energy_converged && force_converged && virial_converged
        return converged
    end

    optimres = solve(problem, solver; maxiters, maxtime, callback=inner_callback, kwargs...)
    converged || @warn "Geometry optimisation not converged."
    (; system=set_dofs(system, dofmgr, optimres.u), converged,
       energy=optimres.objective * u"hartree", geoopt_state.forces, geoopt_state.virial,
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

#
# Some helper functions, which we can hopefully remove at some future point
#

"""
Convert the input system to a kind of system we can work with, i.e. one where
we have the `particles` keyword argument for setting new particles and positions.
"""
function convert_to_updatable(system)
    if length(system) > 0 && !(system[1] isa Atom)
        particles = [Atom(; pairs(atom)...) for atom in system]
        AbstractSystem(system; particles)
    else
        system
    end
end

"""
Clamp given atoms in the system. Clamped atoms are fixed and their positions
will not be optimized. The atoms to be clamped should be given as a list of
indices corresponding to their positions in the system.

!!! warn "Experimental"
    This is this a very experimental interface and will likely change
    in the future.
"""
function clamp_atoms(system, clamped_indexes::Union{AbstractVector{<:Integer},Nothing})
    system = convert_to_updatable(system)
    clamped = falses(length(system))
    clamped[clamped_indexes] .= true
    particles = [Atom(atom; clamped=msk) for (atom, msk) in zip(system, clamped)]
    AbstractSystem(system; particles)
end
