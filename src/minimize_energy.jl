mutable struct GeometryOptimizationState
    start_time::UInt64  # Time when the calculation started from time_ns()
    n_iter::Int
    converged::Bool
    #
    calculator          # Calculator used by the problem
    calc_state          # Current calculator state
    history_energy      # History of all energy values
    forces              # Current force  value
    virial              # Current virial value
    #
    # Cache all objective, gradient, energies, forces, virial
    # These lists get emptied by the callback and are only around
    # to update history_energy, energy, forces and virial properly.
    # The issue here is that potentially multiple objective evaluations
    # per steps are performed and this is a way to make sure
    # what we display is consistent.
    cache_evaluations::Vector{Any}
end
function GeometryOptimizationState(system, calculator; start_time=time_ns())
    GeometryOptimizationState(start_time, 0, false, calculator,
                              AC.get_state(calculator),
                              typeof(AC.zero_energy(system, calculator))[],
                              AC.zero_forces(system, calculator),
                              Array(AC.zero_virial(system, calculator)), [])
end


"""
    GeoOptProblem(system, calculator, dofmgr, geoopt_state; kwargs...)

Turn `system`, `calculator` and `geoopt_state::GeometryOptimizationState`
into an `OptimizationProblem` for `solve_problem`. Note that the `system` is not updated
automatically and that internally atomic units are used.
"""
struct GeoOptProblem{System,Calc,Dof,State}
    system::System
    calculator::Calc
    dofmgr::Dof
    geoopt_state::State
end
function eval_objective_gradient!(G, prob::GeoOptProblem, ps, x)
    geoopt_state = prob.geoopt_state
    res = eval_objective(prob.system, prob.calculator, prob.dofmgr, x, ps, geoopt_state.calc_state)
    objective = res.energy_unitless
    energy = res.energy

    gradnorm = nothing
    forces   = nothing
    virial   = nothing
    if !isnothing(G)
        res = eval_gradient(prob.system, prob.calculator, prob.dofmgr, x, ps, res.state)
        gradnorm = maximum(abs, res.grad)
        haskey(res, :forces) && (forces = res.forces)
        haskey(res, :virial) && (virial = res.virial)
        copy!(G, res.grad)
    end

    # Commit state
    min_energy = minimum(cache.energy for cache in geoopt_state.cache_evaluations;
                         init=abs(energy) + 100u"hartree")
    if energy ≤ min_energy
        geoopt_state.calc_state = res.state
    end
    push!(geoopt_state.cache_evaluations, (; energy, forces, virial, objective, gradnorm))

    objective
end

struct GeoOptConvergence
    tol_energy
    tol_forces
    tol_virial
    check_virial::Bool
end
function is_converged(cvg::GeoOptConvergence, geoopt_state::GeometryOptimizationState)
    is_almost_zero(u::Quantity{T}) where {T} = ustrip(u) < 5eps(T)

    if length(geoopt_state.history_energy) > 1
        ene_diff = abs(geoopt_state.history_energy[end-1] - geoopt_state.history_energy[end])
        energy_converged = austrip(abs(ene_diff)) < austrip(cvg.tol_energy)
    else
        energy_converged = false
    end

    force_converged = austrip(maximum(norm, geoopt_state.forces)) < austrip(cvg.tol_forces)
    force_zero = is_almost_zero(maximum(norm, geoopt_state.forces))

    if cvg.check_virial
        virial_converged = austrip(maximum(abs, geoopt_state.virial)) < austrip(cvg.tol_virial)
        virial_zero = is_almost_zero(maximum(abs, geoopt_state.virial))
    else
        virial_converged = true
        virial_zero = true
    end

    @debug force_zero virial_zero energy_converged force_converged virial_converged

    # If force and virial are zero, nothing can possibly happen in the future
    if force_zero && virial_zero
        return true
    else
        return energy_converged && force_converged && virial_converged
    end
end


"""
Minimise the energy of a `system` using the specified `calculator`. Optimises either only
atomic positions (`variablecell=false`) or both positions and lattice (`variablecell=true`).
Returns a named tuple including the optimised system as first entry. Typical arguments passed
as `solver` are
[`Autoselect()`](@ref) (the default),
[`OptimLBFGS()`](@ref),
[`OptimCG()`](@ref),
[`OptimSD()`](@ref).
These automatically choose some heuristics for setting up the solvers,
which we found to work well in practice.

Beyond that any other first-order `solver`
[compatible with Optimization.jl](https://docs.sciml.ai/Optimization/stable/#Overview-of-the-Optimizers)
can also be employed. Note, that in principle all such solvers should work, but we only
tested a small fraction and you can expect that minor modifications are needed to make
some solvers work (PRs appreciated!). In general only first-order or second-order methods work.

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
                           tol_virial=5e-5u"eV",    # QE default is about 4.6e-5 eV
                           maxstep=0.8u"bohr",
                           verbosity::Integer=0,
                           callback=GeoOptDefaultCallback(verbosity;
                                                          show_virial=variablecell),
                           kwargs...)
    if isempty(system)
        throw(ArgumentError("Cannot optimise a system without atoms."))
    end
    system = convert_to_updatable(system)
    solver = setup_solver(system, calculator, solver; maxstep)

    # TODO Think carefully whether this is the best interface and integration
    clamp = [iatom for (iatom, atom) in enumerate(system) if get(atom, :clamp, false)]
    dofmgr = DofManager(system; variablecell, clamp)
    geoopt_state = GeometryOptimizationState(system, calculator)
    problem = GeoOptProblem(system, calculator, dofmgr, geoopt_state)
    cvg = GeoOptConvergence(tol_energy, tol_forces, tol_virial, variablecell)

    # Run the problem. Note that this mutates geoopt_state
    res = solve_problem(problem, solver, cvg::GeoOptConvergence;
                        callback, maxiters, maxtime, kwargs...)
    geoopt_state.converged || @warn "Geometry optimisation not converged."

    (; system=set_dofs(system, dofmgr, res.minimizer), geoopt_state.converged,
       energy=res.minimum * u"hartree", geoopt_state.forces, geoopt_state.virial,
       state=geoopt_state.calc_state, geoopt_state.history_energy, geoopt_state.n_iter,
       res.optimres)
end

# Default setup_solver function just passes things through
setup_solver(system, calculator, solver::Any; kwargs...) = solver

"""Use a heuristic to automatically select the minimisation algorithm
(Currently [`OptimCG`](@ref), but this may change silently)"""
struct Autoselect end
function setup_solver(system, calculator, ::Autoselect; kwargs...)
    setup_solver(system, calculator, OptimCG(); kwargs...)
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
