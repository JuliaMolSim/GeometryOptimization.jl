#
# Solvers with sane defaults
#

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
# Solve problem
#
function solve_problem(prob::GeoOptProblem, solver::Optim.AbstractOptimizer, cvg::GeoOptConvergence;
                       callback, maxiters, maxtime, kwargs...)
    ps = AC.get_parameters(prob.calculator)
    fg! = function(F, G, x)
        objective = eval_objective_gradient!(G, prob, ps, x)
        if isnothing(F)
            return nothing
        else
            return objective
        end
    end

    geoopt_state = prob.geoopt_state
    inner_callback = function(ts)
        cache_evaluations = geoopt_state.cache_evaluations

        geoopt_state.n_iter = ts.iteration
        if isempty(cache_evaluations)
            # Find out if we already added the current state (if optim cannot
            # make progress it keeps printing iterations, but does not run further
            # function evaluations ... in this case we have no new forces and virials).
            # Also it sometimes does an extra call to the callback even though
            # convergence has already been flagged.
            tol = 10eps(typeof(ts.value))
            is_match = abs(austrip(geoopt_state.history_energy[end]) - ts.value) < tol
            if !geoopt_state.converged && !is_match
                @warn "Discarding optimisation step of iteration $(ts.iteration)"
            end
        else
            # Find position in the cache matching Optim's current state
            i_match = findlast(cache_evaluations) do eval
                isnothing(eval.gradnorm) && return false
                tol = 10eps(typeof(ts.value))
                (   abs(eval.objective - ts.value)  < tol
                 && abs(eval.gradnorm  - ts.g_norm) < tol)
            end
            i_match = @something i_match length(cache_evaluations)

            # Commit data from state and discard the rest
            push!(geoopt_state.history_energy, cache_evaluations[i_match].energy)
            if !isnothing(cache_evaluations[i_match].forces)
                geoopt_state.forces .= cache_evaluations[i_match].forces
            end
            if !isnothing(cache_evaluations[i_match].virial)
                geoopt_state.virial .= cache_evaluations[i_match].virial
            end
            empty!(cache_evaluations)

            # Check for convergence
            geoopt_state.converged = is_converged(cvg, geoopt_state)
        end

        # Callback and possible abortion
        halt = callback(ts, geoopt_state)
        halt && return true

        geoopt_state.converged
    end

    x0 = get_dofs(prob.system, prob.dofmgr)
    T = eltype(x0)

    options = Optim.Options(;
        allow_f_increases=true,
        successive_f_tol=2,
        callback=inner_callback,
        x_abstol=-1, f_abstol=-1, g_tol=10eps(T),
        x_reltol=-1, f_reltol=-1,
        iterations=maxiters,
        time_limit=maxtime,
        kwargs...
    )
    optimres = Optim.optimize(Optim.only_fg!(fg!), x0, solver, options)

    (; minimizer=Optim.minimizer(optimres), minimum=Optim.minimum(optimres), optimres)
end

function solve_problem(prob, solver::Optim.ZerothOrderOptimizer, cvg;
                       callback, maxiters, maxtime, kwargs...)
    # TODO Supporting this needs more fiddeling with the callbacks and convergence checks
    throw(ArgumentError("Zeroth-order optimizers are currently not supported."))
end
