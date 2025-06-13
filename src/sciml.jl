function Optimization.OptimizationProblem(prob::GeoOptProblem; kwargs...)
    system = problem.system
    calculator = problem.calculator
    dofmgr = problem.dofmgr
    geoopt_state = problem.geoopt_state

    f = function(x::AbstractVector{<:Real}, ps)
        eval_objective_gradient!(nothing, prob, ps, x), prob.geoopt_state
    end
    g! = function(G::AbstractVector{<:Real}, x::AbstractVector{<:Real}, ps)
        eval_objective_gradient!(G, prob, ps, x), prob.geoopt_state
        G
    end
    f_opt = OptimizationFunction(f; grad=g!)

    # Note: Some optimisers modify Dofs x0 in-place, so x0 needs to be mutable type.
    x0 = get_dofs(prob.system, prob.dofmgr)
    OptimizationProblem(f_opt, x0, AC.get_parameters(prob.calculator);
                        sense=Optimization.MinSense, kwargs...)
end


function solve_problem(prob::GeoOptProblem, solver, cvg::GeoOptConvergence;
                       callback, maxiters, maxtime, kwargs...)

    function inner_callback(optim_state, ::Any, geoopt_state)
        cache_evaluations = geoopt_state.cache_evaluations

        # Find position in the cache matching the current state
        i_match = findlast(cache_evaluations) do eval
            isnothing(eval.gradnorm) && return false
            tol = 10eps(typeof(optim_state.objective))
            g_norm = maximum(abs, optim_state.gradient)
            (   abs(eval.objective - optim_state.objective)  < tol
             && abs(eval.gradnorm  - g_norm) < tol)
        end
        i_match = @something i_match length(cache_evaluations)

        # Commit data from state and discard the rest
        geoopt_state.n_iter = optim_state.iter
        push!(geoopt_state.history_energy, cache_evaluations[i_match].energy)
        if !isnothing(cache_evaluations[i_match].forces)
            geoopt_state.forces .= cache_evaluations[i_match].forces
        end
        if !isnothing(cache_evaluations[i_match].virial)
            geoopt_state.virial .= cache_evaluations[i_match].virial
        end
        empty!(cache_evaluations)

        # Callback and possible abortion
        halt = callback(geoopt_state, ts)
        halt && return true

        # Check for convergence
        geoopt_state.converged = is_converged(cvg, geoopt_state)
        return geoopt_state.converged
    end

    optimres = solve(Optimization.OptimizationProblem(prob), solver;
                     maxiters, maxtime, callback=inner_callback, kwargs...)
    (; minimizer=optimres.u, optimres.objective, optimres)
end
