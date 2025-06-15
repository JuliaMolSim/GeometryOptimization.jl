module GeometryOptimizationOptimizationExt
using AtomsCalculators
using Optimization
using GeometryOptimization
import GeometryOptimization: GeoOptProblem, GeoOptConvergence
const GO = GeometryOptimization

function GeometryOptimization.solve_problem(
    prob::GeoOptProblem, solver, cvg::GeoOptConvergence;
    callback, maxiters, maxtime, kwargs...)

    function inner_callback(optim_state, ::Any, geoopt_state)
        cache_evaluations = geoopt_state.cache_evaluations

        # Find position in the cache matching the current state
        i_match = findlast(cache_evaluations) do eval
            isnothing(eval.gradnorm) && return false

            tol = 10eps(typeof(optim_state.objective))
            obj_matches = abs(eval.objective - optim_state.objective)  < tol

            if isnothing(optim_state.grad)
                # Nothing we can do, let's just hope it's ok
                grad_matches = true
            else
                g_norm = maximum(abs, optim_state.grad)
                grad_matches = abs(eval.gradnorm  - g_norm) < tol
            end

            obj_matches && grad_matches
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

        # Check for convergence
        geoopt_state.converged = GO.is_converged(cvg, geoopt_state)

        # Callback and possible abortion
        halt = callback(optim_state, geoopt_state)
        halt && return true

        geoopt_state.converged
    end

    optimres = solve(Optimization.OptimizationProblem(prob), solver;
                     maxiters, maxtime, callback=inner_callback, kwargs...)
    (; minimizer=optimres.u, minimum=optimres.objective, optimres)
end


function Optimization.OptimizationProblem(prob::GeoOptProblem; kwargs...)
    f = function(x::AbstractVector{<:Real}, ps)
        GO.eval_objective_gradient!(nothing, prob, ps, x), prob.geoopt_state
    end
    g! = function(G::AbstractVector{<:Real}, x::AbstractVector{<:Real}, ps)
        GO.eval_objective_gradient!(G, prob, ps, x), prob.geoopt_state
        G
    end
    f_opt = OptimizationFunction(f; grad=g!)

    # Note: Some optimisers modify Dofs x0 in-place, so x0 needs to be mutable type.
    x0 = GO.get_dofs(prob.system, prob.dofmgr)
    OptimizationProblem(f_opt, x0, AtomsCalculators.get_parameters(prob.calculator);
                        sense=Optimization.MinSense, kwargs...)
end
end
