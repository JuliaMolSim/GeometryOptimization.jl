import TimerOutputs
using PrettyTables
using Printf

"""
Callback producing a convergence table summarising the geometry
optimisation convergence. If `always_show_header=true` the
header is shown in each iteration. This is helpful if the calculator
produces output as well.
"""
struct GeoOptDefaultCallback
    always_show_header::Bool
    prev_time::Ref{UInt64}
end
function GeoOptDefaultCallback(; always_show_header=false)
    GeoOptDefaultCallback(always_show_header, Ref{UInt64}(0))
end
function GeoOptDefaultCallback(verbosity::Integer=1; kwargs...)
    if verbosity <= 0
        return (os, gs) -> false  # No printing => no callback
    else
        return GeoOptDefaultCallback(; always_show_header=verbosity > 1; kwargs...)
    end
end

format_log8(value) = (value < 0 ? " " : "+") * (@sprintf "%8.2f" log10(abs(value)))

function (cb::GeoOptDefaultCallback)(optim_state, geoopt_state)
    # If first iteration clear a potentially cached previous time
    optim_state.iter ≤ 0 && (cb.prev_time[] = 0)
    runtime_ns = time_ns() - geoopt_state.start_time
    tstr = @sprintf "% 6s" TimerOutputs.prettytime(runtime_ns - cb.prev_time[])
    cb.prev_time[] = runtime_ns

    Estr  = (@sprintf "%+15.12f" round(austrip(geoopt_state.energy), sigdigits=13))[1:15]
    logΔE = optim_state.iter < 1 ? "" : format_log8(austrip(geoopt_state.energy_change))

    maxforce = austrip(maximum(norm, geoopt_state.forces))
    fstr = iszero(maxforce) ? "" : round(maxforce, sigdigits=8)

    fields = [  # Header, width, value
        ("n",           3, optim_state.iter),
        ("Energy",     15, Estr),
        ("log10(ΔE)",   9, logΔE),
        ("max(Force)", 10, fstr),
        # TODO Maximal atomic displacement
        # TODO Current virial, trace of lattice deformation matrix
        ("Δtime",       6, tstr),
        # TODO Would be nice to have some simple way to add in
        #      a few calculator-specific things (e.g. total number of SCF iterations)
    ]

    if cb.always_show_header
        hlines = :all
        show_header = true
    elseif optim_state.iter == 0
        hlines = [0, 1]
        show_header = true
    else
        hlines = :none
        show_header = false
    end
    title = iszero(optim_state.iter) ? "Geometry optimisation convergence (in atomic units)" : ""
    pretty_table(reshape(getindex.(fields, 3), 1, length(fields));
                 header=Symbol.(getindex.(fields, 1)), columns_width=getindex.(fields, 2),
                 title, show_header, hlines)

    flush(stdout)
    return false
end
