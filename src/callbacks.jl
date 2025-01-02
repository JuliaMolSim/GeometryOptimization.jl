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
    verbosity::Int
    always_show_header::Bool
    show_virial::Bool
    prev_time::Ref{UInt64}
end
function GeoOptDefaultCallback(verbosity=1;
                               show_virial=true, always_show_header=verbosity > 1)
    GeoOptDefaultCallback(verbosity, always_show_header, show_virial, Ref{UInt64}(0))
end

format_log8(value) = (value < 0 ? " " : "+") * (@sprintf "%8.2f" log10(abs(value)))

function (cb::GeoOptDefaultCallback)(optim_state, geoopt_state)
    cb.verbosity ≤ 0 && return false  # No printing, just continue iterations

    # If first iteration clear a potentially cached previous time
    optim_state.iter ≤ 0 && (cb.prev_time[] = 0)
    runtime_ns = time_ns() - geoopt_state.start_time
    tstr = @sprintf "% 6s" TimerOutputs.prettytime(runtime_ns - cb.prev_time[])
    cb.prev_time[] = runtime_ns

    Estr  = (@sprintf "%+15.12f" round(austrip(geoopt_state.energy), sigdigits=13))[1:15]
    if iszero(geoopt_state.energy_change) && optim_state.iter < 1
        logΔE = ""
    else
        logΔE = format_log8(austrip(geoopt_state.energy_change))
    end

    maxforce = austrip(maximum(norm, geoopt_state.forces))
    fstr = iszero(maxforce) ? "" : round(maxforce, sigdigits=8)

    fields = [  # Header, width, value
        ("n",           3, optim_state.iter),
        ("Energy",     15, Estr),
        ("log10(ΔE)",   9, logΔE),
        ("max(Force)", 11, fstr),
        # TODO Maximal atomic displacement
    ]
    if cb.show_virial
        maxvirial = austrip(maximum(abs, geoopt_state.virial))
        pressure  = -austrip(tr(geoopt_state.virial)) / 3
        vstr = iszero(maxvirial) ? "" : round(maxvirial, sigdigits=8)
        pstr = iszero(pressure)  ? "" : round(pressure,  sigdigits=2)
        push!(fields, ("max(Virial)", 11, vstr))
        push!(fields, ("Pressure",     8, pstr))
    end
    push!(fields, ("Δtime", 6, tstr))
    # TODO Would be nice to have some simple way to add in
    #      a few calculator-specific things (e.g. total number of SCF iterations)

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

    cb.always_show_header && println(stdout)
    pretty_table(stdout, reshape(getindex.(fields, 3), 1, length(fields));
                 header=Symbol.(getindex.(fields, 1)), columns_width=getindex.(fields, 2),
                 title, show_header, hlines)
    cb.always_show_header && println(stdout)

    flush(stdout)
    return false
end
