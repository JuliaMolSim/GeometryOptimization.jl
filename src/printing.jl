import TimerOutputs
using PrettyTables
using Printf

"""
Callback producing a convergence table summarising the geometry
optimisation convergence. If `always_show_header=true` the
header is shown in each iteration. This is helpful if the calculator
produces output as well.
"""
struct GeoOptDefaultPrint{F}
    always_show_header::Bool
    callback_extra_fields::F
    prev_time::Ref{UInt64}
end
function GeoOptDefaultPrint(; always_show_header=false, callback_extra_fields=nothing)
    GeoOptDefaultPrint(always_show_header, callback_extra_fields, Ref{UInt64}(0))
end

format_log8(value) = (value < 0 ? " " : "+") * (@sprintf "%8.2f" log10(abs(value)))

function (cb::GeoOptDefaultPrint)(optim_state, geoopt_state)
    # If first iteration clear a potentially cached previous time
    optim_state.iter ≤ 0 && (cb.prev_time[] = 0)

    runtime_ns = time_ns() - geoopt_state.start_time
    tstr = @sprintf "% 6s" TimerOutputs.prettytime(runtime_ns - cb.prev_time[])
    cb.prev_time[] = runtime_ns

    Estr  = (@sprintf "%+15.12f" round(austrip(geoopt_state.energy), sigdigits=13))[1:15]
    logΔE = optim_state.iter < 1 ? "" : format_log8(austrip(geoopt_state.energy_change))

    maxforce = austrip(maximum(norm, geoopt_state.forces))
    fstr = iszero(maxforce) ? "" : round(maxforce, sigdigits=8)

    data = [  # Header, width, value
        ("n",           3, optim_state.iter),
        ("Energy",     15, Estr),
        ("log10(ΔE)",   9, logΔE),
        ("max(Force)", 10, fstr),
        # TODO Maximal atomic displacement
        # TODO Current virial, trace of lattice deformation matrix
        ("Δtime",       6, tstr),
    ]
    if !isnothing(cb.callback_extra_fields)
        append!(data, cb.additional_fields(optim_state, geoopt_state))
    end

    title = iszero(optim_state.iter) ? "Geometry optimisation convergence (in atomic units)" : ""
    show_header = optim_state.iter == 0 || cb.always_show_header
    pretty_table(reshape(getindex.(data, 3), 1, length(data));
                 header=Symbol.(getindex.(data, 1)), title,
                 columns_width=getindex.(data, 2),
                 show_header, hlines=show_header ? [0, 1] : :none)

    flush(stdout)
    return false
end
