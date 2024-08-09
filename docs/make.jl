# Setup julia dependencies for docs generation if not yet done
import Pkg
Pkg.activate(@__DIR__)
if !isfile(joinpath(@__DIR__, "Manifest.toml"))
    Pkg.develop(Pkg.PackageSpec(path=joinpath(@__DIR__, "..")))
    Pkg.instantiate()
end

using GeometryOptimization
using Documenter

DocMeta.setdocmeta!(GeometryOptimization, :DocTestSetup, :(using GeometryOptimization); recursive=true)

makedocs(;
    modules=[GeometryOptimization],
    authors="JuliaMolSim community",
    repo="https://github.com/JuliaMolSim/GeometryOptimization.jl/blob/{commit}{path}#{line}",
    sitename="GeometryOptimization.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaMolSim.github.io/GeometryOptimization.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => [
            "examples/aluminium_dftk.md",
            "examples/other_solvers.md",
            "examples/tial_lj.md",
        ],
        "apireference.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaMolSim/GeometryOptimization.jl",
    devbranch="main",
)
