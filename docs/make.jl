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
    ],
)

deploydocs(;
    repo="github.com/JuliaMolSim/GeometryOptimization.jl",
    devbranch="main",
)
