using GeometryOptimization
using Documenter

DocMeta.setdocmeta!(GeometryOptimization, :DocTestSetup, :(using GeometryOptimization); recursive=true)

makedocs(;
    modules=[GeometryOptimization],
    authors="Christoph Ortner <christohortner@gmail.com> and contributors",
    repo="https://github.com/ortner/GeometryOptimization.jl/blob/{commit}{path}#{line}",
    sitename="GeometryOptimization.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ortner.github.io/GeometryOptimization.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ortner/GeometryOptimization.jl",
    devbranch="main",
)
