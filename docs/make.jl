using JADE191
using Documenter

DocMeta.setdocmeta!(JADE191, :DocTestSetup, :(using JADE191); recursive=true)

makedocs(;
    modules=[JADE191],
    authors="kevkevs <22803504+kevkevs@users.noreply.github.com> and contributors",
    repo="https://github.com/kevkevs/JADE191.jl/blob/{commit}{path}#{line}",
    sitename="JADE191.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kevkevs.github.io/JADE191.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kevkevs/JADE191.jl",
    devbranch="main",
)
