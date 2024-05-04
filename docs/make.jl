using CFPlanets
using Documenter

DocMeta.setdocmeta!(CFPlanets, :DocTestSetup, :(using CFPlanets); recursive=true)

makedocs(;
    modules=[CFPlanets],
    authors="The ClimFlows contributors",
    sitename="CFPlanets.jl",
    format=Documenter.HTML(;
        canonical="https://ClimFlows.github.io/CFPlanets.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ClimFlows/CFPlanets.jl",
    devbranch="main",
)
