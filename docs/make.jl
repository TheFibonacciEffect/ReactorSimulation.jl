using ReactorSimulation
using Documenter

DocMeta.setdocmeta!(ReactorSimulation, :DocTestSetup, :(using ReactorSimulation); recursive=true)

makedocs(;
    modules=[ReactorSimulation],
    authors="Caspar Gutsche",
    sitename="ReactorSimulation.jl",
    format=Documenter.HTML(;
        canonical="https://TheFibonacciEffect.github.io/ReactorSimulation.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/TheFibonacciEffect/ReactorSimulation.jl",
    devbranch="main",
)
