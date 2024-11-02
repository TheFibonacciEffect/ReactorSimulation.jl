# ReactorSimulation

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://TheFibonacciEffect.github.io/ReactorSimulation.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://TheFibonacciEffect.github.io/ReactorSimulation.jl/dev/)
[![Build Status](https://github.com/TheFibonacciEffect/ReactorSimulation.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/TheFibonacciEffect/ReactorSimulation.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/TheFibonacciEffect/ReactorSimulation.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/TheFibonacciEffect/ReactorSimulation.jl)

Van neumann boundary conditions: https://math.stackexchange.com/questions/2706701/neumann-boundary-conditions-in-finite-difference

Download Julia using 
`curl -fsSL https://install.julialang.org | sh`
(on Linux), otherwise download julia from https://julialang.org/downloads/

Then start julia in the same directory as the `Project.toml` is located.
Start julia using `julia`.
Go into package mode by pressing the `]` key. Activate the local environment by typing `activate .`, now install all required dependancies by typing `instantiate`. 
Afterwards you can run the files using for example `julia --project=. ./src/ex1.jl`.