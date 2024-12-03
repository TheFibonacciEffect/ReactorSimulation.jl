# ReactorSimulation
Van neumann boundary conditions: https://math.stackexchange.com/questions/2706701/neumann-boundary-conditions-in-finite-difference

Download Julia using 
`curl -fsSL https://install.julialang.org | sh`
(on Linux), otherwise download julia from https://julialang.org/downloads/

Then start julia in the same directory as the `Project.toml` is located.
Start install the dependencies by using using `julia --project=. -e "using Pkg;Pkg.instantiate()"` and waiting a few seconds.
Afterwards you can run all the exercises using for example `julia --project=. src/ReactorSimulation.jl` or `julia --project=. .\src\ReactorSimulation.jl` (on Windows).

If you are on Windows: **Use powershell**, CMD does not have support for unicode output.