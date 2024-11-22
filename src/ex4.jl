using DifferentialEquations
using Plots
# all calculations are done in days and barns!

# Cross sections in barns
sig_c = Dict("U-235" => 12, "U-238" => 4, "Pu-239" => 81, "X" => 5e6, "Y" => 50)
sig_s = Dict("U-235" => 10, "U-238" => 10, "Pu-239" => 10, "X" => 10, "Y" => 10)
sig_f = Dict("U-235" => 56, "U-238" => 1, "Pu-239" => 144, "X" => 0.0, "Y" => 0.0)

# Average number of neutrons per fission
nu = Dict("U-235" => 2.44, "U-238" => 2.79, "Pu-239" => 2.87, "X" => 0.0, "Y" => 0.0)

# Energy release per fission (kappa in MeV)
kappa = Dict("U-235" => 3.24e-11, "U-238" => 3.32e-11, "Pu-239" => 3.33e-11, "X" => 0.0, "Y" => 0.0)

# Half-life (T1/2)
t_half = Dict("U-235" => "stable", "U-238" => "stable", "Pu-239" => "stable", "X" => "9h", "Y" => "stable")
λ =      Dict("U-235" => 0, "U-238" => 0, "Pu-239" => 0, "X" => log(2)/(9), "Y" => 0)
# Fission yields (FY)
fy = Dict("U-235" => 0. , "U-238" => 0. , "Pu-239" => 0. , "X" => 0.06, "Y" => 1.94)

b = 1e-24 #barns in cm^2 
κ = 3.2e-11
println("Cross-section (capture) for Pu-239: ", sig_c["Pu-239"])
function betman_eq!(du,u,p,t)
    Φ = 1e15*b * (60^2*24) # to neutrons per barn per day
    #= 5 =# du[1] = Φ*((- sig_c["U-235"] -  sig_f["U-235"])*u[1] )
    #= 8 =# du[2] = Φ*((- sig_c["U-238"] -  sig_f["U-238"])*u[2] )
    #= 9 =# du[3] = Φ*((- sig_c["Pu-239"] - sig_f["Pu-239"])*u[3] + sig_c["U-238"]*u[2] )
    #= X =# du[4] = -Φ*sig_c["X"]*u[4] + Φ* sig_f["U-235"] * u[1] *fy["X"] - λ["X"]*u[4]
    #= Y =# du[5] = -Φ*sig_c["Y"]*u[5] + Φ* sig_f["U-235"] * u[1] *fy["Y"] - λ["Y"]*u[5]
end

M = 264*2.091044070 #kg
kg_u = 6.0221408E26
N = M * kg_u / ((238*0.97 + 235 * 0.03) + 2*12)

# N = 2.431E24
# TODO
# 20cm x 20cm x 200cm
V = 20^2*400
ρ = N/V
N5 = ρ*0.03
N8 = ρ*0.97
function to_matrix(f,n)
    A = zeros(n,n)
    for i in 1:n
        en = zeros(n)
        en[i] += 1
        v = zeros(n)
        f(v,en,0,0)
        A[:,i] .= v
    end
    return A
end


# u0 = 
prob = ODEProblem(betman_eq!, [N5,N8,0,0,0],(0,365.),dtmax=0.1)
# prob = ODEProblem(betman_eq!, [1,0,0,0,0],(0,1.), dtmax=0.01)
sol = solve(prob)
# plot(sol.t,sol.u', label=["U5" "U8" "P9" "X" "Y"])
plot(sol)
function plotsol(sol)
plot()
for (i, label) in enumerate(["U5" "U8" "P9" "X" "Y"])
    # plot!(sol.t, getindex.(sol.u,i), label=label, yscale=:log)
    plot!(sol.t[2:end], getindex.(sol.u,i)[2:end], label=label, yscale=:log10)
end
plot!()
end
plotsol(sol)
# to_matrix(betman_eq!,5)
