println("Importing Libraries")
using DifferentialEquations
using Plots
using Plots.Measures
using Test
println("Libraries Imported, starting exercise 5")
println("Setting Decault Plotting")
default(left_margin=10mm, yscale=:log10)
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
λ =      Dict("U-235" => 0, "U-238" => 0, "Pu-239" => 0, "X" => log(2)/(9*60^2), "Y" => 0)
# Fission yields (FY)
fy = Dict("U-235" => 0. , "U-238" => 0. , "Pu-239" => 0. , "X" => 0.06, "Y" => 1.94)

Φ = 1e13 # * (60^2*24) # to neutrons per barn per day
b = 1e-24 #barns in cm^2 
κ = 3.2e-11

# V = pi * (0.4cm)^2 * 400cm = 2.010619298E5 mm³
# rho = 10.4g/(cm^3)
# rho × V ≈ 2.091044070 kg
M = 264*2.091044070 #kg
kg_u = 6.0221408E26
N = M * kg_u / ((238*0.97 + 235 * 0.03) + 2*16)

# N = 2.431E24
# 20cm x 20cm x 400cm
V = 20^2*400
ρ = N/V
N5 = ρ*0.03
N8 = ρ*0.97
# The initial concentration of U-235 is 2.31 10^20𝑐𝑚−3 .T
@test N5 ≈ 2.31E20 atol=1E18

function betman_eq!(du,u,p,t)
    #= 5 =# du[1] = Φ*((- sig_c["U-235"]*b -  sig_f["U-235"]*b)*u[1] )
    #= 8 =# du[2] = Φ*((- sig_c["U-238"]*b -  sig_f["U-238"]*b)*u[2] )
    #= 9 =# du[3] = Φ*((- sig_c["Pu-239"]*b - sig_f["Pu-239"]*b)*u[3] + sig_c["U-238"]*b*u[2] )
    #= X =# du[4] = -Φ*sig_c["X"]*b*u[4] + Φ* sig_f["U-235"]*b * u[1] *fy["X"] - λ["X"]*u[4]
    #= Y =# du[5] = -Φ*sig_c["Y"]*b*u[5] + Φ* sig_f["U-235"]*b * u[1] *fy["Y"] - λ["Y"]*u[5]
end

function analytical_solution(t, u, Φ,  fy, λ, b)
    # Extract the initial concentrations
    u1_0, u2_0, u3_0, u4_0, u5_0 = u
    
    # Calculate the decay constants (including the multiplying factor b)
    κ_U235 = Φ * (sig_c["U-235"] + sig_f["U-235"]) * b
    κ_U238 = Φ * (sig_c["U-238"] + sig_f["U-238"]) * b
    κ_Pu239 = Φ * (sig_c["Pu-239"] + sig_f["Pu-239"]) * b
    κ_X = Φ * sig_c["X"] * b + λ["X"]
    κ_Y = Φ * sig_c["Y"] * b + λ["Y"]

    # Solve for each concentration at time t
    
    # u[1](t) = C5(t)
    u1_t = u1_0 * exp(-κ_U235 * t)
    
    # u[2](t) = C8(t)
    u2_t = u2_0 * exp(-κ_U238 * t)
    
    # u[3](t) = C9(t)
    u3_t = u3_0 * exp(-κ_Pu239 * t) +
           (Φ * sig_c["U-238"] * b * u2_0) / (κ_Pu239 - κ_U238) * (exp(-κ_U238 * t) - exp(-κ_Pu239 * t))
    
    # u[4](t) = CX(t)
    u4_t = u4_0 * exp(-κ_X * t) + 
            (Φ * sig_f["U-235"] * b * u1_0 * fy["X"] * b) / (κ_X - κ_U235) * (exp(-κ_U235 * t) - exp(-κ_X * t))
    
    # u[5](t) = CY(t)
    u5_t = u5_0 * exp(-κ_Y * t) + 
            (Φ * sig_f["U-235"] * b * u1_0 * fy["Y"] * b) / (κ_Y - κ_U235) * (exp(-κ_U235 * t) - exp(-κ_Y * t))
    
    return [u1_t, u2_t, u3_t, u4_t, u5_t]
end

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

# testset
@testset "Check Matrix" begin
    A = to_matrix(betman_eq!,5)
    @test  A[1,1] ≈ -6.8e-10 atol=1e-12
    # -7.13930E-05
    @test A[4,1] ≈ 3.36000E-11  atol=1e-12
    @test -7.13930E-05 ≈ A[4,4] atol=1e-7
    @test A[5,5] ≈ 5e-10 atol=1e-7
end

function plotsol(sol)
    plot()
    for (i, label) in enumerate(["U5" "U8" "P9" "X" "Y"])
        # plot!(sol.t, getindex.(sol.u,i), label=label, yscale=:log)
        plot!(sol.t[2:end] ./ seconds_per_day, getindex.(sol.u,i)[2:end], label=label, yscale=:log10)
    end
    plot!()
end


function error(dt)
    prob = ODEProblem(betman_eq!, [N5,N8,0,0,0],(0,1.),dt=dt)
    # prob = ODEProblem(betman_eq!, [1,0,0,0,0],(0,1.), dtmax=0.01)
    sol = solve(prob, Euler)
    return sol.u[end]
end

function plot_analytical(sol)
    plot()
    sola = zero(sol.u)
    for i in eachindex(sol.t)
        sola[i] = analytical_solution(sol.t[i],u0, Φ,  fy, λ, b )
    end
    sola = hcat(sola[2:end]...)
    for (i, label) in enumerate(["U5" "U8" "P9" "X" "Y"])
        # err = getindex.(sol.u,i)[2:end] - sola[i,:]
        err =  sola[i,:]
        # plot!(sol.t, getindex.(sol.u,i), label=label, yscale=:log)
        # plot!(sol.t[2:end], err , label=label, yscale=:log10)
        plot!(sol.t[2:end], err , label=label)
    end
    title!("Analytical Solution")
    plot!()
end

function euler(u0,A,t_end, dt)
    t = 0:dt:t_end
    u = copy(u0)
    us = zeros(length(t),length(u0))
    for i in 1:length(t)
        u += A*u*dt
        us[i,:] .= u
    end
    return t,us
end


function error_euler()
    dt = 60
    t_end = 365*seconds_per_day
    t, u = euler(u0,to_matrix(betman_eq!,5),t_end,dt)
    sola = zeros(length(t),length(u0))
    for i in 1:length(t)
        sola[i,:] = analytical_solution(t[i],u0, Φ,  fy, λ, b )
    end
    err = abs.(u .- sola)
    rel_err = err ./ sola
    return err, rel_err,  0:dt:t_end
end

seconds_per_day  = 24* 60^2

u0 = [N5,N8,0,0,0]
prob = ODEProblem(betman_eq!, u0,(0,365 * seconds_per_day),dtmax=60)
# prob = ODEProblem(betman_eq!, [1,0,0,0,0],(0,1.), dtmax=0.01)
println("Solving ODE")
sol = solve(prob)
# plot(sol.t,sol.u', label=["U5" "U8" "P9" "X" "Y"])
plot(sol)
plotsol(sol)


println("Solving ODE using Euler")
err, rel_err , t = error_euler()
plot(t/seconds_per_day ,err, yscale=:log10)
println(savefig("docs/figs/ex5/err_euler.png"))
plot(t/seconds_per_day ,rel_err, yscale=:log10)
println(savefig("docs/figs/ex5/err_euler_relative.png"))


plot_analytical(sol)
println(savefig("docs/figs/ex5/ode_solution.png"))
t, u = euler(u0,to_matrix(betman_eq!,5),365*seconds_per_day,60)
# getindex.(sol.u,i)[2:end]
plot()
for (i, label) in enumerate(["U5" "U8" "P9" "X" "Y"])
    plot!(t/seconds_per_day, u[:,i], label=label, yscale=:log10)
end
plot!()
# plot(t/seconds_per_day, u, yscale=:log10)
println(savefig("docs/figs/ex5/euler.png"))