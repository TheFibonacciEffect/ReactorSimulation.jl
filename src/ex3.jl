using SparseArrays
using IterativeSolvers
using Unitful
using Plots
using LinearAlgebra
using Statistics
println("Exercise 3")
D_fast = 1.13
D_slow = 0.16
Σa_f = 0.002
Σa_s = 0.06
νΣf_f = 0.001
νΣf_s = 0.069
Σ12 = 0.038
a = 50
b = 10

round5(x) = round(x, digits=5)
function beta(i,j)
    return D
end

function calculate_k(a,d)
    # thermal diffusion length squared, L_th^2
    L_th_sq = D / Σa
    # geometric buckling B_g^2
    @show H = a  # Height of the reactor core
    @show B_g_sq = (π / (H+2d))^2 
    @show η = νΣf / Σa
    # Set f ≈ 1 and p ≈ 1 as approximations
    f = 1
    p = 1
    # non-leakage probability P_TNL
    # TODO this value is very small ≈ 0.001
    P_TNL = 1 / (1 + L_th_sq * B_g_sq)
    
    # Step 7: Calculate k
    k = η * f * p * P_TNL
    return k
end

function apply_boundary_conditions!(A, D, dx, n)
    # nothing to do because boundary is 0
end

get_buckling() = π/(2a + 2*extr_l)

function analytical_reactor_without_reflector(x)
    B = π/(2a + 2*extr_l)
    cos(B*x)
end

# function right_side(B,D,Σa)
#     laplace = spdiagm(-1 => 1* ones(n-1), 0 => -2 * ones(n), 1 => 1* ones(n-1))/dx^2
#     streaming = - D* laplace
#     collision = +  spdiagm(0=> ones(n)) * Σa
#     M = streaming + collision
# end
function left_side(D,Σa,n, dx)
    laplace = spdiagm(-1 => 1* ones(n-1), 0 => -2 * ones(n), 1 => 1* ones(n-1))/dx^2
    streaming = D* laplace
    collision = - Σa * I
    M = streaming + collision
end

function reactor_without_reflector(dx; save = false, do_plot=false, verbose=false, max=false)
    println("------- start run ------")
    # numerical Parameters
    l = 2a + 2extr_l
    n = l ÷ dx -1 |> Int
    x =  range(-l/2, l/2,n)
    B = get_buckling() # assume fast and slow buckling to be the same
    # slow neutrons right hand side
    diffusion_slow = left_side(D_slow, Σa_s, n, dx)
    # fast neutrons right hand side
    diffusion_fast = left_side(D_fast, Σa_f, n, dx)
    # Assume Φf; Φs
    A = [
        diffusion_fast-Σ12 *I   spzeros(n,n)
        + Σ12 * I               diffusion_slow 
    ]
    println("A")
    display(A)
    # fast, slow
    # the minus is important
    F = - [
        νΣf_f * I νΣf_s * I
        spzeros(n,n) spzeros(n,n)
    ]
    println("F")
    display(F)
    # make sparse matrecies dense
    A = Matrix(A)
    F = Matrix(F)
    M = inv(A) * F
    @show k = eigvals(M)[end]
    phi = eigvecs(M)[:,end]
    phi = real.(phi)
    phi = phi ./ phi[n ÷ 2]
    p1 = plot(x,phi[1:n], label="fast neutrons")
    plot!(x,phi[n+1:end], label="slow neutrons")
    
    # reactor_length = n
    # fission = νΣf * spdiagm(0 => ones(reactor_length))
    # F = fission
    # apply_boundary_conditions!(M,D,dx,n)
    # # initial guesses
    # k = 1
    # P = ones(n)
    # P = jacobi_iteration!(M,F,P,k)
    # Pl, Pr = check_boundary(P,dx)
    # p1 = plot(x,P, label="numerical")
    # plot!(x,analytical_reactor_without_reflector.(x), label="analytical")
    # p_err = plot(x, P .- analytical_reactor_without_reflector.(x), label="error")
    # p_boundary_conditions = plot(x[2:end],Pl,title="Boundary Conditions")
    # plot!(x[2:end],Pr)
    # @show P[1], P[end]
    # @show P[Int((1.05 + l/2) ÷ ustrip(dx))] - 9.9723e-1
    # @show analytical_reactor_without_reflector(1.05) - 9.9723e-1
    # JL = (P[2] - P[1])/dx
    # JR = (P[end-1] - P[end])/dx
    # @show calculate_k(a,extr_l) |> round5
    # println("Flux at the Boundary")
    # analytical_reactor_without_reflector(-a/2) |> round5
    # analytical_reactor_without_reflector(a/2) |> round5
    # plot(p1,p_err)
end
# reactor_without_reflector(0.1)

reactor_without_reflector(0.1)
savefig("docs/figs/ex3/bare.png")
