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
b = 0 # Set to 0 for testing

round5(x) = round(x, digits=5)
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

function analytical_reactor_without_reflector(x)
    B = π/(2a)
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
    println("------- start run for bare reactor ------")
    # numerical Parameters
    l = 2a
    n = l ÷ dx -1 |> Int
    x =  range(-l/2, l/2,n)
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
end
reactor_without_reflector(0.1)
savefig("docs/figs/ex3/bare.png")

# functions for reactor with reflector

function beta(i,D)
    # im using a beta without dx here
    # for i+1 it would be
    # return 2 * Di *Dip1/(Dip1+Di)
    return 2*D[i-1] * D[i] / (D[i] + D[i-1])
end


function apply_inner!(A,n,dx,D)
    for i = 2:n-1
        A[i,i] = - beta(i+1,D)/dx^2 + beta(i-1,D)/dx^2
        A[i,i-1] = + beta(i-1,D)/dx^2
        A[i,i+1] = + beta(i+1,D) / dx^2
    end
end

function apply_boundary_conditions!(A, D, dx, n)
    # TODO
    A[1,1] 
    A[end,end]
end

function left_side(D::Array,Σa::Array,n, dx)
    # TODO
    # maybe easiest to make a for loop to fill up the matrix, instead of thinking how to combine the vector and the matrix
    streaming = spzeros(n,n)
    apply_inner!(streaming, n,dx,D)
    collision = - Σa * I
    
    M = streaming + collision
end

function reactor_with_reflector(dx; save = false, do_plot=false, verbose=false, max=false)
    println("------- start run for reactor with reflector ------")
    # numerical Parameters
    l = 2a + 2b
    nc = 2a ÷ dx -1 |> Int
    nr = b ÷ dx -1 |> Int
    x =  range(-l/2, l/2,nc+ 2nr)
    # slow neutrons right hand side
    diffusion_slow = left_side(D_slow, Σa_s, nc, dx)
    # fast neutrons right hand side
    diffusion_fast = left_side(D_fast, Σa_f, nc, dx)
    # Assume Φf; Φs
    rector_reflector_boundary = spzeros(nr, nc)
    rector_reflector_boundary[end,1] = -D/2 # or something like that
    # A = [
    #     diffusion_slow_reflector    reactor_reflector_boundary spzeros(nr,nc)
    #     reactor_reflector_boundary'    diffusion_fast-Σ12 *I   spzeros(nc,nc)
    #     reactor_reflector_boundary'     spzeros(nc,nc)
    #     spzeros(nr,nc)                + Σ12 * I               diffusion_slow 
    # ]
    println("A")
    display(A)
    # fast, slow
    # the minus is important
    F = - [
        νΣf_f * I νΣf_s * I
        spzeros(nc,nc) spzeros(nc,nc)
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
    phi = phi ./ phi[nc ÷ 2]
    p1 = plot(x,phi[1:nc], label="fast neutrons")
    plot!(x,phi[nc+1:end], label="slow neutrons")
end
reactor_with_reflector(0.1)
savefig("docs/figs/ex3/bare.png")
