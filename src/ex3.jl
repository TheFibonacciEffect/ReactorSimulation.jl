using SparseArrays
using IterativeSolvers
using Unitful
using Plots
using LinearAlgebra
using Statistics
using LsqFit

function fit_cos(x,y)
    f(x,p) = @. p[1] * cos(p[2]*x)
    p0 = [1,0.01]
    curve_fit(f,x,y,p0)
end

function get_B(x,y)
    fit_cos(x,y).param[2]#maybe sqrt
end

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
function calculate_k(a,d)
    # thermal diffusion length squared, L_th^2
    L_th_sq = D_slow / Σa_s
    # geometric buckling B_g^2
    @show H = a  # Height of the reactor core
    @show B_g_sq = (π / (H+2d))^2 
    @show η = (νΣf_f + νΣf_s) / (Σa_f + Σa_s)
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

function analytical_reactor_without_reflector(x,A)
    B = π/(2a)
    A*cos(B*x)
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

function compute_flux(D_s, D_f, B, Σ_as, Σ_af, Σ_1_to_2, νΣ_ff, νΣ_fs, k, ϕ_guess)
    # Compute the determinant of A
    det_A = (D_s * B - Σ_af - Σ_1_to_2) * (D_f * B - Σ_as)

    # Compute the inverse of A
    A_inv = [
        (D_s * B - Σ_as) / det_A           0
        -Σ_1_to_2 / det_A                 (D_f * B - Σ_af - Σ_1_to_2) / det_A
    ]

    A = [
        (D_f * B - Σ_af - Σ_1_to_2) 0
        Σ_1_to_2  (D_s * B - Σ_as)
    ]

    # Define matrix F
    F = - [
        νΣ_ff    νΣ_fs
        0        0
    ]

    # Compute -F * ϕ / k
    # Fϕ = -1 / k * F

    # Compute updated flux ϕ = A⁻¹ * (-Fϕ / k)
    @show M = inv(A) * F
    @show k = eigvals(M)[end]
    @show phi = eigvecs(M)[:,end]

    # Return components of ϕ
    return  phi
end



function reactor_without_reflector(dx; save = false, do_plot=false, verbose=false, max=false, error=false)
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
    fast = phi[1:n]
    slow = phi[n+1:end]
    @show B = π/(2a)
    if error
        ϕ = compute_flux(D_slow,D_fast,B,Σa_s, Σa_f,Σ12, νΣf_f, νΣf_s, calculate_k(l,0),[1, 1])
        ϕf, ϕs = ϕ./ϕ[1]
        fa = analytical_reactor_without_reflector.(x,ϕf)
        sl = analytical_reactor_without_reflector.(x,ϕs)
        plot(x,fast - fa, label="error fast")
        return plot!(x, slow - sl, label="error slow") 
        # plot(x,fa, label="error fast")
        # return plot!(x, sl, label="error slow")
    end
    @show get_B(x,fast)
    @show get_B(x,slow)
    p1 = plot(x,fast, label="fast neutrons")
    plot!(x,slow, label="slow neutrons")
end
# functions for reactor with reflector

function beta(i,D)
    # im using a beta without dx here
    # for i+1 it would be
    # return 2 * Di *Dip1/(Dip1+Di)
    return 2*D[i-1] * D[i] / (D[i] + D[i-1])
end

function left_side(D,Σa::Array,n, dx)
    println("left side with arrays")
    # maybe easiest to make a for loop to fill up the matrix, instead of thinking how to combine the vector and the matrix
    # not nessesary, because the diffision length is the same
    streaming = spzeros(n,n)
    # streaming[1,1] 
    # streaming[end,end]
    # for i = 2:n-1
    #     streaming[i,i] = - beta(i+1,D)/dx^2 - beta(i-1,D)/dx^2
    #     streaming[i,i-1] = + beta(i-1,D)/dx^2
    #     streaming[i,i+1] = + beta(i+1,D) / dx^2
    # end
    # apply_inner!(streaming, n,dx,D)
    laplace = spdiagm(-1 => 1* ones(n-1), 0 => -2 * ones(n), 1 => 1* ones(n-1))/dx^2
    streaming = D* laplace
    collision = - spdiagm(0 => Σa)
    M = streaming + collision
end

function reactor_with_reflector(dx; save = false, do_plot=false, verbose=false, max=false)
    println("------- start run for reactor with reflector ------")
    # numerical Parameters
    l = 2a + 2b
    nc = 2a ÷ dx |> Int
    nr = b ÷ dx |> Int
    nt = 2*nr + nc
    x =  range(-l/2, l/2,nt) 
    Σa_s    = [fill(0.012, nr); fill(0.06, nc)   ;   fill(0.012, nr)]
    Σa_f    = [fill(0, nr)    ; fill(0.002, nc)  ;   fill(0.0, nr)]
    # slow neutrons right hand side
    diffusion_slow = left_side(D_slow, Σa_s, nt, dx)
    # fast neutrons right hand side
    diffusion_fast = left_side(D_fast, Σa_f, nt, dx)
    # Assume Φf; Φs
    A = [
        diffusion_fast-Σ12 *I   spzeros(nt,nt)
        + Σ12 * I               diffusion_slow 
    ]
    println("A")
    display(A)
    νΣf_f_array = [zeros(nr) ; νΣf_f * ones(nc) ; zeros(nr)]
    νΣf_s_array = [zeros(nr) ; νΣf_s * ones(nc) ; zeros(nr)]
    # fast, slow
    # the minus is important
    F = - [
        spdiagm(0 => νΣf_f_array)  spdiagm(0 => νΣf_s_array)
        spzeros(nt,nt) spzeros(nt,nt)
    ]
    println("F")
    display(F)
    # make sparse matrecies dense
    A = Matrix(A)
    F = Matrix(F)
    M = inv(A) * F
    @show k = eigvals(M)[end-1]
    phi = eigvecs(M)[:,end-1]
    phi = real.(phi)
    phi = phi ./ phi[nc ÷ 2]
    p1 = plot(x,phi[1:nt], label="fast neutrons")
    plot!(x,phi[nt+1:end], label="slow neutrons")
    println(savefig("docs/figs/ex3/fist_harmonic_reflected.png"))
    @show k = eigvals(M)[end]
    phi = eigvecs(M)[:,end]
    phi = real.(phi)
    phi = phi ./ phi[nc ÷ 2]
    p1 = plot(x,phi[1:nt], label="fast neutrons")
    plot!(x,phi[nt+1:end], label="slow neutrons")
end

# run the code
reactor_without_reflector(0.1)
println(savefig("docs/figs/ex3/bare.png"))
println("file saved to: docs/figs/ex3/bare.png")
reactor_with_reflector(0.1)
println(savefig("docs/figs/ex3/reflected.png"))
println("file saved to: docs/figs/ex3/reflected.png")
reactor_without_reflector(0.1,error=true)
println(savefig("docs/figs/ex3/error.png"))
println("file saved to: docs/figs/ex3/error.png")