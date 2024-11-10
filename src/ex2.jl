using SparseArrays
using IterativeSolvers
using Unitful
using Plots
using LinearAlgebra
using Statistics
# Phyical Parameters
# x0 = 10u"cm"
# S = 1000u"cm^-2*s^-1"
# Σt = Σa + Σs
# L = sqrt(1/(3(Σs + Σa)*Σa))
Σa = 0.01
Σs = 0.3
D = 1/(3 * ( Σs + Σa))
νΣf = 0.015
a = 20
b = 10
λ = 1/(Σa + Σs)
extr_l = 0.71*λ

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
    # first attempt
    #     # Modify left boundary (J^-(0))
    # @show D/dx^2
    # M[1:2,1:2] = [ 1/4-D/dx D/dx
    #             D/dx^2 -2D/dx^2-Σa
    # ]
    # display(M[1:2,1:2])

    # # Modify right boundary (J^+(0))
    # M[n, n] = 1/4 - D / dx
    # M[n, n-1] = D / dx

    # from other solution
    # Left boundary condition at x = 0
    # D * d^2/dx^2 ϕ(0) - Σ_a ϕ(0) = S / 2
    # Modify the first row of A to reflect this boundary condition
    # A[1, 1] = -2 * D / dx^2 - Σa
    # A[1, 2] = D / dx^2
    # A[1,1] = - (1/dx^2)*beta(2,1) - Σa
    # A[1,2] = 1/dx^2 * beta(2,1)

    # boundary conidtion from task 1
    # A[1,1] = 1/(2*(dx/(4D) + 1)) * 1/dx + Σa +  beta(2,1)/dx^2
    # A[1,2] = - beta(1-1,1)/dx^2

    # Right boundary condition at x = n (Last node)
    # J⁺(0) = ϕ(0)/4 - D * dϕ/dx |_{x=0} = 0
    # Set up the right flux condition at the last row in A
    # A[n, n-1] = -D / dx
    # A[n, n] = 1 / 4 + D / dx

    # boundary condition from task 1
    # A[n,n] = 1/(2*(dx/(4D) + 1)) * 1/dx + Σa +  beta(n-1,1)/dx^2
    # A[n,n-1] = - beta(n-1,1)/dx^2
end
    
function check_boundary(P,dx)
    return P[2:end] - D/dx * diff(P) , P[2:end] + D/dx * diff(P)
    # return - D/dx * diff(P) , + D/dx * diff(P)
end

function analytical_reactor_without_reflector(x)
    B = π/(2a + 2*extr_l)
    cos(B*x)
end

function reactor_without_reflector(dx; save = false, do_plot=false, verbose=false, max=false)
    println("------- start run ------")
    # numerical Parameters
    l = 2a + 2extr_l
    n = l ÷ dx -1 |> Int
    x =  range(-l/2, l/2,n)
    # you can include eg. zero boundary conditions by starting and ending the diagonal with zeros
    laplace = D* spdiagm(-1 => 1* ones(n-1), 0 => -2 * ones(n), 1 => 1* ones(n-1))/dx^2
    display(laplace[1:3,1:3])
    @assert laplace[1,2] ≈ D/dx^2
    @show D/dx^2
    streaming =  laplace
    # streaming[1,1] = 
    collision = - spdiagm(0=> ones(n)) * Σa
    M = streaming + collision
    reactor_length = n
    fission = νΣf * spdiagm(0 => ones(reactor_length))
    F = fission
    apply_boundary_conditions!(M,D,dx,n)
    display(M[1:3,1:3])
    # initial guesses
    k = 1
    P = ones(n)
    P = jacobi_iteration!(M,F,P,k)
    Pl, Pr = check_boundary(P,dx)
    p1 = plot(x,P, label="numerical")
    plot!(x,analytical_reactor_without_reflector.(x), label="analytical")
    p_err = plot(x, P .- analytical_reactor_without_reflector.(x), label="error")
    p_boundary_conditions = plot(x[2:end],Pl,title="Boundary Conditions")
    plot!(x[2:end],Pr)
    @show P[1], P[end]
    @show P[Int((1.05 + l/2) ÷ ustrip(dx))] - 9.9723e-1
    @show analytical_reactor_without_reflector(1.05) - 9.9723e-1
    JL = (P[2] - P[1])/dx
    JR = (P[end-1] - P[end])/dx
    @show calculate_k(a,extr_l) |> round5
    println("Flux at the Boundary")
    analytical_reactor_without_reflector(-a/2) |> round5
    analytical_reactor_without_reflector(a/2) |> round5
    plot(p1,p_err)
end
# reactor_without_reflector(0.1)

function reactor_reflector(dx; save = false, do_plot=false, verbose=false, max=false)
    # numerical Parameters
    l = 2a+2b+2extr_l
    n = l ÷ dx - 1 |> Int
    x =  range(-l/2, l/2,n)
    # you can include eg. zero boundary conditions by starting and ending the diagonal with zeros
    laplace = D* spdiagm(-1 => 1* ones(n-1), 0 => -2 * ones(n), 1 => 1* ones(n-1))/dx^2
    display(laplace[1:3,1:3])
    @assert laplace[1,2] ≈ D/dx^2
    @show D/dx^2
    streaming =  laplace
    # streaming[1,1] = 
    collision = - spdiagm(0=> ones(n)) * Σa
    M = streaming + collision
    moderator_length = n * b/l |> round |> Int
    reactor_length = n-2*moderator_length
    fission = νΣf * spdiagm(0 => [zeros(moderator_length) ; ones(reactor_length); zeros(moderator_length)])
    F = fission
    apply_boundary_conditions!(M,D,dx,n)
    display(M[1:3,1:3])
    # initial guesses
    k = 1
    # P = ones(n)
    P = rand(n)
    # P = jacobi_iteration!(M,F,k, P)
    P = jacobi_iteration!(M,F, P,k)
    # @show maximum(eigvals(inv(Matrix(M))*F))
    # @show minimum(eigvals(inv(Matrix(M))*F))
    # @show minimum(eigvals(M ./ F))

    # return plot(eigvecs(inv(Matrix(M))*F)[1,:])
    Pl, Pr = check_boundary(P,dx)
    p1 = plot(x,P)
    p2 = plot(x[2:end],Pl,title="Boundary Conditions")
    plot!(x[2:end],Pr)
    @show P[1], P[end]
    @show P[Int((1.05 + l/2) ÷ ustrip(dx))]
    @show JL = round((P[2] - P[1])/dx, digits=5)
    @show JR = round((P[end-1] - P[end])/dx, digits=5)
    plot(p1)
end

function jacobi_iteration!(M,F,P,k0)
    eps = 0.01
    err = Inf
    maxitter = 2000
    i = 0
    while abs(err) > eps && i < maxitter
        P2 = M \ (F*P)
        # P2 = - F \ (M*P) / k0 # doesnt work for the one with reflector, because F is singular in the reflector
        @show k1 = 1/mean(P2 ./ P) # production / absorbtion
        # err = std(P2 ./ P)
        err = abs(k1-k0)/k1
        P = P2/P2[end÷2]
        k0 = k1
        @show i += 1
    end
    return P
end

function jacobi_iteration_lecture!(M,F,k,P)
    # todo the error osscilattes strangely
    eps = 0.01
    err = Inf
    maxitter = 20
    i = 0
    while abs(err) > eps && i < maxitter
        MP = 1/k * F * P
        P1 = M \ MP
        k1 = k*(norm(F*P1)/norm(F*P))^2
        err = (k1 - k)/k1
        P = P1
        # @show k = k1 + (k1 - k)/1000
        @show k = k1
        i += 1
    end
    @show round(k,digits=5)
    return P/P[end÷2]
end

reactor_without_reflector(0.1)
savefig("docs/figs/ex2/bare.png")
reactor_reflector(0.1)
savefig("docs/figs/ex2/reflector.png")
reactor_reflector(0.01)
