using SparseArrays
using IterativeSolvers
using Unitful
using Plots
using LinearAlgebra
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

function apply_boundary_conditions!(M, D, dx, n)
    # Modify left boundary (J^-(0))
    @show D/dx^2
    M[1:2,1:2] = [ 1/4-D/dx D/dx
                D/dx^2 -2D/dx^2-Σa
    ]
    display(M[1:2,1:2])

    # Modify right boundary (J^+(0))
    M[n, n] = 1/4 - D / dx
    M[n, n-1] = D / dx
end


    
function check_boundary(P,dx)
    return P[2:end] - D/dx * diff(P) , P[2:end] + D/dx * diff(P)
    # return - D/dx * diff(P) , + D/dx * diff(P)
end

function reactor_reflector(dx; save = false, do_plot=false, verbose=false, max=false)
    # numerical Parameters
    l = a+2b
    n = l ÷ dx |> Int
    x =  range(0, l,n)
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
    P = ones(n)
    P = jacobi_iteration!(M,F,k, P)
    Pl, Pr = check_boundary(P,dx)
    p1 = plot(x,P)
    p2 = plot(x[2:end],Pl)
    plot!(x[2:end],Pr)
    plot(p1,p2)
end

function jacobi_iteration!(M,F,k,P)
    eps = 0.01
    err = Inf
    maxitter = 2
    i = 0
    while abs(err) > eps && i < maxitter
        MP = 1/k * F * P
        P1 = M \ MP
        k1 = k*(norm(F*P1)/norm(F*P))^2
        err = (k1 - k)/k1
        P = P1
        k = k1
        i += 1
    end
    return P
end

reactor_reflector(0.1)
reactor_reflector(0.01)
