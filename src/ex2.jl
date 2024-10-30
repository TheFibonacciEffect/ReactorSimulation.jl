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
function reactor_reflector(dx; save = false, do_plot=false, verbose=false, max=false)
    # numerical Parameters
    l = a+2b
    n = l ÷ dx |> Int
    x =  range(0, l,n)
    # you can include eg. zero boundary conditions by starting and ending the diagonal with zeros
    laplace = D* spdiagm(-1 => 1* ones(n-1), 0 => -2 * ones(n), 1 => 1* ones(n-1))/dx^2
    streaming =  laplace
    # streaming[1,1] = 
    collision = - spdiagm(0=> ones(n)) * Σa
    M = streaming + collision
    moderator_length = n * b/l |> round |> Int
    reactor_length = n-2*moderator_length
    fission = νΣf * spdiagm(0 => [zeros(moderator_length) ; ones(reactor_length); zeros(moderator_length)])
    F = fission
    # initial guesses
    k = 1
    P = ones(n)
    jacobi_iteration!(M,F,k, P)
end

function jacobi_iteration!(M,F,k,P)
    eps = 0.01
    err = Inf
    maxitter = 1000
    i = 0
    while abs(err) > eps && i < maxitter
        MP = 1/k * F * P
        P1 = M \ MP
        @infiltrate
        k1 = k*(norm(F*P1)/norm(F*P))^2
        @show err = (k1 - k)/k1
        P = P1
        @show k = k1
        i += 1
    end
    return plot(P)
end

reactor_reflector(0.1)
