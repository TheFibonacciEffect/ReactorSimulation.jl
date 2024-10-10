using Unitful
using Plots

# Phyical Parameters
x0 = 10u"cm"
S = 1000u"cm^-2*s^-1"
Σa = 0.02u"cm^-1"
Σs = 4u"cm^-1"
Σt = Σa + Σs
L = sqrt(1/(3(Σs + Σa)*Σa))
D = 1/(3 * ( Σs + Σa))

# numerical Parameters
dx = 2u"cm"
n = x0/ dx |> Int

# analytical solution
function slab_analytical(x0,S)
    extrapolated_length = 2/3(Σs + Σa)
    a = 2x0  + 2*extrapolated_length
    Φ(x) = (
        S * L / 2 / D * sinh( ( (a-2 * abs(x))/(2L)))/(cosh(a/2/L))
    )
    return Φ
end
x = (0:0.1:10)u"cm"
Φ = slab_analytical(x0,S)
plot(Φ,x, xlabel="l", ylabel="Neutron Flux Density")
Φ(1.05u"cm")
Φ(0u"cm")
Φ(x0)

# numerical solution
# Assuming D is constant, then D=β

using SparseArrays
using IterativeSolvers
function calculate_boundary(dx,n)
    BCp = 1/8 + D/(2dx)
    BCm = 1/8 - D/(2dx)
    S_array = - S/(2BCm*dx^2).*[ 1. ;zeros(n-1)]
    BC_left = -BCm/BCp
    BC_right = -BCp/BCm
    boundary = 1/dx^2 * [BC_left; zeros(n-2); BC_right]
    return S_array,boundary
end

calculate_boundary(dx,2)

# you can include eg. zero boundary conditions by starting and ending the diagonal with zeros
laplace = spdiagm(-1 => ones(n-1), 0 => -2 * ones(n), 1 => ones(n-1))/dx^2
streaming =  laplace
# boundary condition
Q, boundary = calculate_boundary(dx)
collision = - 1/L^2*spdiagm(0=> ones(n))
A = streaming + collision + spdiagm(0 => boundary)
phi = cg(ustrip(A), ustrip(Q))
unit(eltype(Q))/unit(eltype(A)) # the units seem to be correct
plot(phi)

