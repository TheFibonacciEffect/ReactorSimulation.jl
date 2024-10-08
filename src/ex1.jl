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
dx = 0.1u"cm"
n = x0/ dx |> Int
# you can include eg. zero boundary conditions by starting and ending the diagonal with zeros
laplace = spdiagm(-1 => ones(n-1), 0 => -2 * ones(n-1), 1 => ones(n-1))/dx^2
S = zeros(n)
collision = - 1/L*spdiagm(0=> ones(n))
streaming =  D* laplace
C = S/2
BC_right = 1/(2*dx/4D+1)
boundary = [BC_left; zeros(n-2); 1/(2*BC_right-1)]
phi = cg(streaming + collision, S)
plot(phi)

