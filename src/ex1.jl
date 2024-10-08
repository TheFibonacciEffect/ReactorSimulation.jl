using Unitful
using Plots
# using LaTeXStrings
Σa = 0.02u"cm^-1"
Σs = 4u"cm^-1"
Σt = Σa + Σs
# do these units make sense?
S = 1000u"cm^-2*s^-1"
Sleft = S
x0 = 20u"cm"
D = 1/(3 * ( Σs + Σa))
L = sqrt(1/
            (3(Σs + Σa)*Σa)
            )
Φ(x) = (
    Sleft * L / 2 / D * sinh( ( (x0-2 * abs(x))/(2L)))/(cosh(x0/2/L))
)
x = (0:0.1:10)u"cm"
plot(Φ,x, xlabel="l", ylabel="Neutron Flux Density")
Φ(1.05u"cm")
