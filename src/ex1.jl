using Unitful
using Plots
# using LaTeXStrings
Σa = 0.02u"cm^-1"
Σs = 4u"cm^-1"
Σt = Σa + Σs
# do these units make sense?
S = 1000u"cm^-2*s^-1"
Sleft = S/2
x0 = 10u"cm"
n(x) = (Sleft
        *exp(-Σa*x) # not absorbed
        * exp(-0.5*Σs*x) # scattering propability
            )  # half of the particles go to the left
nx0 = n(x0)
println("n(x0) = $nx0")
n0 = n(0u"cm")
println("n(0) = $(n0)")
x = (0:0.1:20)u"cm"
plot(n,x, xlabel="l", ylabel="Neutron Flux Density")

