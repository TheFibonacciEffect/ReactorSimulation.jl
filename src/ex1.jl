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
# numerical solution
# Assuming D is constant, then D=β

using SparseArrays
using IterativeSolvers
function calculate_boundary(S,dx,n)
    BCp = 1/8 + D/(2dx)
    BCm = 1/8 - D/(2dx)
    S_array = -S/2/dx.*[ 1. ;zeros(n-1)]
    BC_left = -BCm/BCp
    BC_right = -BCp/BCm
    BC(i) = 1/(2 * (dx / (4D) + 1))
    # boundary = 1/dx^2 * [D; zeros(n-2); (2*D/((1/4+D/dx)*dx)-1)]
    left = - D + Σa*dx^2
    right = - (D + dx * BC(n) + Σa * dx^2)
    return left/dx^2, right/dx^2, S_array
end

function slab_reactor(n; save = false, do_plot=false, verbose=false, max=false)
    # numerical Parameters
    dx = x0/n
    x =  range(0u"cm", x0,n)
    Φ = slab_analytical(x0,S)
    # you can include eg. zero boundary conditions by starting and ending the diagonal with zeros
    laplace = D* spdiagm(-1 => 1* ones(n-1), 0 => -2 * ones(n), 1 => 1* ones(n-1))/dx^2
    streaming =  laplace
    # boundary condition
    left, right, Q = calculate_boundary(S,dx,n)
    collision = - spdiagm(0=> ones(n)) * Σa
    A = streaming + collision
    A[1,1] = left
    A[end,end] = right
    phi = cg(ustrip(A), ustrip(Q))
    # phi = bicgstabl(ustrip(A), ustrip(Q))
    # TODO the error is much better when using cg, but this shouldnt work, because the matrix is not symmetric
    phi = phi*unit(eltype(Q))/unit(eltype(A))
    if verbose
        @show Φ(1.05u"cm")
        @show Φ(0u"cm")
        @show Φ(x0)
        @show A[5,5]
        @show A[5,6]
        @show A[5,4]
        @show A[4,5]
        @show A[6,5]
        println("boundary conditions")
        @show A[1,1]
        @show A[1,2]
        @show A[2,1]
        @show A[end,end]
        @show A[end-1,end]
        @show A[end, end-1]
        @show unit(eltype(Q))/unit(eltype(A)) # the units seem to be correct
        println("Question 3")
        @show phi[1]
        @show phi[end]
        @show phi[Int(1.05 ÷ ustrip(dx))]
        @show phi[Int(1.05 ÷ ustrip(dx))] - Φ(1.05u"cm")
    end
    if do_plot
        plot(Φ,x, xlabel="l", ylabel="Neutron Flux Density", title="Analytical Solution", legend=false)
        if save savefig("docs/figs/ex1_analytical.png") end
        p_err = plot(x ,phi .- Φ.(x), legend=false)
        xlabel!("l")
        ylabel!("Error")
        # title!("Difference between Numerical and Analytical Solution")
        p_rel = plot(x ,(phi .- Φ.(x)) ./ Φ.(x), legend=false)
        xlabel!("l")
        ylabel!("Relative Error")
        # title!("Difference between Numerical and Analytical Solution")
        if save savefig("docs/figs/ex1_err_$(n).png") end
        return plot(p_err, p_rel)
    end
    if max return maximum(abs.(phi .- Φ.(x))) end
    return abs(sum(phi .- Φ.(x)))
end

function plot_error(n)
    err = slab_reactor.(n; save = false)
    plot(n,err, ylabel="sum err", xlabel="number of grid points")
    savefig("./docs/figs/sum_errors.png")
    err = slab_reactor.(n; save = false, max=true)
    plot(n,err, ylabel="max err", xlabel="number of grid points")
    savefig("./docs/figs/max_errors.png")
end

slab_reactor(100; verbose=true, do_plot=true)
plot_error(100:1000:10000)

calculate_boundary(S,1u"cm", 10)
# @show Q
# @show b
