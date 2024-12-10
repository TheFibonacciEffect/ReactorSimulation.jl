using Unitful
using Plots
using Test
# Phyical Parameters
x0 = 10
S = 1000
Σa = 0.02
Σs = 4
Σt = Σa + Σs
L = sqrt(1/(3(Σs + Σa)*Σa))
D = 1/(3 * ( Σs + Σa))

function beta(i,j)
    return D
end

function round5(x)
    return round(x, digits=5)
end

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
using DelimitedFiles
using SparseArrays

function apply_boundary_conditions!(A, S, dx, n)
    # Left boundary condition at x = 0
    # D * d^2/dx^2 ϕ(0) - Σ_a ϕ(0) = S / 2
    # Modify the first row of A to reflect this boundary condition
    # A[1, 1] = -2 * D / dx^2 - Σa
    # A[1, 2] = D / dx^2
    A[1,1] = - (- (1/dx^2)*beta(2,1) - Σa)
    A[1,2] = -  1/dx^2 * beta(2,1)
    # Right boundary condition at x = n (Last node)
    # J⁺(0) = ϕ(0)/4 - D * dϕ/dx |_{x=0} = 0
    # Set up the right flux condition at the last row in A
    # A[n, n-1] = -D / dx
    # A[n, n] = 1 / 4 + D / dx
    A[n,n] = 1/(2*(dx/(4D) + 1)) * 1/dx + Σa +  beta(n-1,1)/dx^2
    # streaming[end,end] =  -2*D[end]/dx^2 + 2*D[end]/(
    #     (dx/
    #     (2D[end])
    #     + 1) * dx^2
    #     )
    A[n,n-1] = - beta(n-1,1)/dx^2
    # Adjust Q vector to include the source term
    Q = zeros(n)
    Q[1] =  S / 2/dx

    return Q
end

function apply_inner!(A,n,dx)
    for i = 2:n-1
        A[i,i] = beta(i+1,1)/dx^2 + beta(i-1,1)/dx^2 + Σa
        A[i,i-1] = - beta(i-1,1)/dx^2
        A[i,i+1] = -beta(i+1,1) / dx^2
    end
end

function slab_reactor(n; save = false, do_plot=false, verbose=false, max=false, sum=false)
    # numerical Parameters
    dx = x0/n
    x =  dx:dx:x0
    Φ = slab_analytical(x0,S)
    # you can include eg. zero boundary conditions by starting and ending the diagonal with zeros
    # laplace = D* spdiagm(-1 => 1* ones(n-1), 0 => -2 * ones(n), 1 => 1* ones(n-1))/dx^2
    # streaming =  laplace
    # collision = - spdiagm(0=> ones(n)) * Σa
    # A = streaming + collision

    A = spzeros(n,n)
    apply_inner!(A,n,dx)
    # boundary condition
    Q = apply_boundary_conditions!(A,S,dx,n)
    # phi = cg(ustrip(A), ustrip(Q))
    # phi = bicgstabl(ustrip(A), ustrip(Q))
    # phi =  ustrip(A) \ ustrip(Q)
    phi =  A \ ustrip(Q)
    # TODO the error is much better when using cg, but this shouldnt work, because the matrix is not symmetric
    phi = phi*unit(eltype(Q))/unit(eltype(A))
    if verbose
        @test A*phi ≈ Q atol=1e-6
        @test round5(Φ(1.05))  ≈ 7.3299e3 atol=1e-1
        @show round5(Φ(0))
        @show round5(Φ(x0))
        # print("A: ")
        # display(A)
        # writedlm("A.txt", A, " ")
        A_noah = readdlm("src/A.csv", ';', Float64, '\n')
        @test maximum(abs.(-A_noah .- A)) ≈ 0 atol = 1e-3
        display(A[10:11,10:11])
        println("boundary conditions")
        print("left A: ")
        display(A[1:2,1:2])
        print("right: ")
        display(A[end-1:end,end-1:end])
        @show unit(eltype(Q))/unit(eltype(A)) # the units seem to be correct
        println("Exercise 1 Question 3")
        @testset "Exercise 1 Question 3" begin
        @show round5(phi[1])
        @show round5(phi[end])
        include("matA.jl")
        println("Testing if the Matrix that Matieu gave me is the same as the one I calculated")
        @test maximum(abs.(A_matieu .- A)) ≈ 0 atol = 1e-3
        # test middle
        @test round5(A[4,4])  ≈ 16.604 atol = 1e-3
        @test round5(A[4,3]) ≈ -8.2919e0 atol = 1e-3
        @test round5(A[3,4]) ≈ -8.2919e0 atol = 1e-3
        # test left
        @test round5(A[1,1])  ≈ 8.3118740e+00 atol = 1e-3
        @test round5(A[1,2]) ≈  -8.2918740e+00 atol = 1e-3
        @test round5(A[2,1]) ≈  -8.2918740e+00 atol = 1e-3
        # test right
        @test round5(A[end,end])   ≈ 1.2153595e+01 atol = 1e-3
        @test round5(A[end,end-1]) ≈ - 8.2919e0 atol = 1e-3
        @test round5(A[end-1,end]) ≈ - 8.2919e0 atol = 1e-3
        @test Q[1] ≈ 5.0000000e+03 atol=1e-3
        @test_broken round5(phi[Int(1.05 ÷ ustrip(dx))]) ≈ 7.3299e3 atol=1e-3
        @test_broken round5(phi[Int(1.05 ÷ ustrip(dx))] - Φ(1.05)) ≈ 3.7901e-1 atol=1e-3
        @show round5(phi[1] - Φ(0))
        end #testset
    end
    if do_plot
        p_ana =  plot(Φ,x, xlabel="l", ylabel="Neutron Flux Density", title="Analytical Solution", legend=false)
        p_num = plot(x ,phi, legend=false, title="Numerical Solution")
        if save savefig("docs/figs/ex1/ex1_analytical.png") end
        p_err = plot(x ,phi .- Φ.(x), legend=false)
        xlabel!("l")
        ylabel!("Error")
        title!("Difference between Numerical and Analytical Solution")
        if save savefig("docs/figs/ex1/ex1_err_$(n).png") end
        p_rel = plot(x ,(phi .- Φ.(x)) ./ Φ.(x), legend=false)
        xlabel!("l")
        ylabel!("Relative Error")
        # title!("Difference between Numerical and Analytical Solution")
        return plot(p_err, p_rel, p_ana, p_num, layout=(2,2))
    end
    if max return maximum(abs.(phi .- Φ.(x))) end
    if sum return abs(sum(phi .- Φ.(x))) end
    return abs(phi[n] - Φ.(x0))
end

function plot_error(n)
    err = slab_reactor.(n; save = false)./n
    ps = plot(n,err, ylabel="error at x0", xlabel="number of grid points", yscale=:log10, xscale=:log10, marker=:o)
    savefig("./docs/figs/ex1/err_x0.png")
    err = slab_reactor.(n; save = false, max=true)
    pm = plot(n,err, ylabel="max err", xlabel="number of grid points", yscale=:log10, xscale=:log10, marker=:o)
    savefig("./docs/figs/ex1/max_errors.png")
    plot(pm,ps)
end

@time slab_reactor(100; verbose=true, do_plot=true, save=true)
@time slab_reactor(100000; verbose=false, do_plot=false)
plot_error(Int.(round.(10 .^ (1:1:5))))

# @show Q
# @show b
