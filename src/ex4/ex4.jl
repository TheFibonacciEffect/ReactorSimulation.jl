using SparseArrays
using IterativeSolvers
using Unitful
using Plots
using LinearAlgebra
using Statistics
using LsqFit
using Optim

println("Im using Optim for the optimization, I am sorry that it takes so long to compile ...")
println("I know its a heavy Dependency, but its also really good")

function fit_cos(x,y)
    f(x,p) = @. p[1] * cos(p[2]*x)
    p0 = [1,0.01]
    curve_fit(f,x,y,p0)
end

function get_B(x,y)
    fit_cos(x,y).param[2]#maybe sqrt
end

println("Exercise 4")

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


function beta(i,D,n)
    # im using a beta without dx here
    # for i+1 it would be
    # return 2 * Di *Dip1/(Dip1+Di)
    if i == 1 || i == 0
        return 2*D[1] * D[1] / (D[1] + D[1])
    elseif  i == n || i == n
        return 2*D[end] * D[end] / (D[end] + D[end])
    end
    return 2*D[i-1] * D[i] / (D[i] + D[i-1])
end

function left_side(D,Σa::Array,n, dx, half_core)
    println("left side with arrays")
    open_boundary = false # TODO: true leads to negative values for phi[nt]
    # maybe easiest to make a for loop to fill up the matrix, instead of thinking how to combine the vector and the matrix
    # not nessesary, because the diffision length is the same
    streaming = spzeros(n,n)
    
    for i = 2:n-1
        streaming[i,i] = - beta(i+1,D,n)/dx^2 - beta(i-1,D,n)/dx^2
        streaming[i,i-1] = + beta(i-1,D,n)/dx^2
        streaming[i-1,i] = + beta(i+1,D,n) / dx^2
    end
    streaming[end,end-1] = + beta(n-1,D,n)/dx^2
    streaming[end-1,end] = + beta(n,D,n) / dx^2
    if half_core
        streaming[1,1] = -1*D[1]/dx^2
        if open_boundary
            streaming[end,end] =  -2*D[end]/dx^2 + 2*D[end]/(
                (dx/
                (2D[end])
                + 1) * dx^2
                )
        else
            streaming[end,end] = -2*D[end]/dx^2
        end
    else
        println("not half core")
        streaming[1,1] = -2*D[1]/dx^2
        streaming[end,end] = -2*D[end]/dx^2
        # streaming[1,1] = -2*D[1]/dx^2 + 2*D[1]/(
        #     (dx/
        #     (4D[1])
        #     + 1) * dx^2
        #     )
        # streaming[end,end] = -2*D[end]/dx^2 + 2*D[end]/(
        #     (dx/
        #     (4D[end])
        #     + 1) * dx^2
        #     )
    end
    # display(streaming[1:5,1:5])
    # display(streaming[end-5:end,end-5:end])
    # apply_inner!(streaming, n,dx,D)
    # laplace = spdiagm(-1 => 1* ones(n-1), 0 => -2 * ones(n), 1 => 1* ones(n-1))/dx^2
    # streaming = D* laplace
    collision = - spdiagm(0 => Σa)
    M = streaming + collision
end

function reactor_with_reflector(dx, assemblies, half_core; optimization = false, save = false, do_plot=false, verbose=false, max=false)
    println("------- start run for reactor with reflector ------")
    # physical Parameters
    a = 120
    b = 20
    Ds = [1.4824e+00  3.8138e-01  1.4854e+00  3.7045e-01  1.4850e+00  3.6760e-01  1.2000e+00  4.0000e-01]
    D1 = Ds[1:2:end]
    D2 = Ds[2:2:end]
    nusig_f  = [ 7.1695e-03  1.4038e-01  6.0022e-03  1.4267e-01  5.1128e-03  1.2765e-01  0.0000e+00  0.0000e+00]
    nusig_f_1  = nusig_f[1:2:end]
    nusig_f_2  = nusig_f[2:2:end]
    sig_a = [ 9.6159e-03  8.2153e-02  1.0577e-02  9.5616e-02  1.1109e-02  9.3004e-02  1.0000e-03  2.0000e-02]
    sig_a_1  = sig_a[1:2:end]
    sig_a_2  = sig_a[2:2:end]
    from_1 =  [1.9788e-01  1.7369e-02  1.9748e-01  1.6350e-02  1.9754e-01  1.5815e-02  2.5178e-01  2.5000e-02]
    from_2 = [1.6271e-03  7.9024e-01  1.8467e-03  8.0234e-01  1.8112e-03  8.1197e-01  0.0000e+00  8.1333e-01]
    sig_12 = from_1[2:2:end]
    sig_21 = from_2[1:2:end]
    ksig_f = [ 9.1723e-14  1.8644e-12  7.4496e-14  1.7873e-12  6.1958e-14  1.5497e-12  0.0000e+00  0.0000e+00]
    ksig_f_1 = ksig_f[1:2:end]
    ksig_f_2 = ksig_f[2:2:end]

    # numerical Parameters
    λ_fast = 1/(sig_a_1[end] )
    λ_slow = 1/(sig_a_2[end] )
    extrapolatedLength = 0.71 * λ_fast
    l = a + b + extrapolatedLength
    nc = (a+b) ÷ dx |> Int
    nt = l ÷ dx |> Int
    ass_length = nc ÷ length(assemblies)
    x =  range(0, l,nt) 
    Σa_f    = fill(NaN,nt)
    Σa_s    = fill(NaN,nt)
    D_slow = fill(NaN,nt)
    D_fast = fill(NaN,nt)
    νΣf_f_array = fill(NaN,nt)
    νΣf_s_array = fill(NaN,nt)
    ksig_f_1_array = fill(NaN,nt)
    ksig_f_2_array = fill(NaN,nt)
    Σ12 = fill(NaN,nt)
    Σ21 = fill(NaN,nt)
    @show assemblies
    for (i, ia) in enumerate(assemblies)
        for j in ((i-1)*ass_length +1):((i)*ass_length)
            D_fast[j] = D1[ia]
            D_slow[j] = D2[ia]
            νΣf_f_array[j] = nusig_f_1[ia]
            νΣf_s_array[j] = nusig_f_2[ia]
            Σa_f[j]    = sig_a_1[ia]
            Σa_s[j]    = sig_a_2[ia]
            Σ12[j] = sig_12[ia]
            Σ21[j] = sig_21[ia]
            ksig_f_1_array[j] = ksig_f_1[ia]
            ksig_f_2_array[j] = ksig_f_2[ia]
        end
    end
    # the extrapolated length is filled with reflector
    for i in nc:nt
        D_fast[i] = D1[4]
        D_slow[i] = D2[4]
        νΣf_f_array[i] = nusig_f_1[4]
        νΣf_s_array[i] = nusig_f_2[4]
        Σa_f[i]    = sig_a_1[4]
        Σa_s[i]    = sig_a_2[4]
        Σ12[i] = sig_12[4]
        Σ21[i] = sig_21[4]
        ksig_f_1_array[i] = ksig_f_1[4]
        ksig_f_2_array[i] = ksig_f_2[4]
    end
    # slow neutrons right hand side
    diffusion_slow = left_side(D_slow, Σa_s, nt, dx, half_core)
    # fast neutrons right hand side
    diffusion_fast = left_side(D_fast, Σa_f, nt, dx, half_core)
    # Assume Φf; Φs
    A = [
        diffusion_fast- spdiagm(Σ12)    spdiagm(Σ21)
        + spdiagm(Σ12)                  diffusion_slow  - spdiagm(Σ21)
    ]
    # fast, slow
    # the minus is important
    F = - [
        spdiagm(0 => νΣf_f_array)  spdiagm(0 => νΣf_s_array)
        spzeros(nt,nt) spzeros(nt,nt)
    ]
    # make sparse matrecies dense
    A = Matrix(A)
    F = Matrix(F)
    M = inv(A) * F
    eigv = eigvals(M)
    eigenvectors = eigvecs(M)

    for i in 0:2
        phi = eigenvectors[:,end-i]
        phi = real.(phi)
        p_nom_core = 90e6 # MW
        h = 400 # cm
        P_fiss = sum(h*dx* ksig_f_1_array .* phi[1:nt] + h*dx* ksig_f_2_array .* phi[nt+1:end])
        F = p_nom_core / P_fiss
        phi = F * phi
        @show k = eigv[end-i]
        if optimization
            boundary_index = nt - ass_length
            @show phi[nt]
            @show phi[end]
            # fast group at core + reflector boundary, slow group
            # nc is core + reflector
            # nt is total: core + reflector + exrapolated length
            return phi[nc], phi[nt + nc]
        end
        p1 = plot(x,phi[1:nt], label="fast neutrons")
        plot!(x,phi[nt+1:end], label="slow neutrons")
        p2 = twinx(p1)
        plot!(p2,x, νΣf_f_array,label="νΣfission for fast neutrons", legend=:bottomleft, color=:red, linestyle=:dash)
        plot!(p2,x, νΣf_s_array,label="νΣfission for slow neutrons", legend=:bottomleft, color=:green, linestyle=:dash)
        xlims!(0,140)
        println(savefig("docs/figs/ex4/neutrons_core_$i.png"))
        println("saved docs/figs/ex4/neutrons_core_$i.png")
    end
end


# assemblies = [4 1 1 2 2 3 3 3 3 2 2 1 1 4] # fresh fuel outside
# assemblies = fill(1,14)
# assemblies = [4 1 1 2 2 3 3 3 3 2 2 1 1 4] # fresh fuel outside
# assemblies = [4 3 3 2 2 1 1 1 1 2 2 3 3 4] # fresh fuel inside
# assemblies = [4 1 2 3 1 2 3 1 2 3 1 2 3 4] # checkerbord
# for the half core
assemblies = [3 3 2 2 1 1 4] # fresh fuel outside
assemblies = [1 1 2 2 3 3 4] # fresh fuel inside
@show reactor_with_reflector(1, assemblies, true; optimization = false)
@time reactor_with_reflector(1, assemblies, true; optimization = true)


