using LinearAlgebra
using JuMP
using Ipopt
using MathOptInterface
using SparseArrays
using Plots

include("dynamics.jl")



"""
Create the g and H matrix for the following problem with r as the target
trajectory

    J       = -g X + 0.5 X' H X
    X_dot   = A * X + B * U
    Y       = C * X
"""
function build_allAtOnce_system(traj::trajectory)
    h = traj.dt
    r::Array{Array{Float64,1},1} = traj.r
    n = length(r)
    n_tot = (n - 1) * (nx + nu) # x1 vector is the initial condition

    A, B, C = buildLinearSystem()
    Q = Diagonal([100,50,10,5,0,0,100,1,100,1,1000,0.1])
    # Q = Diagonal([0,0,0,0,0,0,100,0,100,0,1000,0])
    R = Diagonal([10,0,0,0,0])
    # R = Diagonal([0,0,0,0,0])

    g = spzeros(n_tot)
    H = spzeros(n_tot, n_tot)
    A_aao = spzeros(n_tot, n_tot)
    B_aao = spzeros(n_tot, n_tot)
    C_aao = spzeros(nx * n, nx * n)
    D_aao = spzeros(n_tot)
    E_aao = spzeros(n_tot, n_tot)
    F_aao = spzeros(n_tot)

    ####### Objective ########

    # populate g vector
    for i in 2:n
        subrange    = nu + (i - 2) * (nx + nu) + 1:nu + (i - 2) * (nx + nu) + nx
        g[subrange] = (r[i]' * Q' * C' + r[i]' * Q * C) / 2
    end

    # populate H matrix
    H_sub = [R                  zeros(nu, nx);
             zeros(nx, nu)      Q               ]
    for i in 1:(n - 1)
        subrange = (i - 1) * (nx + nu) + 1:i * (nx + nu)
        H[subrange,subrange] =  H_sub
    end

    ####### Dynamics ########

    F_aao[nu + 1:nx + nu] = r[1]

    # populate A, B matrices
    submatrix = [A B]
    A_aao[nu + 1:nu + nx, 1:nu] = B
    D_aao[nu + 1:nu + nx] = A * r[1]
    for i in 2:(n - 1)
        rows = (i - 1) * (nx + nu) + nu + 1:i * (nx + nu)
        cols = (i - 2) * (nx + nu) + nu + 1:(i - 1) * (nx + nu) + nu
        A_aao[rows, cols ] = submatrix
    end

    # populate E matrice
    E_aao[1:nu,1:nu] = Diagonal(ones(nu))
    for i in 2:(n - 1)
        rows = (i - 1) * (nx + nu) + 1:((i - 1) * (nx + nu) + nu)
        cols = rows
        E_aao[rows,cols] = Diagonal(ones(nu))
        rows = i * (nx + nu) - nx + 1:i * (nx + nu)
        cols = (i - 1) * (nx + nu) - nx + 1:(i - 1) * (nx + nu)
        E_aao[rows,cols] = Diagonal(ones(nx))
    end
    return g, H, A_aao, D_aao, E_aao, F_aao
end

function main()
    # trajectoire stationnaire
    n = 100
    pts = [[0.,0.,1.],[1.,0.,1.]]
    t = [0.,5.]

    traj = make_linear_trajectory(pts, t, n)
    h = traj.dt
    println("Target trajectory")
    println(traj)
    g, H, A_aao, D_aao, E_aao, F_aao = build_allAtOnce_system(traj)
    A, B, C = buildLinearSystem()

    model = Model(with_optimizer(Ipopt.Optimizer))
    @variable(model, s[1:(n - 1) * (nx + nu)])

    x_id = []
    for i in 1:(n - 1)
        id = collect(1:nx) .+ nu .+ (i - 1) * (nx + nu)
        append!(x_id, id)
    end

    # @constraint(model, con[i=x_id], s[i] == h * (dot(A_aao[i,:], s) + D_aao[i]) + dot(E_aao[i,:], s) + F_aao[i])
    
    for i in 1:(n - 1) * (nx + nu)
        if in(i, x_id)
            @constraint(model, s[i] == h * (dot(A_aao[i,:], s) + D_aao[i]) + dot(E_aao[i,:], s) + F_aao[i])
        end
    end
    @constraint(model, u5[i=1:(n - 1)], s[(i - 1) * (nx + nu) + nu] == 1)

    @objective(model, Min, -dot(g, s) + 0.5 * s' * H * s)
    JuMP.optimize!(model)

    traj_sim, traj_opt = sol_2_trajectory(traj.r[1], value.(s), h)

    println("Optimized trajectory")
    println(traj_opt)

    println("simluated trajectory with non linear dynamics")
    println(traj_sim)

    t = 0:traj.dt:(traj.dt * n)
    n = traj.n

    px = scatter(t, [traj.r[i][xi] for i in 1:n], label="target")
    scatter!(t, [traj_sim.r[i][xi] for i in 1:n], label="simulated")
    xlabel!("time [s]")
    ylabel!("x [m]")
    savefig("myplot.png")
end
main()