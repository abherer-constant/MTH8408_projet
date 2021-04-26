using LinearAlgebra
using JuMP
using Ipopt
using MathOptInterface
using SparseArrays
using Plots
using RipQP
using QuadraticModels
using SolverCore
using NLPModelsJuMP
using NLPModelsIpopt
using ..model

"""
Create all the matrices for the all at once formulation
"""
function build_allAtOnce_system(traj::trajectory)

    h = traj.dt
    r::Array{Array{Float64,1},1} = traj.r
    n = length(r)
    n_tot = (n - 1) * (nx + nu) # x1 vector is the initial condition

    A, B, C = buildLinearSystem()
    Q = Diagonal([100,50,10,5,0,0,100,1,100,1,1000,0.1])
    # Q = Diagonal([0,0,0,0,0,0,100,0,100,0,1000,0])
    R = Diagonal([10,0,0,0])
    # R = Diagonal([0,0,0,0])

    g_vec = zeros(n_tot)
    H = spzeros(n_tot, n_tot)
    A_aao = spzeros(n_tot, n_tot)
    B_aao = spzeros(n_tot, n_tot)
    C_aao = spzeros(nx * n, nx * n)
    D_aao = spzeros(n_tot)
    E_aao = spzeros(n_tot, n_tot)
    F_aao = spzeros(n_tot)
    G_aao = spzeros(n_tot)

    ####### Objective ########

    # populate g vector
    for i in 2:n
        subrange    = nu + (i - 2) * (nx + nu) + 1:nu + (i - 2) * (nx + nu) + nx
        g_vec[subrange] = (r[i]' * Q' * C' + r[i]' * Q * C) / 2
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

    # populate G_aao vector, used to include gravity and disturbance
    disturbance = traj.disturbance
    for i in 1:(n - 1)
        G_aao[(i) * (nu + nx)] = -g + disturbance[i,3] / m
        G_aao[(i) * (nu + nx) - 2] = disturbance[i,2] / m
        G_aao[(i) * (nu + nx) - 4] = disturbance[i,1] / m
    end

    return g_vec, H, A_aao, D_aao, E_aao, F_aao, G_aao
end  


"""
Transform a solution array from [u1 x2 u2 x2 ...] to [x2 x3 x4 ... u1 u2 u3]
"""
function gen_sol(sol)
    n = length(sol)
    x = []
    u = []
    for i in 1:(nx + nu):n
        x = vcat(x, sol[i + nu:i + nu + nx - 1])
        u = vcat(u, sol[i:i + nu - 1])
    end
    return vcat(x, u)
end

function all_at_once_ipopt(traj) 
    n = traj.n
    h = traj.dt

    g_vec, H, A_aao, D_aao, E_aao, F_aao, G_aao = build_allAtOnce_system(traj)
    A, B, C = buildLinearSystem()

    model = Model(Ipopt.Optimizer)
    @variable(model, s[1:(n - 1) * (nx + nu)])

    x_id = []
    for i in 1:(n - 1)
        id = collect(1:nx) .+ nu .+ (i - 1) * (nx + nu)
        append!(x_id, id)
    end
  
    for i in 1:(n - 1) * (nx + nu)
        @constraint(model, s[i] == h * (dot(A_aao[i,:], s) + D_aao[i] + G_aao[i]) + dot(E_aao[i,:], s) + F_aao[i])
    end

    @objective(model, Min, -dot(g_vec, s) + 0.5 * s' * H * s)
        
    nlp = MathOptNLPModel(model)
    stats = ipopt(nlp)
    stats.solution = gen_sol(stats.solution)
    stats.solver_specific[:dt] = traj.dt
    stats.solver_specific[:n] = traj.n 
    return stats
end

function all_at_once_RipQP(traj) 
    n = traj.n
    h = traj.dt
    g_vec, H, A_aao, D_aao, E_aao, F_aao, G_aao = build_allAtOnce_system(traj)
    A, B, C = buildLinearSystem()

    A_qm = I - h * A_aao - E_aao   
    L_qm = h * D_aao + h * G_aao + F_aao 

    QM = QuadraticModel(-g_vec, H, A=A_qm, lcon=L_qm, ucon=L_qm, name="QM")
    
    stats = ripqp(QM)
    stats.solution = gen_sol(stats.solution)
    stats.solver_specific[:dt] = traj.dt
    stats.solver_specific[:n] = traj.n 
    return stats
end