using Plots
using ..model
using SolverCore
using LinearAlgebra

"""
Plot the position in function of time of all trajectories contained in the trajectories dictionnary
"""
function plot_trajectory(trajectories::Dict{String,trajectory})
    
    px = plot([], [], ylabel="x [m]", label="")
    py = plot([], [], ylabel="y [m]", legend=false, label="")
    pz = plot([], [], ylabel="z [m]", xlabel="temps [s]", legend=false, label="")

    for (name, traj) in trajectories
        n = traj.n
        dt = traj.dt
        t = range(0, dt * n; length=n)
        
        x_sol = [traj.r[i][xi] for i in 1:n]
        y_sol = [traj.r[i][yi] for i in 1:n]
        z_sol = [traj.r[i][zi] for i in 1:n]
        
        plot!(px, t, x_sol, label=name)
        plot!(py, t, y_sol, label=name)
        plot!(pz, t, z_sol, label=name)
    end
    
    l = @layout [a ; b; c]
    plt = plot(px, py, pz, layout=l)
    plot!(plt, size=(640, 480))
    
    return plt
end

"""
Plot the u₁ command in function of time of all solutions contained in the stats_dict dictionnary
"""
function plot_u1(stats_dict::Dict{String,GenericExecutionStats})
    
    plt = plot([], [], ylabel="u₁ [N]", xlabel="temps [s]", label="")

    for (name, stats) in stats_dict
        dt = stats.solver_specific[:dt]
        n = stats.solver_specific[:n]
        t = range(0, dt * (n - 1); length=n - 1)
        
        u_sol = stats.solution[nx * (n - 1) + 1:nu:end]
        
        plot!(plt, t, u_sol, label=name)

    end
    
    return plt
end


"""
plot error on position as a function of time and calculate the L2 norm of the error for
    each trajectory in the trajectories dictionnary
"""
function solve_L2_error(trajectories::Dict{String,trajectory}, target_traj::trajectory)

    px = plot([], [], ylabel="erreur sur x [m]", label="")
    py = plot([], [], ylabel="erreur sur y [m]", legend=false, label="")
    pz = plot([], [], ylabel="erreur sur z [m]", xlabel="temps [s]", legend=false, label="")

    n = target_traj.n
    x_target = [target_traj.r[i][xi] for i in 1:n]
    y_target = [target_traj.r[i][yi] for i in 1:n]
    z_target = [target_traj.r[i][zi] for i in 1:n]

    errors = Dict{String,Real}()

    for (name, traj) in trajectories
        n = traj.n
        dt = traj.dt
        t = range(0, dt * n; length=n)
        
        x_sol = [traj.r[i][xi] for i in 1:n]
        y_sol = [traj.r[i][yi] for i in 1:n]
        z_sol = [traj.r[i][zi] for i in 1:n]
        
        x_error = x_sol - x_target
        y_error = y_sol - y_target
        z_error = z_sol - z_target

        plot!(px, t, x_error, label=name)
        plot!(py, t, y_error, label=name)
        plot!(pz, t, z_error, label=name)

        errors[name] = norm(vcat(x_error, y_error, z_error))
    end
    
    l = @layout [a ; b; c]
    plt = plot(px, py, pz, layout=l)
    plot!(plt, size=(640, 480))
    
    return plt, errors
end


function plot_disturbance(disturbance, dt)
    
    disturbance = vcat([0 0 0], disturbance) # add initial condition
    n = size(disturbance, 1)
    t = range(0, dt * n; length=n)
    
    px = plot(t, vec(disturbance[:,1]), ylabel="perturbation\nx [N]", label="")
    py = plot(t, vec(disturbance[:,2]), ylabel="perturbation\ny [N]", legend=false, label="")
    pz = plot(t, vec(disturbance[:,3]), ylabel="perturbation\nz [N]", xlabel="temps [s]", legend=false, label="")

    l = @layout [a ; b; c]
    plt = plot(px, py, pz, layout=l)
    plot!(plt, size=(640, 480))
    
    return plt
end