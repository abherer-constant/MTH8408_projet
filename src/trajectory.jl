using Printf,
    SolverCore
using ..model

"""
Hold information about a drone trajectory

parameters
dt : timestep value
r  : array of array containing the system state a each timestep
n  : number of points
"""
struct trajectory
    dt
    r
    n
end

function Base.show(io::IO, traj::trajectory)
    @printf "dt= %f\n" traj.dt
    @printf "n=  %i\n" traj.n
    @printf "%6s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n" "i" "ϕ" "dϕ" "θ" "dθ" "ψ" "dψ" "x" "dx" "y" "dy" "z" "dz"
    for i in 1:traj.n
        ϕ, dϕ, θ, dθ, ψ, dψ, x, dx, y, dy, z, dz = traj.r[i]
        @printf "%6s  %10f %10f %10f %10f %10f %10f %10f %10f %10f %10f %10f %10f\n" i ϕ dϕ θ dθ ψ dψ x dx y dy z dz
    end
end

"""
Create a trajectory structure from an array of points the drone must pass through

parameters
points  : array of array representing the coordinates the trajectory must go through.
        exemple: [[x1, y1, z1], [x2, y2, z2]], the trajectory will go straight from point 1
        to point 2.
time    : array of timestamps corresponding to each point in points. This will define the 
        velocity of a drone between every 2 points.
n       : number of system states in the final discretized trajectory. Note: the final number
        of system states might be slighty different than the intputed value.
"""
function make_linear_trajectory(points::Array{Array{Float64,1},1}, 
                                time::Array{Float64,1}, 
                                n::Integer)
    for point in points
        @assert(length(point) == 3)
    end
    @assert(length(time) == length(points))
    @assert(time == sort(time))

    @assert(n ≥ length(time))

    # define time allocated for each segment
    time_total = last(time)
    dt = time_total / (n - 1)
    segment_n = Array{Integer,1}(undef, length(time) - 1)
    for i in 2:length(time)
        segment_time = time[i] - time[i - 1]
        segment_n[i - 1] = round(segment_time / time_total * (n - 1))
    end

    n = sum(segment_n) + 1
    r = Array{Array{Float64,1},1}(undef, n)
    ri = 1
    for i in 2:length(points)
        p0 = points[i - 1]
        p1 = points[i]
        n_local = segment_n[i - 1]
        t_local = 0
        f(t) = p0 + (p1 - p0) / (n_local * dt) * t
        for j in 1:n_local
            r[ri] = zeros(nx)
            r[ri][[xi,yi,zi]] = f(t_local)
            t_local += dt
            ri += 1
        end
    end
    r[n] = zeros(nx)
    r[n][[xi,yi,zi]] = last(points)

    # # add gravity
    # for state in r
    #     state[dzi] += g
    # end

    return trajectory(dt, r, n)
end

function sol_2_trajectory(sol::GenericExecutionStats, target_trajectory::trajectory)
    n = target_trajectory.n
    dt = target_trajectory.dt
    r1 = target_trajectory.r[1]

    x_sol = sol.solution[1:nx  * (n - 1)]
    u_sol = sol.solution[nx * (n - 1) + 1:end]

    # 
    r_opt = Array{Array{Float64,1},1}(undef, n) # linear solution
    r_opt[1] = r1
    for i in 2:(n)
        r_opt[i] = x_sol[(i - 2) * nx + 1:(i - 1) * nx ]
    end
    traj_opt = trajectory(dt, r_opt, n)

    r_sim = Array{Array{Float64,1},1}(undef, n) # non-linear solution
    r_sim[1] = r1
    for i in 2:(n)
        r_sim[i] = NonLinearOutput(u_sol[(i - 2) * nu + 1:(i - 1) * nu ], r_sim[i - 1], dt)
    end
    traj_sim = trajectory(dt, r_sim, n)

    return traj_opt, traj_sim
end
