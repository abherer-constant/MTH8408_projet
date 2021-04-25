using ..model
using Distributions

function generate_disturbance(n::Int)
    n -= 1 # remove disturbance on initial condition
    Amplitude = 5 # N
    n_disturb_xy = Int(ceil(0.02 / 30 * n))
    # strong impulsion at 1/6 of the time
    dx_start = Int(round(n / 6))
    dx = zeros(n)
    dx[dx_start:dx_start + n_disturb_xy] .= Amplitude

    # strong impulsion at 1/3 of the time
    dy_start = Int(round(n / 3))
    dy = zeros(n)
    dy[dy_start:dy_start + n_disturb_xy] .= Amplitude

    # gaussian perturbation from 3/6 to 5/6
    dz_start = Int(round(n * 1 / 2))
    dz_end = Int(round(n * 5 / 6))
    dz = zeros(n)

    d = Normal(0.0, 1.0)
    dz[dz_start:dz_end] = rand(d, dz_end - dz_start + 1)

    return hcat(dx, dy, dz)
end

"""
Generate a problem from one of the templates in the reference paper

traj_ID: A or B, the trajectory identifier from the reference paper
n: approximate number of discretisation points
disturbed: add the disturbance described in the reference paper to the trajectory
time_multiplier: multiplier for the total trajectory time. Having a time_multiplier of 2 reduces the velocity by a factor of 2
"""
function generate_problem(traj_ID::String; n::Int=1000, disturbed::Bool=false, time_multiplier::Real=1.)

    if disturbed
        dist = generate_disturbance(n)
    else
        dist = zeros(n, 3)
    end
    
    if traj_ID == "A"
        traj = generate_A(n, time_multiplier, dist)
    elseif traj_ID == "B"
        traj = generate_B(n, time_multiplier, dist)
    end



    return traj
end

"""
Generate the trajectory A from the reference paper
"""
function generate_A(n::Int, time_multiplier::Real, dist)
    pts = [[0.      ,0.     ,0.],
            [15.     ,0.     ,0.],
            [30.     ,0.     ,15.],
            [30.     ,5.     ,15.],
            [15.     ,5.     ,0.],
            [0.      ,5.     ,0.]]
    t = [0.,15.,30.,30.,45.,60.] * time_multiplier
    return make_linear_trajectory(pts, t, n, disturbance=dist)
end

"""
Generate the trajectory B from the reference paper
"""
function generate_B(n::Int, time_multiplier::Real, dist)
    pts = [[0., 0., 1.],[80., 80., 1.]]
    t = [0., 60.] * time_multiplier
    period = Int(round(20 / 60 * n))
    
    B_traj = make_linear_trajectory(pts, t, n, disturbance=dist)
    for i in 1:n
        B_traj.r[i][zi] = 1 + sin(i / period * 2 * π)
    end
    return B_traj
end