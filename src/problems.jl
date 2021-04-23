using ..model
using Distributions

drone_trajectories = Dict{String,trajectory}()

# nombre de points de discrétisation
n = 1000

# points de la trajectoire et le moment correspondant
pts = [[0.      ,0.     ,0.],
    [15.     ,0.     ,0.],
    [30.     ,0.     ,15.],
    [30.     ,5.     ,15.],
    [15.     ,5.     ,0.],
    [0.      ,5.     ,0.]]
t = [0.,15.,30.,30.,45.,60.]

drone_trajectories["A_fast"] = make_linear_trajectory(pts, t, n)
drone_trajectories["A_slow"] = make_linear_trajectory(pts, t * 10, n)
drone_trajectories["A_sample"] = make_linear_trajectory(pts, t, 50)


# trajectoire B de l'article de référence (sinus)
n = 1000
pts = [[0., 0., 1.],[80., 80., 1.]]
t = [0., 60.]
period = Int(round(20 / 60 * n))

B_traj = make_linear_trajectory(pts, t, n)
for i in 1:n
    B_traj.r[i][zi] = 1 + sin(i / period * 2 * π)
end
drone_trajectories["B_fast"] = B_traj
t *= 10
B_traj = make_linear_trajectory(pts, t, n)
for i in 1:n
    B_traj.r[i][zi] = 1 + sin(i / period * 2 * π)
end
drone_trajectories["B_slow"] = B_traj

# définition de la matrice perturbation

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