using ..model

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

drone_trajectories["no_disturbance_fast"] = make_linear_trajectory(pts, t, n)
drone_trajectories["no_disturbance_slow"] = make_linear_trajectory(pts, t * 10, n)
drone_trajectories["no_disturbance_sample"] = make_linear_trajectory(pts, t, 50)

