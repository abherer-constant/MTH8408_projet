using Plots
using LinearAlgebra
using Printf

include("LQT.jl")

n = 100

# points de la trajectoire et le moment correspondant
pts = [[0.      ,0.     ,0.],
       [15.     ,0.     ,0.],
       [30.     ,0.     ,15.],
       [30.     ,5.     ,15.],
       [15.     ,5.     ,0.],
       [0.      ,5.     ,0.]]
t = [0.,15.,30.,30.,45.,60.]

traj_paper = make_linear_trajectory(pts, t, n)

sol = OptimalLQT(traj_paper)


t = range(0, stop=  traj_paper.dt*(length(traj_paper.r)-1), step = traj_paper.dt)*I

plot(t,sol[1][7,:])
