module model

# using Revise

export buildLinearSystem

# trajectory functions
export make_linear_trajectory, trajectory, sol_2_trajectory

# constants
export m, d, g, Ix, Iy, Iz, Kt, nx, nu, ϕi, dϕi, θi, dθi, ψi, dψi, xi, dxi, yi, dyi, zi, dzi

include("trajectory.jl")
include("model.jl")

end

module all_at_once

# using Revise

# all at once functions
export all_at_once_RipQP, all_at_once_ipopt

include("AllAtOnce/AllAtOnce.jl")

end

module LQT

# using Revise

export LQT_D

include("DiscreteTime/LQT_D.jl")

end


module post_processing

# using Revise

export plot_trajectory, plot_u1, solve_L2_error, plot_disturbance

include("display.jl")

end

module problems

# using Revise

export generate_problem

include("problems.jl")

end

