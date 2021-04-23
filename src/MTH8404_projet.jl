module model

using Revise

export buildLinearSystem

# trajectory functions
export make_linear_trajectory, trajectory, sol_2_trajectory

# constants
export m, d, g, Ix, Iy, Iz, Kt, nx, nu, ϕi, dϕi, θi, dθi, ψi, dψi, xi, dxi, yi, dyi, zi, dzi

includet("trajectory.jl")
includet("model.jl")

end

module all_at_once

using Revise

# all at once functions
export all_at_once_RipQP, all_at_once_ipopt

includet("AllAtOnce/AllAtOnce.jl")

end


module post_processing

using Revise

export plot_trajectory, plot_u1, solve_L2_error

includet("display.jl")

end

module problems

using Revise

export drone_trajectories

includet("problems.jl")

end

