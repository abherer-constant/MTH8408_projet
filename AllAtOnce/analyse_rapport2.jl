using Pkg
Pkg.activate(".")
Pkg.instantiate()

include("AllAtOnce.jl")

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

traj_paper = make_linear_trajectory(pts, t, n)

all_at_once_ipopt(traj_paper,savePath="ipopt_paper_comp")
all_at_once_RipQP(traj_paper,savePath="ripqp_paper_comp")

# on ralenti la vitesse de déplacement du drone par un facteur de 4
t = t * 4

traj_paper_slow = make_linear_trajectory(pts, t, n)

all_at_once_ipopt(traj_paper_slow,savePath="ipopt_paper_comp_slow")
all_at_once_RipQP(traj_paper_slow,savePath="ripqp_paper_comp_slow")