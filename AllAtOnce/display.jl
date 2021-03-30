using Plots

function plot_trajectory(t_last, traj_sim, traj_planned; savePath="")
    t = range(0, t_last; length=traj_planned.n)
    n = traj_planned.n
    l = @layout [a ; b; c]


    px = plot(t, [traj_planned.r[i][xi] for i in 1:n], label="cible", ylabel="x [m]")
    plot!(t, [traj_sim.r[i][xi] for i in 1:n], label="simulation")

    py = plot(t, [traj_planned.r[i][yi] for i in 1:n], legend=false, ylabel="y [m]")
    plot!(t, [traj_sim.r[i][yi] for i in 1:n])

    pz = plot(t, [traj_planned.r[i][zi] for i in 1:n], legend=false, ylabel="z [m]", xlabel="time [s]")
    plot!(t, [traj_sim.r[i][zi] for i in 1:n])

    
    plt = plot(px, py, pz, layout=l)
    plot!(plt, size=(600, 800))

    if (savePath != "")
        savefig(savePath * ".png")
    else
        display(plt)
    end
end