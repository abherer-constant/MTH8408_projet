using LinearAlgebra
using Printf

################################################################################
#                           C O N S T A N T S
################################################################################

m = 0.485      # [kg]       
d = 0.35       # [m]       lever length
g = 9.81       # [m s-2]
Ix = 3.4e-3    # [kg m2]
Iy = 3.4e-3    # [kg m2] 
Iz = 4.7e-3    # [kg m2]

Kt = Diagonal([0.1, 0.1, 0.1])

# number of variables in a single time step
nx = 12
nu = 5

# state variable indexes
ϕi      = 1
dϕi     = 2
θi      = 3
dθi     = 4     
ψi      = 5
dψi     = 6
xi      = 7
dxi     = 8
yi      = 9
dyi     = 10
zi      = 11
dzi     = 12

function getLEB(euler)
    Φ, θ, ψ = euler
    LEB = [cos(θ) * cos(ψ)      sin(Φ) * sin(θ) * cos(ψ) - cos(Φ) * sin(ψ)        cos(Φ) * sin(θ) * cos(ψ) + sin(Φ) * sin(ψ)
           cos(θ) * sin(ψ)      sin(Φ) * sin(θ) * sin(ψ) + cos(Φ) * cos(ψ)        cos(Φ) * sin(θ) * sin(ψ) - sin(Φ) * cos(ψ) 
           -sin(θ)              sin(Φ) * cos(θ)                                   cos(ψ) * cos(θ)             ]
    return LEB
end

function NonLinearOutput(U, x0, h)
    LEB = getLEB(x0[[ϕi, θi, ψi]])
    ddx, ddy, ddz = -[0; 0; g] + LEB * [0; 0; U[1] / m] - Kt / m * x0[[dxi, dyi, dzi]]
    ddϕ = ([(Iy - Iz) * x0[dθi] * x0[dψi] / Ix] + [U[2] * d / Ix])[1]
    ddθ = ([(Iz - Ix) * x0[dϕi] * x0[dψi] / Iy] + [U[3] * d / Iy])[1]
    ddψ = ([(Ix - Iy) * x0[dϕi] * x0[dθi] / Iz] + [U[4]  / Iz])[1]
    dx, dy, dz, dϕ, dθ, dψ = x0[[dxi, dyi, dzi, dϕi, dθi, dψi]] .+ h * [ddx, ddy, ddz, ddϕ, ddθ, ddψ]
    x, y, z, ϕ, θ, ψ = x0[[xi, yi, zi, ϕi, θi, ψi]] .+ h * [dx, dy, dz, dϕ, dθ, dψ]
    x1 = [ϕ; dϕ; θ; dθ; ψ; dψ; x; dx; y; dy; z; dz]
    return x1
end

function buildLinearSystem()
    LEB = getLEB([0,0,0])

    n = nx      # number of state variables
    r = nu      # number of inputs -> U1, U2, U3, U4, U5=1
    o = 12      # number of outputs -> x, y, z

    A = zeros(n, n)
    B = zeros(n, r)
    C = zeros(o, n)

    # construction of the matrix

    A[[dxi,dyi,dzi],[dxi,dyi,dzi]] = -Kt / m
    A[xi,dxi] = 1
    A[yi,dyi] = 1
    A[zi,dzi] = 1 
    A[ψi,dψi] = 1
    A[θi,dθi] = 1
    A[ϕi,dϕi] = 1 
    # A[dxi, θi] = 1 / m
    # A[dyi, ϕi] = -1 / m 
    A[dxi, θi] = g
    A[dyi, ϕi] = -g 

    B[[dxi,dyi,dzi],1] = (LEB / m * [0,0,1]) 
    B[dzi,5] = -g
    B[[dϕi,dθi,dψi], [2,3,4]] = Diagonal([d / Ix, d / Iy, Iz])

    C = 1 * Matrix(I, 12, 12)
    return A, B, C
end


struct trajectory
    dt
    r
    n
end

function Base.show(io::IO, traj::trajectory)
    @printf "dt= %f\n" traj.dt
    @printf "n=  %i\n" traj.n
    @printf "%10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n" "ϕ" "dϕ" "θ" "dθ" "ψ" "dψ" "x" "dx" "y" "dy" "z" "dz"
    for i in 1:traj.n
        ϕ, dϕ, θ, dθ, ψ, dψ, x, dx, y, dy, z, dz = traj.r[i]
        @printf "%10f %10f %10f %10f %10f %10f %10f %10f %10f %10f %10f %10f\n" ϕ dϕ θ dθ ψ dψ x dx y dy z dz
    end
end


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
    time_total = sum(time)
    dt = time_total / (n - 1)
    segment_n = Array{Integer,1}(undef, length(time) - 1)
    for i in 2:length(time)
        segment_time = time[i] - time[i - 1]
        segment_n[i - 1] = round(segment_time / time_total) * (n - 1)
    end

    @assert(n - 1 == sum(segment_n))

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

    return trajectory(dt, r, n)
end

function sol_2_trajectory(x1, sol, dt)
    n = Int(length(sol) / (nx + nu))
    r_sim = Array{Array{Float64,1},1}(undef, n + 1)
    r_opt = Array{Array{Float64,1},1}(undef, n + 1)
    r_sim[1] = x1
    r_opt[1] = x1
    for i in 2:(n + 1)
        u = sol[(nx + nu) * (i - 2) + 1:(nx + nu) * (i - 2) + nu]
        r_sim[i] = NonLinearOutput(u, r_sim[i - 1], dt)
        r_opt[i] = sol[(nx + nu) * (i - 2) + nu + 1:(nx + nu) * (i - 1)]
    end
    return trajectory(dt, r_sim, n + 1), trajectory(dt, r_opt, n + 1)
end
