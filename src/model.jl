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
nu = 4

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

"""
Get the non-linearised LEB matrix from a set of euler angles Φ, θ, ψ
"""
function getLEB(euler)
    Φ, θ, ψ = euler
    LEB = [cos(θ) * cos(ψ)      sin(Φ) * sin(θ) * cos(ψ) - cos(Φ) * sin(ψ)        cos(Φ) * sin(θ) * cos(ψ) + sin(Φ) * sin(ψ)
           cos(θ) * sin(ψ)      sin(Φ) * sin(θ) * sin(ψ) + cos(Φ) * cos(ψ)        cos(Φ) * sin(θ) * sin(ψ) - sin(Φ) * cos(ψ) 
           -sin(θ)              sin(Φ) * cos(θ)                                   cos(ψ) * cos(θ)             ]
    return LEB
end

"""
Solve the system step at the next timestep

parameters
U: array of length 4 containing the commands for the last time step
x0: array of length 12 containing the state for the last time step
h: time step length
"""
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

"""
create A, B anc C matrix from the system's dynamic
"""
function buildLinearSystem()
    LEB = getLEB([0,0,0])

    n = nx      # number of state variables
    r = nu      # number of inputs -> U1, U2, U3, U4
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
    # B[dzi,nu] = -g
    B[[dϕi,dθi,dψi], [2,3,4]] = Diagonal([d / Ix, d / Iy, Iz])
    C = 1 * Matrix(I, 12, 12)
    return A, B, C
end


