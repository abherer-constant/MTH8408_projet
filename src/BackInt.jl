using LinearAlgebra
include("linearSys.jl")

struct ControlGains
    n::int32
    L::Array{Array{float64,2},1}
    Lg::Array{Array{float64,2},1}
    g::Array{Array{float64,2},1}
end

function BackInt(n, r::Array{Array{Float64,1},1})
    A, B, C = buildLinearSystem()
    Q = Diagonal([100,50,10,5,0,0,100,1,100,1,1000,0.1])
    R = Diagonal([10,0,0,0])
    I = 1*Matrix(I,12,12)
    V = C'*Q*C
    E = B*inv(R)*B'

    # control gains
    L  = Array{Array{Float64,2},1}(undef, n)#(zeros(4,12),n)
    Lg = Array{Array{Float64,2},1}(undef, n)#zeros(4,12,n)
    g  = Array{Array{Float64,2},1}(undef, n)#zeros(12,1,n)

    # initial conditions
    P = C'*F*C
    g[n] = C'*F*r[n]
    # equation set (36)
    for k in (n-1):1
        g[k] = (A'- A'*P*inv(I+E*P)*E)*g[k+1]+C'*Q*r[k]
        Xstar
        P[k] = A'*P[k]*inv(I+E*P)*A+V
    
    end

end