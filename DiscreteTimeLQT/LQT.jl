using Plots
using LinearAlgebra, SolverTools
using Printf

include("model.jl")

function OptimalLQT(traj::trajectory)

    h = traj.dt
    r_ini = traj.r
    kf = length(r_ini)

    r = zeros(12,kf)
    for i in 1:kf
        r[:, i] = r_ini[i]'
    end

    

    A, B, C = buildLinearSystem()

    Ad = I + A * h 

    Bd = B * h 

    Cd = C 

    ### Calcul des paramètres de la méthode LQT

    # Matrice de poids d'état
    Q = Diagonal([100, 50, 10, 5, 0, 0, 100, 1, 100, 1, 1000, 0.1]);

    # Matrice de poids de contrôle
    R = I

    # Calcul des conditions finales et des solutions au temps final (kf) :
    F = zeros(12,12) ;

    Pk = Cd' * F * Cd 
    gk = Cd' * F * r[1:12, kf] 

    U_star = zeros(4, size(r,2)) 
    X_star = zeros(size(r)) 
    X_star[1:12, 1] = r[1:12, 1] 

    # Paramètres invariants dans le temps :
    V = Cd' * Q * Cd 

    E = Bd * inv(R) * Bd' 

    L = [];
    Lg = [] ;
    P = [] ;
    g = [] ;
    push!(P, Pk)
    push!(g, gk)

    k = kf

    while k !== 1
        k -= 1 ;
    
        Lk = (R + Bd' * Pk * Bd) \ Bd' * Pk * Ad 
        Lgk = (R + Bd' * Pk * Bd) \ Bd' 
        push!(L, Lk)
        push!(Lg, Lgk)
    
        gk = (Ad' - Ad' * Pk / (I + E * Pk) * E) * gk + Cd' * Q * r[1:12, k] ;
        Pk = Ad' * Pk / (I + E * Pk) * Ad + V;
        push!(P, Pk)
        push!(g, gk)
    end

    L = L[kf-1:-1:1]
    Lg = Lg[kf-1:-1:1]

    g = g[kf:-1:1]
    P = P[kf:-1:1]

    for k in 1:kf - 1
        X_star[1:12 , k+1] = (Ad - Bd * L[k]) * X_star[1:12, k] + Bd * Lg[k] * g[k+1] ; # - gravity * h);
        U_star[1:4,  k] = - L[k] * X_star[1:12, k] + Lg[k] * g[k+1] ;
    end

    status = :first_order
    model = Model(Ipopt.Optimizer)
    sol = [[X_star];[U_star]]

    return GenericExecutionStats(status, model, solution = sol)
end