using LinearAlgebra
# Unoptimized implementation of PIPG for MPC with LTI dynamics and static state and input cost
# Can support box and ball constraints on state and input

struct MPC
    #Problem parameters
    A::Array{Float64,2}
    B::Array{Float64,2}
    Q::Array{Float64,2}
    R::Array{Float64,2}

    #Optimization variables
    x::Array{Float64,2}
    u::Array{Float64,2}

    #Proportional and Integral Terms
    v::Array{Float64,2} #Proportional plus Integral term
    w::Array{Float64,2} #Integral Term

    T::Int64

    #Input constraint
    umax::Float64
    function MPC(A,B,Q,R,umax,x0)
        T = size(A)[0]
        new(A,B,Q,R,T,umax)
    end
end
function proj_box(x::Array{Float64,1},l::Array{Float64,1},u::Array{Float64,1})
    #assert(l less than u)
    return min.(max.(x,l),u)
end
function proj_ball(x::Array{Float64,1},ρ::Int64)
    return (ρ/norm(x))*x
end
function solve!(p::MPC)
    for k=1:kmax
        α = 2/((k+1)*p.η_1+2*p.η_2)
        β = (k+1)*p.η_1/(2*p.η_3)
        for t=1:p.T
            p.v[t,:] = p.w[t,:] + β * (p.x[t,:] - p.A * p.x[t-1,:] - p.B * p.u[t-1,:])
            p.u[t-1,:] = proj_ball(p.u[t-1,:] - α * (p.R * u[t-1,:] - B' * p.v[t,:]), umax)
            p.x[t,:] = p.x[t,:] - α * (p.Q*(p.x[t,:]) + p.v[t,:] - p.A' * v[t+1,:])
            p.w[t,:] = p.w[t,:] + β * (p.x[t,:] - p.A * p.x[t-t,:] - p.B * u[t-1,:])
        end
    end
end
let
    T = 10   #time horizon
    N = 2    #number of states

    # Need T-1 A and B matrices for LTV
    # Make Q_packed and R_packed vertically concatanated Q_t and R_t
    Q = I(N)
    R = I(N)
end