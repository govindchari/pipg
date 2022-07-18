using LinearAlgebra
using SparseArrays
using IterativeSolvers
using Printf
using PyPlot

# Unoptimized implementation of PIPG for MPC with LTI dynamics and static state and input cost
# Can support box and ball constraints on state and input

struct MPC
    #Problem parameters
    A::Array{Float64,2}
    B::Array{Float64,2}
    Q::Array{Float64,2}
    R::Array{Float64,2}

    #Time horizon
    T::Int64

    #Optimization variables
    X::Array{Float64,2}
    U::Array{Float64,2}

    #Proportional and Integral Terms
    V::Array{Float64,2} #Proportional plus Integral term
    W::Array{Float64,2} #Integral Term

    #Variable sizes
    nx::Int64
    nu::Int64

    #Eta parameters
    η_1::Float64
    η_2::Float64
    η_3::Float64

    #Input constraint
    umax::Float64

    #Add in support for xT
    function MPC(A,B,Q,R,T,x0,umax)
        nx = size(Q)[1]
        nu = size(R)[1]

        P = sparse(R)
        for i=1:T-1
            P = blockdiag(P,sparse(Q),sparse(R))
        end
        P = blockdiag(P,sparse(Q))

        H = kron(I(T), [-B I(nx)])
        for k = 1:T-1
            H[(k*nx).+(1:nx), (k*(nx+nu)-nx).+(1:nx)] .= -A
        end
        HtH = sparse(H'*H)

        η_1 = real(invpowm(P)[1])
        η_2 = real(powm(P)[1])
        η_3 = real(powm(HtH)[1])

        X = zeros(T+1,nx)
        U = zeros(T,nu)
        V = zeros(T+1,nx)
        W = zeros(T+1,nx)

        X[1,:] = x0

        new(A,B,Q,R,T,X,U,V,W,nx,nu,η_1,η_2,η_3,umax)
    end
end
function proj_box(x::Array{Float64,1},l::Array{Float64,1},u::Array{Float64,1})
    #assert(l less than u)
    return min.(max.(x,l),u)
end
function proj_ball(x::Array{Float64,1},ρ::Float64)
    #assert ρ is positive
    if norm(x) > ρ
        return (ρ/norm(x))*x
    else
        return x
    end
end
function solve!(p::MPC)
    kmax = 10000
    e = zeros(p.T+1)
    println("iter     objv       |Gx-g|\n")
    println("-----------------------------\n")

    for k=1:kmax
        α = 2/((k+1)*p.η_1+2*p.η_2)
        β = (k+1)*p.η_1/(2*p.η_3)
        obj = 0
        for t=2:p.T
            p.V[t,:] = p.W[t,:] + β * (p.X[t,:] - p.A * p.X[t-1,:] - p.B * p.U[t-1,:])
            p.U[t-1,:] = proj_ball(p.U[t-1,:] - α * (p.R * p.U[t-1,:] - p.B' * p.V[t,:]), p.umax)
            p.X[t,:] = proj_box(p.X[t,:] - α * (p.Q*(p.X[t,:]) + p.V[t,:] - p.A' * p.V[t+1,:]),[0.0;-Inf],[Inf,Inf])
            p.W[t,:] = p.W[t,:] + β * (p.X[t,:] - p.A * p.X[t-1,:] - p.B * p.U[t-1,:])
            e[t] = norm(p.X[t,:] - p.A * p.X[t-1,:] - p.B * p.U[t-1,:])^2
            obj += p.X[t,:]'*p.Q*p.X[t,:] + p.U[t-1,:]'*p.R*p.U[t-1,:]
        end
        if (mod(k,1)==0)
            @printf("%3d   %10.3e  %9.2e\n",
            k, obj, sqrt(sum(e)))
        end
    end
end
let
    T = 100  #time horizon
    N = 2    #number of states
    dt = 0.1
    x0 = [10;0]
    umax = 0.1

    A = [1 dt; 0 1]
    B = [0.5*dt^2;dt]
    B = reshape(B,length(B),1)
    Q = Diagonal([1.0,1.0])
    R = Diagonal([10.0])

    # Need T-1 A and B matrices for LTV
    # Make Q_packed and R_packed vertically concatanated Q_t and R_t

    mpc = MPC(A,B,Q,R,T,x0,umax)

    solve!(mpc)
    pygui(true)

    plt.figure()
    plt.plot(mpc.X[:,1])
    plt.title("State Trajectory")
    plt.grid(true)

    plt.figure()
    plt.plot(mpc.U)
    plt.title("Input Trajectory")
    plt.grid(true)
    plt.show()

end