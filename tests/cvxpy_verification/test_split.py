from cvxpy import *
from numpy import array, inf
from scipy import sparse
import pprint
import time

T = 50  #time horizon
nx = 2
nu = 1

dt = 0.1

x0 = array([10, 0])
umax = 0.1

A = sparse.csc_matrix([
    [1, dt],
    [0, 1]
])
B = sparse.csc_matrix([
    [0.5*dt**2],
    [dt]
])
Q = sparse.eye(nx)
R = sparse.eye(nu)
u = Variable((nu, T))
x = Variable((nx, T+1))

times = []
for _ in range(10):
    objective = 0
    constraints = [x[:,0] == x0]
    for k in range(T):
        objective += quad_form(x[:,k], Q) + quad_form(u[:,k], R)
        constraints += [x[:,k+1] == A@x[:,k] + B@u[:,k]]
        constraints += [0.0 <= x[0,k], x[0,k] <= inf]
        constraints += [-0.01 <= x[1,k], x[1,k] <= 0.01]

    objective += quad_form(x[:,T], Q)
    prob = Problem(Minimize(objective), constraints)

    prob.solve(solver=OSQP, warm_start=False)
    times.append(int((prob.solver_stats.solve_time)*1e6))

print("Solve Time: " + str(int(sum(times)/len(times))) + " us")
# pprint.pprint(x[0,:].value)
# pprint.pprint(x[1,:].value)
# pprint.pprint(u[0,:].value)