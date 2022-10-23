#include <iostream>
#include "Eigen/Dense"
#include "pipg.hpp"

int main()
{
    size_t T = 50;
    size_t nx = 2;
    size_t nu = 1;
    double dt = 0.1;
    double umax = 0.1;
    MatrixXd A(nx, nx);
    MatrixXd B(nx, 1);
    VectorXd Q(nx);
    VectorXd R(nu);
    VectorXd x0(nx);

    A << 1, dt, 0, 1;
    B << 0.5 * dt * dt, dt;
    Q << 1, 1;
    R << 1;
    x0 << 10, 0;

    MPC p(T, nx, nu);

    for (size_t i = 0; i < T; i++)
    {
        if (i > 0)
            p.addQ(i, Q);
        p.addR(i, R);
        p.addA(i, A);
        p.addB(i, B);
        p.addBallConstraint(i, 'x', umax);
    }
    p.addQ(T, Q);
    p.addIC(x0);

    p.solve();

    auto X = p.getState();

    std::cout << "[";
    for (auto x : X)
    {
        std::cout << x[0] << ", ";
    }
    std::cout << "]";
}