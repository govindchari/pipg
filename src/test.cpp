#include <iostream>
#include "Eigen/Dense"
#include "pipg.hpp"

int main()
{
    int T = 10;
    size_t nx = 2;
    size_t nu = 2;
    MatrixXd Q0(nx, nx);
    MatrixXd R0(nu, nu);
    Q0 << 1, 2, 3, 4;
    R0 << 5, 6, 7, 8;

    MPC p(T, nx, nu);

    for (size_t i = 0; i < T; i++)
    {
        p.addQ(i, Q0);
        p.addR(i, R0);
    }
    p.updateEta1();
}