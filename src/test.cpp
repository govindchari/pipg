#include <iostream>
#include "Eigen/Dense"
#include "pipg.hpp"

int main()
{
    size_t T = 2;
    size_t nx = 2;
    size_t nu = 2;
    MatrixXd Q1(nx, nx);
    MatrixXd Q2(nx, nx);
    MatrixXd R0(nu, nu);
    MatrixXd R1(nu, nu);
    Q1 << 2, 0, 0, 0;
    Q2 << 4, 0, 0, 1;
    R0 << 0.5, 0, 0, 8;
    R1 << 6, 0, 0, 0.5;

    MPC p(T, nx, nu);

    p.addQ(0, Q1);
    p.addQ(1, Q2);
    p.addR(0, R0);
    p.addR(1, R1);

    // for (size_t i = 0; i < T; i++)
    // {
    //     p.addQ(i, Q0);
    //     p.addR(i, R0);
    // }
    p.updateEta2();
}