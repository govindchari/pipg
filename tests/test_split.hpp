#pragma once
#include "catch.hpp"
#include "pipg.hpp"
#include <iostream>
#include <chrono>

using namespace std::chrono;

TEST_CASE("Split Test")
{
    std::cout << "" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "====================================================================" << std::endl;
    std::cout << "Split Test" << std::endl;
    std::cout << "====================================================================" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "" << std::endl;

    // Setup problem data
    size_t T = 50;
    size_t nx = 2;
    size_t nu = 1;
    double dt = 0.1;
    double umax = 0.1;

    VectorXd rl(1);
    VectorXd ru(1);
    VectorXd vl(1);
    VectorXd vu(1);

    rl << 0.0;
    ru << std::numeric_limits<double>::max();
    vl << -0.01;
    vu << 0.01;

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
    p.splitState('r', 0, 0);
    p.splitState('v', 1, 1);

    for (size_t i = 0; i < T; i++)
    {
        if (i > 0)
            p.addQ(i, Q);
        p.addR(i, R);
        p.addA(i, A);
        p.addB(i, B);
        p.addBoxConstraint(i, 'x', 'r', rl, ru);
        p.addBoxConstraint(i, 'x', 'v', vl, vu);
    }
    p.addQ(T, Q);
    p.addIC(x0);

    auto start = high_resolution_clock::now();
    p.solve(true);
    auto stop = high_resolution_clock::now();
    auto time = duration_cast<microseconds>(stop - start);
    auto X = p.getState();
    auto U = p.getControl();

    VectorXd X0_osqp(T + 1);
    VectorXd X1_osqp(T + 1);
    VectorXd U_osqp(T);
    X0_osqp << 10., 9.9995, 9.9985, 9.9975, 9.9965,
        9.9955, 9.9945, 9.9935, 9.9925, 9.9915,
        9.9905, 9.9895, 9.9885, 9.9875, 9.9865,
        9.9855, 9.9845, 9.9835, 9.9825, 9.9815,
        9.9805, 9.9795, 9.9785, 9.9775, 9.9765,
        9.9755, 9.9745, 9.9735, 9.9725, 9.9715,
        9.9705, 9.9695, 9.9685, 9.9675, 9.9665,
        9.9655, 9.9645, 9.9635, 9.9625, 9.9615,
        9.9605, 9.9595, 9.9585, 9.9575, 9.9565,
        9.9555, 9.9545, 9.9535, 9.9525, 9.9515,
        9.95025866;
    X1_osqp << -3.23357661e-13, -1.00000000e-02, -1.00000000e-02, -1.00000000e-02,
        -1.00000000e-02, -1.00000000e-02, -1.00000000e-02, -1.00000000e-02,
        -1.00000000e-02, -1.00000000e-02, -1.00000000e-02, -1.00000000e-02,
        -1.00000000e-02, -1.00000000e-02, -1.00000000e-02, -1.00000000e-02,
        -1.00000000e-02, -1.00000000e-02, -1.00000000e-02, -1.00000000e-02,
        -1.00000000e-02, -1.00000000e-02, -1.00000000e-02, -1.00000000e-02,
        -1.00000000e-02, -1.00000000e-02, -1.00000000e-02, -1.00000000e-02,
        -1.00000000e-02, -1.00000000e-02, -1.00000000e-02, -1.00000000e-02,
        -1.00000000e-02, -1.00000000e-02, -1.00000000e-02, -1.00000000e-02,
        -1.00000000e-02, -1.00000000e-02, -1.00000000e-02, -1.00000000e-02,
        -1.00000000e-02, -1.00000000e-02, -1.00000000e-02, -1.00000000e-02,
        -1.00000000e-02, -1.00000000e-02, -1.00000000e-02, -1.00000000e-02,
        -1.00000000e-02, -1.00000000e-02, -1.48268607e-02;
    U_osqp << -1.00000000e-01, -3.87883160e-12, -4.69486317e-12, -4.59917482e-12,
        -4.57065894e-12, -4.53612035e-12, -4.49719466e-12, -4.45407860e-12,
        -4.40681439e-12, -4.35544660e-12, -4.30002477e-12, -4.24060124e-12,
        -4.17723358e-12, -4.10998021e-12, -4.03890583e-12, -3.96407787e-12,
        -3.88556639e-12, -3.80344556e-12, -3.71779303e-12, -3.62868882e-12,
        -3.53621725e-12, -3.44046455e-12, -3.34152020e-12, -3.23947795e-12,
        -3.13443197e-12, -3.02648097e-12, -2.91572641e-12, -2.80226996e-12,
        -2.68621865e-12, -2.56767993e-12, -2.44676462e-12, -2.32358438e-12,
        -2.19825413e-12, -2.07089024e-12, -1.94161091e-12, -1.81053561e-12,
        -1.67778649e-12, -1.54348595e-12, -1.40775870e-12, -1.27073016e-12,
        -1.13252746e-12, -9.93278011e-13, -8.53110886e-13, -7.12155445e-13,
        -5.70542100e-13, -4.28326090e-13, -2.87527224e-13, -1.24318297e-13,
        -1.60585141e-13, -4.82686072e-02;

    // Speed Report
    int osqp_time = 135263;
    std::cout << "===========================Speed Results=============================" << std::endl;
    std::cout << "OSQP is " << time.count() / osqp_time << " times faster" << std::endl;
    std::cout << "PIPG Run Time (us):" << std::endl;
    std::cout << time.count() << std::endl;
    std::cout << "OSQP Run Time (us):" << std::endl;
    std::cout << osqp_time << std::endl;

    VectorXd X0_pipg(T + 1);
    VectorXd X1_pipg(T + 1);
    VectorXd U_pipg(T);

    for (int i = 0; i < T + 1; i++)
    {
        X0_pipg[i] = X[i][0];
        X1_pipg[i] = X[i][1];
        if (i != T)
            U_pipg[i] = U[i][0];
    }

    double tol = 1e-3;
    REQUIRE((X0_pipg - X0_osqp).lpNorm<Eigen::Infinity>() <= tol);
    REQUIRE((X1_pipg - X1_osqp).lpNorm<Eigen::Infinity>() <= tol);
    REQUIRE((U_pipg - U_osqp).lpNorm<Eigen::Infinity>() <= tol);
}