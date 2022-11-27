#pragma once
#include "catch.hpp"
#include "pipg.hpp"
#include <iostream>
#include <chrono>

using namespace std::chrono;

TEST_CASE("1D Point Mass MPC")
{
    std::cout << "" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "====================================================================" << std::endl;
    std::cout << "1D Point Mass MPC" << std::endl;
    std::cout << "====================================================================" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "" << std::endl;

    // Setup problem data
    size_t T = 50;
    size_t nx = 2;
    size_t nu = 1;
    double dt = 0.1;
    double umax = 0.1;
    VectorXd uvec(1);
    uvec << umax;
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
        p.addBoxConstraint(i, 'u', -uvec, uvec);
    }
    p.addQ(T, Q);
    p.addIC(x0);

    auto start = high_resolution_clock::now();
    p.solve(false);
    auto stop = high_resolution_clock::now();
    auto time = duration_cast<microseconds>(stop - start);
    auto X = p.getState();
    auto U = p.getControl();

    VectorXd X0_osqp(T + 1);
    VectorXd X1_osqp(T + 1);
    VectorXd U_osqp(T);
    X0_osqp << 10., 9.9995, 9.99799999, 9.99549999, 9.99199999,
        9.98749998, 9.98199998, 9.97549997, 9.96799996, 9.95949996,
        9.94999995, 9.93949994, 9.92799993, 9.91549992, 9.90199992,
        9.88749991, 9.8719999, 9.85549989, 9.83799987, 9.81949986,
        9.79999985, 9.77949984, 9.75799983, 9.73549982, 9.7119998,
        9.68749979, 9.66199978, 9.63549976, 9.60799975, 9.57949974,
        9.54999972, 9.51949971, 9.4879997, 9.45549968, 9.42199967,
        9.38749965, 9.35199964, 9.31549963, 9.27799961, 9.2394996,
        9.19999958, 9.15949957, 9.11799955, 9.07549954, 9.03199953,
        8.98749951, 8.9419995, 8.89549948, 8.84799947, 8.79961068,
        8.75085798;
    X1_osqp << -6.46887489e-09, -1.00000128e-02, -2.00000190e-02, -3.00000251e-02,
        -4.00000310e-02, -5.00000367e-02, -6.00000422e-02, -7.00000476e-02,
        -8.00000528e-02, -9.00000578e-02, -1.00000063e-01, -1.10000067e-01,
        -1.20000072e-01, -1.30000076e-01, -1.40000081e-01, -1.50000085e-01,
        -1.60000088e-01, -1.70000092e-01, -1.80000096e-01, -1.90000099e-01,
        -2.00000102e-01, -2.10000105e-01, -2.20000108e-01, -2.30000111e-01,
        -2.40000114e-01, -2.50000116e-01, -2.60000118e-01, -2.70000121e-01,
        -2.80000123e-01, -2.90000124e-01, -3.00000126e-01, -3.10000128e-01,
        -3.20000129e-01, -3.30000131e-01, -3.40000132e-01, -3.50000133e-01,
        -3.60000134e-01, -3.70000135e-01, -3.80000136e-01, -3.90000136e-01,
        -4.00000137e-01, -4.10000137e-01, -4.20000138e-01, -4.30000138e-01,
        -4.40000138e-01, -4.50000139e-01, -4.60000139e-01, -4.70000139e-01,
        -4.80000139e-01, -4.87775691e-01, -4.87278337e-01;
    U_osqp << -0.1, -0.1, -0.1, -0.1, -0.1,
        -0.1, -0.1, -0.1, -0.1, -0.1,
        -0.1, -0.1, -0.1, -0.1, -0.1,
        -0.1, -0.1, -0.1, -0.1, -0.1,
        -0.1, -0.1, -0.1, -0.1, -0.1,
        -0.1, -0.1, -0.1, -0.1, -0.1,
        -0.1, -0.1, -0.1, -0.1, -0.1,
        -0.1, -0.1, -0.1, -0.1, -0.1,
        -0.1, -0.1, -0.1, -0.1, -0.1,
        -0.1, -0.1, -0.1, -0.07775552, 0.00497354;

    // Speed Report
    int osqp_time = 1961;
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