#pragma once
#include "Eigen/Dense"
#define EIGEN_RUNTIME_NO_MALLOC

using namespace Eigen;

struct TOLERANCE
{
    double pow_tol; // Relative error tolerance for power method iteration
    size_t pow_max_iter;
    TOLERANCE()
    {
        pow_tol = 1e-12;
        pow_max_iter = 100;
    };
};
struct WORKSPACE
{
    VectorXd vec_nx;
    VectorXd vec_nu;
    VectorXd vec_T1;
    VectorXd vec_T2;
    double d;
    WORKSPACE(){};
    WORKSPACE(size_t nx, size_t nu, size_t T)
    {
        vec_nx = VectorXd::Zero(nx);
        vec_nu = VectorXd::Zero(nu);
        vec_T1 = VectorXd::Zero(T);
        vec_T2 = VectorXd::Zero(T);
        d = 0.0;
    }
};