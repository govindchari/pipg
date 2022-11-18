#pragma once
#include "Eigen/Dense"
#define EIGEN_RUNTIME_NO_MALLOC

using namespace Eigen;

struct TOLERANCE
{
    double pow_tol; // Relative error tolerance for power method iteration
    size_t pow_max_iter;
    size_t max_iter;
    size_t k;
    double eq_tol;
    bool stop;
    TOLERANCE()
    {
        // Angle between two iterations of power method must be within acos(1-pow_tol) degrees
        pow_tol = 1e-10;
        pow_max_iter = 100;
        max_iter = 20000;
        eq_tol = 1e-6;
        stop = false;
        k = 1;
    };
};
struct WORKSPACE
{
    VectorXd vec_nx;
    VectorXd vec_nu;
    VectorXd vec_T1;
    VectorXd vec_T2;
    VectorXd vec_T3;
    VectorXd vec_T4;
    WORKSPACE(){};
    WORKSPACE(size_t nx, size_t nu, size_t T)
    {
        vec_nx = VectorXd::Zero(nx);
        vec_nu = VectorXd::Zero(nu);
        vec_T1 = VectorXd::Zero(T);
        vec_T2 = VectorXd::Zero(T);
        vec_T3 = VectorXd::Zero(T);
        vec_T4 = VectorXd::Zero(T);
    }
};