#pragma once
#include "Eigen/Dense"
#define EIGEN_RUNTIME_NO_MALLOC

using namespace Eigen;

struct TOLERANCE
{
    TOLERANCE();
};
struct WORKSPACE
{
    VectorXd vec_nx;
    VectorXd vec_nu;
    double d;
    WORKSPACE(){};
    WORKSPACE(size_t nx, size_t nu)
    {
        vec_nx = VectorXd::Zero(nx);
        vec_nu = VectorXd::Zero(nu);
        d = 0.0;
    }
};