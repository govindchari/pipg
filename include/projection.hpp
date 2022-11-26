#pragma once
#include "Eigen/Dense"
#define EIGEN_RUNTIME_NO_MALLOC

using namespace Eigen;

inline void project_box(VectorXd &x, const VectorXd l, const VectorXd u, size_t startIndex, size_t endIndex)
{
    auto idx = seq(startIndex, endIndex);
    x(idx) = x(idx).cwiseMin(u).cwiseMax(l);
}

inline void project_ball(VectorXd &x, const double r, size_t startIndex, size_t endIndex)
{
    auto idx = seq(startIndex, endIndex);
    if (x(idx).norm() > r)
    {
        x(idx) = (r / x(idx).norm()) * x(idx);
    }
}

inline void project_halfspace(VectorXd &x, const VectorXd c, const double a, size_t startIndex, size_t endIndex)
{
    auto idx = seq(startIndex, endIndex);
    if (c.dot(x(idx)) > a)
    {
        x(idx) = x(idx) - (c.dot(x(idx)) - a) * (c / (c.norm() * c.norm()));
    }
}