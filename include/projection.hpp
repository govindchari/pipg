#pragma once
#include "Eigen/Dense"
#define EIGEN_RUNTIME_NO_MALLOC

using namespace Eigen;

inline void project_box(VectorXd &x, const VectorXd l, const VectorXd u)
{
    x = x.cwiseMin(u).cwiseMax(l);
}

inline void project_ball(VectorXd &x, const double r)
{
    if (x.norm() > r)
    {
        x = (r / x.norm()) * x;
    }
}

inline void project_halfspace(VectorXd &x, const VectorXd c, const double a)
{
    if (c.dot(x) > a)
    {
        x = x - (c.dot(x) - a) * (c / (c.norm() * c.norm()));
    }
}