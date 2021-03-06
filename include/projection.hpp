#pragma once
#include "Eigen/Dense"
#define EIGEN_RUNTIME_NO_MALLOC

using namespace Eigen;

inline void project_box(VectorXd &x, const VectorXd l, const VectorXd u)
{
    assert(l.rows() == u.rows());
    assert(x.rows() == l.rows());
    assert(x.rows() == u.rows());
    assert(l.cwiseMin(u) == l);
    x = x.cwiseMin(u).cwiseMax(l);
}

inline void project_ball(VectorXd &x, const double r)
{
    assert(r > 0);
    if (x.norm() > r)
    {
        x = (r / x.norm()) * x;
    }
}

inline void project_halfspace(VectorXd &x, VectorXd c, double a)
{
    assert(x.rows() == c.rows());
    assert(a != 0);
    if (c.dot(x) > a)
    {
        x = x - (c.dot(x) - a) * (c / (c.norm() * c.norm()));
    }
}