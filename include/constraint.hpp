#pragma once
#include <vector>

#include "Eigen/Dense"
#define EIGEN_RUNTIME_NO_MALLOC

using namespace Eigen;

namespace Constraint
{
    struct Box
    {
        VectorXd _l;
        VectorXd _u;
        char _var;
        Box(const VectorXd l, const VectorXd u, const char var) : _l(l), _u(u), _var(var) {}
    };

    struct Ball
    {
        double _r;
        char _var;
        Ball(const double r, const char var) : _r(r), _var(var) {}
    };

    struct Halfspace
    {
        VectorXd _c;
        double _a;
        char _var;
        Halfspace(const VectorXd c, const double a, const char var) : _c(c), _a(a), _var(var) {}
    };

}
