#pragma once

#include "Eigen/Dense"
#define EIGEN_RUNTIME_NO_MALLOC

using namespace Eigen;

namespace Constraint
{
    struct Box
    {
        size_t _startIndex;
        size_t _endIndex;
        VectorXd _l;
        VectorXd _u;
        char _var;
        Box(size_t startIndex, size_t endIndex, const VectorXd l, const VectorXd u, const char var) : _startIndex(startIndex), _endIndex(endIndex), _l(l), _u(u), _var(var) {}
    };

    struct Ball
    {
        size_t _startIndex;
        size_t _endIndex;
        double _r;
        char _var;
        Ball(size_t startIndex, size_t endIndex, const double r, const char var) : _startIndex(startIndex), _endIndex(endIndex), _r(r), _var(var) {}
    };

    struct Halfspace
    {
        size_t _startIndex;
        size_t _endIndex;
        VectorXd _c;
        double _a;
        char _var;
        Halfspace(size_t startIndex, size_t endIndex, const VectorXd c, const double a, const char var) : _startIndex(startIndex), _endIndex(endIndex), _c(c), _a(a), _var(var) {}
    };

}
