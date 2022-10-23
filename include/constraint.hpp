#pragma once
#include <vector>

#include "Eigen/Dense"
#define EIGEN_RUNTIME_NO_MALLOC

using namespace Eigen;

namespace Constraint
{
    struct Box
    {
        VectorXd l;
        VectorXd u;
        size_t t;
        char var;
        Box(const VectorXd lin, const VectorXd uin, const size_t tin, const char varin) : l(lin), u(uin), t(tin), var(varin) {}
    };

    struct Ball
    {
        double r;
        size_t t;
        char var;
        Ball(const double rin, const size_t tin, const char varin) : r(rin), t(tin), var(varin) {}
    };

    struct Halfspace
    {
        VectorXd c;
        double a;
        size_t t;
        char var;
        Halfspace(const VectorXd cin, const double ain, const size_t tin, const char varin) : c(cin), a(ain), t(tin), var(varin) {}
    };

}
