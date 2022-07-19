#include "Eigen/Dense"

using namespace Eigen;

void project_box(VectorXd &x, const VectorXd l, const VectorXd u){
    // assert(l < u);
    x = x.cwiseMin(u).cwiseMax(l); 
}

void project_ball(VectorXd &x, const double r){
    assert(r > 0);
    if (x.norm() > r){
        x = (r/x.norm())*x;
    }
}

void project_halfspace(VectorXd x, VectorXd c, double a){
    // static_assert(dim(x) == dim(c));
}