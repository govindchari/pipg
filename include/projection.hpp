#include "Eigen/Dense"

using namespace Eigen;

void project_box(VectorXd &x, const VectorXd l, const VectorXd u){
    assert(l.rows() == u.rows());
    assert(x.rows() == l.rows());
    assert(x.rows() == u.rows());
    assert(l.cwiseMin(u) == l);
    x = x.cwiseMin(u).cwiseMax(l); 
}

void project_ball(VectorXd &x, const double r){
    assert(r > 0);
    if (x.norm() > r){
        x = (r/x.norm())*x;
    }
}

void project_halfspace(VectorXd &x, VectorXd c, double a){
    assert(x.rows() == c.rows());
    assert(a != 0);
    if (c.dot(x) > a){
        x = x - (c.dot(x) - a) * (c / (c.norm()*c.norm()));
    }
}