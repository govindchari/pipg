#include <iostream>

#include "projection.hpp"
#include "Eigen/Dense"

int main(){
    double r = 1;
    VectorXd x(2);
    VectorXd l(2);
    VectorXd u(2);
   
    l << 0, 0;
    u << 1, 1;
    x << -1, 2;

    project_box(x, l, u);
    std::cout << x << std::endl;
}