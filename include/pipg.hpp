#pragma once
#include <vector>

#include "Eigen/Dense"
#define EIGEN_RUNTIME_NO_MALLOC

#include "projection.hpp"

using namespace Eigen;
class MPC{
    private:
        // Time Horizon
        int T;

        // Problem Data
        std::vector<MatrixXd> A;
        std::vector<MatrixXd> B;
        std::vector<MatrixXd> C;
        std::vector<MatrixXd> D;

        // Optimization Variables
        std::vector<MatrixXd> X;
        std::vector<MatrixXd> U;

        // Proportional and Integral Terms
        std::vector<MatrixXd> V;
        std::vector<MatrixXd> W;

        // Eta Parameters
        double eta1;
        double eta2;
        double eta3;

    public:
        void addA(const size_t t, const MatrixXd A);
        void addB(const size_t t, const MatrixXd B);
        void addQ(const size_t t, const MatrixXd Q);
        void addR(const size_t t, const MatrixXd R);
        void addBallConstraint(const size_t t, const enum variable, const double r);
        void addBoxConstraint(const size_t t, const enum variable, const VectorXd l, const VectorXd u);
        void addHalfspaceConstraint(const size_t t, const enum variable, const VectorXd c, const double a);
        void solve();
        void solve(bool verbose);

        MPC(const size_t T);

};