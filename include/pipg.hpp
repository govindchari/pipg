#pragma once
#include <vector>

#include "Eigen/Dense"
#define EIGEN_RUNTIME_NO_MALLOC

#include "projection.hpp"
#include "utilities.hpp"

using namespace Eigen;
class MPC
{
private:
    // Time Horizon
    size_t T;
    size_t nx;
    size_t nu;

    // Problem Data
    std::vector<MatrixXd> A;
    std::vector<MatrixXd> B;
    std::vector<MatrixXd> Q;
    std::vector<MatrixXd> R;

    // Optimization Variables
    std::vector<VectorXd> X;
    std::vector<VectorXd> U;

    // Proportional and Integral Terms
    std::vector<VectorXd> V;
    std::vector<VectorXd> W;

    // Eta Parameters
    double eta1;
    double eta2;
    double eta3;
    bool eta1_outdated;
    bool eta2_outdated;
    bool eta3_outdated;

    // Workspace
    WORKSPACE ws;

    // void updateEta1();
    void updateEta2();
    void updateEta3();

public:
    void updateEta1();
    void printQR();
    void addA(const size_t t, const MatrixXd Ain);
    void addB(const size_t t, const MatrixXd Bin);
    void addQ(const size_t t, const MatrixXd Qin);
    void addR(const size_t t, const MatrixXd Rin);
    // void addBallConstraint(const size_t t, const enum variable, const double r);
    // void addBoxConstraint(const size_t t, const enum variable, const VectorXd l, const VectorXd u);
    // void addHalfspaceConstraint(const size_t t, const enum variable, const VectorXd c, const double a);
    void solve();
    void solve(bool verbose);

    MPC(const size_t T_horizon, const size_t nxin, const size_t nuin);
};