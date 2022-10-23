#pragma once
#include <vector>

#include "Eigen/Dense"
#include <Eigen/Eigenvalues>
#define EIGEN_RUNTIME_NO_MALLOC

#include "constraint.hpp"
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
    std::vector<VectorXd> Q;
    std::vector<VectorXd> R;

    MatrixXd H;

    // Constraints
    std::vector<Constraint::Box> box_constraints;
    std::vector<Constraint::Ball> ball_constraints;
    std::vector<Constraint::Halfspace> halfspace_constraints;

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
    TOLERANCE tol;

    void updateEta1();
    void updateEta2();
    void updateEta3();

    void projectAll();

public:
    void printQR();
    void addA(const size_t t, const MatrixXd Ain);
    void addB(const size_t t, const MatrixXd Bin);
    void addQ(const size_t t, const VectorXd Qin);
    void addR(const size_t t, const VectorXd Rin);
    void addIC(const VectorXd x0);
    void addBallConstraint(const size_t t, const char variable, const double r);
    void addBoxConstraint(const size_t t, const char variable, const VectorXd l, const VectorXd u);
    void addHalfspaceConstraint(const size_t t, const char variable, const VectorXd c, const double a);
    void solve();
    void solve(bool verbose);
    std::vector<VectorXd> getState();
    std::vector<VectorXd> getControl();

    MPC(const size_t T_horizon, const size_t nxin, const size_t nuin);
};