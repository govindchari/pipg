#pragma once
#include <vector>
#include <unordered_map>

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
    size_t _T;
    size_t _nx;
    size_t _nu;

    // Problem Data
    std::vector<MatrixXd> _A;
    std::vector<MatrixXd> _B;
    std::vector<VectorXd> _Q;
    std::vector<VectorXd> _R;

    // Equality constraint matrix
    MatrixXd _H;

    // Constraints
    std::vector<std::vector<Constraint::Box>> _box_constraints;
    std::vector<std::vector<Constraint::Ball>> _ball_constraints;
    std::vector<std::vector<Constraint::Halfspace>> _halfspace_constraints;

    // State and control slices
    std::unordered_map<char, std::pair<size_t, size_t>> _state_slice;
    std::unordered_map<char, std::pair<size_t, size_t>> _control_slice;

    // Optimization Variables
    std::vector<VectorXd> _X;
    std::vector<VectorXd> _U;

    // Integral Term and Dual Variable
    std::vector<VectorXd> _V;
    std::vector<VectorXd> _W;

    // Previous Primal and Dual iterates
    std::vector<VectorXd> _Xprev;
    std::vector<VectorXd> _Uprev;
    std::vector<VectorXd> _Wprev;

    // Eta Parameters
    double _eta1;
    double _eta2;
    double _eta3;
    bool _eta1_outdated;
    bool _eta2_outdated;
    bool _eta3_outdated;

    // Workspace
    WORKSPACE ws;
    TOLERANCE _tol;

    void updateEta1();
    void updateEta2();
    void updateEta3();
    bool nonOverlappingSlice(const char variable, const size_t startIndex, const size_t endIndex);
    bool validSplit();
    void project(const size_t t);

public:
    void printQR();
    void addA(const size_t t, const MatrixXd Ain);
    void addB(const size_t t, const MatrixXd Bin);
    void addQ(const size_t t, const VectorXd Qin);
    void addR(const size_t t, const VectorXd Rin);
    void addIC(const VectorXd x0);
    void splitState(const char label, const size_t startIndex, const size_t endIndex);
    void splitControl(const char label, const size_t startIndex, const size_t endIndex);
    void addBallConstraint(const size_t t, const char variable, const char label, const double r);
    void addBallConstraint(const size_t t, const char variable, const double r);
    void addBoxConstraint(const size_t t, const char variable, const char label, const VectorXd l, const VectorXd u);
    void addBoxConstraint(const size_t t, const char variable, const VectorXd l, const VectorXd u);
    void addHalfspaceConstraint(const size_t t, const char variable, const char label, const VectorXd c, const double a);
    void addHalfspaceConstraint(const size_t t, const char variable, const VectorXd c, const double a);
    void solve();
    void solve(bool verbose);
    std::vector<VectorXd> getState();
    std::vector<VectorXd> getControl();

    MPC(const size_t T_horizon, const size_t nxin, const size_t nuin);
};