#include "pipg.hpp"
#include <iostream>

MPC::MPC(const size_t T_horizon, const size_t nxin, const size_t nuin)
{
    T = T_horizon;
    nx = nxin;
    nu = nuin;
    A.resize(T);
    B.resize(T);
    Q.resize(T);
    R.resize(T);
    X.resize(T + 1);
    U.resize(T);
    V.resize(T + 2);
    W.resize(T + 1);

    std::fill(A.begin(), A.end(), MatrixXd::Zero(nx, nx));
    std::fill(B.begin(), B.end(), MatrixXd::Zero(nx, nu));
    std::fill(Q.begin(), Q.end(), MatrixXd::Zero(nx, nx));
    std::fill(R.begin(), R.end(), MatrixXd::Zero(nu, nu));
    std::fill(X.begin(), X.end(), VectorXd::Zero(nx));
    std::fill(U.begin(), U.end(), VectorXd::Zero(nu));
    std::fill(V.begin(), V.end(), VectorXd::Zero(nx));
    std::fill(W.begin(), W.end(), VectorXd::Zero(nx));

    eta1 = 0.0;
    eta2 = 0.0;
    eta3 = 0.0;
    eta1_outdated = false;
    eta2_outdated = false;
    eta3_outdated = false;

    WORKSPACE ws(nx, nu);
}
void MPC::printQR()
{
    std::cout << "Q: " << std::endl;
    for (auto x : Q)
    {
        std::cout << x << std::endl;
    }
    std::cout << "R: " << std::endl;
    for (auto x : R)
    {
        std::cout << x << std::endl;
    }
}
void MPC::updateEta1()
{
    std::fill(X.begin(), X.end(), VectorXd::Random(nx));
    std::fill(U.begin(), U.end(), VectorXd::Random(nu));
    // Add convergence criteria
    for (size_t i = 0; i < 10; i++)
    {
        for (size_t t = 0; t < T; t++)
        {
            ws.vec_nx = Q[t] * X[t + 1];
            ws.d = ws.vec_nx.norm() / X[t + 1].norm();
            eta1 = std::max(ws.d, eta1);
            X[t + 1] = ws.vec_nx / ws.vec_nx.norm();
            ws.vec_nu = R[t] * U[t];
            ws.d = ws.vec_nu.norm() / U[t].norm();
            eta1 = std::max(ws.d, eta1);
            U[t] = ws.vec_nu / ws.vec_nu.norm();
        }
    }
    eta1_outdated = false;
}
void MPC::updateEta2()
{
    eta2_outdated = false;
}
void MPC::updateEta3()
{
    eta3_outdated = false;
}
void MPC::addA(const size_t t, const MatrixXd Ain)
{
    A[t] = Ain;
    eta3_outdated = true;
}
void MPC::addB(const size_t t, const MatrixXd Bin)
{
    B[t] = Bin;
    eta3_outdated = true;
}
void MPC::addQ(const size_t t, const MatrixXd Qin)
{
    Q[t] = Qin;
    eta1_outdated = true;
    eta2_outdated = true;
}
void MPC::addR(const size_t t, const MatrixXd Rin)
{
    R[t] = Rin;
    eta1_outdated = true;
    eta2_outdated = true;
}
// void MPC::addBallConstraint(const size_t t, const enum variable, const double r){

// }
// void MPC::addBoxConstraint(const size_t t, const enum variable, const VectorXd l, const VectorXd u){

// }
// void MPC::addHalfspaceConstraint(const size_t t, const enum variable, const VectorXd c, const double a){

// }
void MPC::solve()
{
    solve(true);
}
void MPC::solve(bool verbose)
{
    if (eta1_outdated)
    {
        updateEta1();
    }
    if (eta2_outdated)
    {
        updateEta2();
    }
    if (eta3_outdated)
    {
        updateEta3();
    }
    for (size_t i = 0; i < 10000; i++)
    {
    }
}