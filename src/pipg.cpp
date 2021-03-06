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

    WORKSPACE ws(nx, nu, T);
    TOLERANCE tol();
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
    auto Xc = X;
    auto Uc = U;
    std::fill(Xc.begin(), Xc.end(), VectorXd::Random(nx).normalized());
    std::fill(Uc.begin(), Uc.end(), VectorXd::Random(nu).normalized());

    // Stores max eigenvalue estimates of Q_t and R_t
    ws.vec_T1 = VectorXd::Zero(T);
    ws.vec_T2 = VectorXd::Zero(T);

    size_t iter = 1;
    bool converged = false;
    double error = std::numeric_limits<double>::max();
    while (!converged && iter < tol.pow_max_iter)
    {
        error = std::numeric_limits<double>::max();
        for (size_t t = 0; t < T; t++)
        {
            ws.vec_nx = Q[t] * Xc[t + 1]; // Corresponds to Q_t * X_t for t > 0 (no Q_0)
            ws.vec_T1[t] = ws.vec_nx.norm() / Xc[t + 1].norm();
            error = std::min(error, Xc[t + 1].dot(ws.vec_nx.normalized()));
            Xc[t + 1] = ws.vec_nx.normalized();
            ws.vec_nu = R[t] * Uc[t];
            ws.vec_T2[t] = ws.vec_nu.norm() / Uc[t].norm();
            error = std::min(error, U[t].dot(ws.vec_nu.normalized()));
            Uc[t] = ws.vec_nu.normalized();
        }
        converged = error > 1 - tol.pow_tol;
        iter++;
    }
    eta1 = std::max(ws.vec_T1.maxCoeff(), ws.vec_T2.maxCoeff());
    eta1_outdated = false;
}
void MPC::updateEta2()
{
    // Uses power method shifted by lambda_max, so needs eta_1
    if (eta1_outdated)
    {
        updateEta1();
        updateEta2();
    }
    else
    {
        auto Xc = X;
        auto Uc = U;
        std::fill(Xc.begin(), Xc.end(), VectorXd::Random(nx).normalized());
        std::fill(Uc.begin(), Uc.end(), VectorXd::Random(nu).normalized());

        // Stores max eigenvalue estimates of Q_t-lambda_max*I and R_t-lambda_max
        ws.vec_T3 = VectorXd::Zero(T);
        ws.vec_T4 = VectorXd::Zero(T);

        size_t iter = 1;
        bool converged = false;
        double error = std::numeric_limits<double>::max();
        auto Inx = MatrixXd::Identity(nx, nx);
        auto Inu = MatrixXd::Identity(nu, nu);
        while (!converged && iter < tol.pow_max_iter)
        {
            error = std::numeric_limits<double>::max();
            for (size_t t = 0; t < T; t++)
            {
                ws.vec_nx = (Q[t] - ws.vec_T1[t] * Inx) * Xc[t + 1]; // Corresponds to Q_t * X_t for t > 0 (no Q_0)
                ws.vec_T3[t] = ws.vec_nx.norm() / Xc[t + 1].norm();
                error = std::min(error, Xc[t + 1].dot(ws.vec_nx.normalized()));
                Xc[t + 1] = ws.vec_nx.normalized();
                ws.vec_nu = (R[t] - ws.vec_T2[t] * Inu) * Uc[t];
                ws.vec_T4[t] = ws.vec_nu.norm() / Uc[t].norm();
                error = std::min(error, U[t].dot(ws.vec_nu.normalized()));
                Uc[t] = ws.vec_nu.normalized();
            }
            converged = error > 1 - tol.pow_tol;
            iter++;
        }
        ws.vec_T3 = -ws.vec_T3 + ws.vec_T1;
        ws.vec_T4 = -ws.vec_T4 + ws.vec_T2;
        eta2 = std::min(ws.vec_T3.minCoeff(), ws.vec_T4.minCoeff());
        eta2_outdated = false;
    }
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
void MPC::addIC(const VectorXd x0)
{
    X[0] = x0;
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
    for (size_t k = 0; k < 10000; k++)
    {
        auto a = 2 / ((k + 1) * eta1 + 2 * eta2);
        auto b = (k + 1) * eta1 / (2 * eta3);
        for (size_t t = 1; t < T + 1; t++)
        {
            V[t] = W[t] + b * (X[t] - A[t - 1] * X[t - 1] - B[t - 1] * U[t - 1]);
            U[t - 1] = U[t - 1] - a * (R[t - 1] * U[t - 1] - B[t - 1].transpose() * V[t]);
            if (t != T + 1)
            {
                X[t] = X[t] - a * (Q[t] * X[t] + V[t] - A[t].transpose() * V[t + 1]);
            }
            W[t] = W[t] + b * (X[t] - A[t - 1] * X[t - 1] - B[t - 1] * U[t - 1]);
        }
    }
}