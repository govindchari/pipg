#include "pipg.hpp"
#include <iostream>

MPC::MPC(const size_t T_horizon, const size_t nxin, const size_t nuin)
{
    T = T_horizon;
    nx = nxin;
    nu = nuin;
    A.resize(T + 1); // Of length T+1 because A_T = 0 (dummy variable for alg to work)
    B.resize(T);
    Q.resize(T); // Of length T+1 since Q_0 is not of meaning (waste of space tbh)
    R.resize(T);
    X.resize(T + 1);
    U.resize(T);
    V.resize(T + 2); // Of length T+2 because V_{T+1} = 0 (dummy variable for alg to work)
    W.resize(T + 1);

    std::fill(A.begin(), A.end(), MatrixXd::Zero(nx, nx));
    std::fill(B.begin(), B.end(), MatrixXd::Zero(nx, nu));
    std::fill(Q.begin(), Q.end(), VectorXd::Zero(nx));
    std::fill(R.begin(), R.end(), VectorXd::Zero(nu));
    std::fill(X.begin(), X.end(), VectorXd::Zero(nx));
    std::fill(U.begin(), U.end(), VectorXd::Zero(nu));
    std::fill(V.begin(), V.end(), VectorXd::Zero(nx));
    std::fill(W.begin(), W.end(), VectorXd::Zero(nx));

    H = MatrixXd::Zero(T * nx, T * (nx + nu));

    for (size_t i = 0; i < T; i++)
    {
        H.block(i * nx, nu + i * (nx + nu), nx, nx) = MatrixXd::Identity(nx, nx);
    }

    eta1 = 0.0;
    eta2 = 0.0;
    eta3 = 0.0;
    eta1_outdated = false;
    eta2_outdated = false;
    eta3_outdated = false;

    WORKSPACE ws(nx, nu, T);
    TOLERANCE tol;
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
    std::cout << eta1 << std::endl;
    std::cout << eta2 << std::endl;
}
void MPC::updateEta1()
{
    // auto Xc = X;
    // auto Uc = U;
    // std::fill(Xc.begin(), Xc.end(), VectorXd::Random(nx).normalized());
    // std::fill(Uc.begin(), Uc.end(), VectorXd::Random(nu).normalized());

    // // Stores max eigenvalue estimates of Q_t and R_t
    // ws.vec_T1 = VectorXd::Zero(T);
    // ws.vec_T2 = VectorXd::Zero(T);

    // size_t iter = 1;
    // bool converged = false;
    // double error = std::numeric_limits<double>::max();
    // while (!converged && iter < tol.pow_max_iter)
    // {
    //     error = std::numeric_limits<double>::max();
    //     for (size_t t = 0; t < T; t++)
    //     {
    //         ws.vec_nx = Q[t] * Xc[t + 1]; // Corresponds to Q_t * X_t for t > 0 (no Q_0)
    //         ws.vec_T1[t] = ws.vec_nx.norm() / Xc[t + 1].norm();
    //         error = std::min(error, Xc[t + 1].dot(ws.vec_nx.normalized()));
    //         Xc[t + 1] = ws.vec_nx.normalized();
    //         ws.vec_nu = R[t] * Uc[t];
    //         ws.vec_T2[t] = ws.vec_nu.norm() / Uc[t].norm();
    //         error = std::min(error, U[t].dot(ws.vec_nu.normalized()));
    //         Uc[t] = ws.vec_nu.normalized();
    //     }
    //     converged = error > 1 - tol.pow_tol;
    //     iter++;
    // }
    // eta1 = std::max(ws.vec_T1.maxCoeff(), ws.vec_T2.maxCoeff());

    // Find min value in Q/R
    double minElement = std::numeric_limits<double>::max();
    for (auto x : Q)
    {
        if (x.minCoeff() < minElement)
        {
            minElement = x.minCoeff();
        }
    }
    for (auto x : R)
    {
        if (x.minCoeff() < minElement)
        {
            minElement = x.minCoeff();
        }
    }
    eta1 = minElement;
    eta1_outdated = false;
}
void MPC::updateEta2()
{
    // Uses power method shifted by lambda_max, so needs eta_1
    // if (eta1_outdated)
    // {
    //     updateEta1();
    //     updateEta2();
    // }
    // else
    // {
    //     auto Xc = X;
    //     auto Uc = U;
    //     std::fill(Xc.begin(), Xc.end(), VectorXd::Random(nx).normalized());
    //     std::fill(Uc.begin(), Uc.end(), VectorXd::Random(nu).normalized());

    //     // Stores max eigenvalue estimates of Q_t-lambda_max*I and R_t-lambda_max
    //     ws.vec_T3 = VectorXd::Zero(T);
    //     ws.vec_T4 = VectorXd::Zero(T);

    //     size_t iter = 1;
    //     bool converged = false;
    //     double error = std::numeric_limits<double>::max();
    //     auto Inx = MatrixXd::Identity(nx, nx);
    //     auto Inu = MatrixXd::Identity(nu, nu);
    //     while (!converged && iter < tol.pow_max_iter)
    //     {
    //         error = std::numeric_limits<double>::max();
    //         for (size_t t = 0; t < T; t++)
    //         {
    //             ws.vec_nx = (Q[t] - ws.vec_T1[t] * Inx) * Xc[t + 1]; // Corresponds to Q_t * X_t for t > 0 (no Q_0)
    //             ws.vec_T3[t] = ws.vec_nx.norm() / Xc[t + 1].norm();
    //             error = std::min(error, Xc[t + 1].dot(ws.vec_nx.normalized()));
    //             Xc[t + 1] = ws.vec_nx.normalized();
    //             ws.vec_nu = (R[t] - ws.vec_T2[t] * Inu) * Uc[t];
    //             ws.vec_T4[t] = ws.vec_nu.norm() / Uc[t].norm();
    //             error = std::min(error, U[t].dot(ws.vec_nu.normalized()));
    //             Uc[t] = ws.vec_nu.normalized();
    //         }
    //         converged = error > 1 - tol.pow_tol;
    //         iter++;
    //     }
    //     ws.vec_T3 = -ws.vec_T3 + ws.vec_T1;
    //     ws.vec_T4 = -ws.vec_T4 + ws.vec_T2;
    //     eta2 = std::min(ws.vec_T3.minCoeff(), ws.vec_T4.minCoeff());
    //     eta2_outdated = false;
    // }
    // Find max element of Q/R
    double maxElement = std::numeric_limits<double>::min();
    for (auto x : Q)
    {
        if (x.maxCoeff() > maxElement)
        {
            maxElement = x.maxCoeff();
        }
    }
    for (auto x : R)
    {
        if (x.maxCoeff() > maxElement)
        {
            maxElement = x.maxCoeff();
        }
    }
    eta2 = maxElement;
    eta2_outdated = false;
}
void MPC::updateEta3()
{
    // Compute HtH
    // Use power method to get max eig
    auto HtH = (H.transpose() * H);
    auto eig = HtH.eigenvalues();
    std::cout << eig[0].real() << std::endl;
    eta3 = eig[0].real();
    eta3_outdated = false;
}
void MPC::addA(const size_t t, const MatrixXd Ain)
{
    assert(t >= 0 && t < T);
    A[t] = Ain;
    if (t != 0)
    {
        H.block(t * nx, nu + (t - 1) * (nu + nx), nx, nx) = -Ain;
    }
    eta3_outdated = true;
}
void MPC::addB(const size_t t, const MatrixXd Bin)
{
    assert(t >= 0 && t < T);
    B[t] = Bin;
    H.block(t * nx, t * (nu + nx), nx, nu) = -Bin;
    eta3_outdated = true;
}
void MPC::addQ(const size_t t, const VectorXd Qin)
{
    assert(t > 0 && t <= T);
    Q[t - 1] = Qin;
    eta1_outdated = true;
    eta2_outdated = true;
}
void MPC::addR(const size_t t, const VectorXd Rin)
{
    assert(t >= 0 && t < T);
    R[t] = Rin;
    eta1_outdated = true;
    eta2_outdated = true;
}
void MPC::addIC(const VectorXd x0)
{
    assert(MPC::nx == std::make_unsigned_t<int>(x0.rows()));
    X[0] = x0;
}
void MPC::addBoxConstraint(const size_t t, const unsigned char variable, const VectorXd l, const VectorXd u)
{
    assert(l.rows() == u.rows());
    assert(MPC::nx == std::make_unsigned_t<int>(l.rows()));
    assert(MPC::nx == std::make_unsigned_t<int>(u.rows()));
    assert(l.cwiseMin(u) == l);

    assert(variable == 'x' || variable == 'u');

    Constraint::Box con(l, u, variable, t);
    MPC::box_constraints.push_back(con);
}
void MPC::addBallConstraint(const size_t t, const unsigned char variable, const double r)
{
    assert(r > 0);

    assert(variable == 'x' || variable == 'u');

    Constraint::Ball con(r, variable, t);
    MPC::ball_constraints.push_back(con);
}
void MPC::addHalfspaceConstraint(const size_t t, unsigned char variable, const VectorXd c, const double a)
{
    assert(MPC::nx == std::make_unsigned_t<int>(c.rows()));
    assert(a != 0);

    assert(variable == 'x' || variable == 'u');

    Constraint::Halfspace con(c, a, variable, t);
    MPC::halfspace_constraints.push_back(con);
}
void MPC::projectAll()
{
    for (auto con : box_constraints)
    {
        if (con.var == 'x')
        {
            project_box(X[con.t], con.l, con.u);
        }
        else if (con.var == 'u')
        {
            project_box(U[con.t], con.l, con.u);
        }
    }
    for (auto con : ball_constraints)
    {
        if (con.var == 'x')
        {
            project_ball(X[con.t], con.r);
        }
        else if (con.var == 'u')
        {
            project_ball(U[con.t], con.r);
        }
    }
    for (auto con : halfspace_constraints)
    {
        if (con.var == 'x')
        {
            project_halfspace(X[con.t], con.c, con.a);
        }
        else if (con.var == 'u')
        {
            project_halfspace(U[con.t], con.c, con.a);
        }
    }
}
std::vector<VectorXd> MPC::getState()
{
    return X;
}
std::vector<VectorXd> MPC::getControl()
{
    return U;
}
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

    VectorXd e = VectorXd::Zero(MPC::T + 1);
    printf("iter     objv       |Gx-g|\n");
    printf("-----------------------------\n");
    for (size_t k = 1; k < 10001; k++)
    {
        auto a = 2 / ((k + 1) * eta1 + 2 * eta2);
        auto b = (k + 1) * eta1 / (2 * eta3);
        double obj = 0;
        for (size_t t = 1; t < T + 1; t++)
        {
            // V[t] = W[t] + b * (X[t] - A[t - 1] * X[t - 1] - B[t - 1] * U[t - 1]);
            // U[t - 1] = U[t - 1] - a * (R[t - 1].asDiagonal() * U[t - 1] - B[t - 1].transpose() * V[t]);
            // X[t] = X[t] - a * (Q[t - 1].asDiagonal() * X[t] + V[t] - A[t].transpose() * V[t + 1]);
            // projectAll();
            // // Only for TC
            // // if (t != T + 1)
            // // {
            // //     X[t] = X[t] - a * (Q[t] * X[t] + V[t] - A[t].transpose() * V[t + 1]);
            // // }
            // W[t] = W[t] + b * (X[t] - A[t - 1] * X[t - 1] - B[t - 1] * U[t - 1]);
            // e[t] = (X[t] - A[t - 1] * X[t - 1] - B[t - 1] * U[t - 1]).squaredNorm();
            // obj += ((X[t].array() * Q[t - 1].array() * X[t].array()).sum() + (U[t - 1].array() * R[t - 1].array() * U[t - 1].array()).sum());

            // LTI
            V[t] = W[t] + b * (X[t] - A[1] * X[t - 1] - B[1] * U[t - 1]);
            U[t - 1] = U[t - 1] - a * (R[1].asDiagonal() * U[t - 1] - B[1].transpose() * V[t]);
            X[t] = X[t] - a * (Q[1].asDiagonal() * X[t] + V[t] - A[1].transpose() * V[t + 1]);
            projectAll();
            W[t] = W[t] + b * (X[t] - A[1] * X[t - 1] - B[1] * U[t - 1]);
            e[t] = (X[t] - A[1] * X[t - 1] - B[1] * U[t - 1]).squaredNorm();
            obj += ((X[t].array() * Q[1].array() * X[t].array()).sum() + (U[t - 1].array() * R[1].array() * U[t - 1].array()).sum());
        }
        if (k % 50 == 0)
        {
            printf("%zu   %10.3e  %9.2e\n", k, obj, sqrt(e.sum()));
        }
    }
}