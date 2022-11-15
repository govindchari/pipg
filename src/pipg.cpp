#include "pipg.hpp"
#include <iostream>

MPC::MPC(const size_t T, const size_t nx, const size_t nu)
{
    _T = T;
    _nx = nx;
    _nu = nu;
    _A.resize(_T + 1); // Of length T+1 because A_T = 0 (dummy variable for alg to work)
    _B.resize(_T);
    _Q.resize(_T);
    _R.resize(_T);
    _X.resize(_T + 1);
    _U.resize(_T);
    _V.resize(_T + 2); // Of length T+2 because V_{T+1} = 0 (dummy variable for alg to work)
    _W.resize(_T + 1);

    _box_constraints.resize(_T + 1);
    _ball_constraints.resize(_T + 1);
    _halfspace_constraints.resize(_T + 1);

    std::fill(_A.begin(), _A.end(), MatrixXd::Zero(_nx, _nx));
    std::fill(_B.begin(), _B.end(), MatrixXd::Zero(_nx, _nu));
    std::fill(_Q.begin(), _Q.end(), VectorXd::Zero(_nx));
    std::fill(_R.begin(), _R.end(), VectorXd::Zero(_nu));
    std::fill(_X.begin(), _X.end(), VectorXd::Zero(_nx));
    std::fill(_U.begin(), _U.end(), VectorXd::Zero(_nu));
    std::fill(_V.begin(), _V.end(), VectorXd::Zero(_nx));
    std::fill(_W.begin(), _W.end(), VectorXd::Zero(_nx));

    _H = MatrixXd::Zero(_T * _nx, _T * (_nx + _nu));

    for (size_t i = 0; i < _T; i++)
    {
        _H.block(i * _nx, _nu + i * (_nx + _nu), _nx, _nx) = MatrixXd::Identity(_nx, _nx);
    }

    _eta1 = 0.0;
    _eta2 = 0.0;
    _eta3 = 0.0;
    _eta1_outdated = false;
    _eta2_outdated = false;
    _eta3_outdated = false;

    WORKSPACE ws(nx, nu, T);
    TOLERANCE tol;
}
void MPC::printQR()
{
    std::cout << "Q: " << std::endl;
    for (auto x : _Q)
    {
        std::cout << x << std::endl;
    }
    std::cout << "R: " << std::endl;
    for (auto x : _R)
    {
        std::cout << x << std::endl;
    }
    std::cout << _eta1 << std::endl;
    std::cout << _eta2 << std::endl;
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
    for (auto x : _Q)
    {
        if (x.minCoeff() < minElement)
        {
            minElement = x.minCoeff();
        }
    }
    for (auto x : _R)
    {
        if (x.minCoeff() < minElement)
        {
            minElement = x.minCoeff();
        }
    }
    _eta1 = minElement;
    _eta1_outdated = false;
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
    for (auto x : _Q)
    {
        if (x.maxCoeff() > maxElement)
        {
            maxElement = x.maxCoeff();
        }
    }
    for (auto x : _R)
    {
        if (x.maxCoeff() > maxElement)
        {
            maxElement = x.maxCoeff();
        }
    }
    _eta2 = maxElement;
    _eta2_outdated = false;
}
void MPC::updateEta3()
{
    // Compute HtH
    // Use power method to get max eig
    auto HtH = (_H.transpose() * _H);
    auto eig = HtH.eigenvalues();
    _eta3 = eig[0].real();
    _eta3_outdated = false;
}
void MPC::addA(const size_t t, const MatrixXd A)
{
    assert(t >= 0 && t < _T);
    assert(std::make_unsigned_t<int>(A.rows()) == _nx && std::make_unsigned_t<int>(A.cols()) == _nx);
    _A[t] = A;
    if (t != 0)
    {
        _H.block(t * _nx, _nu + (t - 1) * (_nu + _nx), _nx, _nx) = -A;
    }
    _eta3_outdated = true;
}
void MPC::addB(const size_t t, const MatrixXd B)
{
    assert(t >= 0 && t < _T);
    assert(std::make_unsigned_t<int>(B.rows()) == _nx && std::make_unsigned_t<int>(B.cols()) == _nu);
    _B[t] = B;
    _H.block(t * _nx, t * (_nu + _nx), _nx, _nu) = -B;
    _eta3_outdated = true;
}
void MPC::addQ(const size_t t, const VectorXd Q)
{
    assert(t > 0 && t <= _T);
    _Q[t - 1] = Q;
    _eta1_outdated = true;
    _eta2_outdated = true;
}
void MPC::addR(const size_t t, const VectorXd R)
{
    assert(t >= 0 && t < _T);
    assert(std::make_unsigned_t<int>(R.rows()) == _nu && std::make_unsigned_t<int>(R.cols()) == _nu);
    _R[t] = R;
    _eta1_outdated = true;
    _eta2_outdated = true;
}
void MPC::addIC(const VectorXd x0)
{
    assert(_nx == std::make_unsigned_t<int>(x0.rows()));
    _X[0] = x0;
}
void MPC::addBoxConstraint(const size_t t, const char variable, const VectorXd l, const VectorXd u)
{
    assert(variable == 'x' || variable == 'u');
    assert(l.rows() == u.rows());
    assert(l.cwiseMin(u) == l);
    if (variable == 'x')
    {
        assert(_nx == std::make_unsigned_t<int>(l.rows()));
        assert(_nx == std::make_unsigned_t<int>(u.rows()));
        Constraint::Box con(l, u, variable);
        _box_constraints[t].push_back(con);
    }
    if (variable == 'u')
    {
        assert(_nu == std::make_unsigned_t<int>(l.rows()));
        assert(_nu == std::make_unsigned_t<int>(u.rows()));
        Constraint::Box con(l, u, variable);
        _box_constraints[t + 1].push_back(con);
    }
}
void MPC::addBallConstraint(const size_t t, const char variable, const double r)
{
    assert(r > 0);

    assert(variable == 'x' || variable == 'u');
    if (variable == 'x')
    {
        Constraint::Ball con(r, variable);
        _ball_constraints[t].push_back(con);
    }
    if (variable == 'u')
    {
        Constraint::Ball con(r, variable);
        _ball_constraints[t + 1].push_back(con);
    }
}
void MPC::addHalfspaceConstraint(const size_t t, const char variable, const VectorXd c, const double a)
{
    assert(variable == 'x' || variable == 'u');
    assert(a != 0);
    if (variable == 'x')
    {
        assert(_nx == std::make_unsigned_t<int>(c.rows()));
        Constraint::Halfspace con(c, a, variable);
        _halfspace_constraints[t].push_back(con);
    }
    if (variable == 'u')
    {
        assert(_nu == std::make_unsigned_t<int>(c.rows()));
        Constraint::Halfspace con(c, a, variable);
        _halfspace_constraints[t + 1].push_back(con);
    }
}
void MPC::project(const size_t t)
{
    for (auto con : _box_constraints[t])
    {
        if (con._var == 'x')
        {
            project_box(_X[t], con._l, con._u);
        }
        else if (con._var == 'u')
        {
            project_box(_U[t - 1], con._l, con._u);
        }
    }
    for (auto con : _ball_constraints[t])
    {
        if (con._var == 'x')
        {
            project_ball(_X[t], con._r);
        }
        else if (con._var == 'u')
        {
            project_ball(_U[t - 1], con._r);
        }
    }
    for (auto con : _halfspace_constraints[t])
    {
        if (con._var == 'x')
        {
            project_halfspace(_X[t], con._c, con._a);
        }
        else if (con._var == 'u')
        {
            project_halfspace(_U[t - 1], con._c, con._a);
        }
    }
}
std::vector<VectorXd> MPC::getState()
{
    return _X;
}
std::vector<VectorXd> MPC::getControl()
{
    return _U;
}
void MPC::solve()
{
    solve(true);
}
void MPC::solve(bool verbose)
{
    if (_eta1_outdated)
    {
        updateEta1();
    }
    if (_eta2_outdated)
    {
        updateEta2();
    }
    if (_eta3_outdated)
    {
        updateEta3();
    }

    VectorXd e = VectorXd::Zero(_T + 1);
    if (verbose)
    {
        printf("iter     objv       |Gx-g|\n");
        printf("-----------------------------\n");
    }
    while (!_tol.stop)
    {
        auto a = 2 / ((_tol.k + 1) * _eta1 + 2 * _eta2);
        auto b = (_tol.k + 1) * _eta1 / (2 * _eta3);
        double obj = 0;
        for (size_t t = 1; t < _T + 1; t++)
        {
            _V[t] = _W[t] + b * (_X[t] - _A[t - 1] * _X[t - 1] - _B[t - 1] * _U[t - 1]);
            _U[t - 1] = _U[t - 1] - a * (_R[t - 1].asDiagonal() * _U[t - 1] - _B[t - 1].transpose() * _V[t]);
            _X[t] = _X[t] - a * (_Q[t - 1].asDiagonal() * _X[t] + _V[t] - _A[t].transpose() * _V[t + 1]);
            project(t);
            // Only for TC
            // if (t != T + 1)
            // {
            //     X[t] = X[t] - a * (Q[t] * X[t] + V[t] - A[t].transpose() * V[t + 1]);
            // }
            _W[t] = _W[t] + b * (_X[t] - _A[t - 1] * _X[t - 1] - _B[t - 1] * _U[t - 1]);
            if (_tol.k % 50 == 0)
            {
                e[t] = (_X[t] - _A[t - 1] * _X[t - 1] - _B[t - 1] * _U[t - 1]).lpNorm<Eigen::Infinity>();
                obj += ((_X[t].array() * _Q[t - 1].array() * _X[t].array()).sum() + (_U[t - 1].array() * _R[t - 1].array() * _U[t - 1].array()).sum());
            }
        }
        if (_tol.k % 50 == 0)
        {
            if (verbose)
                printf("%zu   %10.3e  %9.2e\n", _tol.k, obj, e.lpNorm<Eigen::Infinity>());
            _tol.stop = e.lpNorm<Eigen::Infinity>() < _tol.eq_tol || _tol.k > _tol.max_iter;
        }
        _tol.k++;
    }
}