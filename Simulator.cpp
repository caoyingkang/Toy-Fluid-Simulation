//
// Created by cyk on 18-12-23.
//

#include "Simulator.h"
#include "dCmpUtil.h"
#include <iostream>

using namespace Eigen;
using namespace std;


using Trid = Eigen::Triplet<double>;

void Simulator::Vstep() {
    for (int i = 0; i < 2; ++i)
        addForce(U0[i], F[i]);
    advect(U1, U0, U0);
    for (int i = 0; i < 2; ++i)
        diffuse(U0[i], U1[i], visc);
    project();
}

void Simulator::addForce(MatrixXd &m,
                         const MatrixXd &f) {
    m += dt * f;
}

void Simulator::advect(VecMatXd &vecm1,
                       const VecMatXd &vecm0,
                       const VecMatXd &u) {
    int curr_i, curr_j;
    int pre_i, pre_j;
    double frac_x, frac_y;
    for (curr_i = 0; curr_i < N[0]; ++curr_i) {
        for (curr_j = 0; curr_j < N[1]; ++curr_j) {
            TraceParticle(u, curr_i, curr_j, pre_i, pre_j, frac_x, frac_y);
            if (isZero(frac_x) && isZero(frac_y)) {
                for (int i = 0; i < 2; ++i)
                    vecm1[i](curr_i, curr_j) = vecm0[i](pre_i, pre_j);
            } else if (isZero(frac_x)) {
                for (int i = 0; i < 2; ++i)
                    vecm1[i](curr_i, curr_j) = LinInterp(vecm0[i], pre_i, pre_j, frac_y, false);
            } else if (isZero(frac_y)) {
                for (int i = 0; i < 2; ++i)
                    vecm1[i](curr_i, curr_j) = LinInterp(vecm0[i], pre_i, pre_j, frac_x, true);
            } else {
                for (int i = 0; i < 2; ++i)
                    vecm1[i](curr_i, curr_j) = biLinInterp(vecm0[i], pre_i, pre_j, frac_x, frac_y);
            }
        }
    }
}

void Simulator::advect(MatrixXd &m1,
                       const MatrixXd &m0,
                       const VecMatXd &u) {
    int curr_i, curr_j;
    int pre_i, pre_j;
    double frac_x, frac_y;
    for (curr_i = 0; curr_i < N[0]; ++curr_i) {
        for (curr_j = 0; curr_j < N[1]; ++curr_j) {
            TraceParticle(u, curr_i, curr_j, pre_i, pre_j, frac_x, frac_y);
            if (isZero(frac_x) && isZero(frac_y)) {
                for (int i = 0; i < 2; ++i)
                    m1(curr_i, curr_j) = m0(pre_i, pre_j);
            } else if (isZero(frac_x)) {
                for (int i = 0; i < 2; ++i)
                    m1(curr_i, curr_j) = LinInterp(m0, pre_i, pre_j, frac_y, false);
            } else if (isZero(frac_y)) {
                for (int i = 0; i < 2; ++i)
                    m1(curr_i, curr_j) = LinInterp(m0, pre_i, pre_j, frac_x, true);
            } else {
                for (int i = 0; i < 2; ++i)
                    m1(curr_i, curr_j) = biLinInterp(m0, pre_i, pre_j, frac_x, frac_y);
            }
        }
    }
}

void Simulator::TraceParticle(const VecMatXd &u,
                              int curr_i, int curr_j,
                              int &pre_i, int &pre_j,
                              double &frac_x, double &frac_y) {
    double Di = -dt * u[0](curr_i, curr_j) / dx;
    double Dj = -dt * u[1](curr_i, curr_j) / dx;
    int ni, nj;
    if (isZero(Di)) { // zero x-velocity at this point
        pre_i = curr_i;
        frac_x = 0;
    } else {
        ni = int(Di) - (Di > 0 ? 0 : 1);
        pre_i = curr_i + ni;
        if (pre_i >= 0 && pre_i <= N[0] - 2) {
            frac_x = Di - ni;
        } else { // out-of-domain in x-direction
            pre_i = pre_i < 0 ? 0 : (N[0] - 1);
            frac_x = 0;
        }
    }
    if (isZero(Dj)) { // zero y-velocity at this point
        pre_j = curr_j;
        frac_y = 0;
    } else {
        nj = int(Dj) - (Dj > 0 ? 0 : 1);
        pre_j = curr_j + nj;
        if (pre_j >= 0 && pre_j <= N[1] - 2) {
            frac_y = Dj - nj;
        } else { // out-of-domain in y-direction
            pre_j = pre_j < 0 ? 0 : (N[1] - 1);
            frac_y = 0;
        }
    }
    assert(-EPS < frac_x < 1 + EPS && -EPS < frac_y < 1 + EPS);
    assert(pre_i >= 0 && pre_i <= N[0] - 2);
    assert(pre_j >= 0 && pre_j <= N[1] - 2);
}


double Simulator::LinInterp(const MatrixXd &m,
                            int i, int j,
                            double frac, bool Along_X) {
    assert(Along_X ? (i < N[0] - 1) : (j < N[1] - 1));
    if (Along_X) {
        return (1 - frac) * m(i, j) + frac * m(i + 1, j);
    } else {
        return (1 - frac) * m(i, j) + frac * m(i, j + 1);
    }
}

double Simulator::biLinInterp(const MatrixXd &m,
                              int i, int j,
                              double frac_x, double frac_y) {
    // bi-linear interpolation between grid points (i,j), (i+1,j), (i,j+1), (i+1,j+1)
    // the point to be estimated is (i+frac_x,j+frac_y)
    assert((i < N[0] - 1) && (j < N[1] - 1));
    double tmp1 = (1 - frac_x) * m(i, j) + frac_x * m(i + 1, j);
    double tmp2 = (1 - frac_x) * m(i, j + 1) + frac_x * m(i + 1, j + 1);
    return (1 - frac_y) * tmp1 + frac_y * tmp2;
}

void Simulator::diffuse(MatrixXd &m1,
                        const MatrixXd &m0,
                        double k) {
    int mx = N[0] - 2, my = N[1] - 2;
    int sz = mx * my;
    double alpha = -dt * k / (dx * dx);
    double beta = 1 - 4 * alpha;
    // ****************
    // build A, b
    // ****************
    SparseMatrix<double> A(sz, sz);
    VectorXd x(sz), b(sz);
    b.setZero();
    vector<Trid> coeffs;
    int i, j;
    int idx = 0; // idx = (i-1) + (j-1) * my
    for (j = 1; j <= N[1] - 2; ++j) {
        for (i = 1; i <= N[0] - 2; ++i) {
            coeffs.emplace_back(Trid(idx, idx, beta));
            b(idx) += m0(i, j);
            if (i == 1)
                b(idx) -= alpha * m0(0, j);
            else
                coeffs.emplace_back(Trid(idx, idx - 1, alpha));
            if (i == N[0] - 2)
                b(idx) -= alpha * m0(N[0] - 1, j);
            else
                coeffs.emplace_back(Trid(idx, idx + 1, alpha));
            if (j == 1)
                b(idx) -= alpha * m0(i, 0);
            else
                coeffs.emplace_back(Trid(idx, idx - my, alpha));
            if (j == N[1] - 2)
                b(idx) -= alpha * m0(i, N[1] - 1);
            else
                coeffs.emplace_back(Trid(idx, idx + my, alpha));
            ++idx;
        }
    }
    A.setFromTriplets(coeffs.begin(), coeffs.end());
    // ****************
    // solve
    // ****************
    // todo: analyzePattern and factorize
    ConjugateGradient<SparseMatrix<double>, Lower | Upper, IncompleteLUT<double>> cg;
    // cg.setMaxIterations(maxIt);
    // cg.setTolerance(tol);
    cg.compute(A);
    if (cg.info() != Success) {
        cerr << "decomposition failed!" << endl;
        return;
    }
    x = cg.solve(b);
    if (cg.info() != Success) {
        cerr << "solving failed!" << endl;
        return;
    }
    m1.block(1, 1, mx, my) = Map<MatrixXd>(x.data(), mx, my);
    for (i = 0; i < N[0]; ++i)
        m1(i, 0) = m0(i, 0);
    for (i = 0; i < N[0]; ++i)
        m1(i, N[1] - 1) = m0(i, N[1] - 1);
    for (j = 1; j <= N[1] - 2; ++j) {
        m1(0, j) = m0(0, j);
        m1(N[0] - 1, j) = m0(N[0] - 1, j);
    }
}

void Simulator::project() {

}

void Simulator::Sstep() {

}