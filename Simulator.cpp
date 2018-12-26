//
// Created by cyk on 18-12-23.
//

#include "Simulator.h"
#include "dCmpUtil.h"
#include <iostream>
#include <cmath>
#include <utility>
#include <fstream>

using namespace Eigen;
using namespace std;

using Trid = Eigen::Triplet<double>;

void Simulator::setInlet(int radius_blue,
                         initializer_list<int> center_blue,
                         initializer_list<double> v_blue,
                         int radius_red,
                         initializer_list<int> center_red,
                         initializer_list<double> v_red) {
    v_in1 = vector<double>(v_blue);
    v_in2 = vector<double>(v_red);
    vector<int> c1 = center_blue, c2 = center_red;
    if ((c1[0] - c2[0]) * (c1[0] - c2[0]) +
        (c1[1] - c2[1]) * (c1[1] - c2[1]) <=
        (radius_blue + radius_red) *
        (radius_blue + radius_red)) {
        cerr << "invalid inlet: two inlets intersects!" << endl;
        return;
    }
    if (c1[0] - radius_blue < 0 || c1[0] + radius_blue >= N[0] ||
        c1[1] - radius_blue < 0 || c1[1] + radius_blue >= N[1] ||
        c2[0] - radius_red < 0 || c2[0] + radius_red >= N[0] ||
        c2[1] - radius_red < 0 || c2[1] + radius_red >= N[1]) {
        cerr << "invalid inlet: inlet reaches out-of-domain!" << endl;
        return;
    }
    int i, j, tmp;
    for (j = c1[1] - radius_blue; j <= c1[1] + radius_blue; ++j) {
        tmp = int(sqrt(radius_blue * radius_blue -
                       (j - c1[1]) * (j - c1[1]) + 1e-6));
        for (i = c1[0] - tmp; i <= c1[0] + tmp; ++i) {
            inlet(i, j) = C_BLUE;
        }
    }
    for (j = c2[1] - radius_red; j <= c2[1] + radius_red; ++j) {
        tmp = int(sqrt(radius_red * radius_red -
                       (j - c2[1]) * (j - c2[1]) + 1e-6));
        for (i = c2[0] - tmp; i <= c2[0] + tmp; ++i) {
            inlet(i, j) = C_RED;
        }
    }
}

void Simulator::setDefaultInLet() {
    if (N[0] < 100 || N[1] < 100) {
        cerr << "invalid inlet: domain is too small for default inlet setting!" << endl;
        return;
    }
    v_in1.resize(2);
    v_in1[0] = v_in1[1] = 3 * dx / dt;
    v_in2.resize(2);
    v_in2[0] = -v_in1[0];
    v_in2[1] = v_in1[1];
    int r = N[0] / 70;
    int cx = N[0] / 6, cy = N[1] / 3;
    int i, j, tmp;
    for (j = cy - r; j <= cy + r; ++j) {
        tmp = int(sqrt(r * r - (j - cy) * (j - cy) + 1e-6));
        for (i = cx - tmp; i <= cx + tmp; ++i) {
            inlet(i, j) = C_BLUE;
        }
    }
    cy = 2 * N[1] / 3;
    for (j = cy - r; j <= cy + r; ++j) {
        tmp = int(sqrt(r * r - (j - cy) + 1e-6));
        for (i = cx - tmp; i <= cx + tmp; ++i) {
            inlet(i, j) = C_RED;
        }
    }
    use_default_inlet = true;
}


void Simulator::setForce(const VecMatXd &f) {
    assert(N[0] == f[0].rows() &&
           N[1] == f[0].cols() &&
           N[0] == f[1].rows() &&
           N[1] == f[1].cols());
    F[0] = f[0];
    F[1] = f[0];
    force_setted = true;
}

void Simulator::setForce(double fx, double fy) {
    F[0] = Eigen::MatrixXd::Constant(N[0], N[1], fx);
    F[1] = Eigen::MatrixXd::Constant(N[0], N[1], fy);
    force_setted = true;
}

void Simulator::Forward() {
    Vstep();
    Sstep();
}

void Simulator::Vstep() {
    // newest value is stored in X0
    if (force_setted) {
        for (int i = 0; i < 2; ++i)
            addForce(U0[i], F[i]);
    }
    // newest value is stored in X0
    advect(U1, U0, U0);
    // newest value is stored in X1
    if (!isZero(visc)) {
        for (int i = 0; i < 2; ++i)
            diffuse(U1[i], visc);
    }
    cout << "before project, vx=\n" << U1[0] << endl;
    // newest value is stored in X1
    project(U0, U1);
    // newest value is stored in X0
    cout << "after project, vx=\n" << U0[0] << endl;
}

void Simulator::Sstep() {
    // newest value is stored in X0
    advect(Blue1, Blue0, U0, C_BLUE);
    advect(Red1, Red0, U0, C_RED);
    // newest value is stored in X1
    swap(Blue1, Blue0);
    swap(Red1, Red0);
    // newest value is stored in X0
    if (!isZero(diff)) {
        diffuse(Blue0, diff);
        diffuse(Red0, diff);
    }
    // newest value is stored in X0
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
    for (curr_j = 0; curr_j < N[1]; ++curr_j) {
        for (curr_i = 0; curr_i < N[0]; ++curr_i) {
            if (inlet(curr_i, curr_j) == 1) {
                for (int i = 0; i < 2; ++i)
                    vecm1[i](curr_i, curr_j) = v_in1[i];
            } else if (inlet(curr_i, curr_j) == 2) {
                for (int i = 0; i < 2; ++i)
                    vecm1[i](curr_i, curr_j) = v_in2[i];
            } else { // not belongs to any inlet
                TraceParticle(u, curr_i, curr_j, pre_i, pre_j, frac_x, frac_y);
                if (isZero(frac_x) && isZero(frac_y)) {
                    for (int i = 0; i < 2; ++i)
                        vecm1[i](curr_i, curr_j) = vecm0[i](pre_i, pre_j);
                } else if (isZero(frac_x)) {
                    vecm1[0](curr_i, curr_j) = vecm0[0](pre_i, pre_j);
                    vecm1[1](curr_i, curr_j) = LinInterp(vecm0[1], pre_i, pre_j, frac_y, false);
                } else if (isZero(frac_y)) {
                    vecm1[1](curr_i, curr_j) = vecm0[1](pre_i, pre_j);
                    vecm1[0](curr_i, curr_j) = LinInterp(vecm0[0], pre_i, pre_j, frac_x, true);
                } else {
                    vecm1[0](curr_i, curr_j) = biLinInterp(vecm0[0], pre_i, pre_j, frac_x, frac_y);
                    vecm1[1](curr_i, curr_j) = biLinInterp(vecm0[1], pre_i, pre_j, frac_x, frac_y);
                }
                // force boundary condition
                if (curr_i == 0 || curr_i == N[0] - 1) {
                    vecm1[0](curr_i, curr_j) = 0;
                }
                if (curr_j == 0 || curr_j == N[1] - 1) {
                    vecm1[1](curr_i, curr_j) = 0;
                }
            }
        }
    }
}

void Simulator::advect(MatrixXd &m1,
                       const MatrixXd &m0,
                       const VecMatXd &u,
                       int color) {
    // if not advecting color, set color=0
    int curr_i, curr_j;
    int pre_i, pre_j;
    double frac_x, frac_y;
    for (curr_j = 0; curr_j < N[1]; ++curr_j) {
        for (curr_i = 0; curr_i < N[0]; ++curr_i) {
            if (color != 0 && inlet(curr_i, curr_j) == color) {
                m1(curr_i, curr_j) = 1.0;
            } else {
                TraceParticle(u, curr_i, curr_j, pre_i, pre_j, frac_x, frac_y);
                if (isZero(frac_x) && isZero(frac_y)) {
                    m1(curr_i, curr_j) = m0(pre_i, pre_j);
                } else if (isZero(frac_x)) {
                    m1(curr_i, curr_j) = LinInterp(m0, pre_i, pre_j, frac_y, false);
                } else if (isZero(frac_y)) {
                    m1(curr_i, curr_j) = LinInterp(m0, pre_i, pre_j, frac_x, true);
                } else {
                    m1(curr_i, curr_j) = biLinInterp(m0, pre_i, pre_j, frac_x, frac_y);
                }
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
            assert(0 <= frac_x && frac_x <= 1);
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
            assert(0 <= frac_y && frac_y <= 1);
        } else { // out-of-domain in y-direction
            pre_j = pre_j < 0 ? 0 : (N[1] - 1);
            frac_y = 0;
        }
    }
    assert(pre_i >= 0 && pre_i <= N[0] - 1);
    assert(pre_j >= 0 && pre_j <= N[1] - 1);
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

void Simulator::diffuse(MatrixXd &m,
                        double k) {
    int nx = N[0] - 2, ny = N[1] - 2;
    int sz = nx * ny;
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
    int idx = 0; // idx = (i-1) + (j-1) * ny
    for (j = 1; j <= N[1] - 2; ++j) {
        for (i = 1; i <= N[0] - 2; ++i) {
            coeffs.emplace_back(Trid(idx, idx, beta));
            b(idx) += m(i, j);
            if (i == 1)
                b(idx) -= alpha * m(0, j);
            else
                coeffs.emplace_back(Trid(idx, idx - 1, alpha));
            if (i == N[0] - 2)
                b(idx) -= alpha * m(N[0] - 1, j);
            else
                coeffs.emplace_back(Trid(idx, idx + 1, alpha));
            if (j == 1)
                b(idx) -= alpha * m(i, 0);
            else
                coeffs.emplace_back(Trid(idx, idx - ny, alpha));
            if (j == N[1] - 2)
                b(idx) -= alpha * m(i, N[1] - 1);
            else
                coeffs.emplace_back(Trid(idx, idx + ny, alpha));
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
    m.block(1, 1, nx, ny) = Map<MatrixXd>(x.data(), nx, ny);
}

void Simulator::project(VecMatXd &u1, const VecMatXd &u0) {
    MatrixXd div(N[0], N[1]);
    calcDiv(div, u0);
    // ****************
    // first step: solve for q: div(grad(q))=div(u0)
    // here we solve for q/dx instead
    // here MatrixXd div is actually dx*div(u0)
    // ****************
    // build A, b
    // ****************
    int sz = N[0] * N[1];
    SparseMatrix<double> A(sz, sz);
    VectorXd b = Map<VectorXd>(div.data(), sz);
    b(0) = 0;
    VectorXd x(sz);
    vector<Trid> coeffs;
    int i, j;
    int idx = 0; // idx = i + j * N[1]
    for (j = 0; j < N[1]; ++j) {
        for (i = 0; i < N[0]; ++i) {
            if (i == 0 && j == 0) {
                coeffs.emplace_back(idx, idx, 1);
            } else {
                coeffs.emplace_back(Trid(idx, idx, -4));
                if (i == 0) {
                    coeffs.emplace_back(Trid(idx, idx + 1, 1));
                } else {
                    coeffs.emplace_back(Trid(idx, idx - 1, 1));
                }
                if (i == N[0] - 1) {
                    coeffs.emplace_back(Trid(idx, idx - 1, 1));
                } else {
                    coeffs.emplace_back(Trid(idx, idx + 1, 1));
                }
                if (j == 0) {
                    coeffs.emplace_back(Trid(idx, idx + N[1], 1));
                } else {
                    coeffs.emplace_back(Trid(idx, idx - N[1], 1));
                }
                if (j == N[1] - 1) {
                    coeffs.emplace_back(Trid(idx, idx - N[1], 1));
                } else {
                    coeffs.emplace_back(Trid(idx, idx + N[1], 1));
                }
            }
            ++idx;
        }
    }
    A.setFromTriplets(coeffs.begin(), coeffs.end());
    // ****************
    // solve
    // ****************
    // todo: analyzePattern and factorize
    ConjugateGradient<SparseMatrix<double >> cg;
    cg.compute(A);
    if (cg.info() != Success) {
        cerr << "decomposition in \"project\" failed!" << endl;
        return;
    }
    x = cg.solve(b);
    if (cg.info() != Success) {
        cerr << "solving in \"project\" failed!" << endl;
        return;
    }
    div = Map<MatrixXd>(x.data(), N[0], N[1]);
    // ****************
    // from now, MatrixXd div stores the value for q/dx
    // final step: obtain divergence-free u1
    // ****************
    for (j = 0; j < N[1]; ++j) {
        for (i = 0; i < N[0]; ++i) {
            if (i == 0 || i == N[0] - 1) {
                u1[0](i, j) = 0;
            } else {
                u1[0](i, j) = u0[0](i, j) - 0.5 * (div(i + 1, j) - div(i - 1, j));
            }
            if (j == 0 || j == N[1] - 1) {
                u1[1](i, j) = 0;
            } else {
                u1[1](i, j) = u0[1](i, j) - 0.5 * (div(i, j + 1) - div(i, j - 1));
            }
        }
    }
}


void Simulator::calcDiv(MatrixXd &m, const VecMatXd &u) {
    // store result in m:
    // m = dx * div(u)
    assert(m.rows() == N[0] && m.cols() == N[1]);
    m.setZero();
    int i, j;
    for (j = 0; j < N[1]; ++j) {
        for (i = 0; i < N[0]; ++i) {
            if (i == 0) {
                m(i, j) += -1.5 * u[0](i, j) + 2 * u[0](i + 1, j) - 0.5 * u[0](i + 2, j);
            } else if (i == N[0] - 1) {
                m(i, j) += 1.5 * u[0](i, j) - 2 * u[0](i - 1, j) + 0.5 * u[0](i - 2, j);
            } else {
                m(i, j) += 0.5 * (u[0](i + 1, j) - u[0](i - 1, j));
            }
            if (j == 0) {
                m(i, j) += -1.5 * u[1](i, j) + 2 * u[1](i, j + 1) - 0.5 * u[1](i, j + 2);
            } else if (j == N[1] - 1) {
                m(i, j) += 1.5 * u[1](i, j) - 2 * u[1](i, j - 1) + 0.5 * u[1](i, j - 2);
            } else {
                m(i, j) += 0.5 * (u[1](i, j + 1) - u[1](i, j - 1));
            }
        }
    }
}

void Simulator::printBlue() {
    cout << "blue=\n" << Blue0 << endl;
}

void Simulator::printRed() {
    cout << "red=\n" << Red0 << endl;
}

void Simulator::printVx() {
    cout << "vx=\n" << U0[0] << endl;
}

void Simulator::printVy() {
    cout << "vy=\n" << U0[1] << endl;
}