//
// Created by cyk on 18-12-23.
//

#include "Simulator.h"
#include "cmpUtil.h"
#include "cgUtil.h"
#include <iostream>
#include <cmath>
#include <utility>
#include <fstream>
#include <algorithm>

using namespace Eigen;
using namespace std;


void Simulator::setInlet(int radius_blue,
                         initializer_list<int> center_blue,
                         initializer_list<float> v_blue,
                         int radius_red,
                         initializer_list<int> center_red,
                         initializer_list<float> v_red) {
    v_in1 = vector<float>(v_blue);
    v_in2 = vector<float>(v_red);
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
            Blue0(i, j) = 1;
        }
    }
    for (j = c2[1] - radius_red; j <= c2[1] + radius_red; ++j) {
        tmp = int(sqrt(radius_red * radius_red -
                       (j - c2[1]) * (j - c2[1]) + 1e-6));
        for (i = c2[0] - tmp; i <= c2[0] + tmp; ++i) {
            inlet(i, j) = C_RED;
            Red0(i, j) = 1;
        }
    }
}


void Simulator::setForce(const VecMatXf &f) {
    assert(N[0] == f[0].rows() &&
           N[1] == f[0].cols() &&
           N[0] == f[1].rows() &&
           N[1] == f[1].cols());
    F[0] = f[0];
    F[1] = f[0];
    force_setted = true;
}

void Simulator::setForce(float fx, float fy) {
    F[0] = Eigen::MatrixXf::Constant(N[0], N[1], fx);
    F[1] = Eigen::MatrixXf::Constant(N[0], N[1], fy);
    force_setted = true;
}

void Simulator::setVisc(float val) {
    visc = val;
    visc_setted = true;
}

void Simulator::setDiff(float val) {
    diff = val;
    diff_setted = true;
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
    if (visc_setted && !isZero(visc)) {
        for (int i = 0; i < 2; ++i)
            diffuse(U1[i], visc);
    }

    // newest value is stored in X1
    project(U0, U1);

    // newest value is stored in X0
}

void Simulator::Sstep() {
    // newest value is stored in X0
    advect(Blue1, Blue0, U0, C_BLUE);
    advect(Red1, Red0, U0, C_RED);

    // newest value is stored in X1
    swap(Blue1, Blue0);
    swap(Red1, Red0);

    // newest value is stored in X0
    if (diff_setted && !isZero(diff)) {
        diffuse(Blue0, diff);
        diffuse(Red0, diff);
    }

    // newest value is stored in X0
}

void Simulator::addForce(MatrixXf &m,
                         const MatrixXf &f) {
    m += dt * f;
}

void Simulator::advect(VecMatXf &vecm1,
                       const VecMatXf &vecm0,
                       const VecMatXf &u) {
    int i, j;
    int pre_i, pre_j;
    float frac_x, frac_y;
    for (j = 0; j < N[1]; ++j) {
        for (i = 0; i < N[0]; ++i) {
            if (inlet(i, j) == 1) {
                for (int d = 0; d < 2; ++d)
                    vecm1[d](i, j) = v_in1[d];
            } else if (inlet(i, j) == 2) {
                for (int d = 0; d < 2; ++d)
                    vecm1[d](i, j) = v_in2[d];
            } else { // not belongs to any inlet
                TraceParticle(u, i, j, pre_i, pre_j, frac_x, frac_y);
                if (isZero(frac_x) && isZero(frac_y)) {
                    for (int d = 0; d < 2; ++d)
                        vecm1[d](i, j) = vecm0[d](pre_i, pre_j);
                } else if (isZero(frac_x)) {
                    vecm1[0](i, j) = vecm0[0](pre_i, pre_j);
                    vecm1[1](i, j) = LinInterp(vecm0[1], pre_i, pre_j, frac_y, false);
                } else if (isZero(frac_y)) {
                    vecm1[1](i, j) = vecm0[1](pre_i, pre_j);
                    vecm1[0](i, j) = LinInterp(vecm0[0], pre_i, pre_j, frac_x, true);
                } else {
                    vecm1[0](i, j) = biLinInterp(vecm0[0], pre_i, pre_j, frac_x, frac_y);
                    vecm1[1](i, j) = biLinInterp(vecm0[1], pre_i, pre_j, frac_x, frac_y);
                }
            }
        }
    }
}

void Simulator::advect(MatrixXf &m1,
                       const MatrixXf &m0,
                       const VecMatXf &u,
                       int color) {
    // if not advecting color, set color=0
    int i, j;
    int pre_i, pre_j;
    float frac_x, frac_y;
    for (j = 0; j < N[1]; ++j) {
        for (i = 0; i < N[0]; ++i) {
            if (color != 0 && inlet(i, j) == color) {
                m1(i, j) = 1.0;
            } else {
                TraceParticle(u, i, j, pre_i, pre_j, frac_x, frac_y);
                if (isZero(frac_x) && isZero(frac_y)) {
                    m1(i, j) = m0(pre_i, pre_j);
                } else if (isZero(frac_x)) {
                    m1(i, j) = LinInterp(m0, pre_i, pre_j, frac_y, false);
                } else if (isZero(frac_y)) {
                    m1(i, j) = LinInterp(m0, pre_i, pre_j, frac_x, true);
                } else {
                    m1(i, j) = biLinInterp(m0, pre_i, pre_j, frac_x, frac_y);
                }
            }
        }
    }
}


void Simulator::TraceParticle(const VecMatXf &u,
                              int curr_i, int curr_j,
                              int &pre_i, int &pre_j,
                              float &frac_x, float &frac_y) {
    float Di = -dt * u[0](curr_i, curr_j) / dx;
    float Dj = -dt * u[1](curr_i, curr_j) / dx;
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


float Simulator::LinInterp(const MatrixXf &m,
                           int i, int j,
                           float frac, bool Along_X) {
    assert(Along_X ? (i < N[0] - 1) : (j < N[1] - 1));
    if (Along_X) {
        return (1 - frac) * m(i, j) + frac * m(i + 1, j);
    } else {
        return (1 - frac) * m(i, j) + frac * m(i, j + 1);
    }
}

float Simulator::biLinInterp(const MatrixXf &m,
                             int i, int j,
                             float frac_x, float frac_y) {
    // bi-linear interpolation between grid points (i,j), (i+1,j), (i,j+1), (i+1,j+1)
    // the point to be estimated is (i+frac_x,j+frac_y)
    assert((i < N[0] - 1) && (j < N[1] - 1));
    float tmp1 = (1 - frac_x) * m(i, j) + frac_x * m(i + 1, j);
    float tmp2 = (1 - frac_x) * m(i, j + 1) + frac_x * m(i + 1, j + 1);
    return (1 - frac_y) * tmp1 + frac_y * tmp2;
}


void Simulator::diffuse(MatrixXf &m,
                        float k) {
    int nx = N[0] - 2, ny = N[1] - 2;
    int sz = nx * ny;
    float alpha = -dt * k / (dx * dx);
    float beta = 1 - 4 * alpha;
    // ****************
    // build A, b
    // ****************
    cgMat A(nx, ny);
    VectorXf x(sz), b(sz);
    b.setZero();
    int i, j, ai, aj;
    int idx = 0; // idx = (i-1) + (j-1) * nx
    for (j = 1; j <= N[1] - 2; ++j) {
        for (i = 1; i <= N[0] - 2; ++i) {
            ai = i - 1;
            aj = j - 1;
            A.Adiag(ai, aj) = beta;
            b(idx) += m(i, j);
            if (i == 1)
                b(idx) -= alpha * m(0, j);
            if (i == N[0] - 2)
                b(idx) -= alpha * m(N[0] - 1, j);
            else
                A.Aplusi(ai, aj) = alpha;
            if (j == 1)
                b(idx) -= alpha * m(i, 0);
            if (j == N[1] - 2)
                b(idx) -= alpha * m(i, N[1] - 1);
            else
                A.Aplusj(ai, aj) = alpha;
            ++idx;
        }
    }
    // ****************
    // solve
    // ****************
    int iters;
    cgSolve(x, A, b, 0.01, 100, CG::MIC, &iters);
    //cgSolve(x, A, b, 0.01, 100, CG::Identity, &iters);
    cout << "diffuse phase: iters in CG::cgSolve:    " << iters << endl;
    m.block(1, 1, nx, ny) = Map<MatrixXf>(x.data(), nx, ny);
}


void Simulator::project(VecMatXf &u1, const VecMatXf &u0) {
    MatrixXf div(N[0], N[1]);
    calcDiv(div, u0);

    // ****************
    // first step: solve for q: div(grad(q))=div(u0)
    // here we solve for q/dx instead
    // here MatrixXf div is actually dx*div(u0)
    // ****************
    // build A, b
    // ****************
    int sz = N[0] * N[1];
    VectorXf x(sz);
    VectorXf b = Map<VectorXf>(div.data(), sz);
    b = -b;

    cgMat A(N[0], N[1]);
    int i, j;
    int idx = 0; // idx = i + j * N[0]
    for (j = 0; j < N[1]; ++j) {
        for (i = 0; i < N[0]; ++i) {
            if (isEmpty(i, j)) {
                A.Adiag(i, j) = 1;
                b(idx) = 0;
            } else {
                A.Adiag(i, j) += 4;
                if (i == 0 || i == N[0] - 1) {
                    A.Adiag(i, j) += -1;
                    b(idx) += (i == 0 ? -u0[0](0, j) : u0[0](N[0] - 1, j));
                }
                if (i != N[0] - 1 && !isEmpty(i + 1, j)) {
                    A.Aplusi(i, j) += -1;
                }
                if (j == 0 || j == N[1] - 1) {
                    A.Adiag(i, j) += -1;
                    b(idx) += (j == 0 ? -u0[1](i, 0) : u0[1](i, N[1] - 1));
                }
                if (j != N[1] - 1 && !isEmpty(i, j + 1)) {
                    A.Aplusj(i, j) += -1;
                }
            }
            ++idx;
        }
    }

    // ****************
    // solve
    // ****************

    int iters;
    cgSolve(x, A, b, 0.01, 100, CG::MIC, &iters);
    //cgSolve(x, A, b, 0.1, 100, CG::Identity, &iters);
    cout << "project phase: iters in CG::cgSolve:    " << iters << endl;

    div = Map<MatrixXf>(x.data(), N[0], N[1]);

    // ****************
    // from now, MatrixXf div stores the value for q/dx
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

void Simulator::calcDiv(MatrixXf &m, const VecMatXf &u) {
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


void Simulator::getRenderData(float *vertices) {
    int i, j;
    int idx = 2; // index of vertices
    for (j = 0; j < N[1] - 1; ++j) {
        for (i = 0; i < N[0] - 1; ++i) {
            vertices[idx] = Red0(i, j);
            vertices[idx + 1] = Blue0(i, j);
            idx += 4;
            vertices[idx] = Red0(i + 1, j);
            vertices[idx + 1] = Blue0(i + 1, j);
            idx += 4;
            vertices[idx] = Red0(i, j + 1);
            vertices[idx + 1] = Blue0(i, j + 1);
            idx += 4;
            vertices[idx] = Red0(i + 1, j);
            vertices[idx + 1] = Blue0(i + 1, j);
            idx += 4;
            vertices[idx] = Red0(i, j + 1);
            vertices[idx + 1] = Blue0(i, j + 1);
            idx += 4;
            vertices[idx] = Red0(i + 1, j + 1);
            vertices[idx + 1] = Blue0(i + 1, j + 1);
            idx += 4;
        }
    }
}


bool Simulator::isEmpty(int i, int j) const {
    // strategy 1
    // ----------------------------
    return i == 0 && j == 0;
    // strategy 2
    // ----------------------------
    //return isZero(Blue0(i, j)) && isZero(Red0(i, j));
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

