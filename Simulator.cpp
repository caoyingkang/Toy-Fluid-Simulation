//
// Created by cyk on 18-12-23.
//

#include "Simulator.h"
#include "dCmpUtil.h"
#include "cgUtil.h"
#include <iostream>
#include <cmath>
#include <utility>
#include <fstream>
#include <algorithm>
#include <ctime>

using namespace Eigen;
using namespace std;

using Trid = Eigen::Triplet<float>;


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

void Simulator::setDefaultInLet() {
    if (N[0] < 100 || N[1] < 100) {
        cerr << "invalid inlet: domain is too small for default inlet setting!" << endl;
        return;
    }
    v_in1.resize(2);
    v_in1[0] = v_in1[1] = 10 * dx / dt;
    v_in2.resize(2);
    v_in2[0] = v_in1[0];
    v_in2[1] = +v_in1[1];
    int r = N[0] / 70;
    int cx = N[0] / 6, cy = N[1] / 3;
    int i, j, tmp;
    for (j = cy - r; j <= cy + r; ++j) {
        tmp = int(sqrt(r * r - (j - cy) * (j - cy) + 1e-6));
        for (i = cx - tmp; i <= cx + tmp; ++i) {
            inlet(i, j) = C_BLUE;
            Blue0(i, j) = 1;
        }
    }
    cy = 2 * N[1] / 3;
    for (j = cy - r; j <= cy + r; ++j) {
        tmp = int(sqrt(r * r - (j - cy) + 1e-6));
        for (i = cx - tmp; i <= cx + tmp; ++i) {
            inlet(i, j) = C_RED;
            Red0(i, j) = 1;
        }
    }
    use_default_inlet = true;
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
    if (!isZero(diff)) {
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
    int curr_i, curr_j;
    int pre_i, pre_j;
    float frac_x, frac_y;
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
//                // force boundary condition
//                if (curr_i == 0 || curr_i == N[0] - 1) {
//                    vecm1[0](curr_i, curr_j) = 0;
//                }
//                if (curr_j == 0 || curr_j == N[1] - 1) {
//                    vecm1[1](curr_i, curr_j) = 0;
//                }
            }
        }
    }
}

void Simulator::advect(MatrixXf &m1,
                       const MatrixXf &m0,
                       const VecMatXf &u,
                       int color) {
    // if not advecting color, set color=0
    int curr_i, curr_j;
    int pre_i, pre_j;
    float frac_x, frac_y;
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

//void Simulator::diffuse(MatrixXf &m,
//                        float k) {
//    int nx = N[0] - 2, ny = N[1] - 2;
//    int sz = nx * ny;
//    float alpha = -dt * k / (dx * dx);
//    float beta = 1 - 4 * alpha;
//    // ****************
//    // build A, b
//    // ****************
//    SparseMatrix<float> A(sz, sz);
//    VectorXf
//            x(sz), b(sz);
//    b.setZero();
//    vector<Trid> coeffs;
//    int i, j;
//    int idx = 0; // idx = (i-1) + (j-1) * nx
//    for (j = 1; j <= N[1] - 2; ++j) {
//        for (i = 1; i <= N[0] - 2; ++i) {
//            coeffs.emplace_back(Trid(idx, idx, beta));
//            b(idx) += m(i, j);
//            if (i == 1)
//                b(idx) -= alpha * m(0, j);
//            else
//                coeffs.emplace_back(Trid(idx, idx - 1, alpha));
//            if (i == N[0] - 2)
//                b(idx) -= alpha * m(N[0] - 1, j);
//            else
//                coeffs.emplace_back(Trid(idx, idx + 1, alpha));
//            if (j == 1)
//                b(idx) -= alpha * m(i, 0);
//            else
//                coeffs.emplace_back(Trid(idx, idx - nx, alpha));
//            if (j == N[1] - 2)
//                b(idx) -= alpha * m(i, N[1] - 1);
//            else
//                coeffs.emplace_back(Trid(idx, idx + nx, alpha));
//            ++idx;
//        }
//    }
//    A.setFromTriplets(coeffs.begin(), coeffs.end());
//    // ****************
//    // solve
//    // ****************
//    // todo: analyzePattern and factorize
//    ConjugateGradient<SparseMatrix<float>, Lower | Upper,
//            IncompleteCholesky<float >> cg;
//    // cg.setMaxIterations(maxIt);
//    // cg.setTolerance(tol);
//    cg.compute(A);
//    if (cg.info() != Success) {
//        cerr << "decomposition failed!" << endl;
//        return;
//    }
//    x = cg.solve(b);
//    if (cg.info() != Success) {
//        cerr << "solving failed!" << endl;
//        return;
//    }
//    std::cout << "#iterations:     " << cg.iterations() << std::endl;
//    m.block(1, 1, nx, ny) = Map<MatrixXf>(x.data(), nx, ny);
//}


void Simulator::diffuse(MatrixXf &m,
                        float k) {
    int nx = N[0] - 2, ny = N[1] - 2;
    int sz = nx * ny;
    float alpha = -dt * k / (dx * dx);
    float beta = 1 - 4 * alpha;
    // ****************
    // build A, b
    // ****************
    CG::cgMat A(nx, ny);
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
    //cgSolve(x, A, b, 1e-6f, 100, CG::MIC, &iters);
    cout << "diffuse phase: iters in CG::cgSolve:    " << iters << endl;
    m.block(1, 1, nx, ny) = Map<MatrixXf>(x.data(), nx, ny);
}


//void Simulator::project(VecMatXf &u1, const VecMatXf &u0) {
//    MatrixXf div(N[0], N[1]);
//    calcDiv(div, u0);
//    // ****************
//    // first step: solve for q: div(grad(q))=div(u0)
//    // here we solve for q/dx instead
//    // here MatrixXf div is actually dx*div(u0)
//    // ****************
//    // build A, b
//    // ****************
//    int sz = N[0] * N[1];
//    SparseMatrix<float> A(sz, sz);
//    VectorXf b = Map<VectorXf>(div.data(), sz);
//    b(0) = 0;
//    VectorXf x(sz);
//    vector <Trid> coeffs;
//    int i, j;
//    int idx = 0; // idx = i + j * N[0]
//    for (j = 0; j < N[1]; ++j) {
//        for (i = 0; i < N[0]; ++i) {
//            if (i == 0 && j == 0) {
//                coeffs.emplace_back(idx, idx, 1);
//            } else {
//                coeffs.emplace_back(Trid(idx, idx, -4));
//                if (i == 0) {
//                    coeffs.emplace_back(Trid(idx, idx + 1, 1));
//                } else {
//                    coeffs.emplace_back(Trid(idx, idx - 1, 1));
//                }
//                if (i == N[0] - 1) {
//                    coeffs.emplace_back(Trid(idx, idx - 1, 1));
//                } else {
//                    coeffs.emplace_back(Trid(idx, idx + 1, 1));
//                }
//                if (j == 0) {
//                    coeffs.emplace_back(Trid(idx, idx + N[0], 1));
//                } else {
//                    coeffs.emplace_back(Trid(idx, idx - N[0], 1));
//                }
//                if//void Simulator::project(VecMatXf &u1, const VecMatXf &u0) {
//    MatrixXf div(N[0], N[1]);
//    calcDiv(div, u0);
//    // ****************
//    // first step: solve for q: div(grad(q))=div(u0)
//    // here we solve for q/dx instead
//    // here MatrixXf div is actually dx*div(u0)
//    // ****************
//    // build A, b
//    // ****************
//    int sz = N[0] * N[1];
//    SparseMatrix<float> A(sz, sz);
//    VectorXf b = Map<VectorXf>(div.data(), sz);
//    b(0) = 0;
//    VectorXf x(sz);
//    vector <Trid> coeffs;
//    int i, j;
//    int idx = 0; // idx = i + j * N[0]
//    for (j = 0; j < N[1]; ++j) {
//        for (i = 0; i < N[0]; ++i) {
//            if (i == 0 && j == 0) {
//                coeffs.emplace_back(idx, idx, 1);
//            } else {
//                coeffs.emplace_back(Trid(idx, idx, -4));
//                if (i == 0) {
//                    coeffs.emplace_back(Trid(idx, idx + 1, 1));
//                } else {
//                    coeffs.emplace_back(Trid(idx, idx - 1, 1));
//                }
//                if (i == N[0] - 1) {
//                    coeffs.emplace_back(Trid(idx, idx - 1, 1));
//                } else {
//                    coeffs.emplace_back(Trid(idx, idx + 1, 1));
//                }
//                if (j == 0) {
//                    coeffs.emplace_back(Trid(idx, idx + N[0], 1));
//                } else {
//                    coeffs.emplace_back(Trid(idx, idx - N[0], 1));
//                }
//                if (j == N[1] - 1) {
//                    coeffs.emplace_back(Trid(idx, idx - N[0], 1));
//                } else {
//                    coeffs.emplace_back(Trid(idx, idx + N[0], 1));
//                }
//            }
//            ++idx;
//        }
//    }
//    A.setFromTriplets(coeffs.begin(), coeffs.end());
//    // ****************
//    // solve
//    // ****************
//    // todo: analyzePattern and factorize
//    // ConjugateGradient<SparseMatrix<float >> cg;
//    SparseLU <SparseMatrix<float>> cg;
//    cg.compute(A);
//    if (cg.info() != Success) {
//        cerr << "decomposition in \"project\" failed!" << endl;
//        return;
//    }
//    x = cg.solve(b);
//    if (cg.info() != Success) {
//        cerr << "solving in \"project\" failed!" << endl;
//        return;
//    }
//    div = Map<MatrixXf>(x.data(), N[0], N[1]);
//    // ****************
//    // from now, MatrixXf div stores the value for q/dx
//    // final step: obtain divergence-free u1
//    // ****************
//    for (j = 0; j < N[1]; ++j) {
//        for (i = 0; i < N[0]; ++i) {
//            if (i == 0 || i == N[0] - 1) {
//                u1[0](i, j) = 0;
//            } else {
//                u1[0](i, j) = u0[0](i, j) - 0.5 * (div(i + 1, j) - div(i - 1, j));
//            }
//            if (j == 0 || j == N[1] - 1) {
//                u1[1](i, j) = 0;
//            } else {
//                u1[1](i, j) = u0[1](i, j) - 0.5 * (div(i, j + 1) - div(i, j - 1));
//            }
//        }
//    }
//} (j == N[1] - 1) {
//                    coeffs.emplace_back(Trid(idx, idx - N[0], 1));
//                } else {
//                    coeffs.emplace_back(Trid(idx, idx + N[0], 1));
//                }
//            }
//            ++idx;
//        }
//    }
//    A.setFromTriplets(coeffs.begin(), coeffs.end());
//    // ****************
//    // solve
//    // ****************
//    // todo: analyzePattern and factorize
//    // ConjugateGradient<SparseMatrix<float >> cg;
//    SparseLU <SparseMatrix<float>> cg;
//    cg.compute(A);
//    if (cg.info() != Success) {
//        cerr << "decomposition in \"project\" failed!" << endl;
//        return;
//    }
//    x = cg.solve(b);
//    if (cg.info() != Success) {
//        cerr << "solving in \"project\" failed!" << endl;
//        return;
//    }
//    div = Map<MatrixXf>(x.data(), N[0], N[1]);
//    // ****************
//    // from now, MatrixXf div stores the value for q/dx
//    // final step: obtain divergence-free u1
//    // ****************
//    for (j = 0; j < N[1]; ++j) {
//        for (i = 0; i < N[0]; ++i) {
//            if (i == 0 || i == N[0] - 1) {
//                u1[0](i, j) = 0;
//            } else {
//                u1[0](i, j) = u0[0](i, j) - 0.5 * (div(i + 1, j) - div(i - 1, j));
//            }
//            if (j == 0 || j == N[1] - 1) {
//                u1[1](i, j) = 0;
//            } else {
//                u1[1](i, j) = u0[1](i, j) - 0.5 * (div(i, j + 1) - div(i, j - 1));
//            }
//        }
//    }
//}


//void Simulator::project(VecMatXf &u1, const VecMatXf &u0) {
//    MatrixXf div(N[0], N[1]);
//    calcDiv(div, u0);
//    // ****************
//    // first step: solve for q: div(grad(q))=div(u0)
//    // here we solve for q/dx instead
//    // here MatrixXf div is actually dx*div(u0)
//    // ****************
//    // build A, b
//    // ****************
//    int sz = N[0] * N[1];
//    SparseMatrix<float> A(sz, sz);
//    VectorXf b = Map<VectorXf>(div.data(), sz);
//    b = -b;
//    VectorXf x(sz);
//    vector<Trid> coeffs;
//    int i, j;
//    int idx = 0; // idx = i + j * N[0]
//    for (j = 0; j < N[1]; ++j) {
//        for (i = 0; i < N[0]; ++i) {
//            if (isZero(Blue0(i, j)) && isZero(Red0(i, j))) {
//                coeffs.emplace_back(Trid(idx, idx, 1));
//                b(idx) = 0;
//            } else {
//                coeffs.emplace_back(Trid(idx, idx, 4));
//                if (i == 0) {
//                    coeffs.emplace_back(Trid(idx, idx, -1));
//                } else {
//                    coeffs.emplace_back(Trid(idx, idx - 1, -1));
//                }
//                if (i == N[0] - 1) {
//                    coeffs.emplace_back(Trid(idx, idx, -1));
//                } else {
//                    coeffs.emplace_back(Trid(idx, idx + 1, -1));
//                }
//                if (j == 0) {
//                    coeffs.emplace_back(Trid(idx, idx, -1));
//                } else {
//                    coeffs.emplace_back(Trid(idx, idx - N[0], -1));
//                }
//                if (j == N[1] - 1) {
//                    coeffs.emplace_back(Trid(idx, idx, -1));
//                } else {
//                    coeffs.emplace_back(Trid(idx, idx + N[0], -1));
//                }
//            }
//            ++idx;
//        }
//    }
//    A.setFromTriplets(coeffs.begin(), coeffs.end());
//    // ****************
//    // solve
//    // ****************
//
//    // use eigen conjugate gradient -------------------------------
////    // todo: analyzePattern and factorize
////    ConjugateGradient<SparseMatrix<float >> cg;
////    cg.compute(A);
////    if (cg.info() != Success) {
////        cerr << "decomposition in \"project\" failed!" << endl;
////        return;
////    }
////    cg.setMaxIterations(200);
////    x = cg.solve(b);
////    cout << "# cg iterations:     " << cg.iterations() << endl;
////    if (cg.info() != Success) {
////        cerr << "solving in \"project\" failed!" << endl;
////        return;
////    }
//
//    // use sparse lu ----------------------
//    SparseLU<SparseMatrix<float >> lu;
//    lu.compute(A);
//    if (lu.info() != Success) {
//        cerr << "decomposition in \"project\" failed!" << endl;
//        return;
//    }
//    x = lu.solve(b);
//    if (lu.info() != Success) {
//        cerr << "solving in \"project\" failed!" << endl;
//        return;
//    }
//
//    div = Map<MatrixXf>(x.data(), N[0], N[1]);
//    // ****************
//    // from now, MatrixXf div stores the value for q/dx
//    // final step: obtain divergence-free u1
//    // ****************
//    for (j = 0; j < N[1]; ++j) {
//        for (i = 0; i < N[0]; ++i) {
//            if (i == 0 || i == N[0] - 1) {
//                u1[0](i, j) = 0;
//            } else {
//                u1[0](i, j) = u0[0](i, j) - 0.5 * (div(i + 1, j) - div(i - 1, j));
//            }
//            if (j == 0 || j == N[1] - 1) {
//                u1[1](i, j) = 0;
//            } else {
//                u1[1](i, j) = u0[1](i, j) - 0.5 * (div(i, j + 1) - div(i, j - 1));
//            }
//        }
//    }
//}


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

    CG::cgMat A(N[0], N[1]);
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
    //cgSolve(x, A, b, 0.1, 100, CG::Identity, &iters);
    cgSolve(x, A, b, 0.001, 100, CG::MIC, &iters);
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
    return isZero(Blue0(i, j)) && isZero(Red0(i, j));
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

