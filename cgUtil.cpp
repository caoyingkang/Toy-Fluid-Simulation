//
// Created by cyk on 18-12-27.
//

#include "cgUtil.h"
#include <cassert>
#include <iostream>
#include <cmath>
#include <ctime>

using namespace std;
using namespace Eigen;
using namespace CG;


// z = A * x
void cgMat::multiply(Eigen::VectorXf &z, const Eigen::VectorXf &x) const {
    assert(z.size() == Nx * Ny);
    int i, j;
    int idx = 0; // idx = i + j * Nx
    float val;
    for (j = 0; j < Ny; ++j) {
        for (i = 0; i < Nx; ++i) {
            val = Adiag(i, j) * x(idx);
            if (i > 0) val += Aplusi(i - 1, j) * x(idx - 1);
            if (i < Nx - 1) val += Aplusi(i, j) * x(idx + 1);
            if (j > 0) val += Aplusj(i, j - 1) * x(idx - Nx);
            if (j < Ny - 1) val += Aplusj(i, j) * x(idx + Nx);
            z(idx) = val;
            ++idx;
        }
    }
};

void cgMat::MIC_factorize() {
    assert(!MIC_cond);
    MIC_cond = make_shared<MatrixXf>(Nx, Ny);
    MatrixXf &precon = *MIC_cond;
    float tau = 0.97;
    int i, j;
    float e, tmp;
    for (j = 0; j < Ny; ++j) {
        for (i = 0; i < Nx; ++i) {
            e = Adiag(i, j);
            if (i > 0) {
                tmp = Aplusi(i - 1, j) * precon(i - 1, j);
                e -= tmp * tmp;
                if (j < Ny - 1)
                    e -= tau * tmp * Aplusj(i - 1, j) * precon(i - 1, j);
            }
            if (j > 0) {
                tmp = Aplusj(i, j - 1) * precon(i, j - 1);
                e -= tmp * tmp;
                if (i < Nx - 1)
                    e -= tau * tmp * Aplusi(i, j - 1) * precon(i, j - 1);
            }
            if (e < -1e-30) {
                cerr << "ERROR::CG::MIC_factorize::sqrt_negative   perhaps this is not a PSD matrix?" << endl;
            }
            precon(i, j) = 1 / sqrt(e + 1e-30);
        }
    }
}


void cgMat::ApplyPrecond(VectorXf &z,
                         const VectorXf &x,
                         PrecondType type) {
    if (type == Identity) {
        z = x;
    } else if (type == MIC) {
        assert(x.size() == Nx * Ny);
        assert(z.size() == Nx * Ny);
        if (!MIC_cond) {
            MIC_factorize();
        }
        // Applying the MIC(0) preconditioner
        // ----------------------------------
        MatrixXf &precon = *MIC_cond;
        // first solve Ly=x
        int i, j;
        int idx = 0; // idx = i + j * Nx
        float t;
        for (j = 0; j < Ny; ++j) {
            for (i = 0; i < Nx; ++i) {
                t = x(idx);
                if (i > 0) t -= Aplusi(i - 1, j) * precon(i - 1, j) * z(idx - 1);
                if (j > 0) t -= Aplusj(i, j - 1) * precon(i, j - 1) * z(idx - Nx);
                z(idx) = t * precon(i, j);
                ++idx;
            }
        }
        // next solve Lz=y
        idx = Nx * Ny - 1;
        for (j = Ny - 1; j >= 0; --j) {
            for (i = Nx - 1; i >= 0; --i) {
                t = z(idx);
                if (i < Nx - 1) t -= Aplusi(i, j) * precon(i, j) * z(idx + 1);
                if (j < Ny - 1) t -= Aplusj(i, j) * precon(i, j) * z(idx + Nx);
                z(idx) = t * precon(i, j);
                --idx;
            }
        }
        // ----------------------------------
    } else {
        cerr << "ERROR::cgApplyPrecond::Invalid Preconditioner Type!" << endl;
    }
}

void cgSolve(VectorXf &x,
             cgMat &A,
             const VectorXf &b,
             float TOL,
             int MAXIT,
             PrecondType type,
             int *iterations) {
    assert(x.size() == A.size());
    assert(b.size() == A.size());
    MAXIT = min(MAXIT, int(A.size()));
    x.setZero();
    VectorXf r = b;
    if (r.lpNorm<Infinity>() < TOL) {
        *iterations = 0;
        return;
    }
    VectorXf aux(A.size());
    A.ApplyPrecond(aux, r, type);
    VectorXf s = aux;
    float tmp = r.dot(aux), tmp_new, sAs;
    int iter = 0;
    while (iter < MAXIT) {
        ++iter;
        A.multiply(aux, s);
        sAs = s.dot(aux);
        if (sAs < -1e-15f) {
            cerr << "ERROR::CG::cgSolve::sAs is negative, perhaps this is not a PSD matrix?" << endl;
        }
        float alpha = tmp / ((sAs < 1e-15f && sAs > -1e-15f) ? 1e-15f : sAs);
        x += alpha * s;
        r -= alpha * aux;
        if (r.lpNorm<Infinity>() < TOL)
            break;
        A.ApplyPrecond(aux, r, type);
        tmp_new = r.dot(aux);
        if (tmp < -1e-15f) {
            cerr << "ERROR::CG::cgSolve::rMr is negative, perhaps this is not a PSD matrix?" << endl;
        }
        float beta = tmp_new / ((tmp < 1e-15f && tmp > -1e-15f) ? 1e-15f : tmp);
        s = aux + beta * s;
        tmp = tmp_new;
    }
    *iterations = iter;
}