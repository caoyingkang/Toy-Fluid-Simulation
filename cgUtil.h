//
// Created by cyk on 18-12-27.
//

#ifndef PROJ_CGUTIL_H
#define PROJ_CGUTIL_H

#include <eigen3/Eigen/Dense>
#include <memory>

namespace CG {
    enum PrecondType {
        Identity, MIC
    };
}

class cgMat {
public:
    int Nx, Ny;
    Eigen::MatrixXf Adiag, Aplusi, Aplusj;

    cgMat(int _Nx, int _Ny) : Nx(_Nx), Ny(_Ny),
                              Adiag(Nx, Ny),
                              Aplusi(Nx - 1, Ny),
                              Aplusj(Nx, Ny - 1) {
        Adiag.setZero();
        Aplusi.setZero();
        Aplusj.setZero();
    };

    ~cgMat() = default;

    int size() const { return Nx * Ny; }

    // z = A * x
    void multiply(Eigen::VectorXf &z, const Eigen::VectorXf &x) const;

    void ApplyPrecond(Eigen::VectorXf &z,
                      const Eigen::VectorXf &x,
                      CG::PrecondType type);

private:
    std::shared_ptr<Eigen::MatrixXf> MIC_cond;

    // calculate MIC_cond
    // this function is called the first time ApplyPrecond using MIC is called
    void MIC_factorize();
};

void cgSolve(Eigen::VectorXf &x,
             cgMat &A,
             const Eigen::VectorXf &b,
             float TOL,
             int MAXIT,
             CG::PrecondType type = CG::Identity,
             int *iterations = nullptr);

#endif //PROJ_CGUTIL_H
