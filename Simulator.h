//
// Created by cyk on 18-12-23.
//

#ifndef PROJ_SIMULATOR_H
#define PROJ_SIMULATOR_H

#include <cstddef>
#include <initializer_list>
#include <vector>
#include <cassert>
#include <eigen3/Eigen/Eigen>

const int C_BLUE = 1;
const int C_RED = 2;

class Simulator {
public:
    Simulator(double _dx, double _dt,
              std::initializer_list<int> ln,
              double _visc)
            : dx(_dx), dt(_dt), N(ln), visc(_visc),
              U0(2, Eigen::MatrixXd::Zero(N[0], N[1])),
              U1(2, Eigen::MatrixXd::Zero(N[0], N[1])),
              F(2, Eigen::MatrixXd(N[0], N[1])),
              C0(Eigen::MatrixXi::Zero(N[0], N[1])),
              C1(Eigen::MatrixXi::Zero(N[0], N[1])),
              Entrance(N[0], N[1]){
        assert(ln.size() == 2);
    }

    ~Simulator() = default;

    void setForce(const std::vector<Eigen::MatrixXd> &newF) {
        assert(newF.size() == 2 && newF[0].size() == F[0].size());
        F[0] = newF[0];
        F[1] = newF[1];
    }

    void setForce(double fx, double fy) {
        F[0] = Eigen::MatrixXd::Constant(N[0], N[1], fx);
        F[1] = Eigen::MatrixXd::Constant(N[0], N[1], fy);
    }

private:
    double dx, dt;
    std::vector<int> N; // N[i]: No. of grids in i-th dimension
    double visc;
    std::vector<Eigen::MatrixXd> U0, U1; // velocity of girds
    Eigen::MatrixXi C0, C1; // color of grids
    Eigen::Matrix<bool,Eigen::Dynamic,Eigen::Dynamic> Entrance;
    std::vector<Eigen::MatrixXd> F; // force

    void Vstep();

    void Sstep();

    void addForce(Eigen::MatrixXd &m,
                  const Eigen::MatrixXd &f);

    void advect(std::vector<Eigen::MatrixXd> &vecm1,
                const std::vector<Eigen::MatrixXd> &vecm0,
                const std::vector<Eigen::MatrixXd> &u);

    void advect(Eigen::MatrixXd &m1,
                const Eigen::MatrixXd &m0,
                const std::vector<Eigen::MatrixXd> &u);

    void diffuse(Eigen::MatrixXd &m1,
                 const Eigen::MatrixXd &m0,
                 double k);

    void project();

    void TraceParticle(const std::vector<Eigen::MatrixXd> &u,
                       int curr_i, int curr_j,
                       int &pre_i, int &pre_j,
                       double &frac_x, double &frac_y);

    double LinInterp(const Eigen::MatrixXd &m,
                     int i, int j,
                     double frac, bool Along_X);

    double biLinInterp(const Eigen::MatrixXd &m,
                       int i,int j,
                        double frac_x, double frac_y);
};

#endif //PROJ_SIMULATOR_H
