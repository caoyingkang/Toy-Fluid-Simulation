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

using VecMatXd = std::vector<Eigen::MatrixXd>;

class Simulator {
public:
    Simulator(double _dx, double _dt,
              std::initializer_list<int> ln,
              double _visc, double _diff)
            : dx(_dx), dt(_dt), N(ln),
              visc(_visc), diff(_diff),
              U0(2, Eigen::MatrixXd::Zero(N[0], N[1])),
              U1(2, Eigen::MatrixXd::Zero(N[0], N[1])),
              F(2, Eigen::MatrixXd(N[0], N[1])),
              Blue0(Eigen::MatrixXd::Zero(N[0], N[1])),
              Blue1(Eigen::MatrixXd::Zero(N[0], N[1])),
              Red0(Eigen::MatrixXd::Zero(N[0], N[1])),
              Red1(Eigen::MatrixXd::Zero(N[0], N[1])),
              inlet(Eigen::MatrixXi::Zero(N[0], N[1])),
              use_default_inlet(false),
              force_setted(false) {
        assert(ln.size() == 2);
    }

    ~Simulator() = default;

    void setForce(const VecMatXd &f);

    void setForce(double fx, double fy);

    void setInlet(int radius_blue,
                  std::initializer_list<int> center_blue,
                  std::initializer_list<double> v_blue,
                  int radius_red,
                  std::initializer_list<int> center_red,
                  std::initializer_list<double> v_red);

    void setDefaultInLet();

    void Forward(); // perform one step

    void printBlue();

    void printRed();

    void printVx();

    void printVy();

private:
    double dx, dt;
    std::vector<int> N; // N[i]: No. of grids in i-th dimension
    double visc; // viscosity for fluid
    double diff; // diffusion coefficient for solids
    VecMatXd U0, U1; // velocity of girds

    // color of grids
    // 0: black background
    // 1: blue
    // 2: red
    Eigen::MatrixXd Blue0, Blue1;
    Eigen::MatrixXd Red0, Red1;

    // entrance for fluids
    // 0: not an entrance
    // 1: entrance for blue fluid
    // 2: entrance for red fluid
    Eigen::MatrixXi inlet;
    std::vector<double> v_in1, v_in2;

    VecMatXd F; // force

    // private indicators
    bool use_default_inlet;
    bool force_setted;

    // private methods
    void Vstep();

    void Sstep();

    void addForce(Eigen::MatrixXd &m,
                  const Eigen::MatrixXd &f);

    void advect(VecMatXd &vecm1,
                const VecMatXd &vecm0,
                const VecMatXd &u);

    void advect(Eigen::MatrixXd &m1,
                const Eigen::MatrixXd &m0,
                const VecMatXd &u,
                int color);

    void diffuse(Eigen::MatrixXd &m,
                 double k);

    void project(VecMatXd &u1, const VecMatXd &u0);

    void TraceParticle(const VecMatXd &u,
                       int curr_i, int curr_j,
                       int &pre_i, int &pre_j,
                       double &frac_x, double &frac_y);

    double LinInterp(const Eigen::MatrixXd &m,
                     int i, int j,
                     double frac, bool Along_X);

    double biLinInterp(const Eigen::MatrixXd &m,
                       int i, int j,
                       double frac_x, double frac_y);

    void calcDiv(Eigen::MatrixXd &m,
                 const VecMatXd &u);
};

#endif //PROJ_SIMULATOR_H
