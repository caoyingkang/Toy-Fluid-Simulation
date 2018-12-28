//
// Created by cyk on 18-12-23.
//

#ifndef PROJ_SIMULATOR_H
#define PROJ_SIMULATOR_H

#include <cstddef>
#include <initializer_list>
#include <vector>
#include <cassert>
#include <eigen3/Eigen/Dense>

const int C_BLUE = 1;
const int C_RED = 2;

using VecMatXf = std::vector<Eigen::MatrixXf>;

class Simulator {
public:
    Simulator(float _dx, float _dt,
              std::initializer_list<int> ln)
            : dx(_dx), dt(_dt), N(ln),
              U0(2, Eigen::MatrixXf::Zero(N[0], N[1])),
              U1(2, Eigen::MatrixXf::Zero(N[0], N[1])),
              F(2, Eigen::MatrixXf(N[0], N[1])),
              Blue0(Eigen::MatrixXf::Zero(N[0], N[1])),
              Blue1(Eigen::MatrixXf::Zero(N[0], N[1])),
              Red0(Eigen::MatrixXf::Zero(N[0], N[1])),
              Red1(Eigen::MatrixXf::Zero(N[0], N[1])),
              inlet(Eigen::MatrixXi::Zero(N[0], N[1])),
              force_setted(false),
              visc_setted(false),
              diff_setted(false) {
        assert(ln.size() == 2);
    }

    ~Simulator() = default;

    void setForce(const VecMatXf &f);

    void setForce(float fx, float fy);

    void setVisc(float val);

    void setDiff(float val);

    void setInlet(int radius_blue,
                  std::initializer_list<int> center_blue,
                  std::initializer_list<float> v_blue,
                  int radius_red,
                  std::initializer_list<int> center_red,
                  std::initializer_list<float> v_red);


    void Forward(); // perform one step

    bool isEmpty(int i, int j) const;

    void getRenderData(float *vertices);

    // printXXX: for debug usage
    void printBlue();

    void printRed();

    void printVx();

    void printVy();


private:
    float dx, dt;
    std::vector<int> N; // N[i]: No. of grids in i-th dimension
    float visc; // viscosity for fluid
    float diff; // diffusion coefficient for solids
    VecMatXf U0, U1; // velocity of girds

    // color of grids
    // 0: black background
    // 1: blue
    // 2: red
    Eigen::MatrixXf Blue0, Blue1;
    Eigen::MatrixXf Red0, Red1;

    // entrance for fluids
    // 0: not an entrance
    // 1: entrance for blue fluid
    // 2: entrance for red fluid
    Eigen::MatrixXi inlet;
    std::vector<float> v_in1, v_in2;

    VecMatXf F; // force

    // private indicators
    bool force_setted;
    bool visc_setted;
    bool diff_setted;

    // private methods
    void Vstep();

    void Sstep();

    void addForce(Eigen::MatrixXf &m,
                  const Eigen::MatrixXf &f);

    void advect(VecMatXf &vecm1,
                const VecMatXf &vecm0,
                const VecMatXf &u);

    void advect(Eigen::MatrixXf &m1,
                const Eigen::MatrixXf &m0,
                const VecMatXf &u,
                int color);

    void diffuse(Eigen::MatrixXf &m,
                 float k);

    void project(VecMatXf &u1, const VecMatXf &u0);

    void TraceParticle(const VecMatXf &u,
                       int curr_i, int curr_j,
                       int &pre_i, int &pre_j,
                       float &frac_x, float &frac_y);

    float LinInterp(const Eigen::MatrixXf &m,
                    int i, int j,
                    float frac, bool Along_X);

    float biLinInterp(const Eigen::MatrixXf &m,
                      int i, int j,
                      float frac_x, float frac_y);

    void calcDiv(Eigen::MatrixXf &m,
                 const VecMatXf &u);
};

#endif //PROJ_SIMULATOR_H
