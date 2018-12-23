//
// Created by cyk on 18-12-23.
//

#include "Simulator.h"
#include "dCmpUtil.h"

using namespace Eigen;
using namespace std;

void Simulator::Vstep() {
    for (int i = 0; i < 2; ++i)
        addForce(U0[i], F[i]);
    advect(U1, U0, U0);
    for (int i = 0; i < 2; ++i)
        diffuse(U0[i],U1[i],visc);
    project();
}

void Simulator::addForce(Eigen::MatrixXd &m,
                         const Eigen::MatrixXd &f) {
    m += dt * f;
}

void Simulator::advect(std::vector<Eigen::MatrixXd> &vecm1,
                       const std::vector<Eigen::MatrixXd> &vecm0,
                       const std::vector<Eigen::MatrixXd> &u) {
    int curr_i, curr_j;
    int pre_i,pre_j;
    double frac_x,frac_y;
    for (curr_i = 0; curr_i < N[0]; ++curr_i) {
        for (curr_j = 0; curr_j < N[1]; ++curr_j) {
            TraceParticle(u,curr_i,curr_j,pre_i,pre_j,frac_x,frac_y);
            if(isZero(frac_x)&&isZero(frac_y)){
                for(int i=0;i<2;++i)
                    vecm1[i](curr_i,curr_j)=vecm0[i](pre_i,pre_j);
            }else if(isZero(frac_x)){
                for(int i=0;i<2;++i)
                    vecm1[i](curr_i,curr_j)=LinInterp(vecm0[i],pre_i,pre_j,frac_y, false);
            }else if(isZero(frac_y)){
                for(int i=0;i<2;++i)
                    vecm1[i](curr_i,curr_j)=LinInterp(vecm0[i],pre_i,pre_j,frac_x, true);
            }else{
                for(int i=0;i<2;++i)
                    vecm1[i](curr_i,curr_j)=biLinInterp(vecm0[i],pre_i,pre_j,frac_x,frac_y);
            }
        }
    }
}

void Simulator::advect(Eigen::MatrixXd &m1,
                       const Eigen::MatrixXd &m0,
                       const std::vector<Eigen::MatrixXd> &u) {
    int curr_i, curr_j;
    int pre_i,pre_j;
    double frac_x,frac_y;
    for (curr_i = 0; curr_i < N[0]; ++curr_i) {
        for (curr_j = 0; curr_j < N[1]; ++curr_j) {
            TraceParticle(u,curr_i,curr_j,pre_i,pre_j,frac_x,frac_y);
            if(isZero(frac_x)&&isZero(frac_y)){
                for(int i=0;i<2;++i)
                    m1(curr_i,curr_j)=m0(pre_i,pre_j);
            }else if(isZero(frac_x)){
                for(int i=0;i<2;++i)
                    m1(curr_i,curr_j)=LinInterp(m0,pre_i,pre_j,frac_y, false);
            }else if(isZero(frac_y)){
                for(int i=0;i<2;++i)
                    m1(curr_i,curr_j)=LinInterp(m0,pre_i,pre_j,frac_x, true);
            }else{
                for(int i=0;i<2;++i)
                    m1(curr_i,curr_j)=biLinInterp(m0,pre_i,pre_j,frac_x,frac_y);
            }
        }
    }
}

void Simulator::TraceParticle(const std::vector<Eigen::MatrixXd> &u,
                              int curr_i, int curr_j,
                              int &pre_i, int &pre_j,
                              double &frac_x, double &frac_y) {
    double Di = - dt * u[0](curr_i,curr_j)/dx;
    double Dj = - dt * u[1](curr_i,curr_j)/dx;
    int ni,nj;
    if(isZero(Di)){ // zero x-velocity at this point
        pre_i=curr_i;
        frac_x=0;
    }else{
        ni = int(Di) - (Di>0? 0: 1);
        pre_i = curr_i + ni;
        if(pre_i>=0 && pre_i<=N[0]-2){
            frac_x=Di-ni;
        }else{ // out-of-domain in x-direction
            pre_i = pre_i<0? 0:(N[0]-1);
            frac_x=0;
        }
    }
    if(isZero(Dj)){ // zero y-velocity at this point
        pre_j=curr_j;
        frac_y=0;
    }else{
        nj = int(Dj) - (Dj>0? 0: 1);
        pre_j = curr_j + nj;
        if(pre_j>=0 && pre_j<=N[1]-2){
            frac_y=Dj-nj;
        }else{ // out-of-domain in y-direction
            pre_j = pre_j<0? 0:(N[1]-1);
            frac_y=0;
        }
    }
    assert(-EPS<frac_x<1+EPS && -EPS<frac_y<1+EPS);
    assert(pre_i>=0 && pre_i<=N[0]-2);
    assert(pre_j>=0 && pre_j<=N[1]-2);
}


double Simulator::LinInterp(const Eigen::MatrixXd &m,
                            int i, int j,
                            double frac, bool Along_X){
    assert(Along_X? (i<N[0]-1):(j<N[1]-1));
    if(Along_X){
        return (1-frac)*m(i,j)+frac*m(i+1,j);
    }else{
        return (1-frac)*m(i,j)+frac*m(i,j+1);
    }
}

double Simulator::biLinInterp(const Eigen::MatrixXd &m,
        int i,int j,
        double frac_x, double frac_y) {
    // bi-linear interpolation between grid points (i,j), (i+1,j), (i,j+1), (i+1,j+1)
    // the point to be estimated is (i+frac_x,j+frac_y)
    assert((i<N[0]-1)&&(j<N[1]-1));
    double tmp1 = (1-frac_x)*m(i,j)+frac_x*m(i+1,j);
    double tmp2 = (1-frac_x)*m(i,j+1)+frac_x*m(i+1,j+1);
    return (1-frac_y) * tmp1 + frac_y * tmp2;
}

void Simulator::diffuse(Eigen::MatrixXd &m1,
                        const Eigen::MatrixXd &m0,
                        double k) {

}

void Simulator::project(){

}

void Simulator::Sstep() {
}