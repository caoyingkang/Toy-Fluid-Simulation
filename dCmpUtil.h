//
// Created by cyk on 18-12-24.
//

#ifndef PROJ_DCMPUTIL_H

const double EPS = 1e-8;

bool isZero(double a){
    return -EPS<a && a<EPS;
}

#define PROJ_DCMPUTIL_H

#endif //PROJ_DCMPUTIL_H
