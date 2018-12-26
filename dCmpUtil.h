//
// Created by cyk on 18-12-24.
//

#ifndef PROJ_DCMPUTIL_H
#define PROJ_DCMPUTIL_H

const float EPS = 1e-6;

bool isZero(float a) {
    return -EPS < a && a < EPS;
}


#endif //PROJ_DCMPUTIL_H
