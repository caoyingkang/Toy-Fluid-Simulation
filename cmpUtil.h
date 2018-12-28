//
// Created by cyk on 18-12-29.
//

#ifndef PROJ_CMPUTIL_H
#define PROJ_CMPUTIL_H

namespace CMP {
    const float EPS = 1e-15;
}

inline bool isZero(float a) {
    return -CMP::EPS < a && a < CMP::EPS;
}

inline bool Less(float a, float b) {
    return a < b - CMP::EPS;
}

inline bool Greater(float a, float b) {
    return a > b + CMP::EPS;
}

#endif //PROJ_CMPUTIL_H
