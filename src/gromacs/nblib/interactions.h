//
// Created by sebkelle on 20.11.19.
//

#ifndef GROMACS_INTERACTIONS_H
#define GROMACS_INTERACTIONS_H

#include "gromacs/math/vectypes.h"

//! Type of interaction
struct HarmonicType
{
    real equiDist;
    real forceConstant;
    std::string atomIDi;
    std::string atomIDj;
    //int atomI;
    //int atomJ;
};


#endif //GROMACS_INTERACTIONS_H
