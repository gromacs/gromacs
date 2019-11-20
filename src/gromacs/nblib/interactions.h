//
// Created by sebkelle on 20.11.19.
//

#ifndef GROMACS_INTERACTIONS_H
#define GROMACS_INTERACTIONS_H


//! Type of interaction
struct HarmonicType
{
    real equiDist;
    real forceConstant;
    int atomI;
    int atomJ;
};


#endif //GROMACS_INTERACTIONS_H
