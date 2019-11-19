//
// Created by sebkelle on 19.11.19.
//

#ifndef GROMACS_BOX_H
#define GROMACS_BOX_H

#include <array>
#include "gromacs/math/vectypes.h"

namespace nblib {

class Box {
public:
    Box(real l);

    Box(real x, real y, real z);

    std::array<real, DIM*DIM> getMatrix();

private:
    std::array<real, DIM*DIM> box_;
};

}

#endif //GROMACS_BOX_H
