//
// Created by sebkelle on 19.11.19.
//

#ifndef GROMACS_BOX_H
#define GROMACS_BOX_H

#include <array>
#include "gromacs/math/vectypes.h"

namespace nblib
{

class Box {
    using data_type = std::array<std::array<real, DIM>, DIM>;
public:
    Box(real l);

    Box(real x, real y, real z);

    data_type getMatrix();

private:
    data_type box_;
};

}

#endif //GROMACS_BOX_H
