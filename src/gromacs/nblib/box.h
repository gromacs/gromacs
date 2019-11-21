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
public:
    using data_type = std::array<std::array<real, DIM>, DIM>;

    Box(real l);

    Box(real x, real y, real z);

    data_type matrix();

private:
    data_type box_;
};

} // namespace nblib

#endif //GROMACS_BOX_H
