//
// Created by sebkelle on 19.11.19.
//

#include "gromacs/utility/exceptions.h"

#include "box.h"

namespace nblib {

Box::Box(real l) {
    if (std::isnan(l) or std::isinf(l))
    {
        GMX_THROW(gmx::InvalidInputError("Cannot have NaN or Inf box length."));
    }
            
    box_[XX] = l;
    box_[YY*DIM + ZZ] = l;
    box_[ZZ*DIM + ZZ] = l;
}

    Box::Box(real x, real y, real z)
{
    box_[XX] = x;
    box_[YY*DIM + ZZ] = y;
    box_[ZZ*DIM + ZZ] = z;
}

std::array<real, DIM*DIM> Box::getMatrix() { return box_; }

}
