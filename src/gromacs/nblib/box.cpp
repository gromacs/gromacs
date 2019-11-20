//
// Created by sebkelle on 19.11.19.
//

#include "gromacs/utility/exceptions.h"

#include "box.h"

namespace nblib {

Box::Box(real l)
{
    if (std::isnan(l) or std::isinf(l))
    {
        GMX_THROW(gmx::InvalidInputError("Cannot have NaN or Inf box length."));
    }
            
    box_[XX][XX] = l;
    box_[YY][YY] = l;
    box_[ZZ][ZZ] = l;
}

Box::Box(real x, real y, real z)
{
    if (std::isnan(x) or std::isinf(x) or
        std::isnan(y) or std::isinf(y) or
        std::isnan(z) or std::isinf(z))
    {
        GMX_THROW(gmx::InvalidInputError("Cannot have NaN or Inf box length."));
    }

    box_[XX][XX] = x;
    box_[YY][YY] = y;
    box_[ZZ][ZZ] = z;
}

Box::data_type Box::getMatrix()
{
    return box_;
}

}
