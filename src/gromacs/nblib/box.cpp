//
// Created by sebkelle on 19.11.19.
//

#include "box.h"


Box::Box(real l)
{
    box_[XX][XX] = l;
    box_[YY][YY] = l;
    box_[ZZ][ZZ] = l;
}

Box::Box(real x, real y, real z);
{
    box_[XX][XX] = x;
    box_[YY][YY] = y;
    box_[ZZ][ZZ] = z;
}

matrix Box::getMatrix() { return box_; }

