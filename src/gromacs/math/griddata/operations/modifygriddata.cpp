/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/********************************************************************
 * GridDataReal3D
 */
#include "modifygriddata.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/griddata/operations/realfieldmeasure.h"

namespace gmx
{

ModifyGridData::ModifyGridData(std::vector<real> &grid)
    : grid_ {grid}
{ };

ModifyGridData &ModifyGridData::multiply(real value)
{
    std::for_each(grid_.begin(), grid_.end(), [value](real &v) { v *= value; });
    return *this;
}

real ModifyGridData::normalize()
{
    real integratedDensity =  DataVectorMeasure(grid_).mean();
    multiply(1/integratedDensity);
    return integratedDensity;
}

ModifyGridData &ModifyGridData::add_offset(real value)
{
    std::for_each(grid_.begin(), grid_.end(), [value](real &datum) { datum += value; });
    return *this;
}

ModifyGridData &ModifyGridData::zero()
{
    std::fill(grid_.begin(), grid_.end(), 0);
    return *this;
};

}
