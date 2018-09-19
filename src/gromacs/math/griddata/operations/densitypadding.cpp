/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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

#include "densitypadding.h"

#include <cmath>

#include <array>

#include "gromacs/math/griddata/griddata.h"
#include "gromacs/utility/exceptions.h"

namespace gmx
{

std::unique_ptr<GridDataFloat3D>
padDensityToNextPower2(const GridDataFloat3D &toPad)
{
    const auto &extend   = toPad.getGrid().lattice();
    real        factorXX = pow(2, ceil(log(extend[XX]) / log(2)));
    real        factorYY = pow(2, ceil(log(extend[YY]) / log(2)));
    real        factorZZ = pow(2, ceil(log(extend[ZZ]) / log(2)));

    return padDensity(toPad,
                      {{factorXX / real(extend[XX]), factorYY / real(extend[YY]),
                        factorZZ / real(extend[ZZ])}});
}

std::unique_ptr<GridDataFloat3D>
padDensity(const GridDataFloat3D &toPad, const MdFloatVector<DIM> &paddingFactor)
{
    bounds<DIM> paddedExtend;
    for (int dimension = 0; dimension < DIM; dimension++)
    {
        if (paddingFactor[dimension] > 1)
        {
            paddedExtend[dimension] =
                std::ceil(paddingFactor[dimension] *
                          toPad.getGrid().lattice()[dimension]);
        }
        else
        {
            GMX_THROW(RangeError(
                              "Density padding only viable with padding factor >1 , here : " +
                              std::to_string(paddingFactor[dimension]) + " ."));
        }
    }
    auto paddedGrid = toPad.getGrid();

    paddedGrid.setLatticeAndRescaleCell(paddedExtend);
    std::unique_ptr<GridDataFloat3D> padded(new GridDataFloat3D(paddedGrid));
    std::fill(std::begin(*padded), std::end(*padded), 0.);

    auto extend = toPad.getGrid().lattice();
    for (int iZZ = 0; iZZ < extend[ZZ]; ++iZZ)
    {
        for (int iYY = 0; iYY < extend[YY]; ++iYY)
        {
            auto rowBegin = toPad.begin() + toPad.memoryOffset({0, iYY, iZZ});
            auto rowEnd   = toPad.begin() + toPad.memoryOffset({extend[XX] - 1, iYY, iZZ}) + 1;
            std::copy(rowBegin, rowEnd, padded->begin() + padded->memoryOffset({0, iYY, iZZ}));
        }
    }
    return padded;
}

std::unique_ptr<GridDataFloat3D>
unpadDensity(const GridDataFloat3D &toPad, const bounds<DIM> &unPadExtend)
{
    auto unpaddedGrid = toPad.getGrid();
    unpaddedGrid.setLatticeAndRescaleCell(unPadExtend);
    std::unique_ptr<GridDataFloat3D> unpadded(new GridDataFloat3D(unpaddedGrid));

    for (int iZZ = 0; iZZ < unPadExtend[ZZ]; ++iZZ)
    {
        for (int iYY = 0; iYY < unPadExtend[YY]; ++iYY)
        {
            auto rowBegin = toPad.begin() + toPad.memoryOffset({0, iYY, iZZ});
            auto rowEnd   = toPad.begin() + toPad.memoryOffset({unPadExtend[XX] - 1, iYY, iZZ}) + 1;
            std::copy(rowBegin, rowEnd, unpadded->begin()+unpadded->memoryOffset({0, iYY, iZZ}));
        }
    }

    return unpadded;
}
}
