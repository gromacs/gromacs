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

#include "densitypadding.h"
#include "gromacs/math/griddata/griddata.h"
#include <cmath>
#include <array>
#include "gromacs/utility/exceptions.h"

namespace gmx
{

std::unique_ptr < GridDataReal3D> DensityPadding::padPower2()
{
    const auto &extend   = toPad_.getGrid().lattice().extend();
    real        factorXX = pow(2, ceil(log(extend[XX])/log(2)));
    real        factorYY = pow(2, ceil(log(extend[YY])/log(2)));
    real        factorZZ = pow(2, ceil(log(extend[ZZ])/log(2)));

    return pad({{factorXX/real(extend[XX]), factorYY/real(extend[YY]), factorZZ/real(extend[ZZ])}} );
}

std::unique_ptr < GridDataReal3D>
DensityPadding::pad(const CanonicalVectorBasis<DIM>::NdVector &paddingFactor)
{
    std::array<int, 3> paddedExtend;
    for (int dimension = 0; dimension < DIM; dimension++)
    {
        if (paddingFactor[dimension] > 1)
        {
            paddedExtend[dimension] = std::ceil(paddingFactor[dimension] * toPad_.getGrid().lattice().extend()[dimension]);
        }
        else
        {
            GMX_THROW(RangeError("Density padding only viable with padding factor >1 , here : "+ std::to_string(paddingFactor[dimension])+ " ."));
        }

    }
    auto paddedGrid = Grid<DIM>(toPad_.getGrid().cell(), toPad_.getGrid().lattice());

    paddedGrid.setLatticeAndRescaleCell(paddedExtend);
    std::unique_ptr < GridDataReal3D> padded(new GridDataReal3D(std::unique_ptr < IGrid < 3>>(new Grid<DIM>(paddedGrid))));
    std::fill(std::begin(*padded), std::end(*padded), 0.);

    auto extend = toPad_.getGrid().lattice().extend();
    for (int iZZ = 0; iZZ < extend[ZZ]; ++iZZ)
    {
        for (int iYY = 0; iYY < extend[YY]; ++iYY)
        {
            auto rowBegin = toPad_.iteratorAtMultiIndex({{0, iYY, iZZ}});
            auto rowEnd   = ++(toPad_.iteratorAtMultiIndex({{extend[XX]-1, iYY, iZZ}}));
            std::copy(rowBegin, rowEnd, padded->iteratorAtMultiIndex({{0, iYY, iZZ}}));
        }
    }
    return padded;
}
DensityPadding::DensityPadding(const GridDataReal3D &toPad) : toPad_ {toPad}
{}

std::unique_ptr < GridDataReal3D>
DensityPadding::unpad(const  ColumnMajorLattice<DIM>::MultiIndex &unPadExtend)
{

    auto unpaddedGrid = Grid<3>(toPad_.getGrid().cell(), toPad_.getGrid().lattice());
    unpaddedGrid.setLatticeAndRescaleCell(unPadExtend);
    std::unique_ptr < GridDataReal3D> unpadded(new GridDataReal3D(std::unique_ptr < IGrid < 3>>(new Grid<DIM>(unpaddedGrid))));

    for (int iZZ = 0; iZZ < unPadExtend[ZZ]; ++iZZ)
    {
        for (int iYY = 0; iYY < unPadExtend[YY]; ++iYY)
        {
            auto rowBegin = toPad_.iteratorAtMultiIndex({{0, iYY, iZZ}});
            auto rowEnd   = ++(toPad_.iteratorAtMultiIndex({{unPadExtend[XX]-1, iYY, iZZ}}));
            std::copy(rowBegin, rowEnd, unpadded->iteratorAtMultiIndex({{0, iYY, iZZ}}));
        }
    }

    return unpadded;
}
}
