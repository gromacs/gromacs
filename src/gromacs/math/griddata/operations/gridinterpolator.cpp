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
#include "gridinterpolator.h"

#include "gromacs/math/vec.h"

namespace gmx
{

/*
 * for each target grid point:
 *      find rational number grid cell index in input grid
 *      use fractional part for weights
 */
GridDataFloat3D
interpolateLinearly(const GridDataFloat3D &other, const Grid < DIM, GridWithTranslation < DIM>> &targetGrid)
{
    GridDataFloat3D interpolatedData(targetGrid);

    for (auto index : targetGrid.lattice())
    {
        auto r = targetGrid.multiIndexToCoordinate(index);
        interpolatedData[index] = getLinearInterpolationAt(other, r);
    }

    return interpolatedData;
}


/*
 * for each target grid point:
 *      find rational number grid cell index in input grid
 *      use fractional part for weights
 */
void
interpolateLinearly(const GridDataFloat3D &other, GridDataFloat3D *targetGrid)
{
    for (auto index : targetGrid->getGrid().lattice())
    {
        auto r = targetGrid->getGrid().multiIndexToCoordinate(index);
        (*targetGrid)[index] = getLinearInterpolationAt(other, r);
    }
}

real getLinearInterpolationAt(const GridDataFloat3D &field, const RVec &r)
{
    auto iIndexInGrid = field.getGrid().coordinateToFloorMultiIndex({{r[0], r[1], r[2]}});
    auto w            = field.getGrid().gridVectorFromGridPointToCoordinate({{r[0], r[1], r[2]}}, iIndexInGrid);

    std::array<std::array<std::array<real, 2>, 2>, 2> cube;
    const auto &lattice = field.getGrid().lattice();
    for (int ii_z = 0; ii_z <= 1; ++ii_z)
    {
        for (int ii_y = 0; ii_y <= 1; ++ii_y)
        {
            for (int ii_x = 0; ii_x <= 1; ++ii_x)
            {
                offset<DIM> cube_index = {iIndexInGrid[XX]+ii_x, iIndexInGrid[YY]+ii_y, iIndexInGrid[ZZ]+ii_z};
                if (lattice.contains(cube_index))
                {
                    cube[ii_x][ii_y][ii_z] = field[cube_index];
                }
                else
                {
                    cube[ii_x][ii_y][ii_z] = 0;
                }
            }
        }
    }

    std::array<std::array<real, 2>, 2> interpolated_x;
    for (int ii_z = 0; ii_z <= 1; ++ii_z)
    {
        for (int ii_y = 0; ii_y <= 1; ++ii_y)
        {
            interpolated_x[ii_y][ii_z] =
                (1 - w[XX]) * cube[0][ii_y][ii_z] + (w[XX])*cube[1][ii_y][ii_z];
        }
    }

    std::array<real, 2> interpolated_xy;
    for (int ii_z = 0; ii_z <= 1; ++ii_z)
    {
        interpolated_xy[ii_z] = (1 - w[YY]) * interpolated_x[0][ii_z] +
            (w[YY])*interpolated_x[1][ii_z];
    }

    return ((1 - w[ZZ]) * interpolated_xy[0] + w[ZZ] * interpolated_xy[1]) / 8.0;
}

real getLinearInterpolationAt(const GridDataFloat3D &field, const MdFloatVector<DIM> &r)
{
    auto iIndexInGrid = field.getGrid().coordinateToFloorMultiIndex(r);
    auto w            = field.getGrid().gridVectorFromGridPointToCoordinate(r, iIndexInGrid);

    std::array<std::array<std::array<real, 2>, 2>, 2> cube;
    const auto &lattice = field.getGrid().lattice();
    for (int ii_z = 0; ii_z <= 1; ++ii_z)
    {
        for (int ii_y = 0; ii_y <= 1; ++ii_y)
        {
            for (int ii_x = 0; ii_x <= 1; ++ii_x)
            {
                offset<DIM> cube_index = {iIndexInGrid[XX]+ii_x, iIndexInGrid[YY]+ii_y, iIndexInGrid[ZZ]+ii_z};
                if (lattice.contains(cube_index))
                {
                    cube[ii_x][ii_y][ii_z] = field[cube_index];
                }
                else
                {
                    cube[ii_x][ii_y][ii_z] = 0;
                }
            }
        }
    }

    std::array<std::array<real, 2>, 2> interpolated_x;
    for (int ii_z = 0; ii_z <= 1; ++ii_z)
    {
        for (int ii_y = 0; ii_y <= 1; ++ii_y)
        {
            interpolated_x[ii_y][ii_z] =
                (1 - w[XX]) * cube[0][ii_y][ii_z] + (w[XX])*cube[1][ii_y][ii_z];
        }
    }

    std::array<real, 2> interpolated_xy;
    for (int ii_z = 0; ii_z <= 1; ++ii_z)
    {
        interpolated_xy[ii_z] = (1 - w[YY]) * interpolated_x[0][ii_z] +
            (w[YY])*interpolated_x[1][ii_z];
    }

    return ((1 - w[ZZ]) * interpolated_xy[0] + w[ZZ] * interpolated_xy[1]) / 8.0;
}


}
