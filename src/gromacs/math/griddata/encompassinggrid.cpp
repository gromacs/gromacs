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

 #include "encompassinggrid.h"
namespace gmx
{

GridWithTranslation<DIM> encompassingGridFromCoordinates(const std::vector<RVec> &coordinates, real spacing, real margin)
{
    RVec realSpaceExtend = {0., 0., 0.};
    RVec translation;

    // find the extend of input coordinates
    for (int i = XX; i <= ZZ; i++)
    {
        auto compareIthComponent = [i](RVec a, RVec b) {
                return a[i] < b[i];
            };
        auto minMaxX             =
            minmax_element(coordinates.begin(), coordinates.end(), compareIthComponent);
        realSpaceExtend[i] =
            2 * margin + (*minMaxX.second)[i] - (*minMaxX.first)[i];
        translation[i] = -margin + (*minMaxX.first)[i];
    }

    GridWithTranslation<DIM>::MultiIndex extend {{
                                                     (int)ceil(realSpaceExtend[XX] / spacing), (int)ceil(realSpaceExtend[YY] / spacing), (int)ceil(realSpaceExtend[ZZ] / spacing)
                                                 }};

    GridWithTranslation<DIM> outputdensitygrid( {{{extend[XX] * spacing, extend[YY] * spacing, extend[ZZ] * spacing}}}, extend, {{roundf(translation[XX] / spacing) * spacing, roundf(translation[YY] / spacing) * spacing, roundf(translation[ZZ] / spacing) * spacing}});

    return outputdensitygrid;
}
}
