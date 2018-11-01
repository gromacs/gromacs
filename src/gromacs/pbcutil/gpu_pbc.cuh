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
#ifndef GMX_PBCUTIL_GPU_PBC_CUH
#define GMX_PBCUTIL_GPU_PBC_CUH

#include "gromacs/pbcutil/ishift.h"

struct PbcAiuc
{
    float invBoxDiagZ;
    float boxZX;
    float boxZY;
    float boxZZ;
    float invBoxDiagY;
    float boxYX;
    float boxYY;
    float invBoxDiagX;
    float boxXX;
};

template <bool returnShift>
static __forceinline__ __device__
int pbcDxAiuc(const PbcAiuc &pbcAiuc,
              const float4  &x1,
              const float4  &x2,
              fvec           dx)
{
    dx[XX] = x1.x - x2.x;
    dx[YY] = x1.y - x2.y;
    dx[ZZ] = x1.z - x2.z;

    float shz  = rintf(dx[ZZ]*pbcAiuc.invBoxDiagZ);
    dx[XX]    -= shz*pbcAiuc.boxZX;
    dx[YY]    -= shz*pbcAiuc.boxZY;
    dx[ZZ]    -= shz*pbcAiuc.boxZZ;

    float shy  = rintf(dx[YY]*pbcAiuc.invBoxDiagY);
    dx[XX]    -= shy*pbcAiuc.boxYX;
    dx[YY]    -= shy*pbcAiuc.boxYY;

    float shx  = rintf(dx[XX]*pbcAiuc.invBoxDiagX);
    dx[XX]    -= shx*pbcAiuc.boxXX;

    if (returnShift)
    {
        ivec ishift;

        ishift[XX] = -__float2int_rn(shx);
        ishift[YY] = -__float2int_rn(shy);
        ishift[ZZ] = -__float2int_rn(shz);

        return IVEC2IS(ishift);
    }
    else
    {
        return 0;
    }
}

#endif
