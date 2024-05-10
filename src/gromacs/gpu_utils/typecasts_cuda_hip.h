/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \libinternal \file
 *  \brief Declare functions to be used to cast CPU types to compatible GPU types.
 *
 *  \author Artem Zhmurov <zhmurov@gmail.com>
 *
 *  \inlibraryapi
 */
#ifndef GMX_GPU_UTILS_TYPECASTS_CUDA_HIP_H
#define GMX_GPU_UTILS_TYPECASTS_CUDA_HIP_H

#include "gmxpre.h"

#include "gromacs/math/vectypes.h"

/*! \brief Cast RVec buffer to float3 buffer.
 *
 * \param[in] in The RVec buffer to cast.
 *
 * \returns Buffer, casted to float3*.
 */
__forceinline__ static __host__ __device__ float3* asFloat3(gmx::RVec* in)
{
    static_assert(sizeof(in[0]) == sizeof(float3),
                  "Size of the host-side data-type is different from the size of the device-side "
                  "counterpart.");
    return reinterpret_cast<float3*>(in);
}

/*! \brief Cast pointer RVec buffer to a pointer to float3 buffer.
 *
 * \param[in] in The Pointer to RVec buffer to cast.
 *
 * \returns Buffer pointer, casted to float3*.
 */
__forceinline__ static __host__ __device__ float3** asFloat3Pointer(gmx::RVec** in)
{
    static_assert(sizeof((*in)[0]) == sizeof(float3),
                  "Size of the host-side data-type is different from the size of the device-side "
                  "counterpart.");
    return reinterpret_cast<float3**>(in);
}
static inline __host__ __device__ const float3* const* asFloat3Pointer(const gmx::RVec* const* in)
{
    static_assert(sizeof((*in)[0]) == sizeof(float3),
                  "Size of the host-side data-type is different from the size of the device-side "
                  "counterpart.");
    return reinterpret_cast<const float3* const*>(in);
}

#endif // GMX_GPU_UTILS_TYPECASTS_CUDA_HIP_H
