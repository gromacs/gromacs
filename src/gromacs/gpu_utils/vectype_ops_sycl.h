/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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

#ifndef GMX_GPU_UTILS_VECTYPE_OPS_SYCL_H
#define GMX_GPU_UTILS_VECTYPE_OPS_SYCL_H

#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/gpu_utils/gputraits_sycl.h"

/* \brief Compute the scalar product of two vectors.
 *
 * \param[in] a  First vector.
 * \param[in] b  Second vector.
 * \returns Scalar product.
 */
static inline float gmxDeviceInternalProd(const Float3 a, const Float3 b)
{
    return a.dot(b);
}

/* \brief Compute the vector product of two vectors.
 *
 * \param[in] a  First vector.
 * \param[in] b  Second vector.
 * \returns Vector product.
 */
static inline Float3 gmxDeviceCrossProd(const Float3 a, const Float3 b)
{
    return a.cross(b);
}


static inline float gmxDeviceSin(const float a)
{
    return sycl::sin(a);
}

static inline float gmxDeviceCos(const float a)
{
    return sycl::cos(a);
}

static float inline gmxDeviceSqrt(float input)
{
    return sycl::sqrt(input);
}

static inline float gmxDeviceRSqrt(float input)
{
    return sycl::rsqrt(input);
}

static inline float gmxDeviceAcos(float input)
{
    return sycl::acos(input);
}

static inline float gmxDeviceNorm(Float3 a)
{
    return a.norm();
}

static inline float gmxDeviceNorm2(Float3 a)
{
    return a.norm2();
}

/*! \brief Cosine of an angle between two vectors.
 *
 * Computes cosine using the following formula:
 *
 *                  ax*bx + ay*by + az*bz
 * cos-vec (a,b) =  ---------------------
 *                      ||a|| * ||b||
 *
 * This function also makes sure that the cosine does not leave the [-1, 1]
 * interval, which can happen due to numerical errors.
 *
 * \param[in] a  First vector.
 * \param[in] b  Second vector.
 * \returns Cosine between a and b.
 */
static inline float gmxDeviceCosAngle(const Float3& a, const Float3& b)
{
    float cosval;

    float ipa  = gmxDeviceNorm2(a);
    float ipb  = gmxDeviceNorm2(b);
    float ip   = gmxDeviceInternalProd(a, b);
    float ipab = ipa * ipb;
    if (ipab > 0.0F)
    {
        cosval = ip * gmxDeviceRSqrt(ipab);
    }
    else
    {
        cosval = 1.0F;
    }
    if (cosval > 1.0F)
    {
        return 1.0F;
    }
    if (cosval < -1.0F)
    {
        return -1.0F;
    }

    return cosval;
}

/*! \brief Compute the angle between two vectors.
 *
 * Uses atan( |axb| / a.b ) formula.
 *
 * \param[in] a  First vector.
 * \param[in] b  Second vector.
 * \returns Angle between vectors in radians.
 */
static inline float gmxDeviceAngle(const Float3& a, const Float3& b)
{
    Float3 w = gmxDeviceCrossProd(a, b);

    float wlen = gmxDeviceSqrt(gmxDeviceNorm2(w));
    float s    = gmxDeviceInternalProd(a, b);

    return sycl::atan2(wlen, s);
}

#endif /* GMX_GPU_UTILS_VECTYPE_OPS_SYCL_H */
