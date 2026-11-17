/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
/*! \internal \file
 *
 * \brief Backend-agnostic device functions for virtual site construction and force spreading.
 *
 * \author Lukas Müllender <lukas.muellender@gmail.com>
 *
 * \ingroup module_mdlib
 */
#ifndef GMX_MDLIB_VSITE_GPU_INTERNAL_SHARED_H
#define GMX_MDLIB_VSITE_GPU_INTERNAL_SHARED_H

#include "config.h"

#if GMX_GPU_CUDA
#    include "gromacs/gpu_utils/cuda_kernel_utils.cuh"
#    include "gromacs/gpu_utils/vectype_ops_cuda.h"
#    include "gromacs/pbcutil/pbc_aiuc_cuda.cuh"
#else
#    error "Building GPU vsite for unsupported backend"
#endif

#include "gromacs/gpu_utils/gputraits.h"

#ifndef DOXYGEN

namespace gmx
{

// ---------------------------------------------------------------------------
// Construction device functions
// ---------------------------------------------------------------------------

/*! \brief VSite1 construction: x_v = x_i.
 *
 * \param[in]     idx        Vsite index (one per thread).
 * \param[in]     gm_iatoms  Flat iatom index array; stride 2: [vsite, i].
 * \param[in,out] gm_x       Coordinate array (read for atom i, written for vsite).
 */
GMX_DEVICE_FUNC_ATTRIBUTE static inline void constructVSite1Gpu(int idx,
                                                                DeviceGlobalPtr<const int> gm_iatoms,
                                                                DeviceGlobalPtr<DeviceFloat3> gm_x)
{
    int vsiteAtom   = gm_iatoms[2 * idx];
    int atomI       = gm_iatoms[2 * idx + 1];
    gm_x[vsiteAtom] = gm_x[atomI];
}

/*! \brief VSite2 construction: x_v = x_i + a*(x_j - x_i)
 *
 * \param[in]     idx        Vsite index (one per thread).
 * \param[in]     gm_iatoms  Flat iatom index array; stride 3: [vsite, i, j].
 * \param[in]     gm_params  Per-vsite weight a; one float per vsite.
 * \param[in,out] gm_x       Coordinate array.
 * \param[in]     pbcAiuc    PBC parameters.
 */
GMX_DEVICE_FUNC_ATTRIBUTE static inline void constructVSite2Gpu(int idx,
                                                                DeviceGlobalPtr<const int> gm_iatoms,
                                                                DeviceGlobalPtr<const float> gm_params,
                                                                DeviceGlobalPtr<DeviceFloat3> gm_x,
                                                                const PbcAiuc& pbcAiuc)
{
    int vsiteAtom = gm_iatoms[3 * idx];
    int atomI     = gm_iatoms[3 * idx + 1];
    int atomJ     = gm_iatoms[3 * idx + 2];

    float        a      = gm_params[idx];
    DeviceFloat3 old_xv = gm_x[vsiteAtom];
    DeviceFloat3 xi     = gm_x[atomI];
    DeviceFloat3 xj     = gm_x[atomJ];

    DeviceFloat3 dij    = pbcDxAiuc(pbcAiuc, xj, xi);
    DeviceFloat3 new_xv = xi + a * dij;

    // Keep vsite in the same periodic image as before
    DeviceFloat3 dx = pbcDxAiuc(pbcAiuc, new_xv, old_xv);
    gm_x[vsiteAtom] = old_xv + dx;
}

/*! \brief VSite3 construction: x_v = x_i + a*(x_j - x_i) + b*(x_k - x_i)
 *
 * \param[in]     idx        Vsite index (one per thread).
 * \param[in]     gm_iatoms  Flat iatom index array; stride 4: [vsite, i, j, k].
 * \param[in]     gm_params  Per-vsite weights (a, b): [a0, b0, a1, b1, ...].
 * \param[in,out] gm_x       Coordinate array.
 * \param[in]     pbcAiuc    PBC parameters.
 */
GMX_DEVICE_FUNC_ATTRIBUTE static inline void constructVSite3Gpu(int idx,
                                                                DeviceGlobalPtr<const int> gm_iatoms,
                                                                DeviceGlobalPtr<const float> gm_params,
                                                                DeviceGlobalPtr<DeviceFloat3> gm_x,
                                                                const PbcAiuc& pbcAiuc)
{
    int vsiteAtom = gm_iatoms[4 * idx];
    int atomI     = gm_iatoms[4 * idx + 1];
    int atomJ     = gm_iatoms[4 * idx + 2];
    int atomK     = gm_iatoms[4 * idx + 3];

    float        a      = gm_params[2 * idx];
    float        b      = gm_params[2 * idx + 1];
    DeviceFloat3 old_xv = gm_x[vsiteAtom];
    DeviceFloat3 xi     = gm_x[atomI];
    DeviceFloat3 xj     = gm_x[atomJ];
    DeviceFloat3 xk     = gm_x[atomK];

    DeviceFloat3 dij    = pbcDxAiuc(pbcAiuc, xj, xi);
    DeviceFloat3 dik    = pbcDxAiuc(pbcAiuc, xk, xi);
    DeviceFloat3 new_xv = xi + a * dij + b * dik;

    // Keep vsite in the same periodic image as before.
    DeviceFloat3 dx = pbcDxAiuc(pbcAiuc, new_xv, old_xv);
    gm_x[vsiteAtom] = old_xv + dx;
}

// ---------------------------------------------------------------------------
// Force-spreading device functions
// ---------------------------------------------------------------------------

/*! \brief VSite1 force spreading: f_i += f_v.
 *
 * \param[in]     idx        Vsite index (one per thread).
 * \param[in]     gm_iatoms  Flat iatom index array; stride 2: [vsite, i].
 * \param[in,out] gm_f       Force array.
 */
GMX_DEVICE_FUNC_ATTRIBUTE static inline void spreadVSite1ForceGpu(int idx,
                                                                  DeviceGlobalPtr<const int> gm_iatoms,
                                                                  DeviceGlobalPtr<DeviceFloat3> gm_f)
{
    int          vsiteAtom = gm_iatoms[2 * idx];
    int          atomI     = gm_iatoms[2 * idx + 1];
    DeviceFloat3 fv        = gm_f[vsiteAtom];
    atomicAdd(&gm_f[atomI].x, fv.x);
    atomicAdd(&gm_f[atomI].y, fv.y);
    atomicAdd(&gm_f[atomI].z, fv.z);
}

/*! \brief VSite2 force spreading: f_i += (1-a)*f_v, f_j += a*f_v.
 *
 * \param[in]     idx        Vsite index (one per thread).
 * \param[in]     gm_iatoms  Flat iatom index array; stride 3: [vsite, i, j].
 * \param[in]     gm_params  Per-vsite weight a; one float per vsite.
 * \param[in,out] gm_f       Force array.
 */
GMX_DEVICE_FUNC_ATTRIBUTE static inline void spreadVSite2ForceGpu(int idx,
                                                                  DeviceGlobalPtr<const int> gm_iatoms,
                                                                  DeviceGlobalPtr<const float> gm_params,
                                                                  DeviceGlobalPtr<DeviceFloat3> gm_f)
{
    int          vsiteAtom = gm_iatoms[3 * idx];
    int          atomI     = gm_iatoms[3 * idx + 1];
    int          atomJ     = gm_iatoms[3 * idx + 2];
    float        a         = gm_params[idx];
    DeviceFloat3 fv        = gm_f[vsiteAtom];
    float        b         = 1.0F - a;
    atomicAdd(&gm_f[atomI].x, b * fv.x);
    atomicAdd(&gm_f[atomI].y, b * fv.y);
    atomicAdd(&gm_f[atomI].z, b * fv.z);
    atomicAdd(&gm_f[atomJ].x, a * fv.x);
    atomicAdd(&gm_f[atomJ].y, a * fv.y);
    atomicAdd(&gm_f[atomJ].z, a * fv.z);
}

/*! \brief VSite3 force spreading: f_i += (1-a-b)*f_v, f_j += a*f_v, f_k += b*f_v.
 *
 * \param[in]     idx        Vsite index (one per thread).
 * \param[in]     gm_iatoms  Flat iatom index array; stride 4: [vsite, i, j, k].
 * \param[in]     gm_params  Per-vsite weights (a, b): [a0, b0, a1, b1, ...].
 * \param[in,out] gm_f       Force array.
 */
GMX_DEVICE_FUNC_ATTRIBUTE static inline void spreadVSite3ForceGpu(int idx,
                                                                  DeviceGlobalPtr<const int> gm_iatoms,
                                                                  DeviceGlobalPtr<const float> gm_params,
                                                                  DeviceGlobalPtr<DeviceFloat3> gm_f)
{
    int          vsiteAtom = gm_iatoms[4 * idx];
    int          atomI     = gm_iatoms[4 * idx + 1];
    int          atomJ     = gm_iatoms[4 * idx + 2];
    int          atomK     = gm_iatoms[4 * idx + 3];
    float        a         = gm_params[2 * idx];
    float        b         = gm_params[2 * idx + 1];
    float        c         = 1.0F - a - b;
    DeviceFloat3 fv        = gm_f[vsiteAtom];
    atomicAdd(&gm_f[atomI].x, c * fv.x);
    atomicAdd(&gm_f[atomI].y, c * fv.y);
    atomicAdd(&gm_f[atomI].z, c * fv.z);
    atomicAdd(&gm_f[atomJ].x, a * fv.x);
    atomicAdd(&gm_f[atomJ].y, a * fv.y);
    atomicAdd(&gm_f[atomJ].z, a * fv.z);
    atomicAdd(&gm_f[atomK].x, b * fv.x);
    atomicAdd(&gm_f[atomK].y, b * fv.y);
    atomicAdd(&gm_f[atomK].z, b * fv.z);
}

} // namespace gmx

#endif // DOXYGEN

#endif // GMX_MDLIB_VSITE_GPU_INTERNAL_SHARED_H
