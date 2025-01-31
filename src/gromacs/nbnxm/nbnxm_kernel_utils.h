/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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

#ifndef GMX_GPU_UTILS_NBNXM_KERNEL_UTILS_H
#define GMX_GPU_UTILS_NBNXM_KERNEL_UTILS_H

/*! \internal \file
 *  \brief
 *  NBNXM GPU kernel utility methods
 *
 *  \ingroup module_gpu_utils
 */

#include "config.h"

#include "gromacs/gpu_utils/gpu_kernel_utils.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/nbnxm/gpu_types_common.h"
#include "gromacs/nbnxm/nbnxm.h"

namespace gmx
{

//! \brief Convert \p sigma and \p epsilon VdW parameters to \c c6,c12 pair.
static inline GMX_ALWAYS_INLINE Float2 convertSigmaEpsilonToC6C12(const float sigma, const float epsilon)
{
    const float sigma2 = sigma * sigma;
    const float sigma6 = sigma2 * sigma2 * sigma2;
    const float c6     = epsilon * sigma6;
    const float c12    = c6 * sigma6;

    return { c6, c12 };
}

//! \brief Calculate force and energy for a pair of atoms, VdW force-switch flavor.
template<bool doCalcEnergies>
static inline GMX_ALWAYS_INLINE void ljForceSwitch(const shift_consts_t dispersionShift,
                                                   const shift_consts_t repulsionShift,
                                                   const float          rVdwSwitch,
                                                   const float          c6,
                                                   const float          c12,
                                                   const float          rInv,
                                                   const float          r2,
                                                   float*               fInvR,
                                                   float*               eLJ)
{
    /* force switch constants */
    const float dispShiftV2 = dispersionShift.c2;
    const float dispShiftV3 = dispersionShift.c3;
    const float repuShiftV2 = repulsionShift.c2;
    const float repuShiftV3 = repulsionShift.c3;

    const float r       = r2 * rInv;
    const float rSwitch = gmxGpuFDim(r, rVdwSwitch); // max(r - rVdwSwitch, 0)

    *fInvR += -c6 * (dispShiftV2 + dispShiftV3 * rSwitch) * rSwitch * rSwitch * rInv
              + c12 * (repuShiftV2 + repuShiftV3 * rSwitch) * rSwitch * rSwitch * rInv;

    if constexpr (doCalcEnergies)
    {
        const float dispShiftF2 = dispShiftV2 / 3.0F;
        const float dispShiftF3 = dispShiftV3 / 4.0F;
        const float repuShiftF2 = repuShiftV2 / 3.0F;
        const float repuShiftF3 = repuShiftV3 / 4.0F;
        *eLJ += c6 * (dispShiftF2 + dispShiftF3 * rSwitch) * rSwitch * rSwitch * rSwitch
                - c12 * (repuShiftF2 + repuShiftF3 * rSwitch) * rSwitch * rSwitch * rSwitch;
    }
}

//! \brief Fetch C6 grid contribution coefficients and return the product of these.
template<gmx::VdwType vdwType>
static inline GMX_ALWAYS_INLINE float calculateLJEwaldC6Grid(const Float2* a_nbfpComb,
                                                             const int     typeI,
                                                             const int     typeJ)
{
    if constexpr (vdwType == gmx::VdwType::EwaldGeom)
    {
        Float2Wrapper wrapI(fastLoad(a_nbfpComb, typeI));
        Float2Wrapper wrapJ(fastLoad(a_nbfpComb, typeJ));
        return wrapI[0] * wrapJ[0];
    }
    else
    {
        static_assert(vdwType == gmx::VdwType::EwaldLB);
        /* sigma and epsilon are scaled to give 6*C6 */
        const Float2Wrapper c6c12_i(fastLoad(a_nbfpComb, typeI));
        const Float2Wrapper c6c12_j(fastLoad(a_nbfpComb, typeJ));

        const float sigma   = c6c12_i[0] + c6c12_j[0];
        const float epsilon = c6c12_i[1] * c6c12_j[1];

        const float sigma2 = sigma * sigma;
        return epsilon * sigma2 * sigma2 * sigma2;
    }
}

//! Calculate LJ-PME grid force contribution with geometric or LB combination rule.
template<bool doCalcEnergies, gmx::VdwType vdwType>
static inline GMX_ALWAYS_INLINE void ljEwaldComb(const Float2* a_nbfpComb,
                                                 const float   sh_lj_ewald,
                                                 const int     typeI,
                                                 const int     typeJ,
                                                 const float   r2,
                                                 const float   r2Inv,
                                                 const float   lje_coeff2,
                                                 const float   lje_coeff6_6,
                                                 const float   int_bit,
                                                 float*        fInvR,
                                                 float*        eLJ)
{
    const float c6grid = calculateLJEwaldC6Grid<vdwType>(a_nbfpComb, typeI, typeJ);

    /* Recalculate inv_r6 without exclusion mask */
    const float inv_r6_nm = r2Inv * r2Inv * r2Inv;
    const float cr2       = lje_coeff2 * r2;
    const float expmcr2   = gmxGpuExp(-cr2);
    const float poly      = 1.0F + cr2 + 0.5F * cr2 * cr2;

    /* Subtract the grid force from the total LJ force */
    *fInvR += c6grid * (inv_r6_nm - expmcr2 * (inv_r6_nm * poly + lje_coeff6_6)) * r2Inv;

    if constexpr (doCalcEnergies)
    {
        /* Shift should be applied only to real LJ pairs */
        const float sh_mask = sh_lj_ewald * int_bit;
        *eLJ += c_oneSixth * c6grid * (inv_r6_nm * (1.0F - expmcr2 * poly) + sh_mask);
    }
}

/*! \brief Apply potential switch. */
template<bool doCalcEnergies>
static inline GMX_ALWAYS_INLINE void ljPotentialSwitch(const switch_consts_t vdwSwitch,
                                                       const float           rVdwSwitch,
                                                       const float           rInv,
                                                       const float           r2,
                                                       float*                fInvR,
                                                       float*                eLJ)
{
    /* potential switch constants */
    const float switchV3 = vdwSwitch.c3;
    const float switchV4 = vdwSwitch.c4;
    const float switchV5 = vdwSwitch.c5;
    const float switchF2 = 3.0F * vdwSwitch.c3;
    const float switchF3 = 4.0F * vdwSwitch.c4;
    const float switchF4 = 5.0F * vdwSwitch.c5;

    const float r       = r2 * rInv;
    const float rSwitch = r - rVdwSwitch;

    if (rSwitch > 0.0F)
    {
        const float sw =
                1.0F + (switchV3 + (switchV4 + switchV5 * rSwitch) * rSwitch) * rSwitch * rSwitch * rSwitch;
        const float dsw = (switchF2 + (switchF3 + switchF4 * rSwitch) * rSwitch) * rSwitch * rSwitch;

        *fInvR = (*fInvR) * sw - rInv * (*eLJ) * dsw;
        if constexpr (doCalcEnergies)
        {
            *eLJ *= sw;
        }
    }
}


/*! \brief Calculate analytical Ewald correction term. */
static inline GMX_ALWAYS_INLINE float pmeCorrF(const float z2)
{
    constexpr float FN6 = -1.7357322914161492954e-8F;
    constexpr float FN5 = 1.4703624142580877519e-6F;
    constexpr float FN4 = -0.000053401640219807709149F;
    constexpr float FN3 = 0.0010054721316683106153F;
    constexpr float FN2 = -0.019278317264888380590F;
    constexpr float FN1 = 0.069670166153766424023F;
    constexpr float FN0 = -0.75225204789749321333F;

    constexpr float FD4 = 0.0011193462567257629232F;
    constexpr float FD3 = 0.014866955030185295499F;
    constexpr float FD2 = 0.11583842382862377919F;
    constexpr float FD1 = 0.50736591960530292870F;
    constexpr float FD0 = 1.0F;

    const float z4 = z2 * z2;

    float       polyFD0 = FD4 * z4 + FD2;
    const float polyFD1 = FD3 * z4 + FD1;
    polyFD0             = polyFD0 * z4 + FD0;
    polyFD0             = polyFD1 * z2 + polyFD0;

    polyFD0 = 1.0F / polyFD0;

    float polyFN0 = FN6 * z4 + FN4;
    float polyFN1 = FN5 * z4 + FN3;
    polyFN0       = polyFN0 * z4 + FN2;
    polyFN1       = polyFN1 * z4 + FN1;
    polyFN0       = polyFN0 * z4 + FN0;
    polyFN0       = polyFN1 * z2 + polyFN0;

    return polyFN0 * polyFD0;
}

/*! \brief Interpolate Ewald coulomb force correction using the F*r table. */
static inline GMX_ALWAYS_INLINE float interpolateCoulombForceR(const float* a_coulombTab,
                                                               const float  coulombTabScale,
                                                               const float  r)
{
    const float normalized = coulombTabScale * r;
    const int   index      = static_cast<int>(normalized);
#if (defined(__SYCL_DEVICE_ONLY__) && defined(__AMDGCN__)) || GMX_GPU_HIP
    // TODO: ROCm compiler at least up to 6.0 does not do this transformation. Remove when this is no longer the case.
    const float fraction = __builtin_amdgcn_fractf(normalized);
    // TODO: On AMD, we want to issue a single GLOBAL_LOAD_DWORDX2 which the optimizer
    // does not always recognize with two consecutive loads (despite the lack of alignment
    // requirement on GCN/CDNA/RDNA). Hence we manually force a 64-bit load using a custom
    // type because the native float2 type has alignment requirements that are too strict.
#    if defined(__SYCL_DEVICE_ONLY__)
    // This improves SYCL kernel performance on gfx90a by up to 10% (in particular for more complex
    // kernels like Ewald&LJ force switch). Remove this when no longer the case.
    struct myFloat2
    {
        float x, y;
    };
    static_assert(alignof(myFloat2) == alignof(decltype(*a_coulombTab)),
                  "Custom type for codegen optimization should have same alignment as plain float");
    const myFloat2 lr   = *(reinterpret_cast<const myFloat2*>(indexedAddress(a_coulombTab, index)));
    const float    left = lr.x;
    const float    right = lr.y;
#    else
    // non-SYCL default
    const float left  = *indexedAddress(a_coulombTab, index);
    const float right = *indexedAddress(a_coulombTab, index + 1);
#    endif
#else
    const float fraction = normalized - index;
    const float left     = a_coulombTab[index];
    const float right    = a_coulombTab[index + 1];
#endif

    return lerp(left, right, fraction);
}

} // namespace gmx

#endif
