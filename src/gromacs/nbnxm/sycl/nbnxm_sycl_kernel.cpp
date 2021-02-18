/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020,2021, by the GROMACS development team, led by
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

/*! \internal \file
 *  \brief
 *  NBNXM SYCL kernels
 *
 *  \ingroup module_nbnxm
 */
#include "gmxpre.h"

#include "nbnxm_sycl_kernel.h"

#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/math/functions.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/utility/template_mp.h"

#include "nbnxm_sycl_kernel_utils.h"
#include "nbnxm_sycl_types.h"

namespace Nbnxm
{

//! \brief Set of boolean constants mimicking preprocessor macros.
template<enum ElecType elecType, enum VdwType vdwType>
struct EnergyFunctionProperties {
    static constexpr bool elecCutoff = (elecType == ElecType::Cut); ///< EL_CUTOFF
    static constexpr bool elecRF     = (elecType == ElecType::RF);  ///< EL_RF
    static constexpr bool elecEwaldAna =
            (elecType == ElecType::EwaldAna || elecType == ElecType::EwaldAnaTwin); ///< EL_EWALD_ANA
    static constexpr bool elecEwaldTab =
            (elecType == ElecType::EwaldTab || elecType == ElecType::EwaldTabTwin); ///< EL_EWALD_TAB
    static constexpr bool elecEwaldTwin =
            (elecType == ElecType::EwaldAnaTwin || elecType == ElecType::EwaldTabTwin);
    static constexpr bool elecEwald        = (elecEwaldAna || elecEwaldTab); ///< EL_EWALD_ANY
    static constexpr bool vdwCombLB        = (vdwType == VdwType::CutCombLB);
    static constexpr bool vdwCombGeom      = (vdwType == VdwType::CutCombGeom); ///< LJ_COMB_GEOM
    static constexpr bool vdwComb          = (vdwCombLB || vdwCombGeom);        ///< LJ_COMB
    static constexpr bool vdwEwaldCombGeom = (vdwType == VdwType::EwaldGeom); ///< LJ_EWALD_COMB_GEOM
    static constexpr bool vdwEwaldCombLB   = (vdwType == VdwType::EwaldLB);   ///< LJ_EWALD_COMB_LB
    static constexpr bool vdwEwald         = (vdwEwaldCombGeom || vdwEwaldCombLB); ///< LJ_EWALD
    static constexpr bool vdwFSwitch       = (vdwType == VdwType::FSwitch); ///< LJ_FORCE_SWITCH
    static constexpr bool vdwPSwitch       = (vdwType == VdwType::PSwitch); ///< LJ_POT_SWITCH
};

//! \brief Templated constants to shorten kernel function declaration.
//@{
template<enum VdwType vdwType>
constexpr bool ljComb = EnergyFunctionProperties<ElecType::Count, vdwType>().vdwComb;

template<enum ElecType elecType> // Yes, ElecType
constexpr bool vdwCutoffCheck = EnergyFunctionProperties<elecType, VdwType::Count>().elecEwaldTwin;

template<enum ElecType elecType>
constexpr bool elecEwald = EnergyFunctionProperties<elecType, VdwType::Count>().elecEwald;

template<enum ElecType elecType>
constexpr bool elecEwaldTab = EnergyFunctionProperties<elecType, VdwType::Count>().elecEwaldTab;

template<enum VdwType vdwType>
constexpr bool ljEwald = EnergyFunctionProperties<ElecType::Count, vdwType>().vdwEwald;
//@}

using cl::sycl::access::fence_space;
using cl::sycl::access::mode;
using cl::sycl::access::target;

static inline void convertSigmaEpsilonToC6C12(const float                  sigma,
                                              const float                  epsilon,
                                              cl::sycl::private_ptr<float> c6,
                                              cl::sycl::private_ptr<float> c12)
{
    const float sigma2 = sigma * sigma;
    const float sigma6 = sigma2 * sigma2 * sigma2;
    *c6                = epsilon * sigma6;
    *c12               = (*c6) * sigma6;
}

template<bool doCalcEnergies>
static inline void ljForceSwitch(const shift_consts_t         dispersionShift,
                                 const shift_consts_t         repulsionShift,
                                 const float                  rVdwSwitch,
                                 const float                  c6,
                                 const float                  c12,
                                 const float                  rInv,
                                 const float                  r2,
                                 cl::sycl::private_ptr<float> fInvR,
                                 cl::sycl::private_ptr<float> eLJ)
{
    /* force switch constants */
    const float dispShiftV2 = dispersionShift.c2;
    const float dispShiftV3 = dispersionShift.c3;
    const float repuShiftV2 = repulsionShift.c2;
    const float repuShiftV3 = repulsionShift.c3;

    const float r       = r2 * rInv;
    const float rSwitch = cl::sycl::fdim(r, rVdwSwitch); // max(r - rVdwSwitch, 0)

    *fInvR += -c6 * (dispShiftV2 + dispShiftV3 * rSwitch) * rSwitch * rSwitch * rInv
              + c12 * (repuShiftV2 + repuShiftV3 * rSwitch) * rSwitch * rSwitch * rInv;

    if constexpr (doCalcEnergies)
    {
        const float dispShiftF2 = dispShiftV2 / 3;
        const float dispShiftF3 = dispShiftV3 / 4;
        const float repuShiftF2 = repuShiftV2 / 3;
        const float repuShiftF3 = repuShiftV3 / 4;
        *eLJ += c6 * (dispShiftF2 + dispShiftF3 * rSwitch) * rSwitch * rSwitch * rSwitch
                - c12 * (repuShiftF2 + repuShiftF3 * rSwitch) * rSwitch * rSwitch * rSwitch;
    }
}

//! \brief Fetch C6 grid contribution coefficients and return the product of these.
template<enum VdwType vdwType>
static inline float calculateLJEwaldC6Grid(const DeviceAccessor<float, mode::read> a_nbfpComb,
                                           const int                               typeI,
                                           const int                               typeJ)
{
    if constexpr (vdwType == VdwType::EwaldGeom)
    {
        return a_nbfpComb[2 * typeI] * a_nbfpComb[2 * typeJ];
    }
    else
    {
        static_assert(vdwType == VdwType::EwaldLB);
        /* sigma and epsilon are scaled to give 6*C6 */
        const float c6_i  = a_nbfpComb[2 * typeI];
        const float c12_i = a_nbfpComb[2 * typeI + 1];
        const float c6_j  = a_nbfpComb[2 * typeJ];
        const float c12_j = a_nbfpComb[2 * typeJ + 1];

        const float sigma   = c6_i + c6_j;
        const float epsilon = c12_i * c12_j;

        const float sigma2 = sigma * sigma;
        return epsilon * sigma2 * sigma2 * sigma2;
    }
}

//! Calculate LJ-PME grid force contribution with geometric or LB combination rule.
template<bool doCalcEnergies, enum VdwType vdwType>
static inline void ljEwaldComb(const DeviceAccessor<float, mode::read> a_nbfpComb,
                               const float                             sh_lj_ewald,
                               const int                               typeI,
                               const int                               typeJ,
                               const float                             r2,
                               const float                             r2Inv,
                               const float                             lje_coeff2,
                               const float                             lje_coeff6_6,
                               const float                             int_bit,
                               cl::sycl::private_ptr<float>            fInvR,
                               cl::sycl::private_ptr<float>            eLJ)
{
    const float c6grid = calculateLJEwaldC6Grid<vdwType>(a_nbfpComb, typeI, typeJ);

    /* Recalculate inv_r6 without exclusion mask */
    const float inv_r6_nm = r2Inv * r2Inv * r2Inv;
    const float cr2       = lje_coeff2 * r2;
    const float expmcr2   = cl::sycl::exp(-cr2);
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
static inline void ljPotentialSwitch(const switch_consts_t        vdwSwitch,
                                     const float                  rVdwSwitch,
                                     const float                  rInv,
                                     const float                  r2,
                                     cl::sycl::private_ptr<float> fInvR,
                                     cl::sycl::private_ptr<float> eLJ)
{
    /* potential switch constants */
    const float switchV3 = vdwSwitch.c3;
    const float switchV4 = vdwSwitch.c4;
    const float switchV5 = vdwSwitch.c5;
    const float switchF2 = 3 * vdwSwitch.c3;
    const float switchF3 = 4 * vdwSwitch.c4;
    const float switchF4 = 5 * vdwSwitch.c5;

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
static inline float pmeCorrF(const float z2)
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

/*! \brief Linear interpolation using exactly two FMA operations.
 *
 *  Implements numeric equivalent of: (1-t)*d0 + t*d1.
 */
template<typename T>
static inline T lerp(T d0, T d1, T t)
{
    return fma(t, d1, fma(-t, d0, d0));
}

/*! \brief Interpolate Ewald coulomb force correction using the F*r table. */
static inline float interpolateCoulombForceR(const DeviceAccessor<float, mode::read> a_coulombTab,
                                             const float coulombTabScale,
                                             const float r)
{
    const float normalized = coulombTabScale * r;
    const int   index      = static_cast<int>(normalized);
    const float fraction   = normalized - index;

    const float left  = a_coulombTab[index];
    const float right = a_coulombTab[index + 1];

    return lerp(left, right, fraction); // TODO: cl::sycl::mix
}

static inline void reduceForceJShuffle(float3                                  f,
                                       const cl::sycl::nd_item<1>              itemIdx,
                                       const int                               tidxi,
                                       const int                               aidx,
                                       DeviceAccessor<float, mode::read_write> a_f)
{
    static_assert(c_clSize == 8 || c_clSize == 4);
    sycl_2020::sub_group sg = itemIdx.get_sub_group();

    f[0] += shuffleDown(f[0], 1, sg);
    f[1] += shuffleUp(f[1], 1, sg);
    f[2] += shuffleDown(f[2], 1, sg);
    if (tidxi & 1)
    {
        f[0] = f[1];
    }

    f[0] += shuffleDown(f[0], 2, sg);
    f[2] += shuffleUp(f[2], 2, sg);
    if (tidxi & 2)
    {
        f[0] = f[2];
    }

    if constexpr (c_clSize == 8)
    {
        f[0] += shuffleDown(f[0], 4, sg);
    }

    if (tidxi < 3)
    {
        atomicFetchAdd(a_f, 3 * aidx + tidxi, f[0]);
    }
}


/*! \brief Final i-force reduction.
 *
 * This implementation works only with power of two array sizes.
 */
static inline void reduceForceIAndFShift(cl::sycl::accessor<float, 1, mode::read_write, target::local> sm_buf,
                                         const float3 fCiBuf[c_nbnxnGpuNumClusterPerSupercluster],
                                         const bool   calcFShift,
                                         const cl::sycl::nd_item<1>              itemIdx,
                                         const int                               tidxi,
                                         const int                               tidxj,
                                         const int                               sci,
                                         const int                               shift,
                                         DeviceAccessor<float, mode::read_write> a_f,
                                         DeviceAccessor<float, mode::read_write> a_fShift)
{
    static constexpr int bufStride  = c_clSize * c_clSize;
    static constexpr int clSizeLog2 = gmx::StaticLog2<c_clSize>::value;
    const int            tidx       = tidxi + tidxj * c_clSize;
    float                fShiftBuf  = 0;
    for (int ciOffset = 0; ciOffset < c_nbnxnGpuNumClusterPerSupercluster; ciOffset++)
    {
        const int aidx = (sci * c_nbnxnGpuNumClusterPerSupercluster + ciOffset) * c_clSize + tidxi;
        /* store i forces in shmem */
        sm_buf[tidx]                 = fCiBuf[ciOffset][0];
        sm_buf[bufStride + tidx]     = fCiBuf[ciOffset][1];
        sm_buf[2 * bufStride + tidx] = fCiBuf[ciOffset][2];
        itemIdx.barrier(fence_space::local_space);

        /* Reduce the initial c_clSize values for each i atom to half
         * every step by using c_clSize * i threads. */
        int i = c_clSize / 2;
        for (int j = clSizeLog2 - 1; j > 0; j--)
        {
            if (tidxj < i)
            {
                sm_buf[tidxj * c_clSize + tidxi] += sm_buf[(tidxj + i) * c_clSize + tidxi];
                sm_buf[bufStride + tidxj * c_clSize + tidxi] +=
                        sm_buf[bufStride + (tidxj + i) * c_clSize + tidxi];
                sm_buf[2 * bufStride + tidxj * c_clSize + tidxi] +=
                        sm_buf[2 * bufStride + (tidxj + i) * c_clSize + tidxi];
            }
            i >>= 1;
            itemIdx.barrier(fence_space::local_space);
        }

        /* i == 1, last reduction step, writing to global mem */
        /* Split the reduction between the first 3 line threads
           Threads with line id 0 will do the reduction for (float3).x components
           Threads with line id 1 will do the reduction for (float3).y components
           Threads with line id 2 will do the reduction for (float3).z components. */
        if (tidxj < 3)
        {
            const float f =
                    sm_buf[tidxj * bufStride + tidxi] + sm_buf[tidxj * bufStride + c_clSize + tidxi];
            atomicFetchAdd(a_f, 3 * aidx + tidxj, f);
            if (calcFShift)
            {
                fShiftBuf += f;
            }
        }
        itemIdx.barrier(fence_space::local_space);
    }
    /* add up local shift forces into global mem */
    if (calcFShift)
    {
        /* Only threads with tidxj < 3 will update fshift.
           The threads performing the update must be the same as the threads
           storing the reduction result above. */
        if (tidxj < 3)
        {
            atomicFetchAdd(a_fShift, 3 * shift + tidxj, fShiftBuf);
        }
    }
}


/*! \brief Main kernel for NBNXM.
 *
 */
template<bool doPruneNBL, bool doCalcEnergies, enum ElecType elecType, enum VdwType vdwType>
auto nbnxmKernel(cl::sycl::handler&                                        cgh,
                 DeviceAccessor<float4, mode::read>                        a_xq,
                 DeviceAccessor<float, mode::read_write>                   a_f,
                 DeviceAccessor<float3, mode::read>                        a_shiftVec,
                 DeviceAccessor<float, mode::read_write>                   a_fShift,
                 OptionalAccessor<float, mode::read_write, doCalcEnergies> a_energyElec,
                 OptionalAccessor<float, mode::read_write, doCalcEnergies> a_energyVdw,
                 DeviceAccessor<nbnxn_cj4_t, doPruneNBL ? mode::read_write : mode::read> a_plistCJ4,
                 DeviceAccessor<nbnxn_sci_t, mode::read>                                 a_plistSci,
                 DeviceAccessor<nbnxn_excl_t, mode::read>                    a_plistExcl,
                 OptionalAccessor<float2, mode::read, ljComb<vdwType>>       a_ljComb,
                 OptionalAccessor<int, mode::read, !ljComb<vdwType>>         a_atomTypes,
                 OptionalAccessor<float, mode::read, !ljComb<vdwType>>       a_nbfp,
                 OptionalAccessor<float, mode::read, ljEwald<vdwType>>       a_nbfpComb,
                 OptionalAccessor<float, mode::read, elecEwaldTab<elecType>> a_coulombTab,
                 const int                                                   numTypes,
                 const float                                                 rCoulombSq,
                 const float                                                 rVdwSq,
                 const float                                                 twoKRf,
                 const float                                                 ewaldBeta,
                 const float                                                 rlistOuterSq,
                 const float                                                 ewaldShift,
                 const float                                                 epsFac,
                 const float                                                 ewaldCoeffLJ,
                 const float                                                 cRF,
                 const shift_consts_t                                        dispersionShift,
                 const shift_consts_t                                        repulsionShift,
                 const switch_consts_t                                       vdwSwitch,
                 const float                                                 rVdwSwitch,
                 const float                                                 ljEwaldShift,
                 const float                                                 coulombTabScale,
                 const bool                                                  calcShift)
{
    static constexpr EnergyFunctionProperties<elecType, vdwType> props;

    cgh.require(a_xq);
    cgh.require(a_f);
    cgh.require(a_shiftVec);
    cgh.require(a_fShift);
    if constexpr (doCalcEnergies)
    {
        cgh.require(a_energyElec);
        cgh.require(a_energyVdw);
    }
    cgh.require(a_plistCJ4);
    cgh.require(a_plistSci);
    cgh.require(a_plistExcl);
    if constexpr (!props.vdwComb)
    {
        cgh.require(a_atomTypes);
        cgh.require(a_nbfp);
    }
    else
    {
        cgh.require(a_ljComb);
    }
    if constexpr (props.vdwEwald)
    {
        cgh.require(a_nbfpComb);
    }
    if constexpr (props.elecEwaldTab)
    {
        cgh.require(a_coulombTab);
    }

    // shmem buffer for i x+q pre-loading
    cl::sycl::accessor<float4, 2, mode::read_write, target::local> sm_xq(
            cl::sycl::range<2>(c_nbnxnGpuNumClusterPerSupercluster, c_clSize), cgh);

    // shmem buffer for force reduction
    // SYCL-TODO: Make into 3D; section 4.7.6.11 of SYCL2020 specs
    cl::sycl::accessor<float, 1, mode::read_write, target::local> sm_reductionBuffer(
            cl::sycl::range<1>(c_clSize * c_clSize * DIM), cgh);

    auto sm_atomTypeI = [&]() {
        if constexpr (!props.vdwComb)
        {
            return cl::sycl::accessor<int, 2, mode::read_write, target::local>(
                    cl::sycl::range<2>(c_nbnxnGpuNumClusterPerSupercluster, c_clSize), cgh);
        }
        else
        {
            return nullptr;
        }
    }();

    auto sm_ljCombI = [&]() {
        if constexpr (props.vdwComb)
        {
            return cl::sycl::accessor<float2, 2, mode::read_write, target::local>(
                    cl::sycl::range<2>(c_nbnxnGpuNumClusterPerSupercluster, c_clSize), cgh);
        }
        else
        {
            return nullptr;
        }
    }();

    /* Flag to control the calculation of exclusion forces in the kernel
     * We do that with Ewald (elec/vdw) and RF. Cut-off only has exclusion
     * energy terms. */
    constexpr bool doExclusionForces =
            (props.elecEwald || props.elecRF || props.vdwEwald || (props.elecCutoff && doCalcEnergies));

    constexpr int subGroupSize = c_clSize * c_clSize / 2;

    return [=](cl::sycl::nd_item<1> itemIdx) [[intel::reqd_sub_group_size(subGroupSize)]]
    {
        /* thread/block/warp id-s */
        const cl::sycl::id<3> localId = unflattenId<c_clSize, c_clSize>(itemIdx.get_local_id());
        const unsigned        tidxi   = localId[0];
        const unsigned        tidxj   = localId[1];
        const unsigned        tidx    = tidxj * c_clSize + tidxi;
        const unsigned        tidxz   = 0;

        // Group indexing was flat originally, no need to unflatten it.
        const unsigned bidx = itemIdx.get_group(0);

        const sycl_2020::sub_group sg = itemIdx.get_sub_group();
        // Better use sg.get_group_range, but too much of the logic relies on it anyway
        const unsigned widx = tidx / subGroupSize;

        float3 fCiBuf[c_nbnxnGpuNumClusterPerSupercluster]; // i force buffer
        for (int i = 0; i < c_nbnxnGpuNumClusterPerSupercluster; i++)
        {
            fCiBuf[i] = float3(0.0F, 0.0F, 0.0F);
        }

        const nbnxn_sci_t nbSci     = a_plistSci[bidx];
        const int         sci       = nbSci.sci;
        const int         cij4Start = nbSci.cj4_ind_start;
        const int         cij4End   = nbSci.cj4_ind_end;

        // Only needed if props.elecEwaldAna
        const float beta2 = ewaldBeta * ewaldBeta;
        const float beta3 = ewaldBeta * ewaldBeta * ewaldBeta;

        for (int i = 0; i < c_nbnxnGpuNumClusterPerSupercluster; i += c_clSize)
        {
            /* Pre-load i-atom x and q into shared memory */
            const int             ci       = sci * c_nbnxnGpuNumClusterPerSupercluster + tidxj + i;
            const int             ai       = ci * c_clSize + tidxi;
            const cl::sycl::id<2> cacheIdx = cl::sycl::id<2>(tidxj + i, tidxi);

            const float3 shift = a_shiftVec[nbSci.shift];
            float4       xqi   = a_xq[ai];
            xqi += float4(shift[0], shift[1], shift[2], 0.0F);
            xqi[3] *= epsFac;
            sm_xq[cacheIdx] = xqi;

            if constexpr (!props.vdwComb)
            {
                // Pre-load the i-atom types into shared memory
                sm_atomTypeI[cacheIdx] = a_atomTypes[ai];
            }
            else
            {
                // Pre-load the LJ combination parameters into shared memory
                sm_ljCombI[cacheIdx] = a_ljComb[ai];
            }
        }
        itemIdx.barrier(fence_space::local_space);

        float ewaldCoeffLJ_2, ewaldCoeffLJ_6_6; // Only needed if (props.vdwEwald)
        if constexpr (props.vdwEwald)
        {
            ewaldCoeffLJ_2   = ewaldCoeffLJ * ewaldCoeffLJ;
            ewaldCoeffLJ_6_6 = ewaldCoeffLJ_2 * ewaldCoeffLJ_2 * ewaldCoeffLJ_2 * c_oneSixth;
        }

        float energyVdw, energyElec; // Only needed if (doCalcEnergies)
        if constexpr (doCalcEnergies)
        {
            energyVdw = energyElec = 0.0F;
        }
        if constexpr (doCalcEnergies && doExclusionForces)
        {
            if (nbSci.shift == CENTRAL && a_plistCJ4[cij4Start].cj[0] == sci * c_nbnxnGpuNumClusterPerSupercluster)
            {
                // we have the diagonal: add the charge and LJ self interaction energy term
                for (int i = 0; i < c_nbnxnGpuNumClusterPerSupercluster; i++)
                {
                    // TODO: Are there other options?
                    if constexpr (props.elecEwald || props.elecRF || props.elecCutoff)
                    {
                        const float qi = sm_xq[i][tidxi][3];
                        energyElec += qi * qi;
                    }
                    if constexpr (props.vdwEwald)
                    {
                        energyVdw +=
                                a_nbfp[a_atomTypes[(sci * c_nbnxnGpuNumClusterPerSupercluster + i) * c_clSize + tidxi]
                                       * (numTypes + 1) * 2];
                    }
                }
                /* divide the self term(s) equally over the j-threads, then multiply with the coefficients. */
                if constexpr (props.vdwEwald)
                {
                    energyVdw /= c_clSize;
                    energyVdw *= 0.5F * c_oneSixth * ewaldCoeffLJ_6_6; // c_OneTwelfth?
                }
                if constexpr (props.elecRF || props.elecCutoff)
                {
                    // Correct for epsfac^2 due to adding qi^2 */
                    energyElec /= epsFac * c_clSize;
                    energyElec *= -0.5F * cRF;
                }
                if constexpr (props.elecEwald)
                {
                    // Correct for epsfac^2 due to adding qi^2 */
                    energyElec /= epsFac * c_clSize;
                    energyElec *= -ewaldBeta * c_OneOverSqrtPi; /* last factor 1/sqrt(pi) */
                }
            } // (nbSci.shift == CENTRAL && a_plistCJ4[cij4Start].cj[0] == sci * c_nbnxnGpuNumClusterPerSupercluster)
        }     // (doCalcEnergies && doExclusionForces)

        // Only needed if (doExclusionForces)
        const bool nonSelfInteraction = !(nbSci.shift == CENTRAL & tidxj <= tidxi);

        // loop over the j clusters = seen by any of the atoms in the current super-cluster
        for (int j4 = cij4Start + tidxz; j4 < cij4End; j4 += 1)
        {
            unsigned imask = a_plistCJ4[j4].imei[widx].imask;
            if (!doPruneNBL && !imask)
            {
                continue;
            }
            const int wexclIdx = a_plistCJ4[j4].imei[widx].excl_ind;
            const unsigned wexcl = a_plistExcl[wexclIdx].pair[tidx & (subGroupSize - 1)]; // sg.get_local_linear_id()
            for (int jm = 0; jm < c_nbnxnGpuJgroupSize; jm++)
            {
                const bool maskSet =
                        imask & (superClInteractionMask << (jm * c_nbnxnGpuNumClusterPerSupercluster));
                if (!maskSet)
                {
                    continue;
                }
                unsigned  maskJI = (1U << (jm * c_nbnxnGpuNumClusterPerSupercluster));
                const int cj     = a_plistCJ4[j4].cj[jm];
                const int aj     = cj * c_clSize + tidxj;

                // load j atom data
                const float4 xqj = a_xq[aj];

                const float3 xj(xqj[0], xqj[1], xqj[2]);
                const float  qj = xqj[3];
                int          atomTypeJ; // Only needed if (!props.vdwComb)
                float2       ljCombJ;   // Only needed if (props.vdwComb)
                if constexpr (props.vdwComb)
                {
                    ljCombJ = a_ljComb[aj];
                }
                else
                {
                    atomTypeJ = a_atomTypes[aj];
                }

                float3 fCjBuf(0.0F, 0.0F, 0.0F);

                for (int i = 0; i < c_nbnxnGpuNumClusterPerSupercluster; i++)
                {
                    if (imask & maskJI)
                    {
                        // i cluster index
                        const int ci = sci * c_nbnxnGpuNumClusterPerSupercluster + i;
                        // all threads load an atom from i cluster ci into shmem!
                        const float4 xqi = sm_xq[i][tidxi];
                        const float3 xi(xqi[0], xqi[1], xqi[2]);

                        // distance between i and j atoms
                        const float3 rv = xi - xj;
                        float        r2 = norm2(rv);

                        if constexpr (doPruneNBL)
                        {
                            /* If _none_ of the atoms pairs are in cutoff range,
                             * the bit corresponding to the current
                             * cluster-pair in imask gets set to 0. */
                            if (!sycl_2020::group_any_of(sg, r2 < rlistOuterSq))
                            {
                                imask &= ~maskJI;
                            }
                        }
                        const float pairExclMask = (wexcl & maskJI) ? 1.0F : 0.0F;

                        // cutoff & exclusion check

                        const bool notExcluded = doExclusionForces ? (nonSelfInteraction | (ci != cj))
                                                                   : (wexcl & maskJI);

                        // SYCL-TODO: Check optimal way of branching here.
                        if ((r2 < rCoulombSq) && notExcluded)
                        {
                            const float qi = xqi[3];
                            int         atomTypeI; // Only needed if (!props.vdwComb)
                            float       c6, c12, sigma, epsilon;

                            if constexpr (!props.vdwComb)
                            {
                                /* LJ 6*C6 and 12*C12 */
                                atomTypeI     = sm_atomTypeI[i][tidxi];
                                const int idx = (numTypes * atomTypeI + atomTypeJ) * 2;
                                c6            = a_nbfp[idx]; // TODO: Make a_nbfm into float2
                                c12           = a_nbfp[idx + 1];
                            }
                            else
                            {
                                const float2 ljCombI = sm_ljCombI[i][tidxi];
                                if constexpr (props.vdwCombGeom)
                                {
                                    c6  = ljCombI[0] * ljCombJ[0];
                                    c12 = ljCombI[1] * ljCombJ[1];
                                }
                                else
                                {
                                    static_assert(props.vdwCombLB);
                                    // LJ 2^(1/6)*sigma and 12*epsilon
                                    sigma   = ljCombI[0] + ljCombJ[0];
                                    epsilon = ljCombI[1] * ljCombJ[1];
                                    if constexpr (doCalcEnergies)
                                    {
                                        convertSigmaEpsilonToC6C12(sigma, epsilon, &c6, &c12);
                                    }
                                } // props.vdwCombGeom
                            }     // !props.vdwComb

                            // Ensure distance do not become so small that r^-12 overflows
                            r2 = std::max(r2, c_nbnxnMinDistanceSquared);
                            // SYCL-TODO: sycl::half_precision::rsqrt?
                            const float rInv  = cl::sycl::native::rsqrt(r2);
                            const float r2Inv = rInv * rInv;
                            float       r6Inv, fInvR, energyLJPair;
                            if constexpr (!props.vdwCombLB || doCalcEnergies)
                            {
                                r6Inv = r2Inv * r2Inv * r2Inv;
                                if constexpr (doExclusionForces)
                                {
                                    // SYCL-TODO: Check if true for SYCL
                                    /* We could mask r2Inv, but with Ewald masking both
                                     * r6Inv and fInvR is faster */
                                    r6Inv *= pairExclMask;
                                }
                                fInvR = r6Inv * (c12 * r6Inv - c6) * r2Inv;
                            }
                            else
                            {
                                float sig_r  = sigma * rInv;
                                float sig_r2 = sig_r * sig_r;
                                float sig_r6 = sig_r2 * sig_r2 * sig_r2;
                                if constexpr (doExclusionForces)
                                {
                                    sig_r6 *= pairExclMask;
                                }
                                fInvR = epsilon * sig_r6 * (sig_r6 - 1.0F) * r2Inv;
                            } // (!props.vdwCombLB || doCalcEnergies)
                            if constexpr (doCalcEnergies || props.vdwPSwitch)
                            {
                                energyLJPair = pairExclMask
                                               * (c12 * (r6Inv * r6Inv + repulsionShift.cpot) * c_oneTwelfth
                                                  - c6 * (r6Inv + dispersionShift.cpot) * c_oneSixth);
                            }
                            if constexpr (props.vdwFSwitch)
                            {
                                ljForceSwitch<doCalcEnergies>(
                                        dispersionShift, repulsionShift, rVdwSwitch, c6, c12, rInv, r2, &fInvR, &energyLJPair);
                            }
                            if constexpr (props.vdwEwald)
                            {
                                ljEwaldComb<doCalcEnergies, vdwType>(a_nbfpComb,
                                                                     ljEwaldShift,
                                                                     atomTypeI,
                                                                     atomTypeJ,
                                                                     r2,
                                                                     r2Inv,
                                                                     ewaldCoeffLJ_2,
                                                                     ewaldCoeffLJ_6_6,
                                                                     pairExclMask,
                                                                     &fInvR,
                                                                     &energyLJPair);
                            } // (props.vdwEwald)
                            if constexpr (props.vdwPSwitch)
                            {
                                ljPotentialSwitch<doCalcEnergies>(
                                        vdwSwitch, rVdwSwitch, rInv, r2, &fInvR, &energyLJPair);
                            }
                            if constexpr (props.elecEwaldTwin)
                            {
                                // Separate VDW cut-off check to enable twin-range cut-offs
                                // (rVdw < rCoulomb <= rList)
                                const float vdwInRange = (r2 < rVdwSq) ? 1.0F : 0.0F;
                                fInvR *= vdwInRange;
                                if constexpr (doCalcEnergies)
                                {
                                    energyLJPair *= vdwInRange;
                                }
                            }
                            if constexpr (doCalcEnergies)
                            {
                                energyVdw += energyLJPair;
                            }

                            if constexpr (props.elecCutoff)
                            {
                                if constexpr (doExclusionForces)
                                {
                                    fInvR += qi * qj * pairExclMask * r2Inv * rInv;
                                }
                                else
                                {
                                    fInvR += qi * qj * r2Inv * rInv;
                                }
                            }
                            if constexpr (props.elecRF)
                            {
                                fInvR += qi * qj * (pairExclMask * r2Inv * rInv - twoKRf);
                            }
                            if constexpr (props.elecEwaldAna)
                            {
                                fInvR += qi * qj
                                         * (pairExclMask * r2Inv * rInv + pmeCorrF(beta2 * r2) * beta3);
                            }
                            if constexpr (props.elecEwaldTab)
                            {
                                fInvR += qi * qj
                                         * (pairExclMask * r2Inv
                                            - interpolateCoulombForceR(
                                                      a_coulombTab, coulombTabScale, r2 * rInv))
                                         * rInv;
                            }

                            if constexpr (doCalcEnergies)
                            {
                                if constexpr (props.elecCutoff)
                                {
                                    energyElec += qi * qj * (pairExclMask * rInv - cRF);
                                }
                                if constexpr (props.elecRF)
                                {
                                    energyElec +=
                                            qi * qj * (pairExclMask * rInv + 0.5f * twoKRf * r2 - cRF);
                                }
                                if constexpr (props.elecEwald)
                                {
                                    energyElec +=
                                            qi * qj
                                            * (rInv * (pairExclMask - cl::sycl::erf(r2 * rInv * ewaldBeta))
                                               - pairExclMask * ewaldShift);
                                }
                            }

                            const float3 forceIJ = rv * fInvR;

                            /* accumulate j forces in registers */
                            fCjBuf -= forceIJ;
                            /* accumulate i forces in registers */
                            fCiBuf[i] += forceIJ;
                        } // (r2 < rCoulombSq) && notExcluded
                    }     // (imask & maskJI)
                    /* shift the mask bit by 1 */
                    maskJI += maskJI;
                } // for (int i = 0; i < c_nbnxnGpuNumClusterPerSupercluster; i++)
                /* reduce j forces */
                reduceForceJShuffle(fCjBuf, itemIdx, tidxi, aj, a_f);
            } // for (int jm = 0; jm < c_nbnxnGpuJgroupSize; jm++)
            if constexpr (doPruneNBL)
            {
                /* Update the imask with the new one which does not contain the
                 * out of range clusters anymore. */
                a_plistCJ4[j4].imei[widx].imask = imask;
            }
        } // for (int j4 = cij4Start; j4 < cij4End; j4 += 1)

        /* skip central shifts when summing shift forces */
        const bool doCalcShift = (calcShift && !(nbSci.shift == CENTRAL));

        reduceForceIAndFShift(
                sm_reductionBuffer, fCiBuf, doCalcShift, itemIdx, tidxi, tidxj, sci, nbSci.shift, a_f, a_fShift);

        if constexpr (doCalcEnergies)
        {
            const float energyVdwGroup = sycl_2020::group_reduce(
                    itemIdx.get_group(), energyVdw, 0.0F, sycl_2020::plus<float>());
            const float energyElecGroup = sycl_2020::group_reduce(
                    itemIdx.get_group(), energyElec, 0.0F, sycl_2020::plus<float>());

            if (tidx == 0)
            {
                atomicFetchAdd(a_energyVdw, 0, energyVdwGroup);
                atomicFetchAdd(a_energyElec, 0, energyElecGroup);
            }
        }
    };
}

// SYCL 1.2.1 requires providing a unique type for a kernel. Should not be needed for SYCL2020.
template<bool doPruneNBL, bool doCalcEnergies, enum ElecType elecType, enum VdwType vdwType>
class NbnxmKernelName;

template<bool doPruneNBL, bool doCalcEnergies, enum ElecType elecType, enum VdwType vdwType, class... Args>
cl::sycl::event launchNbnxmKernel(const DeviceStream& deviceStream, const int numSci, Args&&... args)
{
    // Should not be needed for SYCL2020.
    using kernelNameType = NbnxmKernelName<doPruneNBL, doCalcEnergies, elecType, vdwType>;

    /* Kernel launch config:
     * - The thread block dimensions match the size of i-clusters, j-clusters,
     *   and j-cluster concurrency, in x, y, and z, respectively.
     * - The 1D block-grid contains as many blocks as super-clusters.
     */
    const int                   numBlocks = numSci;
    const cl::sycl::range<3>    blockSize{ c_clSize, c_clSize, 1 };
    const cl::sycl::range<3>    globalSize{ numBlocks * blockSize[0], blockSize[1], blockSize[2] };
    const cl::sycl::nd_range<3> range{ globalSize, blockSize };

    cl::sycl::queue q = deviceStream.stream();

    cl::sycl::event e = q.submit([&](cl::sycl::handler& cgh) {
        auto kernel = nbnxmKernel<doPruneNBL, doCalcEnergies, elecType, vdwType>(
                cgh, std::forward<Args>(args)...);
        cgh.parallel_for<kernelNameType>(flattenNDRange(range), kernel);
    });

    return e;
}

template<class... Args>
cl::sycl::event chooseAndLaunchNbnxmKernel(bool          doPruneNBL,
                                           bool          doCalcEnergies,
                                           enum ElecType elecType,
                                           enum VdwType  vdwType,
                                           Args&&... args)
{
    return gmx::dispatchTemplatedFunction(
            [&](auto doPruneNBL_, auto doCalcEnergies_, auto elecType_, auto vdwType_) {
                return launchNbnxmKernel<doPruneNBL_, doCalcEnergies_, elecType_, vdwType_>(
                        std::forward<Args>(args)...);
            },
            doPruneNBL,
            doCalcEnergies,
            elecType,
            vdwType);
}

void launchNbnxmKernel(NbnxmGpu* nb, const gmx::StepWorkload& stepWork, const InteractionLocality iloc)
{
    sycl_atomdata_t*    adat         = nb->atdat;
    NBParamGpu*         nbp          = nb->nbparam;
    gpu_plist*          plist        = nb->plist[iloc];
    const bool          doPruneNBL   = (plist->haveFreshList && !nb->didPrune[iloc]);
    const DeviceStream& deviceStream = *nb->deviceStreams[iloc];

    // Casting to float simplifies using atomic ops in the kernel
    cl::sycl::buffer<float3, 1> f(*adat->f.buffer_);
    auto                        fAsFloat = f.reinterpret<float, 1>(f.get_count() * DIM);
    cl::sycl::buffer<float3, 1> fShift(*adat->fShift.buffer_);
    auto fShiftAsFloat = fShift.reinterpret<float, 1>(fShift.get_count() * DIM);

    cl::sycl::event e = chooseAndLaunchNbnxmKernel(doPruneNBL,
                                                   stepWork.computeEnergy,
                                                   nbp->elecType,
                                                   nbp->vdwType,
                                                   deviceStream,
                                                   plist->nsci,
                                                   adat->xq,
                                                   fAsFloat,
                                                   adat->shiftVec,
                                                   fShiftAsFloat,
                                                   adat->eElec,
                                                   adat->eLJ,
                                                   plist->cj4,
                                                   plist->sci,
                                                   plist->excl,
                                                   adat->ljComb,
                                                   adat->atomTypes,
                                                   nbp->nbfp,
                                                   nbp->nbfp_comb,
                                                   nbp->coulomb_tab,
                                                   adat->numTypes,
                                                   nbp->rcoulomb_sq,
                                                   nbp->rvdw_sq,
                                                   nbp->two_k_rf,
                                                   nbp->ewald_beta,
                                                   nbp->rlistOuter_sq,
                                                   nbp->sh_ewald,
                                                   nbp->epsfac,
                                                   nbp->ewaldcoeff_lj,
                                                   nbp->c_rf,
                                                   nbp->dispersion_shift,
                                                   nbp->repulsion_shift,
                                                   nbp->vdw_switch,
                                                   nbp->rvdw_switch,
                                                   nbp->sh_lj_ewald,
                                                   nbp->coulomb_tab_scale,
                                                   stepWork.computeVirial);
}

} // namespace Nbnxm
