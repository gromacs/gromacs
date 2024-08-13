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

/*! \internal \file
 *  \brief
 *  NBNXM SYCL kernels
 *
 *  \ingroup module_nbnxm
 */
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/gpu_utils/packed_float.h"
#include "gromacs/math/functions.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/utility/template_mp.h"

#include "nbnxm_sycl_kernel.h"
#include "nbnxm_sycl_kernel_utils.h"
#include "nbnxm_sycl_types.h"

namespace gmx
{

//! \brief Class name for NBNXM kernel
template<bool doPruneNBL, bool doCalcEnergies, enum ElecType elecType, enum VdwType vdwType, int subGroupSize>
class NbnxmKernel;

/*! \brief Macro to control the enablement of manually-packed Float3 structure.
 *
 * If enabled (default), the explicit packed math will be used on devices where
 * it is known to be beneficial (currently, AMD MI250X / gfx90a).
 * If disabled, packed math will never be explicitly  used.
 *
 * This only controls the use of AmdPackedFloat3 datastructure, not the layout
 * of fCi buffer.
 *
 * See issue #4854 */
#define GMX_NBNXM_ENABLE_PACKED_FLOAT3 1


#if (GMX_NBNXM_ENABLE_PACKED_FLOAT3 && defined(__AMDGCN__) && defined(__gfx90a__))
using FCiFloat3 = AmdPackedFloat3;
#else
using FCiFloat3 = Float3;
#endif

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
            (elecType == ElecType::EwaldAnaTwin || elecType == ElecType::EwaldTabTwin); ///< Use twin cut-off.
    static constexpr bool elecEwald = (elecEwaldAna || elecEwaldTab);  ///< EL_EWALD_ANY
    static constexpr bool vdwCombLB = (vdwType == VdwType::CutCombLB); ///< LJ_COMB && !LJ_COMB_GEOM
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

template<enum ElecType elecType>
constexpr bool elecEwald = EnergyFunctionProperties<elecType, VdwType::Count>().elecEwald;

template<enum ElecType elecType>
constexpr bool elecEwaldTab = EnergyFunctionProperties<elecType, VdwType::Count>().elecEwaldTab;

template<enum VdwType vdwType>
constexpr bool ljEwald = EnergyFunctionProperties<ElecType::Count, vdwType>().vdwEwald;
//@}

/*! \brief Should we avoid FP atomics to the same location from the same work-group?
 *
 * Intel GPUs without native floating-point operations emulate them via CAS-loop,
 * which is very, very slow when two threads from the same group write to the same
 * global location. We don't specialize the kernels by vendor, so we use c_clSize == 4
 * as a proxy to detect such devices.
 */
constexpr bool c_avoidFloatingPointAtomics = (c_clSize == 4);

using sycl::access::fence_space;
using mode = sycl::access_mode;

//! \brief Convert \p sigma and \p epsilon VdW parameters to \c c6,c12 pair.
static inline Float2 convertSigmaEpsilonToC6C12(const float sigma, const float epsilon)
{
    const float sigma2 = sigma * sigma;
    const float sigma6 = sigma2 * sigma2 * sigma2;
    const float c6     = epsilon * sigma6;
    const float c12    = c6 * sigma6;

    return { c6, c12 };
}

//! \brief Calculate force and energy for a pair of atoms, VdW force-switch flavor.
template<bool doCalcEnergies>
static inline void ljForceSwitch(const shift_consts_t     dispersionShift,
                                 const shift_consts_t     repulsionShift,
                                 const float              rVdwSwitch,
                                 const float              c6,
                                 const float              c12,
                                 const float              rInv,
                                 const float              r2,
                                 sycl::private_ptr<float> fInvR,
                                 sycl::private_ptr<float> eLJ)
{
    /* force switch constants */
    const float dispShiftV2 = dispersionShift.c2;
    const float dispShiftV3 = dispersionShift.c3;
    const float repuShiftV2 = repulsionShift.c2;
    const float repuShiftV3 = repulsionShift.c3;

    const float r       = r2 * rInv;
    const float rSwitch = sycl::fdim(r, rVdwSwitch); // max(r - rVdwSwitch, 0)

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
template<enum VdwType vdwType>
static inline float calculateLJEwaldC6Grid(const sycl::global_ptr<const Float2> a_nbfpComb,
                                           const int                            typeI,
                                           const int                            typeJ)
{
    if constexpr (vdwType == VdwType::EwaldGeom)
    {
        return a_nbfpComb[typeI][0] * a_nbfpComb[typeJ][0];
    }
    else
    {
        static_assert(vdwType == VdwType::EwaldLB);
        /* sigma and epsilon are scaled to give 6*C6 */
        const Float2 c6c12_i = a_nbfpComb[typeI];
        const Float2 c6c12_j = a_nbfpComb[typeJ];

        const float sigma   = c6c12_i[0] + c6c12_j[0];
        const float epsilon = c6c12_i[1] * c6c12_j[1];

        const float sigma2 = sigma * sigma;
        return epsilon * sigma2 * sigma2 * sigma2;
    }
}

//! Calculate LJ-PME grid force contribution with geometric or LB combination rule.
template<bool doCalcEnergies, enum VdwType vdwType>
static inline void ljEwaldComb(const sycl::global_ptr<const Float2> a_nbfpComb,
                               const float                          sh_lj_ewald,
                               const int                            typeI,
                               const int                            typeJ,
                               const float                          r2,
                               const float                          r2Inv,
                               const float                          lje_coeff2,
                               const float                          lje_coeff6_6,
                               const float                          int_bit,
                               sycl::private_ptr<float>             fInvR,
                               sycl::private_ptr<float>             eLJ)
{
    const float c6grid = calculateLJEwaldC6Grid<vdwType>(a_nbfpComb, typeI, typeJ);

    /* Recalculate inv_r6 without exclusion mask */
    const float inv_r6_nm = r2Inv * r2Inv * r2Inv;
    const float cr2       = lje_coeff2 * r2;
    const float expmcr2   = sycl::exp(-cr2);
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
static inline void ljPotentialSwitch(const switch_consts_t    vdwSwitch,
                                     const float              rVdwSwitch,
                                     const float              rInv,
                                     const float              r2,
                                     sycl::private_ptr<float> fInvR,
                                     sycl::private_ptr<float> eLJ)
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
    return sycl::fma(t, d1, sycl::fma(-t, d0, d0));
}

/*! \brief Interpolate Ewald coulomb force correction using the F*r table. */
static inline float interpolateCoulombForceR(const sycl::global_ptr<const float> a_coulombTab,
                                             const float                         coulombTabScale,
                                             const float                         r)
{
    const float normalized = coulombTabScale * r;
    const int   index      = static_cast<int>(normalized);
#if defined(__SYCL_DEVICE_ONLY__) && defined(__AMDGCN__)
    // TODO: up to ROCm v5.3 compiler does not do this transformation. Remove when this is no longer the case.
    const float fraction = __builtin_amdgcn_fractf(normalized);
#else
    const float fraction = normalized - index;
#endif

    const float left  = a_coulombTab[index];
    const float right = a_coulombTab[index + 1];

    return lerp(left, right, fraction); // TODO: sycl::mix
}

#if defined(__SYCL_DEVICE_ONLY__) && defined(__AMDGCN__)
/*! \brief Reduce c_clSize j-force components using AMD DPP instruction and atomically accumulate into a_f.
 *
 * c_clSize consecutive threads hold the force components of a j-atom which we
 * reduced in log2(cl_Size) steps using shift and atomically accumulate them into \p a_f.
 *
 * Note: This causes massive amount of spills with the tabulated kernel on gfx803 using ROCm 5.3.
 * We don't disable it only for the tabulated kernel as the analytical is the default anyway.
 */
static inline void reduceForceJAmdDpp(Float3 f, const int tidxi, const int aidx, sycl::global_ptr<Float3> a_f)
{
    static_assert(c_clSize == 8);

    f[0] += amdDppUpdateShfl<float, /* row_shl:1 */ 0x101>(f[0]);
    f[1] += amdDppUpdateShfl<float, /* row_shr:1 */ 0x111>(f[1]);
    f[2] += amdDppUpdateShfl<float, /* row_shl:1 */ 0x101>(f[2]);

    if (tidxi & 1)
    {
        f[0] = f[1];
    }

    f[0] += amdDppUpdateShfl<float, /* row_shl:2 */ 0x102>(f[0]);
    f[2] += amdDppUpdateShfl<float, /* row_shr:2 */ 0x112>(f[2]);

    if (tidxi & 2)
    {
        f[0] = f[2];
    }

    f[0] += amdDppUpdateShfl<float, /* row_shl:4 */ 0x104>(f[0]);

    if (tidxi < 3)
    {
        atomicFetchAdd(a_f[aidx][tidxi], f[0]);
    }
}
#endif

/*! \brief Reduce c_clSize j-force components using shifts and atomically accumulate into a_f.
 *
 * c_clSize consecutive threads hold the force components of a j-atom which we
 * reduced in log2(cl_Size) steps using shift and atomically accumulate them into \p a_f.
 */
static inline void reduceForceJShuffle(Float3                   f,
                                       const sycl::nd_item<3>&  itemIdx,
                                       const int                tidxi,
                                       const int                aidx,
                                       sycl::global_ptr<Float3> a_f)
{
    static_assert(c_clSize == 8 || c_clSize == 4);
    sycl::sub_group sg = itemIdx.get_sub_group();

    f[0] += sycl::shift_group_left(sg, f[0], 1);
    f[1] += sycl::shift_group_right(sg, f[1], 1);
    f[2] += sycl::shift_group_left(sg, f[2], 1);
    if (tidxi & 1)
    {
        f[0] = f[1];
    }

    f[0] += sycl::shift_group_left(sg, f[0], 2);
    f[2] += sycl::shift_group_right(sg, f[2], 2);
    if (tidxi & 2)
    {
        f[0] = f[2];
    }

    if constexpr (c_clSize == 8)
    {
        f[0] += sycl::shift_group_left(sg, f[0], 4);
    }

    if (tidxi < 3)
    {
        atomicFetchAdd(a_f[aidx][tidxi], f[0]);
    }
}

/*!
 * \brief Do workgroup-level reduction of a single \c float.
 *
 * While SYCL has \c sycl::reduce_over_group, it currently (oneAPI 2021.3.0) uses a very large
 * shared memory buffer, which leads to a reduced occupancy.
 *
 * \note The caller must make sure there are no races when reusing the \p sm_buf.
 *
 * \tparam subGroupSize Size of a sub-group.
 * \tparam groupSize Size of a work-group.
 * \param itemIdx Current thread's \c sycl::nd_item.
 * \param tidxi Current thread's linearized local index.
 * \param sm_buf Accessor for local reduction buffer.
 * \param valueToReduce Current thread's value. Must have length of at least 1.
 * \return For thread with \p tidxi 0: sum of all \p valueToReduce. Other threads: unspecified.
 */
template<int subGroupSize, int groupSize>
static inline float groupReduce(const sycl::nd_item<3> itemIdx,
                                const unsigned int     tidxi,
                                sycl::local_ptr<float> sm_buf,
                                float                  valueToReduce)
{
    constexpr int numSubGroupsInGroup = groupSize / subGroupSize;
    static_assert(numSubGroupsInGroup == 1 || numSubGroupsInGroup == 2);
    sycl::sub_group sg = itemIdx.get_sub_group();
    valueToReduce      = sycl::reduce_over_group(sg, valueToReduce, sycl::plus<float>());
    // If we have two sub-groups, we should reduce across them.
    if constexpr (numSubGroupsInGroup == 2)
    {
        if (tidxi == subGroupSize)
        {
            sm_buf[0] = valueToReduce;
        }
        itemIdx.barrier(fence_space::local_space);
        if (tidxi == 0)
        {
            valueToReduce += sm_buf[0];
        }
    }
    return valueToReduce;
}

/*! \brief Reduce c_clSize j-force components using local memory and atomically accumulate into a_f.
 *
 * c_clSize consecutive threads hold the force components of a j-atom which we
 * reduced in cl_Size steps using shift and atomically accumulate them into \p a_f.
 *
 * TODO: implement binary reduction flavor for the case where cl_Size is power of two.
 */
static inline void reduceForceJGeneric(sycl::local_ptr<float>   sm_buf,
                                       Float3                   f,
                                       const sycl::nd_item<3>&  itemIdx,
                                       const int                tidxi,
                                       const int                tidxj,
                                       const int                aidx,
                                       sycl::global_ptr<Float3> a_f)
{
    static constexpr int sc_fBufferStride = c_clSizeSq;
    int                  tidx             = tidxi + tidxj * c_clSize;
    sm_buf[0 * sc_fBufferStride + tidx]   = f[0];
    sm_buf[1 * sc_fBufferStride + tidx]   = f[1];
    sm_buf[2 * sc_fBufferStride + tidx]   = f[2];

    subGroupBarrier(itemIdx);

    // reducing data 8-by-by elements on the leader of same threads as those storing above
    SYCL_ASSERT(itemIdx.get_sub_group().get_max_local_range()[0] >= c_clSize);

    if (tidxi < 3)
    {
        float fSum = 0.0F;
        for (int j = tidxj * c_clSize; j < (tidxj + 1) * c_clSize; j++)
        {
            fSum += sm_buf[sc_fBufferStride * tidxi + j];
        }

        atomicFetchAdd(a_f[aidx][tidxi], fSum);
    }
}


/*! \brief Reduce c_clSize j-force components using either shifts or local memory and atomically accumulate into a_f.
 */
template<bool useShuffleReduction>
static inline void reduceForceJ(sycl::local_ptr<float>   sm_buf,
                                Float3                   f,
                                const sycl::nd_item<3>   itemIdx,
                                const int                tidxi,
                                const int                tidxj,
                                const int                aidx,
                                sycl::global_ptr<Float3> a_f)
{
    if constexpr (!useShuffleReduction)
    {
        reduceForceJGeneric(sm_buf, f, itemIdx, tidxi, tidxj, aidx, a_f);
    }
    else
    {
#if defined(__SYCL_DEVICE_ONLY__) && defined(__AMDGCN__)
        reduceForceJAmdDpp(f, tidxi, aidx, a_f);
#else
        reduceForceJShuffle(f, itemIdx, tidxi, aidx, a_f);
#endif
    }
}

/*! \brief Local memory-based i-force reduction.
 *
 * Note that this reduction is unoptimized and some of the barrier synchronization
 * used could be avoided on >=8-wide architectures.
 */
template<typename FCiBufferWrapperX, typename FCiBufferWrapperY, typename FCiBufferWrapperZ>
static inline void reduceForceIAndFShiftGeneric(sycl::local_ptr<float>   sm_buf,
                                                const FCiBufferWrapperX& fCiBufX,
                                                const FCiBufferWrapperY& fCiBufY,
                                                const FCiBufferWrapperZ& fCiBufZ,
                                                const bool               calcFShift,
                                                const sycl::nd_item<3>   itemIdx,
                                                const int                tidxi,
                                                const int                tidxj,
                                                const int                sci,
                                                const int                shift,
                                                sycl::global_ptr<Float3> a_f,
                                                sycl::global_ptr<Float3> a_fShift)
{
    static constexpr int bufStride  = c_clSize * c_clSize;
    static constexpr int clSizeLog2 = gmx::StaticLog2<c_clSize>::value;
    const int            tidx       = tidxi + tidxj * c_clSize;
    float                fShiftBuf  = 0.0F;
#pragma unroll c_nbnxnGpuNumClusterPerSupercluster
    for (int ciOffset = 0; ciOffset < c_nbnxnGpuNumClusterPerSupercluster; ciOffset++)
    {
        const int aidx = (sci * c_nbnxnGpuNumClusterPerSupercluster + ciOffset) * c_clSize + tidxi;
        // Store i-forces in local memory
        sm_buf[tidx]                 = fCiBufX(ciOffset);
        sm_buf[bufStride + tidx]     = fCiBufY(ciOffset);
        sm_buf[2 * bufStride + tidx] = fCiBufZ(ciOffset);
        itemIdx.barrier(fence_space::local_space);

        // Reduce the initial c_clSize values for each i atom to half every step by using c_clSize * i threads.
        int i = c_clSize / 2;
        for (int j = clSizeLog2 - 1; j > 0; j--)
        {
            if (tidxj < i)
            {
                sm_buf[tidx] += sm_buf[tidx + i * c_clSize];
                sm_buf[bufStride + tidx] += sm_buf[bufStride + tidx + i * c_clSize];
                sm_buf[2 * bufStride + tidx] += sm_buf[2 * bufStride + tidx + i * c_clSize];
            }
            i >>= 1;
            itemIdx.barrier(fence_space::local_space);
        }

        /* i == 1, last reduction step, combined with writing to global mem.
         * Split the reduction between the first 3 threads in a warp.
         * Threads with lane id 0 will do the reduction for X components, 1 will do Y etc.
         * */
        if (tidxj < 3)
        {
            const float f =
                    sm_buf[tidxj * bufStride + tidxi] + sm_buf[tidxj * bufStride + c_clSize + tidxi];
            atomicFetchAdd(a_f[aidx][tidxj], f);
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
            if constexpr (c_avoidFloatingPointAtomics)
            {
                /* Intel Xe (Gen12LP) and earlier GPUs implement floating-point atomics via
                 * a compare-and-swap (CAS) loop. It has particularly poor performance when
                 * updating the same memory location from the same work-group.
                 * Such optimization might be slightly beneficial for NVIDIA and AMD as well,
                 * but it is unlikely to make a big difference and thus was not evaluated.
                 */
                auto sg = itemIdx.get_sub_group();
                fShiftBuf += sycl::shift_group_left(sg, fShiftBuf, 1);
                fShiftBuf += sycl::shift_group_left(sg, fShiftBuf, 2);
                if (tidxi == 0)
                {
                    atomicFetchAdd(a_fShift[shift][tidxj], fShiftBuf);
                }
            }
            else
            {
                atomicFetchAdd(a_fShift[shift][tidxj], fShiftBuf);
            }
        }
    }
}

/*! \brief Shuffle-based i-force reduction.
 *
 * We need to reduce c_clSize values spaced c_clSize threads apart (hardware threads are consecutive
 * for \c tidxi, have stride c_clSize for \c tidxj).
 *
 * We can have up to three reduction steps done with shuffles:
 *
 * One step (e.g, Intel iGPU, c_clSize == 4, subGroupSize == 8): handled in a separate
 * specialization.
 * Two steps (e.g., NVIDIA, c_clSize == 8, subGroupSize == 32): after two shuffle reduction steps,
 * we do atomicFetchAdd from each sub-group.
 * Three steps (e.g., AMD CDNA, c_clSize == 8, subGroupSize == 64): similar to the two-step
 * approach, but we have two times less atomicFetchAdd's.
 */
template<int numShuffleReductionSteps, typename FCiBufferWrapperX, typename FCiBufferWrapperY, typename FCiBufferWrapperZ>
typename std::enable_if_t<numShuffleReductionSteps != 1, void> static inline reduceForceIAndFShiftShuffles(
        const FCiBufferWrapperX& fCiBufX,
        const FCiBufferWrapperY& fCiBufY,
        const FCiBufferWrapperZ& fCiBufZ,
        const bool               calcFShift,
        const sycl::nd_item<3>   itemIdx,
        const int                tidxi,
        const int                tidxj,
        const int                sci,
        const int                shift,
        sycl::global_ptr<Float3> a_f,
        sycl::global_ptr<Float3> a_fShift)
{
    const sycl::sub_group sg = itemIdx.get_sub_group();
    static_assert(numShuffleReductionSteps == 2 || numShuffleReductionSteps == 3);
    SYCL_ASSERT(sg.get_max_local_range()[0] >= 4 * c_clSize
                && "Subgroup too small for two-step shuffle reduction, use 1-step");

    float fShiftBuf = 0.0F;

#if defined(__SYCL_DEVICE_ONLY__) && defined(__AMDGCN__) && (__AMDGCN_WAVEFRONT_SIZE == 64)
    // Use AMD's cross-lane DPP reduction only for 64-wide exec
    // can't use static_assert because 32-wide compiler passes will trip on it
    SYCL_ASSERT(numShuffleReductionSteps == 3);
#    pragma unroll c_nbnxnGpuNumClusterPerSupercluster
    for (int ciOffset = 0; ciOffset < c_nbnxnGpuNumClusterPerSupercluster; ciOffset++)
    {
        const int aidx = (sci * c_nbnxnGpuNumClusterPerSupercluster + ciOffset) * c_clSize + tidxj;
        float     fx   = fCiBufX(ciOffset);
        float     fy   = fCiBufY(ciOffset);
        float     fz   = fCiBufZ(ciOffset);

        // Transpose values so DPP-based reduction can be used later
        fx = sycl::select_from_group(sg, fx, tidxi * c_clSize + tidxj);
        fy = sycl::select_from_group(sg, fy, tidxi * c_clSize + tidxj);
        fz = sycl::select_from_group(sg, fz, tidxi * c_clSize + tidxj);

        fx += amdDppUpdateShfl<float, /* row_shl:1 */ 0x101>(fx);
        fy += amdDppUpdateShfl<float, /* row_shr:1 */ 0x111>(fy);
        fz += amdDppUpdateShfl<float, /* row_shl:1 */ 0x101>(fz);

        if (tidxi & 1)
        {
            fx = fy;
        }

        fx += amdDppUpdateShfl<float, /* row_shl:2 */ 0x102>(fx);
        fz += amdDppUpdateShfl<float, /* row_shr:2 */ 0x112>(fz);

        if (tidxi & 2)
        {
            fx = fz;
        }

        fx += amdDppUpdateShfl<float, /* row_shl:4 */ 0x104>(fx);

        // Threads 0,1,2 increment X, Y, Z for their sub-groups
        if (tidxi < 3)
        {
            atomicFetchAdd(a_f[aidx][(tidxi)], fx);

            if (calcFShift)
            {
                fShiftBuf += fx;
            }
        }
    }
    /* add up local shift forces into global mem */
    if (calcFShift)
    {
        if ((tidxi) < 3)
        {
            atomicFetchAdd(a_fShift[shift][(tidxi)], fShiftBuf);
        }
    }

#else // GMX_SYCL_HIPSYCL && HIPSYCL_LIBKERNEL_IS_DEVICE_PASS_HIP

    // Thread mask to use to select first three threads (in tidxj) in each reduction "tree".
    // Two bits for two steps, three bits for three steps.
    constexpr int threadBitMask = (1U << numShuffleReductionSteps) - 1;

#    pragma unroll c_nbnxnGpuNumClusterPerSupercluster
    for (int ciOffset = 0; ciOffset < c_nbnxnGpuNumClusterPerSupercluster; ciOffset++)
    {
        const int aidx = (sci * c_nbnxnGpuNumClusterPerSupercluster + ciOffset) * c_clSize + tidxi;
        float     fx   = fCiBufX(ciOffset);
        float     fy   = fCiBufY(ciOffset);
        float     fz   = fCiBufZ(ciOffset);

        // First reduction step
        fx += sycl::shift_group_left(sg, fx, c_clSize);
        fy += sycl::shift_group_right(sg, fy, c_clSize);
        fz += sycl::shift_group_left(sg, fz, c_clSize);
        if (tidxj & 1)
        {
            fx = fy;
        }
        // Second reduction step
        fx += sycl::shift_group_left(sg, fx, 2 * c_clSize);
        fz += sycl::shift_group_right(sg, fz, 2 * c_clSize);
        if (tidxj & 2)
        {
            fx = fz;
        }
        // Third reduction step if possible
        if constexpr (numShuffleReductionSteps == 3)
        {
            fx += sycl::shift_group_left(sg, fx, 4 * c_clSize);
        }
        // Threads 0,1,2 (and 4,5,6 in case of numShuffleReductionSteps == 2) increment X, Y, Z for their sub-groups
        if ((tidxj & threadBitMask) < 3)
        {
            atomicFetchAdd(a_f[aidx][(tidxj & threadBitMask)], fx);

            if (calcFShift)
            {
                fShiftBuf += fx;
            }
        }
    }
    /* add up local shift forces into global mem */
    if (calcFShift)
    {
        if ((tidxj & threadBitMask) < 3)
        {
            atomicFetchAdd(a_fShift[shift][(tidxj & threadBitMask)], fShiftBuf);
        }
    }
#endif // GMX_SYCL_HIPSYCL && HIPSYCL_LIBKERNEL_IS_DEVICE_PASS_HIP
}

/*! \brief \c reduceForceIAndFShiftShuffles specialization for single-step reduction (e.g., Intel iGPUs).
 *
 * We have three components to reduce, but only one reduction step, so it is not possible
 * to gather the components in \c fx of different threads, like we do with 2 and more reduction steps.
 *
 * Therefore, first even threads handle X and odd threads handle Y; then, even threads additionally
 * handle Z. This also requires the second fShift buffer register.
 *
 * After one step of reduction using shuffles is complete, we use atomicAdd to accumulate the results
 * in the global memory. That causes a lot of atomic operations on a single memory location, which
 * is poorly handled by some hardware (e.g., Intel Gen9-11 and Xe LP). This can be remediated
 * by using local memory reduction after shuffles, but that's a TODO.
 */
template<int numShuffleReductionSteps, typename FCiBufferWrapperX, typename FCiBufferWrapperY, typename FCiBufferWrapperZ>
typename std::enable_if_t<numShuffleReductionSteps == 1, void> static inline reduceForceIAndFShiftShuffles(
        const FCiBufferWrapperX& fCiBufX,
        const FCiBufferWrapperY& fCiBufY,
        const FCiBufferWrapperZ& fCiBufZ,
        const bool               calcFShift,
        const sycl::nd_item<3>   itemIdx,
        const int                tidxi,
        const int                tidxj,
        const int                sci,
        const int                shift,
        sycl::global_ptr<Float3> a_f,
        sycl::global_ptr<Float3> a_fShift)
{
    const sycl::sub_group sg = itemIdx.get_sub_group();
    SYCL_ASSERT(sg.get_max_local_range()[0] >= 2 * c_clSize
                && "Subgroup too small even for 1-step shuffle reduction");
    SYCL_ASSERT(sg.get_max_local_range()[0] < 4 * c_clSize
                && "One-step shuffle reduction inefficient, use two-step version");
    float fShiftBufXY = 0.0F;
    float fShiftBufZ  = 0.0F;
#pragma unroll c_nbnxnGpuNumClusterPerSupercluster
    for (int ciOffset = 0; ciOffset < c_nbnxnGpuNumClusterPerSupercluster; ciOffset++)
    {
        const int aidx = (sci * c_nbnxnGpuNumClusterPerSupercluster + ciOffset) * c_clSize + tidxi;
        float     fx   = fCiBufX(ciOffset);
        float     fy   = fCiBufY(ciOffset);
        float     fz   = fCiBufZ(ciOffset);
        // First reduction step
        fx += sycl::shift_group_left(sg, fx, c_clSize);
        fy += sycl::shift_group_right(sg, fy, c_clSize);
        fz += sycl::shift_group_left(sg, fz, c_clSize);
        if (tidxj & 1)
        {
            fx = fy;
        }
        // Can not use shuffles to reduce further, do global atomics
        // Add the current X and Y values to the global buffer
        atomicFetchAdd(a_f[aidx][(tidxj & 1)], fx);
        if (calcFShift)
        {
            fShiftBufXY += fx;
        }
        // Threads tidxj == 0 or 2 increment Z
        if ((tidxj & 1) == 0)
        {
            atomicFetchAdd(a_f[aidx][2], fz);
            if (calcFShift)
            {
                fShiftBufZ += fz;
            }
        }
        subGroupBarrier(itemIdx);
    }
    /* add up local shift forces into global mem */
    if (calcFShift)
    {
        // Update X and Y by even and odd threads, respectively
        atomicFetchAdd(a_fShift[shift][tidxj & 1], fShiftBufXY);
        if ((tidxj & 1) == 0)
        {
            atomicFetchAdd(a_fShift[shift][2], fShiftBufZ);
        }
    }
}

/*! \brief Final i-force reduction.
 *
 * Reduce c_nbnxnGpuNumClusterPerSupercluster i-force components stored in \p fCiBuf[]
 * accumulating atomically into \p a_f.
 * If \p calcFShift is true, further reduce shift forces and atomically accumulate into \p a_fShift.
 *
 * This implementation works only with power of two array sizes.
 */
template<bool useShuffleReduction, int subGroupSize, typename FCiBufferWrapperX, typename FCiBufferWrapperY, typename FCiBufferWrapperZ>
static inline void reduceForceIAndFShift(sycl::local_ptr<float>   sm_buf,
                                         const FCiBufferWrapperX& fCiBufX,
                                         const FCiBufferWrapperY& fCiBufY,
                                         const FCiBufferWrapperZ& fCiBufZ,
                                         const bool               calcFShift,
                                         const sycl::nd_item<3>   itemIdx,
                                         const int                tidxi,
                                         const int                tidxj,
                                         const int                sci,
                                         const int                shift,
                                         sycl::global_ptr<Float3> a_f,
                                         sycl::global_ptr<Float3> a_fShift)
{
    // must have power of two elements in fCiBuf
    static_assert(gmx::isPowerOfTwo(c_nbnxnGpuNumClusterPerSupercluster));

    if constexpr (useShuffleReduction)
    {
        constexpr int numSteps = gmx::StaticLog2<subGroupSize / c_clSize>::value;
        static_assert(numSteps > 0 && numSteps <= 3,
                      "Invalid combination of sub-group size and cluster size");
        reduceForceIAndFShiftShuffles<numSteps>(
                fCiBufX, fCiBufY, fCiBufZ, calcFShift, itemIdx, tidxi, tidxj, sci, shift, a_f, a_fShift);
    }
    else
    {
        reduceForceIAndFShiftGeneric(
                sm_buf, fCiBufX, fCiBufY, fCiBufZ, calcFShift, itemIdx, tidxi, tidxj, sci, shift, a_f, a_fShift);
    }
}

/*! \brief Main kernel for NBNXM.
 *
 */
template<int subGroupSize, bool doPruneNBL, bool doCalcEnergies, enum ElecType elecType, enum VdwType vdwType>
static auto nbnxmKernel(sycl::handler& cgh,
                        const Float4* __restrict__ gm_xq,
                        Float3* __restrict__ gm_f,
                        const Float3* __restrict__ gm_shiftVec,
                        Float3* __restrict__ gm_fShift,
                        float* __restrict__ gm_energyElec,
                        float* __restrict__ gm_energyVdw,
                        nbnxn_cj_packed_t* __restrict__ gm_plistCJPacked,
                        const nbnxn_sci_t* __restrict__ gm_plistSci,
                        const nbnxn_excl_t* __restrict__ gm_plistExcl,
                        const Float2* __restrict__ gm_ljComb /* used iff ljComb<vdwType> */,
                        const int* __restrict__ gm_atomTypes /* used iff !ljComb<vdwType> */,
                        const Float2* __restrict__ gm_nbfp /* used iff !ljComb<vdwType> */,
                        const Float2* __restrict__ gm_nbfpComb /* used iff ljEwald<vdwType> */,
                        const float* __restrict__ gm_coulombTab /* used iff elecEwaldTab<elecType> */,
                        int* __restrict__ gm_sciHistogram,      /* used iff doPruneNBL */
                        int* __restrict__ gm_sciCount,          /* used iff doPruneNBL */
                        const int             numTypes,
                        const float           rCoulombSq,
                        const float           rVdwSq,
                        const float           twoKRf,
                        const float           ewaldBeta,
                        const float           rlistOuterSq,
                        const float           ewaldShift,
                        const float           epsFac,
                        const float           ewaldCoeffLJ_2,
                        const float           cRF,
                        const shift_consts_t  dispersionShift,
                        const shift_consts_t  repulsionShift,
                        const switch_consts_t vdwSwitch,
                        const float           rVdwSwitch,
                        const float           ljEwaldShift,
                        const float           coulombTabScale,
                        const bool            calcShift)
{
    static constexpr EnergyFunctionProperties<elecType, vdwType> props;

    // The post-prune j-i cluster-pair organization is linked to how exclusion and interaction mask
    // data is stored. Currently, this is ideally suited for 32-wide subgroup size but slightly less
    // so for others, e.g. subGroupSize > prunedClusterPairSize on AMD GCN / CDNA.
    constexpr int prunedClusterPairSize = c_clSize * c_splitClSize;

    constexpr int numReductionSteps = gmx::StaticLog2<subGroupSize / c_clSize>::value;
    /* We use shuffles only if we:
     * - use no more than three reduction steps (we only implement 1-3), and
     * - have acceptable cluster size (only 4x4 and 8x8 supported).

     * However, exceptionally we don't do 1-step reduction on hardware with poor global
     * floating-point atomics because that reduction needs many such atomics.
     * Currently (mid-2022), it disables shuffle reduction on all low-end Intel devices, because
     * it causes up to 20x slowdown compared to generic, local memory-based reduction. */
    constexpr bool useShuffleReductionForceI =
            (numReductionSteps <= 3) && (c_clSize == 8 || c_clSize == 4)
            && !(numReductionSteps == 1 && c_avoidFloatingPointAtomics);
    constexpr bool useShuffleReductionForceJ = gmx::isPowerOfTwo(c_nbnxnGpuNumClusterPerSupercluster);

    // Local memory buffer for i x+q pre-loading
    sycl::local_accessor<Float4, 1> sm_xq(
            sycl::range<1>(c_nbnxnGpuNumClusterPerSupercluster * c_clSize), cgh);

    auto sm_atomTypeI = [&]() {
        if constexpr (!props.vdwComb)
        {
            return sycl::local_accessor<int, 1>(
                    sycl::range<1>(c_nbnxnGpuNumClusterPerSupercluster * c_clSize), cgh);
        }
        else
        {
            return nullptr;
        }
    }();

    auto sm_ljCombI = [&]() {
        if constexpr (props.vdwComb)
        {
            return sycl::local_accessor<Float2, 1>(
                    sycl::range<1>(c_nbnxnGpuNumClusterPerSupercluster * c_clSize), cgh);
        }
        else
        {
            return nullptr;
        }
    }();

    /* Local memory buffer for force and energy reduction.
     * For i- and j-force reduction, we need DIM elements for each thread.
     * For energy reduction, we need only 1 element (or none at all) per work-group.
     * But using nullptr as local buffers here triggers a bug with DPC++/OpenCL
     * (https://github.com/intel/llvm/issues/4969), so we always allocate at least one element;
     * it also simplifies the conditional a bit.
     */
    constexpr bool haveAnyLocalMemoryForceReduction =
            !useShuffleReductionForceI || !useShuffleReductionForceJ;
    // need one element for energy reduction (or dummy for DPCPP, https://github.com/intel/llvm/issues/4969)
    constexpr bool needExtraElementForReduction = doCalcEnergies || (GMX_SYCL_DPCPP != 0);
    constexpr int  sm_reductionBufferSize       = haveAnyLocalMemoryForceReduction
                                                          ? c_clSize * c_clSize * DIM
                                                          : (needExtraElementForReduction ? 1 : 0);
    sycl::local_accessor<float, 1> sm_reductionBuffer(sycl::range<1>(sm_reductionBufferSize), cgh);

    auto sm_prunedPairCount = [&]() {
        if constexpr (doPruneNBL && nbnxmSortListsOnGpu())
        {
            return sycl::local_accessor<int, 1>(sycl::range<1>(1), cgh);
        }
        else
        {
            return nullptr;
        }
    }();

    /* Flag to control the calculation of exclusion forces in the kernel
     * We do that with Ewald (elec/vdw) and RF. Cut-off only has exclusion for energy terms. */
    constexpr bool doExclusionForces =
            (props.elecEwald || props.elecRF || props.vdwEwald || (props.elecCutoff && doCalcEnergies));


    return [=](sycl::nd_item<3> itemIdx) [[intel::reqd_sub_group_size(subGroupSize)]]
    {
        if constexpr (skipKernelCompilation<subGroupSize>())
        {
            return;
        }
        /* thread/block/warp id-s */
        const unsigned tidxi = itemIdx.get_local_id(2);
        const unsigned tidxj = itemIdx.get_local_id(1);
        const unsigned tidx  = tidxj * c_clSize + tidxi;

        const unsigned bidx = itemIdx.get_group(0);

        const sycl::sub_group sg = itemIdx.get_sub_group();
        // Could use sg.get_group_range to compute the imask & exclusion Idx, but too much of the logic relies on it anyway
        // and in cases where prunedClusterPairSize != subGroupSize we can't use it anyway
        const unsigned imeiIdx = tidx / prunedClusterPairSize;

#if defined(__SYCL_DEVICE_ONLY__) && defined(__AMDGCN__)
        FCiFloat3 fCiBuf_[c_nbnxnGpuNumClusterPerSupercluster]; // i force buffer
        for (int i = 0; i < c_nbnxnGpuNumClusterPerSupercluster; i++)
        {
            fCiBuf_[i] = { 0.0F, 0.0F, 0.0F };
        }
        auto fCiBufX = [&](auto i) -> float& { return fCiBuf_[i][XX]; };
        auto fCiBufY = [&](auto i) -> float& { return fCiBuf_[i][YY]; };
        auto fCiBufZ = [&](auto i) -> float& { return fCiBuf_[i][ZZ]; };
#else
        float fCiBufX_[c_nbnxnGpuNumClusterPerSupercluster] = { 0.0F }; // i force buffer
        float fCiBufY_[c_nbnxnGpuNumClusterPerSupercluster] = { 0.0F }; // i force buffer
        float fCiBufZ_[c_nbnxnGpuNumClusterPerSupercluster] = { 0.0F }; // i force buffer
        auto  fCiBufX = [&](auto i) -> float& { return fCiBufX_[i]; };
        auto  fCiBufY = [&](auto i) -> float& { return fCiBufY_[i]; };
        auto  fCiBufZ = [&](auto i) -> float& { return fCiBufZ_[i]; };
#endif

        const nbnxn_sci_t nbSci          = gm_plistSci[bidx];
        const int         sci            = nbSci.sci;
        const int         cijPackedBegin = nbSci.cjPackedBegin;
        const int         cijPackedEnd   = nbSci.cjPackedEnd;

        // Only needed if props.elecEwaldAna
        const float beta2 = ewaldBeta * ewaldBeta;
        const float beta3 = ewaldBeta * ewaldBeta * ewaldBeta;

        // We may need only a subset of threads active for preloading i-atoms
        // depending on the super-cluster and cluster / thread-block size.
        constexpr bool c_loadUsingAllXYThreads = (c_clSize == c_nbnxnGpuNumClusterPerSupercluster);
        if (c_loadUsingAllXYThreads || tidxj < c_nbnxnGpuNumClusterPerSupercluster)
        {
            for (int i = 0; i < c_nbnxnGpuNumClusterPerSupercluster; i += c_clSize)
            {
                /* Pre-load i-atom x and q into shared memory */
                const int ci       = sci * c_nbnxnGpuNumClusterPerSupercluster + tidxj + i;
                const int ai       = ci * c_clSize + tidxi;
                const int cacheIdx = (tidxj + i) * c_clSize + tidxi;

                const Float3 shift = gm_shiftVec[nbSci.shift];
                Float4       xqi   = gm_xq[ai];
                xqi += Float4(shift[0], shift[1], shift[2], 0.0F);
                xqi[3] *= epsFac;
                sm_xq[cacheIdx] = xqi;

                if constexpr (!props.vdwComb)
                {
                    // Pre-load the i-atom types into shared memory
                    sm_atomTypeI[cacheIdx] = gm_atomTypes[ai];
                }
                else
                {
                    // Pre-load the LJ combination parameters into shared memory
                    sm_ljCombI[cacheIdx] = gm_ljComb[ai];
                }
            }
        }

        if constexpr (doPruneNBL && nbnxmSortListsOnGpu())
        {
            /* Initialise one int for reducing prunedPairCount over warps */
            if (tidx == 0)
            {
                sm_prunedPairCount[0] = 0;
            }
        }
        int prunedPairCount = 0;

        itemIdx.barrier(fence_space::local_space);

        float ewaldCoeffLJ_6_6; // Only needed if (props.vdwEwald)
        if constexpr (props.vdwEwald)
        {
            ewaldCoeffLJ_6_6 = ewaldCoeffLJ_2 * ewaldCoeffLJ_2 * ewaldCoeffLJ_2 * c_oneSixth;
        }

        float energyVdw, energyElec; // Only needed if (doCalcEnergies)
        if constexpr (doCalcEnergies)
        {
            energyVdw = energyElec = 0.0F;
        }
        if constexpr (doCalcEnergies && doExclusionForces)
        {
            if (nbSci.shift == gmx::c_centralShiftIndex
                && gm_plistCJPacked[cijPackedBegin].cj[0] == sci * c_nbnxnGpuNumClusterPerSupercluster)
            {
                // we have the diagonal: add the charge and LJ self interaction energy term
                for (int i = 0; i < c_nbnxnGpuNumClusterPerSupercluster; i++)
                {
                    // TODO: Are there other options?
                    if constexpr (props.elecEwald || props.elecRF || props.elecCutoff)
                    {
                        const float qi = sm_xq[i * c_clSize + tidxi][3];
                        energyElec += qi * qi;
                    }
                    if constexpr (props.vdwEwald)
                    {
                        energyVdw +=
                                gm_nbfp[gm_atomTypes[(sci * c_nbnxnGpuNumClusterPerSupercluster + i) * c_clSize + tidxi]
                                        * (numTypes + 1)][0];
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
            } // (nbSci.shift == gmx::c_centralShiftIndex && a_plistCJPacked[cijPackedBegin].cj[0] == sci * c_nbnxnGpuNumClusterPerSupercluster)
        }     // (doCalcEnergies && doExclusionForces)

        // Only needed if (doExclusionForces)
        // Note that we use & instead of && for performance (benchmarked in 2017)
        const bool nonSelfInteraction = !(nbSci.shift == gmx::c_centralShiftIndex & tidxj <= tidxi);

        // loop over the j clusters = seen by any of the atoms in the current super-cluster
        for (int jPacked = cijPackedBegin; jPacked < cijPackedEnd; jPacked += 1)
        {
            unsigned imask = UNIFORM_LOAD_CLUSTER_PAIR_DATA(gm_plistCJPacked[jPacked].imei[imeiIdx].imask);
            if (!doPruneNBL && !imask)
            {
                continue;
            }
            const int wexclIdx =
                    UNIFORM_LOAD_CLUSTER_PAIR_DATA(gm_plistCJPacked[jPacked].imei[imeiIdx].excl_ind);

            static_assert(gmx::isPowerOfTwo(prunedClusterPairSize));
            const unsigned wexcl = gm_plistExcl[wexclIdx].pair[tidx & (prunedClusterPairSize - 1)];
            // Unrolling has been verified to improve performance on AMD and Nvidia
#if defined(__AMDGCN__)
            constexpr int unrollFactor =
                    c_nbnxnGpuJgroupSize; // Unrolling has been verified to improve performance on AMD
#elif defined(__SYCL_CUDA_ARCH__) && __SYCL_CUDA_ARCH__ >= 800
            // Unrolling parameters follow CUDA implementation for Ampere and later.
            constexpr int unrollFactor = [=]() {
                if constexpr (!doCalcEnergies && !doPruneNBL)
                {
                    if constexpr (props.elecCutoff || props.elecRF
                                  || (props.elecEwald && !props.vdwFSwitch && !props.vdwPSwitch
                                      && (props.vdwCombLB || __SYCL_CUDA_ARCH__ == 800)))
                    {
                        return 4;
                    }
                    else
                    {
                        return 2;
                    }
                }
                else
                {
                    if constexpr (props.elecCutoff
                                  || (props.elecRF && !props.vdwFSwitch && !props.vdwPSwitch))
                    {
                        return 2;
                    }
                    else
                    {
                        return 1;
                    }
                }
            }();
#else
            constexpr int unrollFactor = 1; // No unrolling.
#endif

#pragma unroll unrollFactor
            for (int jm = 0; jm < c_nbnxnGpuJgroupSize; jm++)
            {
                const bool maskSet =
                        imask & (superClInteractionMask << (jm * c_nbnxnGpuNumClusterPerSupercluster));
                if (!maskSet)
                {
                    continue;
                }
                unsigned  maskJI = (1U << (jm * c_nbnxnGpuNumClusterPerSupercluster));
                const int cj     = gm_plistCJPacked[jPacked].cj[jm];
                const int aj     = cj * c_clSize + tidxj;

                // load j atom data
                const Float4 xqj = gm_xq[aj];

                const Float3 xj(xqj[0], xqj[1], xqj[2]);
                const float  qj = xqj[3];
                int          atomTypeJ; // Only needed if (!props.vdwComb)
                Float2       ljCombJ;   // Only needed if (props.vdwComb)
                if constexpr (props.vdwComb)
                {
                    ljCombJ = gm_ljComb[aj];
                }
                else
                {
                    atomTypeJ = gm_atomTypes[aj];
                }

                Float3 fCjBuf(0.0F, 0.0F, 0.0F);

#pragma unroll c_nbnxnGpuNumClusterPerSupercluster
                for (int i = 0; i < c_nbnxnGpuNumClusterPerSupercluster; i++)
                {
                    if (imask & maskJI)
                    {
                        // i cluster index
                        const int ci = sci * c_nbnxnGpuNumClusterPerSupercluster + i;
                        // all threads load an atom from i cluster ci into shmem!
                        const Float4 xqi = sm_xq[i * c_clSize + tidxi];
                        const Float3 xi(xqi[0], xqi[1], xqi[2]);

                        // distance between i and j atoms
                        const Float3 rv = xi - xj;
                        float        r2 = norm2(rv);

                        if constexpr (doPruneNBL)
                        {
                            /* If _none_ of the atoms pairs are in cutoff range,
                             * the bit corresponding to the current
                             * cluster-pair in imask gets set to 0. */
                            if (!sycl::any_of_group(sg, r2 < rlistOuterSq))
                            {
                                imask &= ~maskJI;
                            }
                        }
                        const float pairExclMask = (wexcl & maskJI) ? 1.0F : 0.0F;

                        // cutoff & exclusion check

                        const bool notExcluded = doExclusionForces ? (nonSelfInteraction | (ci != cj))
                                                                   : (wexcl & maskJI);

#if defined(__SYCL_CUDA_ARCH__)
                        /* The use of * was benchmarked in 2024 for DPC++ 2024.1 CUDA and
                         * found to be faster than the use of &&, just like for CUDA.
                         * "&&" is still better for Intel and AMD devices. */
                        if ((r2 < rCoulombSq) * notExcluded)
#else // Intel and AMD paths
                        if ((r2 < rCoulombSq) && notExcluded)
#endif
                        {
                            const float qi = xqi[3];
                            int         atomTypeI; // Only needed if (!props.vdwComb)
                            float       sigma, epsilon;
                            Float2      c6c12;

                            if constexpr (!props.vdwComb)
                            {
                                /* LJ 6*C6 and 12*C12 */
                                atomTypeI = sm_atomTypeI[i * c_clSize + tidxi];
                                c6c12     = gm_nbfp[numTypes * atomTypeI + atomTypeJ];
                            }
                            else
                            {
                                const Float2 ljCombI = sm_ljCombI[i * c_clSize + tidxi];
                                if constexpr (props.vdwCombGeom)
                                {
                                    c6c12 = Float2(ljCombI[0] * ljCombJ[0], ljCombI[1] * ljCombJ[1]);
                                }
                                else
                                {
                                    static_assert(props.vdwCombLB);
                                    // LJ 2^(1/6)*sigma and 12*epsilon
                                    sigma   = ljCombI[0] + ljCombJ[0];
                                    epsilon = ljCombI[1] * ljCombJ[1];
                                    if constexpr (doCalcEnergies)
                                    {
                                        c6c12 = convertSigmaEpsilonToC6C12(sigma, epsilon);
                                    }
                                } // props.vdwCombGeom
                            }     // !props.vdwComb

                            // c6 and c12 are unused and garbage iff props.vdwCombLB && !doCalcEnergies
                            const float c6  = c6c12[0];
                            const float c12 = c6c12[1];

                            // Ensure distance do not become so small that r^-12 overflows
                            r2 = sycl::max(r2, c_nbnxnMinDistanceSquared);
#if GMX_SYCL_HIPSYCL
                            // No fast/native functions in some compilation passes
                            const float rInv = sycl::rsqrt(r2);
#else
                            // SYCL-TODO: sycl::half_precision::rsqrt?
                            const float rInv = sycl::native::rsqrt(r2);
#endif
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
                                ljEwaldComb<doCalcEnergies, vdwType>(gm_nbfpComb,
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
                                                    gm_coulombTab, coulombTabScale, r2 * rInv))
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
                                            qi * qj * (pairExclMask * rInv + 0.5F * twoKRf * r2 - cRF);
                                }
                                if constexpr (props.elecEwald)
                                {
                                    energyElec += qi * qj
                                                  * (rInv * (pairExclMask - sycl::erf(r2 * rInv * ewaldBeta))
                                                     - pairExclMask * ewaldShift);
                                }
                            }

                            const Float3 forceIJ = rv * fInvR;

                            /* accumulate j forces in registers */
                            fCjBuf -= forceIJ;
                            /* accumulate i forces in registers */
#if defined(__SYCL_DEVICE_ONLY__) && defined(__AMDGCN__)
                            fCiBuf_[i] += FCiFloat3(forceIJ);
#else
                            fCiBufX(i) += forceIJ[0];
                            fCiBufY(i) += forceIJ[1];
                            fCiBufZ(i) += forceIJ[2];
#endif
                        } // (r2 < rCoulombSq) && notExcluded
                    }     // (imask & maskJI)
                    /* shift the mask bit by 1 */
                    maskJI += maskJI;
                } // for (int i = 0; i < c_nbnxnGpuNumClusterPerSupercluster; i++)
                /* reduce j forces */
                reduceForceJ<useShuffleReductionForceJ>(
                        sm_reductionBuffer, fCjBuf, itemIdx, tidxi, tidxj, aj, gm_f);
            } // for (int jm = 0; jm < c_nbnxnGpuJgroupSize; jm++)
            if constexpr (doPruneNBL)
            {
                /* Update the imask with the new one which does not contain the
                 * out of range clusters anymore. */
                gm_plistCJPacked[jPacked].imei[imeiIdx].imask = imask;
                if constexpr (nbnxmSortListsOnGpu())
                {
                    prunedPairCount += sycl::popcount(imask);
                }
            }
        } // for (int jPacked = cijPackedBegin; jPacked < cijPackedEnd; jPacked += 1)

        /* skip central shifts when summing shift forces */
        const bool doCalcShift = (calcShift && nbSci.shift != gmx::c_centralShiftIndex);

        reduceForceIAndFShift<useShuffleReductionForceI, subGroupSize>(
                sm_reductionBuffer, fCiBufX, fCiBufY, fCiBufZ, doCalcShift, itemIdx, tidxi, tidxj, sci, nbSci.shift, gm_f, gm_fShift);

        if constexpr (doCalcEnergies)
        {
            const float energyVdwGroup =
                    groupReduce<subGroupSize, c_clSizeSq>(itemIdx, tidx, sm_reductionBuffer, energyVdw);
            itemIdx.barrier(fence_space::local_space); // Prevent the race on sm_reductionBuffer.
            const float energyElecGroup = groupReduce<subGroupSize, c_clSizeSq>(
                    itemIdx, tidx, sm_reductionBuffer, energyElec);

            if (tidx == 0)
            {
                atomicFetchAdd(gm_energyVdw[0], energyVdwGroup);
                atomicFetchAdd(gm_energyElec[0], energyElecGroup);
            }
        }
        if constexpr (doPruneNBL && nbnxmSortListsOnGpu())
        {
            /* aggregate neighbour counts, to be used in bucket sci sort */
            /* One thread in each warp contributes the count for that warp as soon as it reaches
             * here. Masks are calculated per warp in a warp synchronising operation, so no
             * syncthreads required here. */
            if (sg.leader())
            {
                atomicFetchAddLocal(sm_prunedPairCount[0], prunedPairCount);
            }
            itemIdx.barrier(fence_space::local_space);
            prunedPairCount = sm_prunedPairCount[0];
            if (tidxi == 0 && tidxj == 0)
            {
                /* one thread in the block writes the final count for this sci */
                int index = sycl::max(c_sciHistogramSize - prunedPairCount - 1, 0);
                atomicFetchAdd(gm_sciHistogram[index], 1);
                gm_sciCount[bidx] = index;
            }
        }
    };
}

//! \brief NBNXM kernel launch code.
template<int subGroupSize, bool doPruneNBL, bool doCalcEnergies, enum ElecType elecType, enum VdwType vdwType, class... Args>
static void launchNbnxmKernel(const DeviceStream& deviceStream, const int numSci, Args&&... args)
{
    using kernelNameType = NbnxmKernel<doPruneNBL, doCalcEnergies, elecType, vdwType, subGroupSize>;

    /* Kernel launch config:
     * - The thread block dimensions match the size of i-clusters, j-clusters,
     *   and j-cluster concurrency, in x, y, and z, respectively.
     * - The 1D block-grid contains as many blocks as super-clusters.
     */
    const int               numBlocks = numSci;
    const sycl::range<3>    blockSize{ 1, c_clSize, c_clSize };
    const sycl::range<3>    globalSize{ numBlocks * blockSize[0], blockSize[1], blockSize[2] };
    const sycl::nd_range<3> range{ globalSize, blockSize };

    sycl::queue q = deviceStream.stream();

    q.submit(GMX_SYCL_DISCARD_EVENT[&](sycl::handler & cgh) {
        auto kernel = nbnxmKernel<subGroupSize, doPruneNBL, doCalcEnergies, elecType, vdwType>(
                cgh, std::forward<Args>(args)...);
        cgh.parallel_for<kernelNameType>(range, kernel);
    });
}

//! \brief Select templated kernel and launch it.
template<int subGroupSize, bool doPruneNBL, bool doCalcEnergies, class... Args>
void chooseAndLaunchNbnxmKernel(enum ElecType elecType, enum VdwType vdwType, Args&&... args)
{
    gmx::dispatchTemplatedFunction(
            [&](auto elecType_, auto vdwType_) {
                return launchNbnxmKernel<subGroupSize, doPruneNBL, doCalcEnergies, elecType_, vdwType_>(
                        std::forward<Args>(args)...);
            },
            elecType,
            vdwType);
}

template<int subGroupSize, bool doPruneNBL, bool doCalcEnergies>
void launchNbnxmKernelHelper(NbnxmGpu* nb, const gmx::StepWorkload& stepWork, const InteractionLocality iloc)
{
    NBAtomDataGpu*      adat         = nb->atdat;
    NBParamGpu*         nbp          = nb->nbparam;
    auto*               plist        = nb->plist[iloc].get();
    const DeviceStream& deviceStream = *nb->deviceStreams[iloc];

    GMX_ASSERT(doPruneNBL == (plist->haveFreshList && !nb->didPrune[iloc]), "Wrong template called");
    GMX_ASSERT(doCalcEnergies == stepWork.computeEnergy, "Wrong template called");

    chooseAndLaunchNbnxmKernel<subGroupSize, doPruneNBL, doCalcEnergies>(
            nbp->elecType,
            nbp->vdwType,
            deviceStream,
            plist->numSci,
            adat->xq.get_pointer(),
            adat->f.get_pointer(),
            adat->shiftVec.get_pointer(),
            adat->fShift.get_pointer(),
            adat->eElec.get_pointer(),
            adat->eLJ.get_pointer(),
            plist->cjPacked.get_pointer(),
            (doPruneNBL || !nbnxmSortListsOnGpu()) ? plist->sci.get_pointer()
                                                   : plist->sorting.sciSorted.get_pointer(),
            plist->excl.get_pointer(),
            adat->ljComb.get_pointer(),
            adat->atomTypes.get_pointer(),
            nbp->nbfp.get_pointer(),
            nbp->nbfp_comb.get_pointer(),
            nbp->coulomb_tab.get_pointer(),
            plist->sorting.sciHistogram.get_pointer(),
            plist->sorting.sciCount.get_pointer(),
            adat->numTypes,
            nbp->rcoulomb_sq,
            nbp->rvdw_sq,
            nbp->two_k_rf,
            nbp->ewald_beta,
            nbp->rlistOuter_sq,
            nbp->sh_ewald,
            nbp->epsfac,
            nbp->ewaldcoeff_lj * nbp->ewaldcoeff_lj,
            nbp->c_rf,
            nbp->dispersion_shift,
            nbp->repulsion_shift,
            nbp->vdw_switch,
            nbp->rvdw_switch,
            nbp->sh_lj_ewald,
            nbp->coulomb_tab_scale,
            stepWork.computeVirial);
}

} // namespace gmx
