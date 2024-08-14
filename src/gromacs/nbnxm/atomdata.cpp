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

#include "gmxpre.h"

#include "atomdata.h"

#include <cmath>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <array>
#include <filesystem>
#include <string>

#include "gromacs/math/functions.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/atominfo.h"
#include "gromacs/mdtypes/forcerec.h" // only for GET_CGINFO_*
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/simd/simd.h"
#include "gromacs/utility/allocator.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/range.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringutil.h"

#include "grid.h"
#include "gridset.h"
#include "nbnxm_enums.h"
#include "nbnxm_geometry.h"
#include "nbnxm_gpu.h"
#include "pairlist.h"
#include "simd_energy_accumulator.h"

namespace gmx
{

const char* enumValueToString(LJCombinationRule enumValue)
{
    static constexpr EnumerationArray<LJCombinationRule, const char*> s_ljCombinationRuleNames = {
        "Geometric", "Lorentz-Berthelot", "None"
    };
    return s_ljCombinationRuleNames[enumValue];
}

void nbnxn_atomdata_t::resizeCoordinateBuffer(const int numAtoms, const int domainDecompositioZone)
{
    numAtoms_ = numAtoms;

    if (domainDecompositioZone == 0)
    {
        numLocalAtoms_ = numAtoms;
    }

    x_.resize(numAtoms * xstride);
}

void nbnxn_atomdata_t::resizeForceBuffers()
{
    /* Force buffers need padding up to a multiple of the buffer flag size */
    const int paddedSize =
            (numAtoms() + NBNXN_BUFFERFLAG_SIZE - 1) / NBNXN_BUFFERFLAG_SIZE * NBNXN_BUFFERFLAG_SIZE;

    /* Should we let each thread allocate it's own data instead? */
    for (nbnxn_atomdata_output_t& outputBuffer : outputBuffers_)
    {
        outputBuffer.f.resize(paddedSize * fstride);
    }
}

/* Initializes an nbnxn_atomdata_output_t data structure */
nbnxn_atomdata_output_t::nbnxn_atomdata_output_t(NbnxmKernelType kernelType,
                                                 int             numEnergyGroups,
                                                 PinningPolicy   pinningPolicy) :
    f({}, { pinningPolicy }),
    fshift({}, { pinningPolicy }),
    Vvdw({}, { pinningPolicy }),
    Vc({}, { pinningPolicy })
{
    fshift.resize(c_numShiftVectors * DIM);
    Vvdw.resize(numEnergyGroups * numEnergyGroups);
    Vc.resize(numEnergyGroups * numEnergyGroups);

    if (kernelTypeIsSimd(kernelType))
    {
        if (numEnergyGroups == 1)
        {
            accumulatorSingleEnergies = std::make_unique<EnergyAccumulator<false, true>>();
        }
        else
        {
            const int c_iClusterSize = sc_iClusterSize(kernelType);
            const int c_jClusterSize = sc_jClusterSize(kernelType);

            accumulatorGroupEnergies = std::make_unique<EnergyAccumulator<true, true>>(
                    numEnergyGroups, c_iClusterSize, c_jClusterSize);
        }
    }
}

nbnxn_atomdata_output_t::nbnxn_atomdata_output_t(nbnxn_atomdata_output_t&&) noexcept = default;

nbnxn_atomdata_output_t::~nbnxn_atomdata_output_t() = default;

static void copy_int_to_nbat_int(const int* a, int na, int na_round, const int* in, int fill, int* innb)
{
    int j = 0;
    for (int i = 0; i < na; i++)
    {
        innb[j++] = in[a[i]];
    }
    /* Complete the partially filled last cell with fill */
    for (int i = na; i < na_round; i++)
    {
        innb[j++] = fill;
    }
}

void copy_rvec_to_nbat_real(const int* a, int na, int na_round, const rvec* x, int nbatFormat, real* xnb, int a0)
{
    /* We complete partially filled cells, can only be the last one in each
     * column, with coordinates farAway. The actual coordinate value does
     * not influence the results, since these filler particles do not interact.
     * Clusters with normal atoms + fillers have a bounding box based only
     * on the coordinates of the atoms. Clusters with only fillers have as
     * the bounding box the coordinates of the first filler. Such clusters
     * are not considered as i-entries, but they are considered as j-entries.
     * So for performance it is better to have their bounding boxes far away,
     * such that filler only clusters don't end up in the pair list.
     */
    const real farAway = -1000000;

    if (nbatFormat == nbatXYZ)
    {
        int i = 0;
        int j = a0 * STRIDE_XYZ;
        for (; i < na; i++)
        {
            xnb[j++] = x[a[i]][XX];
            xnb[j++] = x[a[i]][YY];
            xnb[j++] = x[a[i]][ZZ];
        }
        /* Complete the partially filled last cell with farAway elements */
        for (; i < na_round; i++)
        {
            xnb[j++] = farAway;
            xnb[j++] = farAway;
            xnb[j++] = farAway;
        }
    }
    else if (nbatFormat == nbatXYZQ)
    {
        int i = 0;
        int j = a0 * STRIDE_XYZQ;
        for (; i < na; i++)
        {
            xnb[j++] = x[a[i]][XX];
            xnb[j++] = x[a[i]][YY];
            xnb[j++] = x[a[i]][ZZ];
            j++;
        }
        /* Complete the partially filled last cell with zeros */
        for (; i < na_round; i++)
        {
            xnb[j++] = farAway;
            xnb[j++] = farAway;
            xnb[j++] = farAway;
            j++;
        }
    }
    else if (nbatFormat == nbatX4)
    {
        int i = 0;
        int j = atom_to_x_index<c_packX4>(a0);
        int c = a0 & (c_packX4 - 1);
        for (; i < na; i++)
        {
            xnb[j + XX * c_packX4] = x[a[i]][XX];
            xnb[j + YY * c_packX4] = x[a[i]][YY];
            xnb[j + ZZ * c_packX4] = x[a[i]][ZZ];
            j++;
            c++;
            if (c == c_packX4)
            {
                j += (DIM - 1) * c_packX4;
                c = 0;
            }
        }
        /* Complete the partially filled last cell with zeros */
        for (; i < na_round; i++)
        {
            xnb[j + XX * c_packX4] = farAway;
            xnb[j + YY * c_packX4] = farAway;
            xnb[j + ZZ * c_packX4] = farAway;
            j++;
            c++;
            if (c == c_packX4)
            {
                j += (DIM - 1) * c_packX4;
                c = 0;
            }
        }
    }
    else if (nbatFormat == nbatX8)
    {
        int i = 0;
        int j = atom_to_x_index<c_packX8>(a0);
        int c = a0 & (c_packX8 - 1);
        for (; i < na; i++)
        {
            xnb[j + XX * c_packX8] = x[a[i]][XX];
            xnb[j + YY * c_packX8] = x[a[i]][YY];
            xnb[j + ZZ * c_packX8] = x[a[i]][ZZ];
            j++;
            c++;
            if (c == c_packX8)
            {
                j += (DIM - 1) * c_packX8;
                c = 0;
            }
        }
        /* Complete the partially filled last cell with zeros */
        for (; i < na_round; i++)
        {
            xnb[j + XX * c_packX8] = farAway;
            xnb[j + YY * c_packX8] = farAway;
            xnb[j + ZZ * c_packX8] = farAway;
            j++;
            c++;
            if (c == c_packX8)
            {
                j += (DIM - 1) * c_packX8;
                c = 0;
            }
        }
    }
    else
    {
        gmx_incons("Unsupported nbnxn_atomdata_t format");
    }
}

/* Stores the LJ parameter data in a format convenient for different kernels */
static void set_lj_parameter_data(nbnxn_atomdata_t::Params* params, gmx_bool bSIMD)
{
    int nt = params->numTypes;

    if (bSIMD)
    {
#if GMX_SIMD
        /* nbfp_aligned stores two parameters using the stride most suitable
         * for the present SIMD architecture, as specified by the constant
         * c_simdBestPairAlignment from the SIMD header.
         * There's a slight inefficiency in allocating and initializing nbfp_aligned
         * when it might not be used, but introducing the conditional code is not
         * really worth it.
         */
        params->nbfp_aligned.resize(nt * nt * c_simdBestPairAlignment);

        for (int i = 0; i < nt; i++)
        {
            for (int j = 0; j < nt; j++)
            {
                params->nbfp_aligned[(i * nt + j) * c_simdBestPairAlignment + 0] =
                        params->nbfp[(i * nt + j) * 2 + 0];
                params->nbfp_aligned[(i * nt + j) * c_simdBestPairAlignment + 1] =
                        params->nbfp[(i * nt + j) * 2 + 1];
                if (c_simdBestPairAlignment > 2)
                {
                    params->nbfp_aligned[(i * nt + j) * c_simdBestPairAlignment + 2] = 0;
                    params->nbfp_aligned[(i * nt + j) * c_simdBestPairAlignment + 3] = 0;
                }
            }
        }
#endif
    }

    /* We use combination rule data for SIMD combination rule kernels
     * and with LJ-PME kernels. We then only need parameters per atom type,
     * not per pair of atom types.
     */
    params->nbfp_comb.resize(nt * 2);
    switch (params->ljCombinationRule)
    {
        case LJCombinationRule::Geometric:
            for (int i = 0; i < nt; i++)
            {
                /* Store the sqrt of the diagonal from the nbfp matrix */
                params->nbfp_comb[i * 2]     = std::sqrt(params->nbfp[(i * nt + i) * 2]);
                params->nbfp_comb[i * 2 + 1] = std::sqrt(params->nbfp[(i * nt + i) * 2 + 1]);
            }
            break;
        case LJCombinationRule::LorentzBerthelot:
            for (int i = 0; i < nt; i++)
            {
                /* Get 6*C6 and 12*C12 from the diagonal of the nbfp matrix */
                const real c6  = params->nbfp[(i * nt + i) * 2];
                const real c12 = params->nbfp[(i * nt + i) * 2 + 1];
                if (c6 > 0 && c12 > 0)
                {
                    /* We store 0.5*2^1/6*sigma and sqrt(4*3*eps),
                     * so we get 6*C6 and 12*C12 after combining.
                     */
                    params->nbfp_comb[i * 2]     = 0.5 * sixthroot(c12 / c6);
                    params->nbfp_comb[i * 2 + 1] = std::sqrt(c6 * c6 / c12);
                }
                else
                {
                    params->nbfp_comb[i * 2]     = 0;
                    params->nbfp_comb[i * 2 + 1] = 0;
                }
            }
            break;
        case LJCombinationRule::None:
            /* We always store the full matrix (see code above) */
            break;
        default: gmx_incons("Unknown combination rule");
    }
}

nbnxn_atomdata_t::SimdMasks::SimdMasks(const NbnxmKernelType gmx_unused kernelType)
{
#if GMX_SIMD
    constexpr int simd_width = GMX_SIMD_REAL_WIDTH;
    /* Set the diagonal cluster pair exclusion mask setup data.
     * In the kernel we check 0 < j - i to generate the masks.
     * Here we store j - i for generating the mask for the first i,
     * we subtract 0.5 to avoid rounding issues.
     * In the kernel we can subtract 1 to generate the subsequent mask.
     */
    const int simd_4xn_diag_size = std::max(sc_iClusterSize(kernelType), simd_width);
    diagonal_4xn_j_minus_i.resize(simd_4xn_diag_size);
    for (int j = 0; j < simd_4xn_diag_size; j++)
    {
        diagonal_4xn_j_minus_i[j] = j - 0.5;
    }

    diagonal_2xnn_j_minus_i.resize(simd_width);
    for (int j = 0; j < simd_width / 2; j++)
    {
        /* The j-cluster size is half the SIMD width */
        diagonal_2xnn_j_minus_i[j] = j - 0.5;
        /* The next half of the SIMD width is for i + 1 */
        diagonal_2xnn_j_minus_i[simd_width / 2 + j] = j - 1 - 0.5;
    }

    /* We use up to 32 bits for exclusion masking.
     * The same masks are used for the 4xN and 2x(N+N) kernels.
     * The masks are read either into integer SIMD registers or into
     * real SIMD registers (together with a cast).
     * In single precision this means the real and integer SIMD registers
     * are of equal size.
     */
    const int simd_excl_size = sc_iClusterSize(kernelType) * simd_width;
#    if GMX_DOUBLE && !GMX_SIMD_HAVE_INT32_LOGICAL
    exclusion_filter64.resize(simd_excl_size);
#    else
    exclusion_filter.resize(simd_excl_size);
#    endif

    for (int j = 0; j < simd_excl_size; j++)
    {
        /* Set the consecutive bits for masking pair exclusions */
#    if GMX_DOUBLE && !GMX_SIMD_HAVE_INT32_LOGICAL
        exclusion_filter64[j] = (1U << j);
#    else
        exclusion_filter[j] = (1U << j);
#    endif
    }

#endif // GMX_SIMD
}

nbnxn_atomdata_t::Params::Params(PinningPolicy pinningPolicy) :
    numTypes(0),
    nbfp({}, { pinningPolicy }),
    nbfp_comb({}, { pinningPolicy }),
    type({}, { pinningPolicy }),
    lj_comb({}, { pinningPolicy }),
    q({}, { pinningPolicy }),
    numEnergyGroups(0)
{
}

/* Initializes an nbnxn_atomdata_t::Params data structure */
static void nbnxn_atomdata_params_init(const MDLogger&                         mdlog,
                                       nbnxn_atomdata_t::Params*               params,
                                       const NbnxmKernelType                   kernelType,
                                       const std::optional<LJCombinationRule>& ljCombinationRule,
                                       const LJCombinationRule                 pmeLJCombinationRule,
                                       ArrayRef<const real>                    nbfp,
                                       const bool                              addFillerAtomType,
                                       const int                               numEnergyGroups)
{
    const bool usingLJPme = (pmeLJCombinationRule != LJCombinationRule::None);

    GMX_RELEASE_ASSERT(!usingLJPme || !ljCombinationRule
                               || ljCombinationRule.value() == LJCombinationRule::None,
                       "Only one of ljCombinationRule and pmeLJCombinationRule can be active");

    const int numTypes = std::sqrt(ssize(nbfp) / 2);
    GMX_RELEASE_ASSERT(size_t(numTypes) * size_t(numTypes) * 2 == nbfp.size(),
                       "The size of nbfp should be two times a square number");

    if (!addFillerAtomType)
    {
        for (int i = 0; i < numTypes; i++)
        {
            if (nbfp[(i * numTypes + numTypes - 1) * 2 + 0] != 0
                || nbfp[(i * numTypes + numTypes - 1) * 2 + 1] != 0
                || nbfp[((numTypes - 1) * numTypes + i) * 2 + 0] != 0
                || nbfp[((numTypes - 1) * numTypes + i) * 2 + 1] != 0)
            {
                GMX_RELEASE_ASSERT(
                        false,
                        "With addFillerAtomType=false atom-type numTypes-1 should have params 0");
            }
        }
    }

    if (debug)
    {
        fprintf(debug,
                "There are %d atom types in the system%s\n",
                numTypes,
                addFillerAtomType ? ", adding one for nbnxn_atomdata_t" : "");
    }
    // We add one type for the filler particles
    params->numTypes = numTypes + (addFillerAtomType ? 1 : 0);
    params->nbfp.resize(params->numTypes * params->numTypes * 2);
    params->nbfp_comb.resize(params->numTypes * 2);

    /* A tolerance of 1e-5 seems reasonable for (possibly hand-typed)
     * force-field floating point parameters.
     */
    real        tol               = 1e-5;
    const char* tolOverrideString = getenv("GMX_LJCOMB_TOL");
    if (tolOverrideString != nullptr)
    {
        double tolOverride = std::strtod(tolOverrideString, nullptr);
        tol                = tolOverride;
    }
    bool bCombGeom = true;
    bool bCombLB   = true;

    /* Temporarily fill params->nbfp_comb with sigma and epsilon
     * to check for the LB rule.
     */
    for (int i = 0; i < numTypes; i++)
    {
        const real c6  = nbfp[(i * numTypes + i) * 2] / 6.0;
        const real c12 = nbfp[(i * numTypes + i) * 2 + 1] / 12.0;
        if (c6 > 0 && c12 > 0)
        {
            params->nbfp_comb[i * 2]     = sixthroot(c12 / c6);
            params->nbfp_comb[i * 2 + 1] = 0.25 * c6 * c6 / c12;
        }
        else if (c6 == 0 && c12 == 0)
        {
            params->nbfp_comb[i * 2]     = 0;
            params->nbfp_comb[i * 2 + 1] = 0;
        }
        else
        {
            /* Can not use LB rule with only dispersion or repulsion */
            bCombLB = false;
        }
    }

    for (int i = 0; i < params->numTypes; i++)
    {
        for (int j = 0; j < params->numTypes; j++)
        {
            if (i < numTypes && j < numTypes)
            {
                /* fr->nbfp has been updated, so that array too now stores c6/c12 including
                 * the 6.0/12.0 prefactors to save 2 flops in the most common case (force-only).
                 */
                real c6  = nbfp[(i * numTypes + j) * 2];
                real c12 = nbfp[(i * numTypes + j) * 2 + 1];

                params->nbfp[(i * params->numTypes + j) * 2]     = c6;
                params->nbfp[(i * params->numTypes + j) * 2 + 1] = c12;

                /* Compare 6*C6 and 12*C12 for geometric cobination rule */
                bCombGeom =
                        bCombGeom
                        && gmx_within_tol(
                                c6 * c6, nbfp[(i * numTypes + i) * 2] * nbfp[(j * numTypes + j) * 2], tol)
                        && gmx_within_tol(
                                c12 * c12,
                                nbfp[(i * numTypes + i) * 2 + 1] * nbfp[(j * numTypes + j) * 2 + 1],
                                tol);

                /* Compare C6 and C12 for Lorentz-Berthelot combination rule */
                c6 /= 6.0;
                c12 /= 12.0;
                bCombLB = bCombLB
                          && ((c6 == 0 && c12 == 0
                               && (params->nbfp_comb[i * 2 + 1] == 0 || params->nbfp_comb[j * 2 + 1] == 0))
                              || (c6 > 0 && c12 > 0
                                  && gmx_within_tol(
                                          sixthroot(c12 / c6),
                                          0.5 * (params->nbfp_comb[i * 2] + params->nbfp_comb[j * 2]),
                                          tol)
                                  && gmx_within_tol(0.25 * c6 * c6 / c12,
                                                    std::sqrt(params->nbfp_comb[i * 2 + 1]
                                                              * params->nbfp_comb[j * 2 + 1]),
                                                    tol)));
            }
            else
            {
                /* Add zero parameters for the additional dummy atom type */
                params->nbfp[(i * params->numTypes + j) * 2]     = 0;
                params->nbfp[(i * params->numTypes + j) * 2 + 1] = 0;
            }
        }
    }
    if (debug)
    {
        fprintf(debug,
                "Combination rules: geometric %s Lorentz-Berthelot %s\n",
                boolToString(bCombGeom),
                boolToString(bCombLB));
    }

    const bool simple = kernelTypeUsesSimplePairlist(kernelType);

    if (usingLJPme || ljCombinationRule)
    {
        params->ljCombinationRule = (usingLJPme ? pmeLJCombinationRule : ljCombinationRule.value());
    }
    else
    {
        /* We prefer the geometric combination rule,
         * as that gives a slightly faster kernel than the LB rule.
         */
        if (bCombGeom)
        {
            params->ljCombinationRule = LJCombinationRule::Geometric;
        }
        else if (bCombLB)
        {
            params->ljCombinationRule = LJCombinationRule::LorentzBerthelot;
        }
        else
        {
            params->ljCombinationRule = LJCombinationRule::None;

            params->nbfp_comb.clear();
        }

        {
            std::string mesg;
            if (params->ljCombinationRule == LJCombinationRule::None)
            {
                mesg = "Using full Lennard-Jones parameter combination matrix";
            }
            else
            {
                mesg = formatString("Using %s Lennard-Jones combination rule",
                                    enumValueToString(params->ljCombinationRule));
            }
            GMX_LOG(mdlog.info).asParagraph().appendText(mesg);
        }
    }

    const bool bSIMD = kernelTypeIsSimd(kernelType);

    set_lj_parameter_data(params, bSIMD);

    params->numEnergyGroups = numEnergyGroups;
    if (!simple)
    {
        // We now check for energy groups already when starting mdrun
        GMX_RELEASE_ASSERT(numEnergyGroups == 1, "GPU kernels do not support energy groups");
    }
    /* Temporary storage goes as #grp^3*simd_width^2/2, so limit to 64 */
    if (params->numEnergyGroups > 64)
    {
        gmx_fatal(FARGS, "With NxN kernels not more than 64 energy groups are supported\n");
    }
    if (params->numEnergyGroups > 1)
    {
        params->energyGroupsPerCluster = std::make_unique<EnergyGroupsPerCluster>(
                params->numEnergyGroups, sc_iClusterSize(kernelType));
    }
}

/* Initializes an nbnxn_atomdata_t data structure */
nbnxn_atomdata_t::nbnxn_atomdata_t(PinningPolicy                           pinningPolicy,
                                   const MDLogger&                         mdlog,
                                   const NbnxmKernelType                   kernelType,
                                   const std::optional<LJCombinationRule>& ljCombinationRule,
                                   const LJCombinationRule                 pmeLJCombinationRule,
                                   ArrayRef<const real>                    nbfp,
                                   const bool                              addFillerAtomType,
                                   const int                               numEnergyGroups,
                                   const int                               numOutputBuffers) :
    params_(pinningPolicy),
    numAtoms_(0),
    numLocalAtoms_(0),
    shift_vec({}, { pinningPolicy }),
    x_({}, { pinningPolicy }),
    simdMasks_(kernelType),
    useBufferFlags_(numOutputBuffers > 1)
{
    nbnxn_atomdata_params_init(mdlog,
                               &paramsDeprecated(),
                               kernelType,
                               ljCombinationRule,
                               pmeLJCombinationRule,
                               nbfp,
                               addFillerAtomType,
                               numEnergyGroups);

    const bool simple = kernelTypeUsesSimplePairlist(kernelType);
    const bool bSIMD  = kernelTypeIsSimd(kernelType);

    if (simple)
    {
        if (bSIMD)
        {
            int pack_x = std::max(sc_iClusterSize(kernelType), sc_jClusterSize(kernelType));
            switch (pack_x)
            {
                case 4: XFormat = nbatX4; break;
                case 8: XFormat = nbatX8; break;
                default: gmx_incons("Unsupported packing width");
            }
        }
        else
        {
            XFormat = nbatXYZ;
        }

        FFormat = XFormat;
    }
    else
    {
        XFormat = nbatXYZQ;
        FFormat = nbatXYZ;
    }

    shift_vec.resize(c_numShiftVectors);

    xstride = (XFormat == nbatXYZQ ? STRIDE_XYZQ : DIM);
    fstride = (FFormat == nbatXYZQ ? STRIDE_XYZQ : DIM);

    /* Initialize the output data structures */
    for (int i = 0; i < numOutputBuffers; i++)
    {
        const auto& outputPinningPolicy = params().type.get_allocator().pinningPolicy();
        outputBuffers_.emplace_back(kernelType, params().numEnergyGroups, outputPinningPolicy);
    }

    bufferFlags_.clear();
}

nbnxn_atomdata_t::~nbnxn_atomdata_t() = default;

template<int packSize>
static void copy_lj_to_nbat_lj_comb(ArrayRef<const real> ljparam_type, const int* type, int na, real* ljparam_at)
{
    /* The LJ params follow the combination rule:
     * copy the params for the type array to the atom array.
     */
    for (int is = 0; is < na; is += packSize)
    {
        for (int k = 0; k < packSize; k++)
        {
            int i                             = is + k;
            ljparam_at[is * 2 + k]            = ljparam_type[type[i] * 2];
            ljparam_at[is * 2 + packSize + k] = ljparam_type[type[i] * 2 + 1];
        }
    }
}

/* Sets the atom type in nbnxn_atomdata_t */
static void nbnxn_atomdata_set_atomtypes(nbnxn_atomdata_t::Params* params,
                                         const GridSet&            gridSet,
                                         ArrayRef<const int>       atomTypes)
{
    params->type.resize(gridSet.numGridAtomsTotal());

    for (const Grid& grid : gridSet.grids())
    {
        /* Loop over all columns and copy and fill */
        for (int i = 0; i < grid.numColumns(); i++)
        {
            const int numAtoms   = grid.paddedNumAtomsInColumn(i);
            const int atomOffset = grid.firstAtomInColumn(i);

            copy_int_to_nbat_int(gridSet.atomIndices().data() + atomOffset,
                                 grid.numAtomsInColumn(i),
                                 numAtoms,
                                 atomTypes.data(),
                                 params->numTypes - 1,
                                 params->type.data() + atomOffset);
        }
    }
}

/* Sets the LJ combination rule parameters in nbnxn_atomdata_t */
static void nbnxn_atomdata_set_ljcombparams(nbnxn_atomdata_t::Params* params,
                                            const int                 XFormat,
                                            const GridSet&            gridSet)
{
    params->lj_comb.resize(gridSet.numGridAtomsTotal() * 2);

    if (params->ljCombinationRule != LJCombinationRule::None)
    {
        for (const Grid& grid : gridSet.grids())
        {
            /* Loop over all columns and copy and fill */
            for (int i = 0; i < grid.numColumns(); i++)
            {
                const int numAtoms   = grid.paddedNumAtomsInColumn(i);
                const int atomOffset = grid.firstAtomInColumn(i);

                if (XFormat == nbatX4)
                {
                    copy_lj_to_nbat_lj_comb<c_packX4>(params->nbfp_comb,
                                                      params->type.data() + atomOffset,
                                                      numAtoms,
                                                      params->lj_comb.data() + atomOffset * 2);
                }
                else if (XFormat == nbatX8)
                {
                    copy_lj_to_nbat_lj_comb<c_packX8>(params->nbfp_comb,
                                                      params->type.data() + atomOffset,
                                                      numAtoms,
                                                      params->lj_comb.data() + atomOffset * 2);
                }
                else if (XFormat == nbatXYZQ)
                {
                    copy_lj_to_nbat_lj_comb<1>(params->nbfp_comb,
                                               params->type.data() + atomOffset,
                                               numAtoms,
                                               params->lj_comb.data() + atomOffset * 2);
                }
            }
        }
    }
}

/* Sets the charges in nbnxn_atomdata_t *nbat */
static void nbnxn_atomdata_set_charges(nbnxn_atomdata_t* nbat, const GridSet& gridSet, ArrayRef<const real> charges)
{
    if (nbat->XFormat != nbatXYZQ)
    {
        nbat->paramsDeprecated().q.resize(nbat->numAtoms());
    }

    for (const Grid& grid : gridSet.grids())
    {
        /* Loop over all columns and copy and fill */
        for (int cxy = 0; cxy < grid.numColumns(); cxy++)
        {
            const int atomOffset     = grid.firstAtomInColumn(cxy);
            const int numAtoms       = grid.numAtomsInColumn(cxy);
            const int paddedNumAtoms = grid.paddedNumAtomsInColumn(cxy);

            if (nbat->XFormat == nbatXYZQ)
            {
                real* q = nbat->x().data() + atomOffset * STRIDE_XYZQ + ZZ + 1;
                for (int i = 0; i < numAtoms; i++)
                {
                    *q = charges[gridSet.atomIndices()[atomOffset + i]];
                    q += STRIDE_XYZQ;
                }
                /* Complete the partially filled last cell with zeros */
                for (int i = numAtoms; i < paddedNumAtoms; i++)
                {
                    *q = 0;
                    q += STRIDE_XYZQ;
                }
            }
            else
            {
                real* q = nbat->paramsDeprecated().q.data() + atomOffset;
                for (int i = 0; i < numAtoms; i++)
                {
                    *q = charges[gridSet.atomIndices()[atomOffset + i]];
                    q++;
                }
                /* Complete the partially filled last cell with zeros */
                for (int i = numAtoms; i < paddedNumAtoms; i++)
                {
                    *q = 0;
                    q++;
                }
            }
        }
    }
}

/* Set the charges of perturbed atoms in nbnxn_atomdata_t to 0.
 * This is to automatically remove the RF/PME self term in the nbnxn kernels.
 * Part of the zero interactions are still calculated in the normal kernels.
 * All perturbed interactions are calculated in the free energy kernel,
 * using the original charge and LJ data, not nbnxn_atomdata_t.
 */
static void nbnxn_atomdata_mask_fep(nbnxn_atomdata_t* nbat, const GridSet& gridSet)
{
    nbnxn_atomdata_t::Params& params = nbat->paramsDeprecated();

    const bool formatIsXYZQ = (nbat->XFormat == nbatXYZQ);

    real* q        = formatIsXYZQ ? (nbat->x().data() + ZZ + 1) : params.q.data();
    int   stride_q = formatIsXYZQ ? STRIDE_XYZQ : 1;

    for (const Grid& grid : gridSet.grids())
    {
        const int nsubc = (grid.geometry().isSimple) ? 1 : c_gpuNumClusterPerCell;

        const int c_offset = grid.firstAtomInColumn(0);

        /* Loop over all columns and copy and fill */
        for (int c = 0; c < grid.numCells() * nsubc; c++)
        {
            /* Does this cluster contain perturbed particles? */
            if (grid.clusterIsPerturbed(c))
            {
                const int numAtomsPerCluster = grid.geometry().numAtomsICluster;
                for (int i = 0; i < numAtomsPerCluster; i++)
                {
                    /* Is this a perturbed particle? */
                    if (grid.atomIsPerturbed(c, i))
                    {
                        int ind = c_offset + c * numAtomsPerCluster + i;
                        /* Set atom type and charge to non-interacting */
                        params.type[ind]  = params.numTypes - 1;
                        q[ind * stride_q] = 0;
                    }
                }
            }
        }
    }
}

/* Set the energy group indices for atoms in nbnxn_atomdata_t */
static void nbnxn_atomdata_set_energygroups(const GridSet&          gridSet,
                                            ArrayRef<const int32_t> atomInfo,
                                            EnergyGroupsPerCluster* energyGroupsPerCluster)
{
    for (const Grid& grid : gridSet.grids())
    {
        /* Loop over all columns and copy and fill */
        for (int i = 0; i < grid.numColumns(); i++)
        {
            const int numAtoms   = grid.paddedNumAtomsInColumn(i);
            const int atomOffset = grid.firstAtomInColumn(i);

            energyGroupsPerCluster->setEnergyGroups(gridSet.atomIndices().subArray(atomOffset, numAtoms),
                                                    atomInfo,
                                                    sc_atomInfo_EnergyGroupIdMask,
                                                    grid.atomToCluster(atomOffset));
        }
    }
}

/* Sets all required atom parameter data in nbnxn_atomdata_t */
void nbnxn_atomdata_set(nbnxn_atomdata_t*       nbat,
                        const GridSet&          gridSet,
                        ArrayRef<const int>     atomTypes,
                        ArrayRef<const real>    atomCharges,
                        ArrayRef<const int32_t> atomInfo)
{
    nbnxn_atomdata_t::Params& params = nbat->paramsDeprecated();

    nbnxn_atomdata_set_atomtypes(&params, gridSet, atomTypes);

    nbnxn_atomdata_set_charges(nbat, gridSet, atomCharges);

    if (gridSet.haveFep())
    {
        nbnxn_atomdata_mask_fep(nbat, gridSet);
    }

    /* This must be done after masking types for FEP */
    nbnxn_atomdata_set_ljcombparams(&params, nbat->XFormat, gridSet);

    if (nbat->params().energyGroupsPerCluster)
    {
        nbnxn_atomdata_set_energygroups(gridSet, atomInfo, nbat->params().energyGroupsPerCluster.get());
    }
}

/* Copies the shift vector array to nbnxn_atomdata_t */
void nbnxn_atomdata_copy_shiftvec(gmx_bool bDynamicBox, ArrayRef<RVec> shift_vec, nbnxn_atomdata_t* nbat)
{
    nbat->bDynamicBox = bDynamicBox;
    std::copy(shift_vec.begin(), shift_vec.end(), nbat->shift_vec.begin());
}

// Returns the used range of grids for the given locality
static Range<int> getGridRange(const GridSet& gridSet, const AtomLocality locality)
{
    int gridBegin = 0;
    int gridEnd   = 0;

    switch (locality)
    {
        case AtomLocality::All:
            gridBegin = 0;
            gridEnd   = gridSet.grids().size();
            break;
        case AtomLocality::Local:
            gridBegin = 0;
            gridEnd   = 1;
            break;
        case AtomLocality::NonLocal:
            gridBegin = 1;
            gridEnd   = gridSet.grids().size();
            break;
        default: GMX_ASSERT(false, "Invalid locality specifier"); break;
    }

    return Range<int>(gridBegin, gridEnd);
}

/* Copies (and reorders) the coordinates to nbnxn_atomdata_t */
void nbnxn_atomdata_copy_x_to_nbat_x(const GridSet&     gridSet,
                                     const AtomLocality locality,
                                     const rvec*        coordinates,
                                     nbnxn_atomdata_t*  nbat)
{
    const auto gridRange = getGridRange(gridSet, locality);

    const int nth = gmx_omp_nthreads_get(ModuleMultiThread::Pairsearch);
#pragma omp parallel for num_threads(nth) schedule(static)
    for (int th = 0; th < nth; th++)
    {
        try
        {
            for (int g : gridRange)
            {
                const Grid& grid       = gridSet.grid(g);
                const int   numCellsXY = grid.numColumns();

                const int cxy0 = (numCellsXY * th + nth - 1) / nth;
                const int cxy1 = (numCellsXY * (th + 1) + nth - 1) / nth;

                for (int cxy = cxy0; cxy < cxy1; cxy++)
                {
                    const int na  = grid.numAtomsInColumn(cxy);
                    const int ash = grid.firstAtomInColumn(cxy);

                    copy_rvec_to_nbat_real(gridSet.atomIndices().data() + ash,
                                           na,
                                           na,
                                           coordinates,
                                           nbat->XFormat,
                                           nbat->x().data(),
                                           ash);
                }
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }
}

/* Copies (and reorders) the coordinates to nbnxn_atomdata_t on the GPU*/
void nbnxn_atomdata_x_to_nbat_x_gpu(const GridSet&        gridSet,
                                    const AtomLocality    locality,
                                    NbnxmGpu*             gpu_nbv,
                                    DeviceBuffer<RVec>    d_x,
                                    GpuEventSynchronizer* xReadyOnDevice)
{
    const auto gridRange = getGridRange(gridSet, locality);

    for (int g : gridRange)
    {
        nbnxn_gpu_x_to_nbat_x(gridSet.grids()[g],
                              gpu_nbv,
                              d_x,
                              (g == *gridRange.begin()) ? xReadyOnDevice
                                                        : nullptr, // Sync on first iteration only
                              locality,
                              g,
                              gridSet.numColumnsMax(),
                              (g == *gridRange.end() - 1));
    }
}

static void nbnxn_atomdata_clear_reals(ArrayRef<real> dest, int i0, int i1)
{
    for (int i = i0; i < i1; i++)
    {
        dest[i] = 0;
    }
}

gmx_unused static void nbnxn_atomdata_reduce_reals(real* gmx_restrict        dest,
                                                   gmx_bool                  bDestSet,
                                                   const real** gmx_restrict src,
                                                   int                       nsrc,
                                                   int                       i0,
                                                   int                       i1)
{
    if (bDestSet)
    {
        /* The destination buffer contains data, add to it */
        for (int i = i0; i < i1; i++)
        {
            for (int s = 0; s < nsrc; s++)
            {
                dest[i] += src[s][i];
            }
        }
    }
    else
    {
        /* The destination buffer is uninitialized, set it first */
        for (int i = i0; i < i1; i++)
        {
            dest[i] = src[0][i];
            for (int s = 1; s < nsrc; s++)
            {
                dest[i] += src[s][i];
            }
        }
    }
}

gmx_unused static void nbnxn_atomdata_reduce_reals_simd(real gmx_unused* gmx_restrict dest,
                                                        gmx_bool gmx_unused           bDestSet,
                                                        const gmx_unused real** gmx_restrict src,
                                                        int gmx_unused                       nsrc,
                                                        int gmx_unused                       i0,
                                                        int gmx_unused                       i1)
{
#if GMX_SIMD
    /* The SIMD width here is actually independent of that in the kernels,
     * but we use the same width for simplicity (usually optimal anyhow).
     */
    SimdReal dest_SSE, src_SSE;

    if (bDestSet)
    {
        for (int i = i0; i < i1; i += GMX_SIMD_REAL_WIDTH)
        {
            dest_SSE = load<SimdReal>(dest + i);
            for (int s = 0; s < nsrc; s++)
            {
                src_SSE  = load<SimdReal>(src[s] + i);
                dest_SSE = dest_SSE + src_SSE;
            }
            store(dest + i, dest_SSE);
        }
    }
    else
    {
        for (int i = i0; i < i1; i += GMX_SIMD_REAL_WIDTH)
        {
            dest_SSE = load<SimdReal>(src[0] + i);
            for (int s = 1; s < nsrc; s++)
            {
                src_SSE  = load<SimdReal>(src[s] + i);
                dest_SSE = dest_SSE + src_SSE;
            }
            store(dest + i, dest_SSE);
        }
    }
#endif
}

/* Add part of the force array(s) from nbnxn_atomdata_t to f
 *
 * Note: Adding restrict to f makes this function 50% slower with gcc 7.3
 */
static void nbnxn_atomdata_add_nbat_f_to_f_part(const GridSet&                 gridSet,
                                                const nbnxn_atomdata_t&        nbat,
                                                const nbnxn_atomdata_output_t& out,
                                                const int                      a0,
                                                const int                      a1,
                                                rvec*                          f)
{
    ArrayRef<const int> cell = gridSet.cells();
    // Note: Using ArrayRef instead makes this code 25% slower with gcc 7.3
    const real* fnb = out.f.data();

    /* Loop over all columns and copy and fill */
    switch (nbat.FFormat)
    {
        case nbatXYZ:
        case nbatXYZQ:
            for (int a = a0; a < a1; a++)
            {
                int i = cell[a] * nbat.fstride;

                f[a][XX] += fnb[i];
                f[a][YY] += fnb[i + 1];
                f[a][ZZ] += fnb[i + 2];
            }
            break;
        case nbatX4:
            for (int a = a0; a < a1; a++)
            {
                int i = atom_to_x_index<c_packX4>(cell[a]);

                f[a][XX] += fnb[i + XX * c_packX4];
                f[a][YY] += fnb[i + YY * c_packX4];
                f[a][ZZ] += fnb[i + ZZ * c_packX4];
            }
            break;
        case nbatX8:
            for (int a = a0; a < a1; a++)
            {
                int i = atom_to_x_index<c_packX8>(cell[a]);

                f[a][XX] += fnb[i + XX * c_packX8];
                f[a][YY] += fnb[i + YY * c_packX8];
                f[a][ZZ] += fnb[i + ZZ * c_packX8];
            }
            break;
        default: gmx_incons("Unsupported nbnxn_atomdata_t format");
    }
}

void nbnxn_atomdata_t::reduceForcesOverThreads()
{
    // The number of output buffers should match the number of OpenMP threads
    const int nth = gmx::ssize(outputBuffers_);

#pragma omp parallel for num_threads(nth) schedule(static)
    for (int th = 0; th < nth; th++)
    {
        try
        {
            const real* fptr[NBNXN_BUFFERFLAG_MAX_THREADS];

            ArrayRef<const gmx_bitmask_t> flags = bufferFlags_;

            /* Calculate the cell-block range for our thread */
            const int b0 = (flags.size() * th) / nth;
            const int b1 = (flags.size() * (th + 1)) / nth;

            for (int b = b0; b < b1; b++)
            {
                const int i0 = b * NBNXN_BUFFERFLAG_SIZE * fstride;
                const int i1 = (b + 1) * NBNXN_BUFFERFLAG_SIZE * fstride;

                int nfptr = 0;
                for (Index out = 1; out < gmx::ssize(outputBuffers_); out++)
                {
                    if (bitmask_is_set(flags[b], out))
                    {
                        fptr[nfptr++] = outputBuffers_[out].f.data();
                    }
                }
                if (nfptr > 0)
                {
#if GMX_SIMD
                    nbnxn_atomdata_reduce_reals_simd
#else
                    nbnxn_atomdata_reduce_reals
#endif
                            (outputBuffers_[0].f.data(), bitmask_is_set(flags[b], 0), fptr, nfptr, i0, i1);
                }
                else if (!bitmask_is_set(flags[b], 0))
                {
                    nbnxn_atomdata_clear_reals(outputBuffers_[0].f, i0, i1);
                }
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }
}

// Return the atom range for the given locality
static Range<int> getAtomRange(const AtomLocality locality, const GridSet& gridSet)
{
    int atomStart = 0;
    int atomEnd   = 0;

    switch (locality)
    {
        case AtomLocality::All:
            atomStart = 0;
            atomEnd   = gridSet.numRealAtomsTotal();
            break;
        case AtomLocality::Local:
            atomStart = 0;
            atomEnd   = gridSet.numRealAtomsLocal();
            break;
        case AtomLocality::NonLocal:
            atomStart = gridSet.numRealAtomsLocal();
            atomEnd   = gridSet.numRealAtomsTotal();
            break;
        default: GMX_ASSERT(false, "Invalid locality specifier"); break;
    }

    return Range<int>(atomStart, atomEnd);
}

/* Add the force array(s) from nbnxn_atomdata_t to f */
void nbnxn_atomdata_t::reduceForces(const AtomLocality locality, const GridSet& gridSet, rvec* f)
{
    const auto atomRange = getAtomRange(locality, gridSet);

    if (atomRange.empty())
    {
        /* The are no atoms for this reduction, avoid some overhead */
        return;
    }

    int nth = gmx_omp_nthreads_get(ModuleMultiThread::Nonbonded);

    if (outputBuffers_.size() > 1)
    {
        if (locality != AtomLocality::All)
        {
            gmx_incons("add_f_to_f called with nout>1 and locality!=eatAll");
        }

        /* Reduce the force thread output buffers into buffer 0, before adding
         * them to the, differently ordered, "real" force buffer.
         */
        reduceForcesOverThreads();
    }
#pragma omp parallel for num_threads(nth) schedule(static)
    for (int th = 0; th < nth; th++)
    {
        try
        {
            nbnxn_atomdata_add_nbat_f_to_f_part(gridSet,
                                                *this,
                                                outputBuffers_[0],
                                                *atomRange.begin() + ((th + 0) * atomRange.size()) / nth,
                                                *atomRange.begin() + ((th + 1) * atomRange.size()) / nth,
                                                f);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }
}

void nbnxn_atomdata_add_nbat_fshift_to_fshift(const nbnxn_atomdata_t& nbat, ArrayRef<RVec> fshift)
{
    ArrayRef<const nbnxn_atomdata_output_t> outputBuffers = nbat.outputBuffers();

    for (int s = 0; s < c_numShiftVectors; s++)
    {
        rvec sum;
        clear_rvec(sum);
        for (const nbnxn_atomdata_output_t& out : outputBuffers)
        {
            sum[XX] += out.fshift[s * DIM + XX];
            sum[YY] += out.fshift[s * DIM + YY];
            sum[ZZ] += out.fshift[s * DIM + ZZ];
        }
        fshift[s] += sum;
    }
}

//! Clears all elements of buffer
static void clearBufferAll(ArrayRef<real> buffer)
{
    std::fill(buffer.begin(), buffer.end(), 0.0_real);
}

/*! \brief Clears elements of size and stride \p numComponentsPerElement
 *
 * Only elements with flags in \p nbat set for index \p outputIndex
 * are cleared.
 */
template<int numComponentsPerElement>
static void clearBufferFlagged(const int outputIndex, ArrayRef<const gmx_bitmask_t> flags, ArrayRef<real> buffer)
{
    gmx_bitmask_t our_flag; // NOLINT(cppcoreguidelines-init-variables)
    bitmask_init_bit(&our_flag, outputIndex);

    constexpr size_t numComponentsPerBlock = NBNXN_BUFFERFLAG_SIZE * numComponentsPerElement;

    for (size_t b = 0; b < flags.size(); b++)
    {
        if (!bitmask_is_disjoint(flags[b], our_flag))
        {
            clearBufferAll(buffer.subArray(b * numComponentsPerBlock, numComponentsPerBlock));
        }
    }
}

void nbnxn_atomdata_t::clearForceBuffer(const int outputIndex)
{
    if (useBufferFlags_)
    {
        GMX_ASSERT(fstride == DIM, "Only fstride=3 is currently handled here");

        clearBufferFlagged<DIM>(outputIndex, bufferFlags_, outputBuffers_[outputIndex].f);
    }
    else
    {
        clearBufferAll(outputBuffers_[outputIndex].f);
    }
}

} // namespace gmx
