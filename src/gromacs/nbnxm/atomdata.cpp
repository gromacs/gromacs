/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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

#include "gmxpre.h"

#include "atomdata.h"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <algorithm>

#include "thread_mpi/atomic.h"

#include "gromacs/math/functions.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/forcerec.h" // only for GET_CGINFO_*
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/nbnxm_geometry.h"
#include "gromacs/nbnxm/pairlist.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/simd/simd.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringutil.h"

#include "grid.h"
#include "internal.h"

using namespace gmx; // TODO: Remove when this file is moved into gmx namespace

void nbnxn_atomdata_t::resizeCoordinateBuffer(int numAtoms)
{
    numAtoms_ = numAtoms;

    x_.resize(numAtoms*xstride);
}

void nbnxn_atomdata_t::resizeForceBuffers()
{
    /* Force buffers need padding up to a multiple of the buffer flag size */
    const int paddedSize = (numAtoms() + NBNXN_BUFFERFLAG_SIZE - 1)/NBNXN_BUFFERFLAG_SIZE*NBNXN_BUFFERFLAG_SIZE;

    /* Should we let each thread allocate it's own data instead? */
    for (nbnxn_atomdata_output_t &outBuffer : out)
    {
        outBuffer.f.resize(paddedSize*fstride);
    }
}

/* Initializes an nbnxn_atomdata_output_t data structure */
nbnxn_atomdata_output_t::nbnxn_atomdata_output_t(Nbnxm::KernelType  kernelType,
                                                 int                numEnergyGroups,
                                                 int                simdEnergyBufferStride,
                                                 gmx::PinningPolicy pinningPolicy) :
    f({}, {pinningPolicy}),
    fshift({}, {pinningPolicy}),
    Vvdw({}, {pinningPolicy}),
    Vc({}, {pinningPolicy})
{
    fshift.resize(SHIFTS*DIM);
    Vvdw.resize(numEnergyGroups*numEnergyGroups);
    Vc.resize(numEnergyGroups*numEnergyGroups);

    if (Nbnxm::kernelTypeIsSimd(kernelType))
    {
        int cj_size     = Nbnxm::JClusterSizePerKernelType[kernelType];
        int numElements = numEnergyGroups*numEnergyGroups*simdEnergyBufferStride*(cj_size/2)*cj_size;
        VSvdw.resize(numElements);
        VSc.resize(numElements);
    }
}

static void copy_int_to_nbat_int(const int *a, int na, int na_round,
                                 const int *in, int fill, int *innb)
{
    int i, j;

    j = 0;
    for (i = 0; i < na; i++)
    {
        innb[j++] = in[a[i]];
    }
    /* Complete the partially filled last cell with fill */
    for (; i < na_round; i++)
    {
        innb[j++] = fill;
    }
}

void copy_rvec_to_nbat_real(const int *a, int na, int na_round,
                            const rvec *x, int nbatFormat,
                            real *xnb, int a0)
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

    int        i, j, c;

    switch (nbatFormat)
    {
        case nbatXYZ:
            j = a0*STRIDE_XYZ;
            for (i = 0; i < na; i++)
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
            break;
        case nbatXYZQ:
            j = a0*STRIDE_XYZQ;
            for (i = 0; i < na; i++)
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
            break;
        case nbatX4:
            j = atom_to_x_index<c_packX4>(a0);
            c = a0 & (c_packX4-1);
            for (i = 0; i < na; i++)
            {
                xnb[j+XX*c_packX4] = x[a[i]][XX];
                xnb[j+YY*c_packX4] = x[a[i]][YY];
                xnb[j+ZZ*c_packX4] = x[a[i]][ZZ];
                j++;
                c++;
                if (c == c_packX4)
                {
                    j += (DIM-1)*c_packX4;
                    c  = 0;
                }
            }
            /* Complete the partially filled last cell with zeros */
            for (; i < na_round; i++)
            {
                xnb[j+XX*c_packX4] = farAway;
                xnb[j+YY*c_packX4] = farAway;
                xnb[j+ZZ*c_packX4] = farAway;
                j++;
                c++;
                if (c == c_packX4)
                {
                    j += (DIM-1)*c_packX4;
                    c  = 0;
                }
            }
            break;
        case nbatX8:
            j = atom_to_x_index<c_packX8>(a0);
            c = a0 & (c_packX8 - 1);
            for (i = 0; i < na; i++)
            {
                xnb[j+XX*c_packX8] = x[a[i]][XX];
                xnb[j+YY*c_packX8] = x[a[i]][YY];
                xnb[j+ZZ*c_packX8] = x[a[i]][ZZ];
                j++;
                c++;
                if (c == c_packX8)
                {
                    j += (DIM-1)*c_packX8;
                    c  = 0;
                }
            }
            /* Complete the partially filled last cell with zeros */
            for (; i < na_round; i++)
            {
                xnb[j+XX*c_packX8] = farAway;
                xnb[j+YY*c_packX8] = farAway;
                xnb[j+ZZ*c_packX8] = farAway;
                j++;
                c++;
                if (c == c_packX8)
                {
                    j += (DIM-1)*c_packX8;
                    c  = 0;
                }
            }
            break;
        default:
            gmx_incons("Unsupported nbnxn_atomdata_t format");
    }
}

/* Stores the LJ parameter data in a format convenient for different kernels */
static void set_lj_parameter_data(nbnxn_atomdata_t::Params *params, gmx_bool bSIMD)
{
    int  nt = params->numTypes;

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
        params->nbfp_aligned.resize(nt*nt*c_simdBestPairAlignment);

        for (int i = 0; i < nt; i++)
        {
            for (int j = 0; j < nt; j++)
            {
                params->nbfp_aligned[(i*nt+j)*c_simdBestPairAlignment+0] = params->nbfp[(i*nt+j)*2+0];
                params->nbfp_aligned[(i*nt+j)*c_simdBestPairAlignment+1] = params->nbfp[(i*nt+j)*2+1];
                params->nbfp_aligned[(i*nt+j)*c_simdBestPairAlignment+2] = 0;
                params->nbfp_aligned[(i*nt+j)*c_simdBestPairAlignment+3] = 0;
            }
        }
#endif
    }

    /* We use combination rule data for SIMD combination rule kernels
     * and with LJ-PME kernels. We then only need parameters per atom type,
     * not per pair of atom types.
     */
    params->nbfp_comb.resize(nt*2);
    switch (params->comb_rule)
    {
        case ljcrGEOM:
            params->comb_rule = ljcrGEOM;

            for (int i = 0; i < nt; i++)
            {
                /* Store the sqrt of the diagonal from the nbfp matrix */
                params->nbfp_comb[i*2  ] = std::sqrt(params->nbfp[(i*nt+i)*2  ]);
                params->nbfp_comb[i*2+1] = std::sqrt(params->nbfp[(i*nt+i)*2+1]);
            }
            break;
        case ljcrLB:
            for (int i = 0; i < nt; i++)
            {
                /* Get 6*C6 and 12*C12 from the diagonal of the nbfp matrix */
                const real c6  = params->nbfp[(i*nt+i)*2  ];
                const real c12 = params->nbfp[(i*nt+i)*2+1];
                if (c6 > 0 && c12 > 0)
                {
                    /* We store 0.5*2^1/6*sigma and sqrt(4*3*eps),
                     * so we get 6*C6 and 12*C12 after combining.
                     */
                    params->nbfp_comb[i*2  ] = 0.5*gmx::sixthroot(c12/c6);
                    params->nbfp_comb[i*2+1] = std::sqrt(c6*c6/c12);
                }
                else
                {
                    params->nbfp_comb[i*2  ] = 0;
                    params->nbfp_comb[i*2+1] = 0;
                }
            }
            break;
        case ljcrNONE:
            /* We always store the full matrix (see code above) */
            break;
        default:
            gmx_incons("Unknown combination rule");
    }
}

nbnxn_atomdata_t::SimdMasks::SimdMasks()
{
#if GMX_SIMD
    constexpr int simd_width = GMX_SIMD_REAL_WIDTH;
    /* Set the diagonal cluster pair exclusion mask setup data.
     * In the kernel we check 0 < j - i to generate the masks.
     * Here we store j - i for generating the mask for the first i,
     * we substract 0.5 to avoid rounding issues.
     * In the kernel we can subtract 1 to generate the subsequent mask.
     */
    const int simd_4xn_diag_size = std::max(c_nbnxnCpuIClusterSize, simd_width);
    diagonal_4xn_j_minus_i.resize(simd_4xn_diag_size);
    for (int j = 0; j < simd_4xn_diag_size; j++)
    {
        diagonal_4xn_j_minus_i[j] = j - 0.5;
    }

    diagonal_2xnn_j_minus_i.resize(simd_width);
    for (int j = 0; j < simd_width/2; j++)
    {
        /* The j-cluster size is half the SIMD width */
        diagonal_2xnn_j_minus_i[j]                = j - 0.5;
        /* The next half of the SIMD width is for i + 1 */
        diagonal_2xnn_j_minus_i[simd_width/2 + j] = j - 1 - 0.5;
    }

    /* We use up to 32 bits for exclusion masking.
     * The same masks are used for the 4xN and 2x(N+N) kernels.
     * The masks are read either into integer SIMD registers or into
     * real SIMD registers (together with a cast).
     * In single precision this means the real and integer SIMD registers
     * are of equal size.
     */
    const int simd_excl_size = c_nbnxnCpuIClusterSize*simd_width;
#if GMX_DOUBLE && !GMX_SIMD_HAVE_INT32_LOGICAL
    exclusion_filter64.resize(simd_excl_size);
#else
    exclusion_filter.resize(simd_excl_size);
#endif

    for (int j = 0; j < simd_excl_size; j++)
    {
        /* Set the consecutive bits for masking pair exclusions */
#if GMX_DOUBLE && !GMX_SIMD_HAVE_INT32_LOGICAL
        exclusion_filter64[j] = (1U << j);
#else
        exclusion_filter[j]   = (1U << j);
#endif
    }

    if (!GMX_SIMD_HAVE_LOGICAL && !GMX_SIMD_HAVE_INT32_LOGICAL) // NOLINT(misc-redundant-expression)
    {
        // If the SIMD implementation has no bitwise logical operation support
        // whatsoever we cannot use the normal masking. Instead,
        // we generate a vector of all 2^4 possible ways an i atom
        // interacts with its 4 j atoms. Each array entry contains
        // GMX_SIMD_REAL_WIDTH values that are read with a single aligned SIMD load.
        // Since there is no logical value representation in this case, we use
        // any nonzero value to indicate 'true', while zero mean 'false'.
        // This can then be converted to a SIMD boolean internally in the SIMD
        // module by comparing to zero.
        // Each array entry encodes how this i atom will interact with the 4 j atoms.
        // Matching code exists in set_ci_top_excls() to generate indices into this array.
        // Those indices are used in the kernels.

        const int  simd_excl_size = c_nbnxnCpuIClusterSize*c_nbnxnCpuIClusterSize;
        const real simdFalse      = 0.0;
        const real simdTrue       = 1.0;

        interaction_array.resize(simd_excl_size * GMX_SIMD_REAL_WIDTH);
        for (int j = 0; j < simd_excl_size; j++)
        {
            const int index = j * GMX_SIMD_REAL_WIDTH;
            for (int i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
            {
                interaction_array[index + i] = (j & (1 << i)) ? simdTrue : simdFalse;
            }
        }
    }
#endif
}

nbnxn_atomdata_t::Params::Params(gmx::PinningPolicy pinningPolicy) :
    numTypes(0),
    nbfp({}, {pinningPolicy}),
    nbfp_comb({}, {pinningPolicy}),
    type({}, {pinningPolicy}),
    lj_comb({}, {pinningPolicy}),
    q({}, {pinningPolicy}),
    nenergrp(0),
    neg_2log(0),
    energrp({}, {pinningPolicy})
{
}

nbnxn_atomdata_t::nbnxn_atomdata_t(gmx::PinningPolicy pinningPolicy) :
    params_(pinningPolicy),
    numAtoms_(0),
    natoms_local(0),
    shift_vec({}, {pinningPolicy}),
    x_({}, {pinningPolicy}),
    simdMasks(),
    bUseBufferFlags(FALSE),
    bUseTreeReduce(FALSE)
{
}

/* Initializes an nbnxn_atomdata_t::Params data structure */
static void nbnxn_atomdata_params_init(const gmx::MDLogger &mdlog,
                                       nbnxn_atomdata_t::Params *params,
                                       const Nbnxm::KernelType kernelType,
                                       int enbnxninitcombrule,
                                       int ntype, const real *nbfp,
                                       int n_energygroups)
{
    real     c6, c12, tol;
    char    *ptr;
    gmx_bool simple, bCombGeom, bCombLB, bSIMD;

    if (debug)
    {
        fprintf(debug, "There are %d atom types in the system, adding one for nbnxn_atomdata_t\n", ntype);
    }
    params->numTypes = ntype + 1;
    params->nbfp.resize(params->numTypes*params->numTypes*2);
    params->nbfp_comb.resize(params->numTypes*2);

    /* A tolerance of 1e-5 seems reasonable for (possibly hand-typed)
     * force-field floating point parameters.
     */
    tol = 1e-5;
    ptr = getenv("GMX_LJCOMB_TOL");
    if (ptr != nullptr)
    {
        double dbl;

        sscanf(ptr, "%lf", &dbl);
        tol = dbl;
    }
    bCombGeom = TRUE;
    bCombLB   = TRUE;

    /* Temporarily fill params->nbfp_comb with sigma and epsilon
     * to check for the LB rule.
     */
    for (int i = 0; i < ntype; i++)
    {
        c6  = nbfp[(i*ntype+i)*2    ]/6.0;
        c12 = nbfp[(i*ntype+i)*2 + 1]/12.0;
        if (c6 > 0 && c12 > 0)
        {
            params->nbfp_comb[i*2    ] = gmx::sixthroot(c12/c6);
            params->nbfp_comb[i*2 + 1] = 0.25*c6*c6/c12;
        }
        else if (c6 == 0 && c12 == 0)
        {
            params->nbfp_comb[i*2    ] = 0;
            params->nbfp_comb[i*2 + 1] = 0;
        }
        else
        {
            /* Can not use LB rule with only dispersion or repulsion */
            bCombLB = FALSE;
        }
    }

    for (int i = 0; i < params->numTypes; i++)
    {
        for (int j = 0; j < params->numTypes; j++)
        {
            if (i < ntype && j < ntype)
            {
                /* fr->nbfp has been updated, so that array too now stores c6/c12 including
                 * the 6.0/12.0 prefactors to save 2 flops in the most common case (force-only).
                 */
                c6  = nbfp[(i*ntype+j)*2    ];
                c12 = nbfp[(i*ntype+j)*2 + 1];
                params->nbfp[(i*params->numTypes+j)*2    ] = c6;
                params->nbfp[(i*params->numTypes+j)*2 + 1] = c12;

                /* Compare 6*C6 and 12*C12 for geometric cobination rule */
                bCombGeom = bCombGeom &&
                    gmx_within_tol(c6*c6, nbfp[(i*ntype+i)*2  ]*nbfp[(j*ntype+j)*2  ], tol) &&
                    gmx_within_tol(c12*c12, nbfp[(i*ntype+i)*2 + 1]*nbfp[(j*ntype+j)*2 + 1], tol);

                /* Compare C6 and C12 for Lorentz-Berthelot combination rule */
                c6     /= 6.0;
                c12    /= 12.0;
                bCombLB = bCombLB &&
                    ((c6 == 0 && c12 == 0 &&
                      (params->nbfp_comb[i*2 + 1] == 0 || params->nbfp_comb[j*2 + 1] == 0)) ||
                     (c6 > 0 && c12 > 0 &&
                      gmx_within_tol(gmx::sixthroot(c12/c6),
                                     0.5*(params->nbfp_comb[i*2]+params->nbfp_comb[j*2]), tol) &&
                      gmx_within_tol(0.25*c6*c6/c12, std::sqrt(params->nbfp_comb[i*2 + 1]*params->nbfp_comb[j*2 + 1]), tol)));
            }
            else
            {
                /* Add zero parameters for the additional dummy atom type */
                params->nbfp[(i*params->numTypes + j)*2  ] = 0;
                params->nbfp[(i*params->numTypes + j)*2+1] = 0;
            }
        }
    }
    if (debug)
    {
        fprintf(debug, "Combination rules: geometric %s Lorentz-Berthelot %s\n",
                gmx::boolToString(bCombGeom), gmx::boolToString(bCombLB));
    }

    simple = Nbnxm::kernelTypeUsesSimplePairlist(kernelType);

    switch (enbnxninitcombrule)
    {
        case enbnxninitcombruleDETECT:
            /* We prefer the geometic combination rule,
             * as that gives a slightly faster kernel than the LB rule.
             */
            if (bCombGeom)
            {
                params->comb_rule = ljcrGEOM;
            }
            else if (bCombLB)
            {
                params->comb_rule = ljcrLB;
            }
            else
            {
                params->comb_rule = ljcrNONE;

                params->nbfp_comb.clear();
            }

            {
                std::string mesg;
                if (params->comb_rule == ljcrNONE)
                {
                    mesg = "Using full Lennard-Jones parameter combination matrix";
                }
                else
                {
                    mesg = gmx::formatString("Using %s Lennard-Jones combination rule",
                                             params->comb_rule == ljcrGEOM ? "geometric" : "Lorentz-Berthelot");
                }
                GMX_LOG(mdlog.info).asParagraph().appendText(mesg);
            }
            break;
        case enbnxninitcombruleGEOM:
            params->comb_rule = ljcrGEOM;
            break;
        case enbnxninitcombruleLB:
            params->comb_rule = ljcrLB;
            break;
        case enbnxninitcombruleNONE:
            params->comb_rule = ljcrNONE;

            params->nbfp_comb.clear();
            break;
        default:
            gmx_incons("Unknown enbnxninitcombrule");
    }

    bSIMD = Nbnxm::kernelTypeIsSimd(kernelType);

    set_lj_parameter_data(params, bSIMD);

    params->nenergrp = n_energygroups;
    if (!simple)
    {
        // We now check for energy groups already when starting mdrun
        GMX_RELEASE_ASSERT(n_energygroups == 1, "GPU kernels do not support energy groups");
    }
    /* Temporary storage goes as #grp^3*simd_width^2/2, so limit to 64 */
    if (params->nenergrp > 64)
    {
        gmx_fatal(FARGS, "With NxN kernels not more than 64 energy groups are supported\n");
    }
    params->neg_2log = 1;
    while (params->nenergrp > (1<<params->neg_2log))
    {
        params->neg_2log++;
    }
}

/* Initializes an nbnxn_atomdata_t data structure */
void nbnxn_atomdata_init(const gmx::MDLogger &mdlog,
                         nbnxn_atomdata_t *nbat,
                         const Nbnxm::KernelType kernelType,
                         int enbnxninitcombrule,
                         int ntype, const real *nbfp,
                         int n_energygroups,
                         int nout)
{
    nbnxn_atomdata_params_init(mdlog, &nbat->paramsDeprecated(), kernelType,
                               enbnxninitcombrule, ntype, nbfp, n_energygroups);

    const bool simple = Nbnxm::kernelTypeUsesSimplePairlist(kernelType);
    const bool bSIMD  = Nbnxm::kernelTypeIsSimd(kernelType);

    if (simple)
    {
        int pack_x;

        if (bSIMD)
        {
            pack_x = std::max(c_nbnxnCpuIClusterSize,
                              Nbnxm::JClusterSizePerKernelType[kernelType]);
            switch (pack_x)
            {
                case 4:
                    nbat->XFormat = nbatX4;
                    break;
                case 8:
                    nbat->XFormat = nbatX8;
                    break;
                default:
                    gmx_incons("Unsupported packing width");
            }
        }
        else
        {
            nbat->XFormat = nbatXYZ;
        }

        nbat->FFormat = nbat->XFormat;
    }
    else
    {
        nbat->XFormat = nbatXYZQ;
        nbat->FFormat = nbatXYZ;
    }

    nbat->shift_vec.resize(SHIFTS);

    nbat->xstride = (nbat->XFormat == nbatXYZQ ? STRIDE_XYZQ : DIM);
    nbat->fstride = (nbat->FFormat == nbatXYZQ ? STRIDE_XYZQ : DIM);

    /* Initialize the output data structures */
    for (int i = 0; i < nout; i++)
    {
        const auto &pinningPolicy = nbat->params().type.get_allocator().pinningPolicy();
        nbat->out.emplace_back(kernelType, nbat->params().nenergrp, 1 << nbat->params().neg_2log,
                               pinningPolicy);
    }

    nbat->buffer_flags.flag        = nullptr;
    nbat->buffer_flags.flag_nalloc = 0;

    const int   nth = gmx_omp_nthreads_get(emntNonbonded);

    const char *ptr = getenv("GMX_USE_TREEREDUCE");
    if (ptr != nullptr)
    {
        nbat->bUseTreeReduce = (strtol(ptr, nullptr, 10) != 0);
    }
#if defined __MIC__
    else if (nth > 8) /*on the CPU we currently don't benefit even at 32*/
    {
        nbat->bUseTreeReduce = 1;
    }
#endif
    else
    {
        nbat->bUseTreeReduce = false;
    }
    if (nbat->bUseTreeReduce)
    {
        GMX_LOG(mdlog.info).asParagraph().appendText("Using tree force reduction");

        nbat->syncStep = new tMPI_Atomic[nth];
    }
}

template<int packSize>
static void copy_lj_to_nbat_lj_comb(gmx::ArrayRef<const real> ljparam_type,
                                    const int *type, int na,
                                    real *ljparam_at)
{
    /* The LJ params follow the combination rule:
     * copy the params for the type array to the atom array.
     */
    for (int is = 0; is < na; is += packSize)
    {
        for (int k = 0; k < packSize; k++)
        {
            int i = is + k;
            ljparam_at[is*2            + k] = ljparam_type[type[i]*2    ];
            ljparam_at[is*2 + packSize + k] = ljparam_type[type[i]*2 + 1];
        }
    }
}

static int numAtomsFromGrids(const nbnxn_search &nbs)
{
    return nbs.grid.back().atomIndexEnd();
}

/* Sets the atom type in nbnxn_atomdata_t */
static void nbnxn_atomdata_set_atomtypes(nbnxn_atomdata_t::Params *params,
                                         const nbnxn_search       *nbs,
                                         const int                *type)
{
    params->type.resize(numAtomsFromGrids(*nbs));

    for (const Nbnxm::Grid &grid : nbs->grid)
    {
        /* Loop over all columns and copy and fill */
        for (int i = 0; i < grid.numColumns(); i++)
        {
            const int numAtoms   = grid.paddedNumAtomsInColumn(i);
            const int atomOffset = grid.firstAtomInColumn(i);

            copy_int_to_nbat_int(nbs->a.data() + atomOffset, grid.numAtomsInColumn(i), numAtoms,
                                 type, params->numTypes - 1, params->type.data() + atomOffset);
        }
    }
}

/* Sets the LJ combination rule parameters in nbnxn_atomdata_t */
static void nbnxn_atomdata_set_ljcombparams(nbnxn_atomdata_t::Params *params,
                                            const int                 XFormat,
                                            const nbnxn_search       *nbs)
{
    params->lj_comb.resize(numAtomsFromGrids(*nbs)*2);

    if (params->comb_rule != ljcrNONE)
    {
        for (const Nbnxm::Grid &grid : nbs->grid)
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
                                                      params->lj_comb.data() + atomOffset*2);
                }
                else if (XFormat == nbatX8)
                {
                    copy_lj_to_nbat_lj_comb<c_packX8>(params->nbfp_comb,
                                                      params->type.data() + atomOffset,
                                                      numAtoms,
                                                      params->lj_comb.data() + atomOffset*2);
                }
                else if (XFormat == nbatXYZQ)
                {
                    copy_lj_to_nbat_lj_comb<1>(params->nbfp_comb,
                                               params->type.data() + atomOffset,
                                               numAtoms,
                                               params->lj_comb.data() + atomOffset*2);
                }
            }
        }
    }
}

/* Sets the charges in nbnxn_atomdata_t *nbat */
static void nbnxn_atomdata_set_charges(nbnxn_atomdata_t    *nbat,
                                       const nbnxn_search  *nbs,
                                       const real          *charge)
{
    if (nbat->XFormat != nbatXYZQ)
    {
        nbat->paramsDeprecated().q.resize(nbat->numAtoms());
    }

    for (const Nbnxm::Grid &grid : nbs->grid)
    {
        /* Loop over all columns and copy and fill */
        for (int cxy = 0; cxy < grid.numColumns(); cxy++)
        {
            const int atomOffset     = grid.firstAtomInColumn(cxy);
            const int numAtoms       = grid.numAtomsInColumn(cxy);
            const int paddedNumAtoms = grid.paddedNumAtomsInColumn(cxy);

            if (nbat->XFormat == nbatXYZQ)
            {
                real *q = nbat->x().data() + atomOffset*STRIDE_XYZQ + ZZ + 1;
                int   i;
                for (i = 0; i < numAtoms; i++)
                {
                    *q = charge[nbs->a[atomOffset + i]];
                    q += STRIDE_XYZQ;
                }
                /* Complete the partially filled last cell with zeros */
                for (; i < paddedNumAtoms; i++)
                {
                    *q = 0;
                    q += STRIDE_XYZQ;
                }
            }
            else
            {
                real *q = nbat->paramsDeprecated().q.data() + atomOffset;
                int   i;
                for (i = 0; i < numAtoms; i++)
                {
                    *q = charge[nbs->a[atomOffset + i]];
                    q++;
                }
                /* Complete the partially filled last cell with zeros */
                for (; i < paddedNumAtoms; i++)
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
static void nbnxn_atomdata_mask_fep(nbnxn_atomdata_t    *nbat,
                                    const nbnxn_search  *nbs)
{
    nbnxn_atomdata_t::Params &params = nbat->paramsDeprecated();
    real                     *q;
    int                       stride_q;

    if (nbat->XFormat == nbatXYZQ)
    {
        q        = nbat->x().data() + ZZ + 1;
        stride_q = STRIDE_XYZQ;
    }
    else
    {
        q        = params.q.data();
        stride_q = 1;
    }

    for (const Nbnxm::Grid &grid : nbs->grid)
    {
        int nsubc;
        if (grid.geometry().isSimple)
        {
            nsubc = 1;
        }
        else
        {
            nsubc = c_gpuNumClusterPerCell;
        }

        int c_offset = grid.firstAtomInColumn(0);

        /* Loop over all columns and copy and fill */
        for (int c = 0; c < grid.numCells()*nsubc; c++)
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
                        int ind = c_offset + c*numAtomsPerCluster + i;
                        /* Set atom type and charge to non-interacting */
                        params.type[ind] = params.numTypes - 1;
                        q[ind*stride_q]  = 0;
                    }
                }
            }
        }
    }
}

/* Copies the energy group indices to a reordered and packed array */
static void copy_egp_to_nbat_egps(const int *a, int na, int na_round,
                                  int na_c, int bit_shift,
                                  const int *in, int *innb)
{
    int i;
    int comb;

    int j = 0;
    for (i = 0; i < na; i += na_c)
    {
        /* Store na_c energy group numbers into one int */
        comb = 0;
        for (int sa = 0; sa < na_c; sa++)
        {
            int at = a[i+sa];
            if (at >= 0)
            {
                comb |= (GET_CGINFO_GID(in[at]) << (sa*bit_shift));
            }
        }
        innb[j++] = comb;
    }
    /* Complete the partially filled last cell with fill */
    for (; i < na_round; i += na_c)
    {
        innb[j++] = 0;
    }
}

/* Set the energy group indices for atoms in nbnxn_atomdata_t */
static void nbnxn_atomdata_set_energygroups(nbnxn_atomdata_t::Params *params,
                                            const nbnxn_search       *nbs,
                                            const int                *atinfo)
{
    if (params->nenergrp == 1)
    {
        return;
    }

    params->energrp.resize(numAtomsFromGrids(*nbs));

    for (const Nbnxm::Grid &grid : nbs->grid)
    {
        /* Loop over all columns and copy and fill */
        for (int i = 0; i < grid.numColumns(); i++)
        {
            const int numAtoms   = grid.paddedNumAtomsInColumn(i);
            const int atomOffset = grid.firstAtomInColumn(i);

            copy_egp_to_nbat_egps(nbs->a.data() + atomOffset, grid.numAtomsInColumn(i), numAtoms,
                                  c_nbnxnCpuIClusterSize, params->neg_2log,
                                  atinfo,
                                  params->energrp.data() + grid.atomToCluster(atomOffset));
        }
    }
}

/* Sets all required atom parameter data in nbnxn_atomdata_t */
void nbnxn_atomdata_set(nbnxn_atomdata_t    *nbat,
                        const nbnxn_search  *nbs,
                        const t_mdatoms     *mdatoms,
                        const int           *atinfo)
{
    nbnxn_atomdata_t::Params &params = nbat->paramsDeprecated();

    nbnxn_atomdata_set_atomtypes(&params, nbs, mdatoms->typeA);

    nbnxn_atomdata_set_charges(nbat, nbs, mdatoms->chargeA);

    if (nbs->bFEP)
    {
        nbnxn_atomdata_mask_fep(nbat, nbs);
    }

    /* This must be done after masking types for FEP */
    nbnxn_atomdata_set_ljcombparams(&params, nbat->XFormat, nbs);

    nbnxn_atomdata_set_energygroups(&params, nbs, atinfo);
}

/* Copies the shift vector array to nbnxn_atomdata_t */
void nbnxn_atomdata_copy_shiftvec(gmx_bool          bDynamicBox,
                                  rvec             *shift_vec,
                                  nbnxn_atomdata_t *nbat)
{
    int i;

    nbat->bDynamicBox = bDynamicBox;
    for (i = 0; i < SHIFTS; i++)
    {
        copy_rvec(shift_vec[i], nbat->shift_vec[i]);
    }
}

/* Copies (and reorders) the coordinates to nbnxn_atomdata_t */
void nbnxn_atomdata_copy_x_to_nbat_x(const nbnxn_search       *nbs,
                                     const Nbnxm::AtomLocality locality,
                                     gmx_bool                  FillLocal,
                                     rvec                     *x,
                                     nbnxn_atomdata_t         *nbat,
                                     gmx_wallcycle            *wcycle)
{
    wallcycle_start(wcycle, ewcNB_XF_BUF_OPS);
    wallcycle_sub_start(wcycle, ewcsNB_X_BUF_OPS);

    int g0 = 0, g1 = 0;
    int nth, th;

    switch (locality)
    {
        case Nbnxm::AtomLocality::All:
        case Nbnxm::AtomLocality::Count:
            g0 = 0;
            g1 = nbs->grid.size();
            break;
        case Nbnxm::AtomLocality::Local:
            g0 = 0;
            g1 = 1;
            break;
        case Nbnxm::AtomLocality::NonLocal:
            g0 = 1;
            g1 = nbs->grid.size();
            break;
    }

    if (FillLocal)
    {
        nbat->natoms_local = nbs->grid[0].atomIndexEnd();
    }

    nth = gmx_omp_nthreads_get(emntPairsearch);

#pragma omp parallel for num_threads(nth) schedule(static)
    for (th = 0; th < nth; th++)
    {
        try
        {
            for (int g = g0; g < g1; g++)
            {
                const Nbnxm::Grid  &grid       = nbs->grid[g];
                const int           numCellsXY = grid.numColumns();

                const int           cxy0 = (numCellsXY* th      + nth - 1)/nth;
                const int           cxy1 = (numCellsXY*(th + 1) + nth - 1)/nth;

                for (int cxy = cxy0; cxy < cxy1; cxy++)
                {
                    const int na  = grid.numAtomsInColumn(cxy);
                    const int ash = grid.firstAtomInColumn(cxy);

                    int       na_fill;
                    if (g == 0 && FillLocal)
                    {
                        na_fill = grid.paddedNumAtomsInColumn(cxy);
                    }
                    else
                    {
                        /* We fill only the real particle locations.
                         * We assume the filling entries at the end have been
                         * properly set before during pair-list generation.
                         */
                        na_fill = na;
                    }
                    copy_rvec_to_nbat_real(nbs->a.data() + ash, na, na_fill, x,
                                           nbat->XFormat, nbat->x().data(), ash);
                }
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }

    wallcycle_sub_stop(wcycle, ewcsNB_X_BUF_OPS);
    wallcycle_stop(wcycle, ewcNB_XF_BUF_OPS);
}

static void
nbnxn_atomdata_clear_reals(gmx::ArrayRef<real> dest,
                           int i0, int i1)
{
    for (int i = i0; i < i1; i++)
    {
        dest[i] = 0;
    }
}

gmx_unused static void
nbnxn_atomdata_reduce_reals(real * gmx_restrict dest,
                            gmx_bool bDestSet,
                            const real ** gmx_restrict src,
                            int nsrc,
                            int i0, int i1)
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
        /* The destination buffer is unitialized, set it first */
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

gmx_unused static void
nbnxn_atomdata_reduce_reals_simd(real gmx_unused * gmx_restrict dest,
                                 gmx_bool gmx_unused bDestSet,
                                 const gmx_unused real ** gmx_restrict src,
                                 int gmx_unused nsrc,
                                 int gmx_unused i0, int gmx_unused i1)
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
            dest_SSE = load<SimdReal>(dest+i);
            for (int s = 0; s < nsrc; s++)
            {
                src_SSE  = load<SimdReal>(src[s]+i);
                dest_SSE = dest_SSE + src_SSE;
            }
            store(dest+i, dest_SSE);
        }
    }
    else
    {
        for (int i = i0; i < i1; i += GMX_SIMD_REAL_WIDTH)
        {
            dest_SSE = load<SimdReal>(src[0]+i);
            for (int s = 1; s < nsrc; s++)
            {
                src_SSE  = load<SimdReal>(src[s]+i);
                dest_SSE = dest_SSE + src_SSE;
            }
            store(dest+i, dest_SSE);
        }
    }
#endif
}

/* Add part of the force array(s) from nbnxn_atomdata_t to f */
static void
nbnxn_atomdata_add_nbat_f_to_f_part(const nbnxn_search *nbs,
                                    const nbnxn_atomdata_t *nbat,
                                    gmx::ArrayRef<const nbnxn_atomdata_output_t> out,
                                    int nfa,
                                    int a0, int a1,
                                    rvec *f)
{
    gmx::ArrayRef<const int> cell = nbs->cell;

    /* Loop over all columns and copy and fill */
    switch (nbat->FFormat)
    {
        case nbatXYZ:
        case nbatXYZQ:
            if (nfa == 1)
            {
                const real *fnb = out[0].f.data();

                for (int a = a0; a < a1; a++)
                {
                    int i = cell[a]*nbat->fstride;

                    f[a][XX] += fnb[i];
                    f[a][YY] += fnb[i+1];
                    f[a][ZZ] += fnb[i+2];
                }
            }
            else
            {
                for (int a = a0; a < a1; a++)
                {
                    int i = cell[a]*nbat->fstride;

                    for (int fa = 0; fa < nfa; fa++)
                    {
                        f[a][XX] += out[fa].f[i];
                        f[a][YY] += out[fa].f[i+1];
                        f[a][ZZ] += out[fa].f[i+2];
                    }
                }
            }
            break;
        case nbatX4:
            if (nfa == 1)
            {
                const real *fnb = out[0].f.data();

                for (int a = a0; a < a1; a++)
                {
                    int i = atom_to_x_index<c_packX4>(cell[a]);

                    f[a][XX] += fnb[i+XX*c_packX4];
                    f[a][YY] += fnb[i+YY*c_packX4];
                    f[a][ZZ] += fnb[i+ZZ*c_packX4];
                }
            }
            else
            {
                for (int a = a0; a < a1; a++)
                {
                    int i = atom_to_x_index<c_packX4>(cell[a]);

                    for (int fa = 0; fa < nfa; fa++)
                    {
                        f[a][XX] += out[fa].f[i+XX*c_packX4];
                        f[a][YY] += out[fa].f[i+YY*c_packX4];
                        f[a][ZZ] += out[fa].f[i+ZZ*c_packX4];
                    }
                }
            }
            break;
        case nbatX8:
            if (nfa == 1)
            {
                const real *fnb = out[0].f.data();

                for (int a = a0; a < a1; a++)
                {
                    int i = atom_to_x_index<c_packX8>(cell[a]);

                    f[a][XX] += fnb[i+XX*c_packX8];
                    f[a][YY] += fnb[i+YY*c_packX8];
                    f[a][ZZ] += fnb[i+ZZ*c_packX8];
                }
            }
            else
            {
                for (int a = a0; a < a1; a++)
                {
                    int i = atom_to_x_index<c_packX8>(cell[a]);

                    for (int fa = 0; fa < nfa; fa++)
                    {
                        f[a][XX] += out[fa].f[i+XX*c_packX8];
                        f[a][YY] += out[fa].f[i+YY*c_packX8];
                        f[a][ZZ] += out[fa].f[i+ZZ*c_packX8];
                    }
                }
            }
            break;
        default:
            gmx_incons("Unsupported nbnxn_atomdata_t format");
    }
}

static inline unsigned char reverse_bits(unsigned char b)
{
    /* http://graphics.stanford.edu/~seander/bithacks.html#ReverseByteWith64BitsDiv */
    return (b * 0x0202020202ULL & 0x010884422010ULL) % 1023;
}

static void nbnxn_atomdata_add_nbat_f_to_f_treereduce(nbnxn_atomdata_t *nbat,
                                                      int               nth)
{
    const nbnxn_buffer_flags_t *flags = &nbat->buffer_flags;

    int                         next_pow2 = 1<<(gmx::log2I(nth-1)+1);

    const int                   numOutputBuffers = nbat->out.size();
    GMX_ASSERT(numOutputBuffers == nth,
               "tree-reduce currently only works for numOutputBuffers==nth");

    memset(nbat->syncStep, 0, sizeof(*(nbat->syncStep))*nth);

#pragma omp parallel num_threads(nth)
    {
        try
        {
            int   b0, b1, b;
            int   i0, i1;
            int   group_size, th;

            th = gmx_omp_get_thread_num();

            for (group_size = 2; group_size < 2*next_pow2; group_size *= 2)
            {
                int index[2], group_pos, partner_pos, wu;
                int partner_th = th ^ (group_size/2);

                if (group_size > 2)
                {
#ifdef TMPI_ATOMICS
                    /* wait on partner thread - replaces full barrier */
                    int sync_th, sync_group_size;

                    tMPI_Atomic_memory_barrier();                         /* gurantee data is saved before marking work as done */
                    tMPI_Atomic_set(&(nbat->syncStep[th]), group_size/2); /* mark previous step as completed */

                    /* find thread to sync with. Equal to partner_th unless nth is not a power of two. */
                    for (sync_th = partner_th, sync_group_size = group_size; sync_th >= nth && sync_group_size > 2; sync_group_size /= 2)
                    {
                        sync_th &= ~(sync_group_size/4);
                    }
                    if (sync_th < nth) /* otherwise nothing to sync index[1] will be >=nout */
                    {
                        /* wait on the thread which computed input data in previous step */
                        while (tMPI_Atomic_get(static_cast<volatile tMPI_Atomic_t*>(&(nbat->syncStep[sync_th]))) < group_size/2)
                        {
                            gmx_pause();
                        }
                        /* guarantee that no later load happens before wait loop is finisehd */
                        tMPI_Atomic_memory_barrier();
                    }
#else               /* TMPI_ATOMICS */
#pragma omp barrier
#endif
                }

                /* Calculate buffers to sum (result goes into first buffer) */
                group_pos = th % group_size;
                index[0]  = th - group_pos;
                index[1]  = index[0] + group_size/2;

                /* If no second buffer, nothing to do */
                if (index[1] >= numOutputBuffers && group_size > 2)
                {
                    continue;
                }

#if NBNXN_BUFFERFLAG_MAX_THREADS > 256
#error reverse_bits assumes max 256 threads
#endif
                /* Position is permuted so that one of the 2 vectors being added was computed on the same thread in the previous step.
                   This improves locality and enables to sync with just a single thread between steps (=the levels in the btree).
                   The permutation which allows this corresponds to reversing the bits of the group position.
                 */
                group_pos = reverse_bits(group_pos)/(256/group_size);

                partner_pos = group_pos ^ 1;

                /* loop over two work-units (own and partner) */
                for (wu = 0; wu < 2; wu++)
                {
                    if (wu == 1)
                    {
                        if (partner_th < nth)
                        {
                            break; /* partner exists we don't have to do his work */
                        }
                        else
                        {
                            group_pos = partner_pos;
                        }
                    }

                    /* Calculate the cell-block range for our thread */
                    b0 = (flags->nflag* group_pos   )/group_size;
                    b1 = (flags->nflag*(group_pos+1))/group_size;

                    for (b = b0; b < b1; b++)
                    {
                        i0 =  b   *NBNXN_BUFFERFLAG_SIZE*nbat->fstride;
                        i1 = (b+1)*NBNXN_BUFFERFLAG_SIZE*nbat->fstride;

                        if (bitmask_is_set(flags->flag[b], index[1]) || group_size > 2)
                        {
                            const real *fIndex1 = nbat->out[index[1]].f.data();
#if GMX_SIMD
                            nbnxn_atomdata_reduce_reals_simd
#else
                            nbnxn_atomdata_reduce_reals
#endif
                                (nbat->out[index[0]].f.data(),
                                bitmask_is_set(flags->flag[b], index[0]) || group_size > 2,
                                &fIndex1, 1, i0, i1);

                        }
                        else if (!bitmask_is_set(flags->flag[b], index[0]))
                        {
                            nbnxn_atomdata_clear_reals(nbat->out[index[0]].f,
                                                       i0, i1);
                        }
                    }
                }
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }
}


static void nbnxn_atomdata_add_nbat_f_to_f_stdreduce(nbnxn_atomdata_t *nbat,
                                                     int               nth)
{
#pragma omp parallel for num_threads(nth) schedule(static)
    for (int th = 0; th < nth; th++)
    {
        try
        {
            const nbnxn_buffer_flags_t *flags;
            int                         nfptr;
            const real                 *fptr[NBNXN_BUFFERFLAG_MAX_THREADS];

            flags = &nbat->buffer_flags;

            /* Calculate the cell-block range for our thread */
            int b0 = (flags->nflag* th   )/nth;
            int b1 = (flags->nflag*(th+1))/nth;

            for (int b = b0; b < b1; b++)
            {
                int i0 =  b   *NBNXN_BUFFERFLAG_SIZE*nbat->fstride;
                int i1 = (b+1)*NBNXN_BUFFERFLAG_SIZE*nbat->fstride;

                nfptr = 0;
                for (int out = 1; out < gmx::ssize(nbat->out); out++)
                {
                    if (bitmask_is_set(flags->flag[b], out))
                    {
                        fptr[nfptr++] = nbat->out[out].f.data();
                    }
                }
                if (nfptr > 0)
                {
#if GMX_SIMD
                    nbnxn_atomdata_reduce_reals_simd
#else
                    nbnxn_atomdata_reduce_reals
#endif
                        (nbat->out[0].f.data(),
                        bitmask_is_set(flags->flag[b], 0),
                        fptr, nfptr,
                        i0, i1);
                }
                else if (!bitmask_is_set(flags->flag[b], 0))
                {
                    nbnxn_atomdata_clear_reals(nbat->out[0].f,
                                               i0, i1);
                }
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }
}

/* Add the force array(s) from nbnxn_atomdata_t to f */
void
nonbonded_verlet_t::atomdata_add_nbat_f_to_f(const Nbnxm::AtomLocality  locality,
                                             rvec                      *f,
                                             gmx_wallcycle             *wcycle)
{
    /* Skip the non-local reduction if there was no non-local work to do */
    if (locality == Nbnxm::AtomLocality::NonLocal &&
        pairlistSets().pairlistSet(Nbnxm::InteractionLocality::NonLocal).nblGpu[0]->sci.empty())
    {
        return;
    }

    wallcycle_start(wcycle, ewcNB_XF_BUF_OPS);
    wallcycle_sub_start(wcycle, ewcsNB_F_BUF_OPS);

    int a0 = 0, na = 0;

    nbs_cycle_start(&nbs->cc[enbsCCreducef]);

    switch (locality)
    {
        case Nbnxm::AtomLocality::All:
        case Nbnxm::AtomLocality::Count:
            a0 = 0;
            na = nbs->natoms_nonlocal;
            break;
        case Nbnxm::AtomLocality::Local:
            a0 = 0;
            na = nbs->natoms_local;
            break;
        case Nbnxm::AtomLocality::NonLocal:
            a0 = nbs->natoms_local;
            na = nbs->natoms_nonlocal - nbs->natoms_local;
            break;
    }

    int nth = gmx_omp_nthreads_get(emntNonbonded);

    if (nbat->out.size() > 1)
    {
        if (locality != Nbnxm::AtomLocality::All)
        {
            gmx_incons("add_f_to_f called with nout>1 and locality!=eatAll");
        }

        /* Reduce the force thread output buffers into buffer 0, before adding
         * them to the, differently ordered, "real" force buffer.
         */
        if (nbat->bUseTreeReduce)
        {
            nbnxn_atomdata_add_nbat_f_to_f_treereduce(nbat.get(), nth);
        }
        else
        {
            nbnxn_atomdata_add_nbat_f_to_f_stdreduce(nbat.get(), nth);
        }
    }
#pragma omp parallel for num_threads(nth) schedule(static)
    for (int th = 0; th < nth; th++)
    {
        try
        {
            nbnxn_atomdata_add_nbat_f_to_f_part(nbs.get(), nbat.get(),
                                                nbat->out,
                                                1,
                                                a0+((th+0)*na)/nth,
                                                a0+((th+1)*na)/nth,
                                                f);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }

    nbs_cycle_stop(&nbs->cc[enbsCCreducef]);

    wallcycle_sub_stop(wcycle, ewcsNB_F_BUF_OPS);
    wallcycle_stop(wcycle, ewcNB_XF_BUF_OPS);
}

/* Adds the shift forces from nbnxn_atomdata_t to fshift */
void nbnxn_atomdata_add_nbat_fshift_to_fshift(const nbnxn_atomdata_t *nbat,
                                              rvec                   *fshift)
{
    gmx::ArrayRef<const nbnxn_atomdata_output_t> outputBuffers = nbat->out;

    if (outputBuffers.size() == 1)
    {
        /* When there is a single output object, with CPU or GPU, shift forces
         * have been written directly to the main buffer instead of to the
         * (single) thread local output object. There is nothing to reduce.
         */
        return;
    }

    for (int s = 0; s < SHIFTS; s++)
    {
        rvec sum;
        clear_rvec(sum);
        for (const nbnxn_atomdata_output_t &out : outputBuffers)
        {
            sum[XX] += out.fshift[s*DIM+XX];
            sum[YY] += out.fshift[s*DIM+YY];
            sum[ZZ] += out.fshift[s*DIM+ZZ];
        }
        rvec_inc(fshift[s], sum);
    }
}
