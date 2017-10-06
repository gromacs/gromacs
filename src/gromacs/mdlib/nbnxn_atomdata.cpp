/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include "nbnxn_atomdata.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <cmath>

#include <algorithm>

#include "thread_mpi/atomic.h"

#include "gromacs/math/functions.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/nbnxn_consts.h"
#include "gromacs/mdlib/nbnxn_internal.h"
#include "gromacs/mdlib/nbnxn_search.h"
#include "gromacs/mdlib/nbnxn_util.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/simd/simd.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/smalloc.h"

using namespace gmx; // TODO: Remove when this file is moved into gmx namespace

/* Default nbnxn allocation routine, allocates NBNXN_MEM_ALIGN byte aligned */
void nbnxn_alloc_aligned(void **ptr, size_t nbytes)
{
    *ptr = save_malloc_aligned("ptr", __FILE__, __LINE__, nbytes, 1, NBNXN_MEM_ALIGN);
}

/* Free function for memory allocated with nbnxn_alloc_aligned */
void nbnxn_free_aligned(void *ptr)
{
    sfree_aligned(ptr);
}

/* Reallocation wrapper function for nbnxn data structures */
void nbnxn_realloc_void(void **ptr,
                        int nbytes_copy, int nbytes_new,
                        nbnxn_alloc_t *ma,
                        nbnxn_free_t  *mf)
{
    void *ptr_new;

    ma(&ptr_new, nbytes_new);

    if (nbytes_new > 0 && ptr_new == nullptr)
    {
        gmx_fatal(FARGS, "Allocation of %d bytes failed", nbytes_new);
    }

    if (nbytes_copy > 0)
    {
        if (nbytes_new < nbytes_copy)
        {
            gmx_incons("In nbnxn_realloc_void: new size less than copy size");
        }
        memcpy(ptr_new, *ptr, nbytes_copy);
    }
    if (*ptr != nullptr)
    {
        mf(*ptr);
    }
    *ptr = ptr_new;
}

/* Reallocate the nbnxn_atomdata_t for a size of n atoms */
void nbnxn_atomdata_realloc(nbnxn_atomdata_t *nbat, int n)
{
    int t;

    nbnxn_realloc_void((void **)&nbat->type,
                       nbat->natoms*sizeof(*nbat->type),
                       n*sizeof(*nbat->type),
                       nbat->alloc, nbat->free);
    nbnxn_realloc_void((void **)&nbat->lj_comb,
                       nbat->natoms*2*sizeof(*nbat->lj_comb),
                       n*2*sizeof(*nbat->lj_comb),
                       nbat->alloc, nbat->free);
    if (nbat->XFormat != nbatXYZQ)
    {
        nbnxn_realloc_void((void **)&nbat->q,
                           nbat->natoms*sizeof(*nbat->q),
                           n*sizeof(*nbat->q),
                           nbat->alloc, nbat->free);
    }
    if (nbat->nenergrp > 1)
    {
        nbnxn_realloc_void((void **)&nbat->energrp,
                           nbat->natoms/nbat->na_c*sizeof(*nbat->energrp),
                           n/nbat->na_c*sizeof(*nbat->energrp),
                           nbat->alloc, nbat->free);
    }
    nbnxn_realloc_void((void **)&nbat->x,
                       nbat->natoms*nbat->xstride*sizeof(*nbat->x),
                       n*nbat->xstride*sizeof(*nbat->x),
                       nbat->alloc, nbat->free);
    for (t = 0; t < nbat->nout; t++)
    {
        /* Allocate one element extra for possible signaling with GPUs */
        nbnxn_realloc_void((void **)&nbat->out[t].f,
                           nbat->natoms*nbat->fstride*sizeof(*nbat->out[t].f),
                           n*nbat->fstride*sizeof(*nbat->out[t].f),
                           nbat->alloc, nbat->free);
    }
    nbat->nalloc = n;
}

/* Initializes an nbnxn_atomdata_output_t data structure */
static void nbnxn_atomdata_output_init(nbnxn_atomdata_output_t *out,
                                       int nb_kernel_type,
                                       int nenergrp, int stride,
                                       nbnxn_alloc_t *ma)
{
    out->f = nullptr;
    ma((void **)&out->fshift, SHIFTS*DIM*sizeof(*out->fshift));
    out->nV = nenergrp*nenergrp;
    ma((void **)&out->Vvdw, out->nV*sizeof(*out->Vvdw));
    ma((void **)&out->Vc, out->nV*sizeof(*out->Vc  ));

    if (nb_kernel_type == nbnxnk4xN_SIMD_4xN ||
        nb_kernel_type == nbnxnk4xN_SIMD_2xNN)
    {
        int cj_size  = nbnxn_kernel_to_cluster_j_size(nb_kernel_type);
        out->nVS = nenergrp*nenergrp*stride*(cj_size>>1)*cj_size;
        ma((void **)&out->VSvdw, out->nVS*sizeof(*out->VSvdw));
        ma((void **)&out->VSc, out->nVS*sizeof(*out->VSc  ));
    }
    else
    {
        out->nVS = 0;
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
static void set_lj_parameter_data(nbnxn_atomdata_t *nbat, gmx_bool bSIMD)
{
    real c6, c12;

    int  nt = nbat->ntype;

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
        nbat->alloc((void **)&nbat->nbfp_aligned,
                    nt*nt*c_simdBestPairAlignment*sizeof(*nbat->nbfp_aligned));
        for (int i = 0; i < nt; i++)
        {
            for (int j = 0; j < nt; j++)
            {
                nbat->nbfp_aligned[(i*nt+j)*c_simdBestPairAlignment+0] = nbat->nbfp[(i*nt+j)*2+0];
                nbat->nbfp_aligned[(i*nt+j)*c_simdBestPairAlignment+1] = nbat->nbfp[(i*nt+j)*2+1];
                nbat->nbfp_aligned[(i*nt+j)*c_simdBestPairAlignment+2] = 0;
                nbat->nbfp_aligned[(i*nt+j)*c_simdBestPairAlignment+3] = 0;
            }
        }
#endif
    }

    /* We use combination rule data for SIMD combination rule kernels
     * and with LJ-PME kernels. We then only need parameters per atom type,
     * not per pair of atom types.
     */
    switch (nbat->comb_rule)
    {
        case ljcrGEOM:
            nbat->comb_rule = ljcrGEOM;

            for (int i = 0; i < nt; i++)
            {
                /* Store the sqrt of the diagonal from the nbfp matrix */
                nbat->nbfp_comb[i*2  ] = std::sqrt(nbat->nbfp[(i*nt+i)*2  ]);
                nbat->nbfp_comb[i*2+1] = std::sqrt(nbat->nbfp[(i*nt+i)*2+1]);
            }
            break;
        case ljcrLB:
            for (int i = 0; i < nt; i++)
            {
                /* Get 6*C6 and 12*C12 from the diagonal of the nbfp matrix */
                c6  = nbat->nbfp[(i*nt+i)*2  ];
                c12 = nbat->nbfp[(i*nt+i)*2+1];
                if (c6 > 0 && c12 > 0)
                {
                    /* We store 0.5*2^1/6*sigma and sqrt(4*3*eps),
                     * so we get 6*C6 and 12*C12 after combining.
                     */
                    nbat->nbfp_comb[i*2  ] = 0.5*gmx::sixthroot(c12/c6);
                    nbat->nbfp_comb[i*2+1] = std::sqrt(c6*c6/c12);
                }
                else
                {
                    nbat->nbfp_comb[i*2  ] = 0;
                    nbat->nbfp_comb[i*2+1] = 0;
                }
            }
            break;
        case ljcrNONE:
            /* We always store the full matrix (see code above) */
            break;
        default:
            gmx_incons("Unknown combination rule");
            break;
    }
}

#if GMX_SIMD
static void
nbnxn_atomdata_init_simple_exclusion_masks(nbnxn_atomdata_t *nbat)
{
    const int simd_width = GMX_SIMD_REAL_WIDTH;
    int       simd_excl_size;
    /* Set the diagonal cluster pair exclusion mask setup data.
     * In the kernel we check 0 < j - i to generate the masks.
     * Here we store j - i for generating the mask for the first i,
     * we substract 0.5 to avoid rounding issues.
     * In the kernel we can subtract 1 to generate the subsequent mask.
     */
    int        simd_4xn_diag_size;

    simd_4xn_diag_size = std::max(NBNXN_CPU_CLUSTER_I_SIZE, simd_width);
    snew_aligned(nbat->simd_4xn_diagonal_j_minus_i, simd_4xn_diag_size, NBNXN_MEM_ALIGN);
    for (int j = 0; j < simd_4xn_diag_size; j++)
    {
        nbat->simd_4xn_diagonal_j_minus_i[j] = j - 0.5;
    }

    snew_aligned(nbat->simd_2xnn_diagonal_j_minus_i, simd_width, NBNXN_MEM_ALIGN);
    for (int j = 0; j < simd_width/2; j++)
    {
        /* The j-cluster size is half the SIMD width */
        nbat->simd_2xnn_diagonal_j_minus_i[j]              = j - 0.5;
        /* The next half of the SIMD width is for i + 1 */
        nbat->simd_2xnn_diagonal_j_minus_i[simd_width/2+j] = j - 1 - 0.5;
    }

    /* We use up to 32 bits for exclusion masking.
     * The same masks are used for the 4xN and 2x(N+N) kernels.
     * The masks are read either into integer SIMD registers or into
     * real SIMD registers (together with a cast).
     * In single precision this means the real and integer SIMD registers
     * are of equal size.
     */
    simd_excl_size = NBNXN_CPU_CLUSTER_I_SIZE*simd_width;
#if GMX_DOUBLE && !GMX_SIMD_HAVE_INT32_LOGICAL
    snew_aligned(nbat->simd_exclusion_filter64, simd_excl_size,   NBNXN_MEM_ALIGN);
#else
    snew_aligned(nbat->simd_exclusion_filter, simd_excl_size,   NBNXN_MEM_ALIGN);
#endif

    for (int j = 0; j < simd_excl_size; j++)
    {
        /* Set the consecutive bits for masking pair exclusions */
#if GMX_DOUBLE && !GMX_SIMD_HAVE_INT32_LOGICAL
        nbat->simd_exclusion_filter64[j]     = (1U << j);
#else
        nbat->simd_exclusion_filter[j]       = (1U << j);
#endif
    }

#if !GMX_SIMD_HAVE_LOGICAL && !GMX_SIMD_HAVE_INT32_LOGICAL
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

    simd_excl_size = NBNXN_CPU_CLUSTER_I_SIZE*NBNXN_CPU_CLUSTER_I_SIZE;
    const real simdFalse =  0.0;
    const real simdTrue  =  1.0;
    real      *simd_interaction_array;

    snew_aligned(simd_interaction_array, simd_excl_size * GMX_SIMD_REAL_WIDTH, NBNXN_MEM_ALIGN);
    for (int j = 0; j < simd_excl_size; j++)
    {
        int index = j * GMX_SIMD_REAL_WIDTH;
        for (int i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
        {
            simd_interaction_array[index + i] = (j & (1 << i)) ? simdTrue : simdFalse;
        }
    }
    nbat->simd_interaction_array = simd_interaction_array;
#endif
}
#endif

/* Initializes an nbnxn_atomdata_t data structure */
void nbnxn_atomdata_init(FILE *fp,
                         nbnxn_atomdata_t *nbat,
                         int nb_kernel_type,
                         int enbnxninitcombrule,
                         int ntype, const real *nbfp,
                         int n_energygroups,
                         int nout,
                         nbnxn_alloc_t *alloc,
                         nbnxn_free_t  *free)
{
    int      nth;
    real     c6, c12, tol;
    char    *ptr;
    gmx_bool simple, bCombGeom, bCombLB, bSIMD;

    if (alloc == nullptr)
    {
        nbat->alloc = nbnxn_alloc_aligned;
    }
    else
    {
        nbat->alloc = alloc;
    }
    if (free == nullptr)
    {
        nbat->free = nbnxn_free_aligned;
    }
    else
    {
        nbat->free = free;
    }

    if (debug)
    {
        fprintf(debug, "There are %d atom types in the system, adding one for nbnxn_atomdata_t\n", ntype);
    }
    nbat->ntype = ntype + 1;
    nbat->alloc((void **)&nbat->nbfp,
                nbat->ntype*nbat->ntype*2*sizeof(*nbat->nbfp));
    nbat->alloc((void **)&nbat->nbfp_comb, nbat->ntype*2*sizeof(*nbat->nbfp_comb));

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

    /* Temporarily fill nbat->nbfp_comb with sigma and epsilon
     * to check for the LB rule.
     */
    for (int i = 0; i < ntype; i++)
    {
        c6  = nbfp[(i*ntype+i)*2  ]/6.0;
        c12 = nbfp[(i*ntype+i)*2+1]/12.0;
        if (c6 > 0 && c12 > 0)
        {
            nbat->nbfp_comb[i*2  ] = gmx::sixthroot(c12/c6);
            nbat->nbfp_comb[i*2+1] = 0.25*c6*c6/c12;
        }
        else if (c6 == 0 && c12 == 0)
        {
            nbat->nbfp_comb[i*2  ] = 0;
            nbat->nbfp_comb[i*2+1] = 0;
        }
        else
        {
            /* Can not use LB rule with only dispersion or repulsion */
            bCombLB = FALSE;
        }
    }

    for (int i = 0; i < nbat->ntype; i++)
    {
        for (int j = 0; j < nbat->ntype; j++)
        {
            if (i < ntype && j < ntype)
            {
                /* fr->nbfp has been updated, so that array too now stores c6/c12 including
                 * the 6.0/12.0 prefactors to save 2 flops in the most common case (force-only).
                 */
                c6  = nbfp[(i*ntype+j)*2  ];
                c12 = nbfp[(i*ntype+j)*2+1];
                nbat->nbfp[(i*nbat->ntype+j)*2  ] = c6;
                nbat->nbfp[(i*nbat->ntype+j)*2+1] = c12;

                /* Compare 6*C6 and 12*C12 for geometric cobination rule */
                bCombGeom = bCombGeom &&
                    gmx_within_tol(c6*c6, nbfp[(i*ntype+i)*2  ]*nbfp[(j*ntype+j)*2  ], tol) &&
                    gmx_within_tol(c12*c12, nbfp[(i*ntype+i)*2+1]*nbfp[(j*ntype+j)*2+1], tol);

                /* Compare C6 and C12 for Lorentz-Berthelot combination rule */
                c6     /= 6.0;
                c12    /= 12.0;
                bCombLB = bCombLB &&
                    ((c6 == 0 && c12 == 0 &&
                      (nbat->nbfp_comb[i*2+1] == 0 || nbat->nbfp_comb[j*2+1] == 0)) ||
                     (c6 > 0 && c12 > 0 &&
                      gmx_within_tol(gmx::sixthroot(c12/c6),
                                     0.5*(nbat->nbfp_comb[i*2]+nbat->nbfp_comb[j*2]), tol) &&
                      gmx_within_tol(0.25*c6*c6/c12, std::sqrt(nbat->nbfp_comb[i*2+1]*nbat->nbfp_comb[j*2+1]), tol)));
            }
            else
            {
                /* Add zero parameters for the additional dummy atom type */
                nbat->nbfp[(i*nbat->ntype+j)*2  ] = 0;
                nbat->nbfp[(i*nbat->ntype+j)*2+1] = 0;
            }
        }
    }
    if (debug)
    {
        fprintf(debug, "Combination rules: geometric %d Lorentz-Berthelot %d\n",
                bCombGeom, bCombLB);
    }

    simple = nbnxn_kernel_pairlist_simple(nb_kernel_type);

    switch (enbnxninitcombrule)
    {
        case enbnxninitcombruleDETECT:
            /* We prefer the geometic combination rule,
             * as that gives a slightly faster kernel than the LB rule.
             */
            if (bCombGeom)
            {
                nbat->comb_rule = ljcrGEOM;
            }
            else if (bCombLB)
            {
                nbat->comb_rule = ljcrLB;
            }
            else
            {
                nbat->comb_rule = ljcrNONE;

                nbat->free(nbat->nbfp_comb);
            }

            if (fp)
            {
                if (nbat->comb_rule == ljcrNONE)
                {
                    fprintf(fp, "Using full Lennard-Jones parameter combination matrix\n\n");
                }
                else
                {
                    fprintf(fp, "Using %s Lennard-Jones combination rule\n\n",
                            nbat->comb_rule == ljcrGEOM ? "geometric" : "Lorentz-Berthelot");
                }
            }
            break;
        case enbnxninitcombruleGEOM:
            nbat->comb_rule = ljcrGEOM;
            break;
        case enbnxninitcombruleLB:
            nbat->comb_rule = ljcrLB;
            break;
        case enbnxninitcombruleNONE:
            nbat->comb_rule = ljcrNONE;

            nbat->free(nbat->nbfp_comb);
            break;
        default:
            gmx_incons("Unknown enbnxninitcombrule");
    }

    bSIMD = (nb_kernel_type == nbnxnk4xN_SIMD_4xN ||
             nb_kernel_type == nbnxnk4xN_SIMD_2xNN);

    set_lj_parameter_data(nbat, bSIMD);

    nbat->natoms  = 0;
    nbat->type    = nullptr;
    nbat->lj_comb = nullptr;
    if (simple)
    {
        int pack_x;

        if (bSIMD)
        {
            pack_x = std::max(NBNXN_CPU_CLUSTER_I_SIZE,
                              nbnxn_kernel_to_cluster_j_size(nb_kernel_type));
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
    nbat->q        = nullptr;
    nbat->nenergrp = n_energygroups;
    if (!simple)
    {
        /* Energy groups not supported yet for super-sub lists */
        if (n_energygroups > 1 && fp != nullptr)
        {
            fprintf(fp, "\nNOTE: With GPUs, reporting energy group contributions is not supported\n\n");
        }
        nbat->nenergrp = 1;
    }
    /* Temporary storage goes as #grp^3*simd_width^2/2, so limit to 64 */
    if (nbat->nenergrp > 64)
    {
        gmx_fatal(FARGS, "With NxN kernels not more than 64 energy groups are supported\n");
    }
    nbat->neg_2log = 1;
    while (nbat->nenergrp > (1<<nbat->neg_2log))
    {
        nbat->neg_2log++;
    }
    nbat->energrp = nullptr;
    nbat->alloc((void **)&nbat->shift_vec, SHIFTS*sizeof(*nbat->shift_vec));
    nbat->xstride = (nbat->XFormat == nbatXYZQ ? STRIDE_XYZQ : DIM);
    nbat->fstride = (nbat->FFormat == nbatXYZQ ? STRIDE_XYZQ : DIM);
    nbat->x       = nullptr;

#if GMX_SIMD
    if (simple)
    {
        nbnxn_atomdata_init_simple_exclusion_masks(nbat);
    }
#endif

    /* Initialize the output data structures */
    nbat->nout    = nout;
    snew(nbat->out, nbat->nout);
    nbat->nalloc  = 0;
    for (int i = 0; i < nbat->nout; i++)
    {
        nbnxn_atomdata_output_init(&nbat->out[i],
                                   nb_kernel_type,
                                   nbat->nenergrp, 1<<nbat->neg_2log,
                                   nbat->alloc);
    }
    nbat->buffer_flags.flag        = nullptr;
    nbat->buffer_flags.flag_nalloc = 0;

    nth = gmx_omp_nthreads_get(emntNonbonded);

    ptr = getenv("GMX_USE_TREEREDUCE");
    if (ptr != nullptr)
    {
        nbat->bUseTreeReduce = strtol(ptr, nullptr, 10);
    }
#if defined __MIC__
    else if (nth > 8) /*on the CPU we currently don't benefit even at 32*/
    {
        nbat->bUseTreeReduce = 1;
    }
#endif
    else
    {
        nbat->bUseTreeReduce = 0;
    }
    if (nbat->bUseTreeReduce)
    {
        if (fp)
        {
            fprintf(fp, "Using tree force reduction\n\n");
        }
        snew(nbat->syncStep, nth);
    }
}

template<int packSize>
static void copy_lj_to_nbat_lj_comb(const real *ljparam_type,
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

/* Sets the atom type in nbnxn_atomdata_t */
static void nbnxn_atomdata_set_atomtypes(nbnxn_atomdata_t    *nbat,
                                         const nbnxn_search_t nbs,
                                         const int           *type)
{
    for (int g = 0; g < nbs->ngrid; g++)
    {
        const nbnxn_grid_t * grid = &nbs->grid[g];

        /* Loop over all columns and copy and fill */
        for (int i = 0; i < grid->ncx*grid->ncy; i++)
        {
            int ncz = grid->cxy_ind[i+1] - grid->cxy_ind[i];
            int ash = (grid->cell0 + grid->cxy_ind[i])*grid->na_sc;

            copy_int_to_nbat_int(nbs->a+ash, grid->cxy_na[i], ncz*grid->na_sc,
                                 type, nbat->ntype-1, nbat->type+ash);
        }
    }
}

/* Sets the LJ combination rule parameters in nbnxn_atomdata_t */
static void nbnxn_atomdata_set_ljcombparams(nbnxn_atomdata_t    *nbat,
                                            const nbnxn_search_t nbs)
{
    if (nbat->comb_rule != ljcrNONE)
    {
        for (int g = 0; g < nbs->ngrid; g++)
        {
            const nbnxn_grid_t * grid = &nbs->grid[g];

            /* Loop over all columns and copy and fill */
            for (int i = 0; i < grid->ncx*grid->ncy; i++)
            {
                int ncz = grid->cxy_ind[i+1] - grid->cxy_ind[i];
                int ash = (grid->cell0 + grid->cxy_ind[i])*grid->na_sc;

                if (nbat->XFormat == nbatX4)
                {
                    copy_lj_to_nbat_lj_comb<c_packX4>(nbat->nbfp_comb,
                                                      nbat->type + ash,
                                                      ncz*grid->na_sc,
                                                      nbat->lj_comb + ash*2);
                }
                else if (nbat->XFormat == nbatX8)
                {
                    copy_lj_to_nbat_lj_comb<c_packX8>(nbat->nbfp_comb,
                                                      nbat->type + ash,
                                                      ncz*grid->na_sc,
                                                      nbat->lj_comb + ash*2);
                }
                else if (nbat->XFormat == nbatXYZQ)
                {
                    copy_lj_to_nbat_lj_comb<1>(nbat->nbfp_comb,
                                               nbat->type + ash,
                                               ncz*grid->na_sc,
                                               nbat->lj_comb + ash*2);
                }
            }
        }
    }
}

/* Sets the charges in nbnxn_atomdata_t *nbat */
static void nbnxn_atomdata_set_charges(nbnxn_atomdata_t    *nbat,
                                       const nbnxn_search_t nbs,
                                       const real          *charge)
{
    int                 i;
    real               *q;

    for (int g = 0; g < nbs->ngrid; g++)
    {
        const nbnxn_grid_t * grid = &nbs->grid[g];

        /* Loop over all columns and copy and fill */
        for (int cxy = 0; cxy < grid->ncx*grid->ncy; cxy++)
        {
            int ash      = (grid->cell0 + grid->cxy_ind[cxy])*grid->na_sc;
            int na       = grid->cxy_na[cxy];
            int na_round = (grid->cxy_ind[cxy+1] - grid->cxy_ind[cxy])*grid->na_sc;

            if (nbat->XFormat == nbatXYZQ)
            {
                q = nbat->x + ash*STRIDE_XYZQ + ZZ + 1;
                for (i = 0; i < na; i++)
                {
                    *q = charge[nbs->a[ash+i]];
                    q += STRIDE_XYZQ;
                }
                /* Complete the partially filled last cell with zeros */
                for (; i < na_round; i++)
                {
                    *q = 0;
                    q += STRIDE_XYZQ;
                }
            }
            else
            {
                q = nbat->q + ash;
                for (i = 0; i < na; i++)
                {
                    *q = charge[nbs->a[ash+i]];
                    q++;
                }
                /* Complete the partially filled last cell with zeros */
                for (; i < na_round; i++)
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
                                    const nbnxn_search_t nbs)
{
    real               *q;
    int                 stride_q, nsubc;

    if (nbat->XFormat == nbatXYZQ)
    {
        q        = nbat->x + ZZ + 1;
        stride_q = STRIDE_XYZQ;
    }
    else
    {
        q        = nbat->q;
        stride_q = 1;
    }

    for (int g = 0; g < nbs->ngrid; g++)
    {
        const nbnxn_grid_t * grid = &nbs->grid[g];
        if (grid->bSimple)
        {
            nsubc = 1;
        }
        else
        {
            nsubc = c_gpuNumClusterPerCell;
        }

        int c_offset = grid->cell0*grid->na_sc;

        /* Loop over all columns and copy and fill */
        for (int c = 0; c < grid->nc*nsubc; c++)
        {
            /* Does this cluster contain perturbed particles? */
            if (grid->fep[c] != 0)
            {
                for (int i = 0; i < grid->na_c; i++)
                {
                    /* Is this a perturbed particle? */
                    if (grid->fep[c] & (1 << i))
                    {
                        int ind = c_offset + c*grid->na_c + i;
                        /* Set atom type and charge to non-interacting */
                        nbat->type[ind] = nbat->ntype - 1;
                        q[ind*stride_q] = 0;
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
static void nbnxn_atomdata_set_energygroups(nbnxn_atomdata_t    *nbat,
                                            const nbnxn_search_t nbs,
                                            const int           *atinfo)
{
    if (nbat->nenergrp == 1)
    {
        return;
    }

    for (int g = 0; g < nbs->ngrid; g++)
    {
        const nbnxn_grid_t * grid = &nbs->grid[g];

        /* Loop over all columns and copy and fill */
        for (int i = 0; i < grid->ncx*grid->ncy; i++)
        {
            int ncz = grid->cxy_ind[i+1] - grid->cxy_ind[i];
            int ash = (grid->cell0 + grid->cxy_ind[i])*grid->na_sc;

            copy_egp_to_nbat_egps(nbs->a+ash, grid->cxy_na[i], ncz*grid->na_sc,
                                  nbat->na_c, nbat->neg_2log,
                                  atinfo, nbat->energrp+(ash>>grid->na_c_2log));
        }
    }
}

/* Sets all required atom parameter data in nbnxn_atomdata_t */
void nbnxn_atomdata_set(nbnxn_atomdata_t    *nbat,
                        const nbnxn_search_t nbs,
                        const t_mdatoms     *mdatoms,
                        const int           *atinfo)
{
    nbnxn_atomdata_set_atomtypes(nbat, nbs, mdatoms->typeA);

    nbnxn_atomdata_set_charges(nbat, nbs, mdatoms->chargeA);

    if (nbs->bFEP)
    {
        nbnxn_atomdata_mask_fep(nbat, nbs);
    }

    /* This must be done after masking types for FEP */
    nbnxn_atomdata_set_ljcombparams(nbat, nbs);

    nbnxn_atomdata_set_energygroups(nbat, nbs, atinfo);
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
void nbnxn_atomdata_copy_x_to_nbat_x(const nbnxn_search_t nbs,
                                     int                  locality,
                                     gmx_bool             FillLocal,
                                     rvec                *x,
                                     nbnxn_atomdata_t    *nbat)
{
    int g0 = 0, g1 = 0;
    int nth, th;

    switch (locality)
    {
        case eatAll:
            g0 = 0;
            g1 = nbs->ngrid;
            break;
        case eatLocal:
            g0 = 0;
            g1 = 1;
            break;
        case eatNonlocal:
            g0 = 1;
            g1 = nbs->ngrid;
            break;
    }

    if (FillLocal)
    {
        nbat->natoms_local = nbs->grid[0].nc*nbs->grid[0].na_sc;
    }

    nth = gmx_omp_nthreads_get(emntPairsearch);

#pragma omp parallel for num_threads(nth) schedule(static)
    for (th = 0; th < nth; th++)
    {
        try
        {
            for (int g = g0; g < g1; g++)
            {
                const nbnxn_grid_t *grid;
                int                 cxy0, cxy1;

                grid = &nbs->grid[g];

                cxy0 = (grid->ncx*grid->ncy* th   +nth-1)/nth;
                cxy1 = (grid->ncx*grid->ncy*(th+1)+nth-1)/nth;

                for (int cxy = cxy0; cxy < cxy1; cxy++)
                {
                    int na, ash, na_fill;

                    na  = grid->cxy_na[cxy];
                    ash = (grid->cell0 + grid->cxy_ind[cxy])*grid->na_sc;

                    if (g == 0 && FillLocal)
                    {
                        na_fill =
                            (grid->cxy_ind[cxy+1] - grid->cxy_ind[cxy])*grid->na_sc;
                    }
                    else
                    {
                        /* We fill only the real particle locations.
                         * We assume the filling entries at the end have been
                         * properly set before during pair-list generation.
                         */
                        na_fill = na;
                    }
                    copy_rvec_to_nbat_real(nbs->a+ash, na, na_fill, x,
                                           nbat->XFormat, nbat->x, ash);
                }
            }
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }
}

static void
nbnxn_atomdata_clear_reals(real * gmx_restrict dest,
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
                            real ** gmx_restrict src,
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
                                 real gmx_unused ** gmx_restrict src,
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
nbnxn_atomdata_add_nbat_f_to_f_part(const nbnxn_search_t nbs,
                                    const nbnxn_atomdata_t *nbat,
                                    nbnxn_atomdata_output_t *out,
                                    int nfa,
                                    int a0, int a1,
                                    rvec *f)
{
    const int  *cell;
    const real *fnb;

    cell = nbs->cell;

    /* Loop over all columns and copy and fill */
    switch (nbat->FFormat)
    {
        case nbatXYZ:
        case nbatXYZQ:
            if (nfa == 1)
            {
                fnb = out[0].f;

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
                fnb = out[0].f;

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
                fnb = out[0].f;

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

static gmx_inline unsigned char reverse_bits(unsigned char b)
{
    /* http://graphics.stanford.edu/~seander/bithacks.html#ReverseByteWith64BitsDiv */
    return (b * 0x0202020202ULL & 0x010884422010ULL) % 1023;
}

static void nbnxn_atomdata_add_nbat_f_to_f_treereduce(const nbnxn_atomdata_t *nbat,
                                                      int                     nth)
{
    const nbnxn_buffer_flags_t *flags = &nbat->buffer_flags;

    int next_pow2 = 1<<(gmx::log2I(nth-1)+1);

    assert(nbat->nout == nth); /* tree-reduce currently only works for nout==nth */

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
                        while (tMPI_Atomic_get((volatile tMPI_Atomic_t*)&(nbat->syncStep[sync_th])) < group_size/2)
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
                if (index[1] >= nbat->nout && group_size > 2)
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
#if GMX_SIMD
                            nbnxn_atomdata_reduce_reals_simd
#else
                            nbnxn_atomdata_reduce_reals
#endif
                                (nbat->out[index[0]].f,
                                bitmask_is_set(flags->flag[b], index[0]) || group_size > 2,
                                &(nbat->out[index[1]].f), 1, i0, i1);

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


static void nbnxn_atomdata_add_nbat_f_to_f_stdreduce(const nbnxn_atomdata_t *nbat,
                                                     int                     nth)
{
#pragma omp parallel for num_threads(nth) schedule(static)
    for (int th = 0; th < nth; th++)
    {
        try
        {
            const nbnxn_buffer_flags_t *flags;
            int   nfptr;
            real *fptr[NBNXN_BUFFERFLAG_MAX_THREADS];

            flags = &nbat->buffer_flags;

            /* Calculate the cell-block range for our thread */
            int b0 = (flags->nflag* th   )/nth;
            int b1 = (flags->nflag*(th+1))/nth;

            for (int b = b0; b < b1; b++)
            {
                int i0 =  b   *NBNXN_BUFFERFLAG_SIZE*nbat->fstride;
                int i1 = (b+1)*NBNXN_BUFFERFLAG_SIZE*nbat->fstride;

                nfptr = 0;
                for (int out = 1; out < nbat->nout; out++)
                {
                    if (bitmask_is_set(flags->flag[b], out))
                    {
                        fptr[nfptr++] = nbat->out[out].f;
                    }
                }
                if (nfptr > 0)
                {
#if GMX_SIMD
                    nbnxn_atomdata_reduce_reals_simd
#else
                    nbnxn_atomdata_reduce_reals
#endif
                        (nbat->out[0].f,
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
void nbnxn_atomdata_add_nbat_f_to_f(const nbnxn_search_t    nbs,
                                    int                     locality,
                                    const nbnxn_atomdata_t *nbat,
                                    rvec                   *f)
{
    int a0 = 0, na = 0;

    nbs_cycle_start(&nbs->cc[enbsCCreducef]);

    switch (locality)
    {
        case eatAll:
            a0 = 0;
            na = nbs->natoms_nonlocal;
            break;
        case eatLocal:
            a0 = 0;
            na = nbs->natoms_local;
            break;
        case eatNonlocal:
            a0 = nbs->natoms_local;
            na = nbs->natoms_nonlocal - nbs->natoms_local;
            break;
    }

    int nth = gmx_omp_nthreads_get(emntNonbonded);

    if (nbat->nout > 1)
    {
        if (locality != eatAll)
        {
            gmx_incons("add_f_to_f called with nout>1 and locality!=eatAll");
        }

        /* Reduce the force thread output buffers into buffer 0, before adding
         * them to the, differently ordered, "real" force buffer.
         */
        if (nbat->bUseTreeReduce)
        {
            nbnxn_atomdata_add_nbat_f_to_f_treereduce(nbat, nth);
        }
        else
        {
            nbnxn_atomdata_add_nbat_f_to_f_stdreduce(nbat, nth);
        }
    }
#pragma omp parallel for num_threads(nth) schedule(static)
    for (int th = 0; th < nth; th++)
    {
        try
        {
            nbnxn_atomdata_add_nbat_f_to_f_part(nbs, nbat,
                                                nbat->out,
                                                1,
                                                a0+((th+0)*na)/nth,
                                                a0+((th+1)*na)/nth,
                                                f);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }

    nbs_cycle_stop(&nbs->cc[enbsCCreducef]);
}

/* Adds the shift forces from nbnxn_atomdata_t to fshift */
void nbnxn_atomdata_add_nbat_fshift_to_fshift(const nbnxn_atomdata_t *nbat,
                                              rvec                   *fshift)
{
    const nbnxn_atomdata_output_t * out = nbat->out;

    for (int s = 0; s < SHIFTS; s++)
    {
        rvec sum;
        clear_rvec(sum);
        for (int th = 0; th < nbat->nout; th++)
        {
            sum[XX] += out[th].fshift[s*DIM+XX];
            sum[YY] += out[th].fshift[s*DIM+YY];
            sum[ZZ] += out[th].fshift[s*DIM+ZZ];
        }
        rvec_inc(fshift[s], sum);
    }
}
