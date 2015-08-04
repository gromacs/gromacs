/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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

#include "config.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "thread_mpi/atomic.h"

#include "gromacs/legacyheaders/gmx_omp_nthreads.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/nbnxn_consts.h"
#include "gromacs/mdlib/nbnxn_internal.h"
#include "gromacs/mdlib/nbnxn_search.h"
#include "gromacs/mdlib/nbnxn_simd.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/simd/simd.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/smalloc.h"

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

    if (nbytes_new > 0 && ptr_new == NULL)
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
    if (*ptr != NULL)
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
    int cj_size;

    out->f = NULL;
    ma((void **)&out->fshift, SHIFTS*DIM*sizeof(*out->fshift));
    out->nV = nenergrp*nenergrp;
    ma((void **)&out->Vvdw, out->nV*sizeof(*out->Vvdw));
    ma((void **)&out->Vc, out->nV*sizeof(*out->Vc  ));

    if (nb_kernel_type == nbnxnk4xN_SIMD_4xN ||
        nb_kernel_type == nbnxnk4xN_SIMD_2xNN)
    {
        cj_size  = nbnxn_kernel_to_cj_size(nb_kernel_type);
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

static void clear_nbat_real(int na, int nbatFormat, real *xnb, int a0)
{
    int a, d, j, c;

    switch (nbatFormat)
    {
        case nbatXYZ:
            for (a = 0; a < na; a++)
            {
                for (d = 0; d < DIM; d++)
                {
                    xnb[(a0+a)*STRIDE_XYZ+d] = 0;
                }
            }
            break;
        case nbatXYZQ:
            for (a = 0; a < na; a++)
            {
                for (d = 0; d < DIM; d++)
                {
                    xnb[(a0+a)*STRIDE_XYZQ+d] = 0;
                }
            }
            break;
        case nbatX4:
            j = X4_IND_A(a0);
            c = a0 & (PACK_X4-1);
            for (a = 0; a < na; a++)
            {
                xnb[j+XX*PACK_X4] = 0;
                xnb[j+YY*PACK_X4] = 0;
                xnb[j+ZZ*PACK_X4] = 0;
                j++;
                c++;
                if (c == PACK_X4)
                {
                    j += (DIM-1)*PACK_X4;
                    c  = 0;
                }
            }
            break;
        case nbatX8:
            j = X8_IND_A(a0);
            c = a0 & (PACK_X8-1);
            for (a = 0; a < na; a++)
            {
                xnb[j+XX*PACK_X8] = 0;
                xnb[j+YY*PACK_X8] = 0;
                xnb[j+ZZ*PACK_X8] = 0;
                j++;
                c++;
                if (c == PACK_X8)
                {
                    j += (DIM-1)*PACK_X8;
                    c  = 0;
                }
            }
            break;
    }
}

void copy_rvec_to_nbat_real(const int *a, int na, int na_round,
                            rvec *x, int nbatFormat, real *xnb, int a0,
                            int cx, int cy, int cz)
{
    int i, j, c;

/* We might need to place filler particles to fill up the cell to na_round.
 * The coefficients (LJ and q) for such particles are zero.
 * But we might still get NaN as 0*NaN when distances are too small.
 * We hope that -107 nm is far away enough from to zero
 * to avoid accidental short distances to particles shifted down for pbc.
 */
#define NBAT_FAR_AWAY 107

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
            /* Complete the partially filled last cell with copies of the last element.
             * This simplifies the bounding box calculation and avoid
             * numerical issues with atoms that are coincidentally close.
             */
            for (; i < na_round; i++)
            {
                xnb[j++] = -NBAT_FAR_AWAY*(1 + cx);
                xnb[j++] = -NBAT_FAR_AWAY*(1 + cy);
                xnb[j++] = -NBAT_FAR_AWAY*(1 + cz + i);
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
            /* Complete the partially filled last cell with particles far apart */
            for (; i < na_round; i++)
            {
                xnb[j++] = -NBAT_FAR_AWAY*(1 + cx);
                xnb[j++] = -NBAT_FAR_AWAY*(1 + cy);
                xnb[j++] = -NBAT_FAR_AWAY*(1 + cz + i);
                j++;
            }
            break;
        case nbatX4:
            j = X4_IND_A(a0);
            c = a0 & (PACK_X4-1);
            for (i = 0; i < na; i++)
            {
                xnb[j+XX*PACK_X4] = x[a[i]][XX];
                xnb[j+YY*PACK_X4] = x[a[i]][YY];
                xnb[j+ZZ*PACK_X4] = x[a[i]][ZZ];
                j++;
                c++;
                if (c == PACK_X4)
                {
                    j += (DIM-1)*PACK_X4;
                    c  = 0;
                }
            }
            /* Complete the partially filled last cell with particles far apart */
            for (; i < na_round; i++)
            {
                xnb[j+XX*PACK_X4] = -NBAT_FAR_AWAY*(1 + cx);
                xnb[j+YY*PACK_X4] = -NBAT_FAR_AWAY*(1 + cy);
                xnb[j+ZZ*PACK_X4] = -NBAT_FAR_AWAY*(1 + cz + i);
                j++;
                c++;
                if (c == PACK_X4)
                {
                    j += (DIM-1)*PACK_X4;
                    c  = 0;
                }
            }
            break;
        case nbatX8:
            j = X8_IND_A(a0);
            c = a0 & (PACK_X8 - 1);
            for (i = 0; i < na; i++)
            {
                xnb[j+XX*PACK_X8] = x[a[i]][XX];
                xnb[j+YY*PACK_X8] = x[a[i]][YY];
                xnb[j+ZZ*PACK_X8] = x[a[i]][ZZ];
                j++;
                c++;
                if (c == PACK_X8)
                {
                    j += (DIM-1)*PACK_X8;
                    c  = 0;
                }
            }
            /* Complete the partially filled last cell with particles far apart */
            for (; i < na_round; i++)
            {
                xnb[j+XX*PACK_X8] = -NBAT_FAR_AWAY*(1 + cx);
                xnb[j+YY*PACK_X8] = -NBAT_FAR_AWAY*(1 + cy);
                xnb[j+ZZ*PACK_X8] = -NBAT_FAR_AWAY*(1 + cz + i);
                j++;
                c++;
                if (c == PACK_X8)
                {
                    j += (DIM-1)*PACK_X8;
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
    int  nt, i, j;
    real c6, c12;

    nt = nbat->ntype;

    if (bSIMD)
    {
        /* nbfp_s4 stores two parameters using a stride of 4,
         * because this would suit x86 SIMD single-precision
         * quad-load intrinsics. There's a slight inefficiency in
         * allocating and initializing nbfp_s4 when it might not
         * be used, but introducing the conditional code is not
         * really worth it. */
        nbat->alloc((void **)&nbat->nbfp_s4, nt*nt*4*sizeof(*nbat->nbfp_s4));
        for (i = 0; i < nt; i++)
        {
            for (j = 0; j < nt; j++)
            {
                nbat->nbfp_s4[(i*nt+j)*4+0] = nbat->nbfp[(i*nt+j)*2+0];
                nbat->nbfp_s4[(i*nt+j)*4+1] = nbat->nbfp[(i*nt+j)*2+1];
                nbat->nbfp_s4[(i*nt+j)*4+2] = 0;
                nbat->nbfp_s4[(i*nt+j)*4+3] = 0;
            }
        }
    }

    /* We use combination rule data for SIMD combination rule kernels
     * and with LJ-PME kernels. We then only need parameters per atom type,
     * not per pair of atom types.
     */
    switch (nbat->comb_rule)
    {
        case ljcrGEOM:
            nbat->comb_rule = ljcrGEOM;

            for (i = 0; i < nt; i++)
            {
                /* Store the sqrt of the diagonal from the nbfp matrix */
                nbat->nbfp_comb[i*2  ] = sqrt(nbat->nbfp[(i*nt+i)*2  ]);
                nbat->nbfp_comb[i*2+1] = sqrt(nbat->nbfp[(i*nt+i)*2+1]);
            }
            break;
        case ljcrLB:
            for (i = 0; i < nt; i++)
            {
                /* Get 6*C6 and 12*C12 from the diagonal of the nbfp matrix */
                c6  = nbat->nbfp[(i*nt+i)*2  ];
                c12 = nbat->nbfp[(i*nt+i)*2+1];
                if (c6 > 0 && c12 > 0)
                {
                    /* We store 0.5*2^1/6*sigma and sqrt(4*3*eps),
                     * so we get 6*C6 and 12*C12 after combining.
                     */
                    nbat->nbfp_comb[i*2  ] = 0.5*pow(c12/c6, 1.0/6.0);
                    nbat->nbfp_comb[i*2+1] = sqrt(c6*c6/c12);
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

#ifdef GMX_NBNXN_SIMD
static void
nbnxn_atomdata_init_simple_exclusion_masks(nbnxn_atomdata_t *nbat)
{
    int       i, j;
    const int simd_width = GMX_SIMD_REAL_WIDTH;
    int       simd_excl_size;
    /* Set the diagonal cluster pair exclusion mask setup data.
     * In the kernel we check 0 < j - i to generate the masks.
     * Here we store j - i for generating the mask for the first i,
     * we substract 0.5 to avoid rounding issues.
     * In the kernel we can subtract 1 to generate the subsequent mask.
     */
    int        simd_4xn_diag_size;
    const real simdFalse = -1, simdTrue = 1;
    real      *simd_interaction_array;

    simd_4xn_diag_size = max(NBNXN_CPU_CLUSTER_I_SIZE, simd_width);
    snew_aligned(nbat->simd_4xn_diagonal_j_minus_i, simd_4xn_diag_size, NBNXN_MEM_ALIGN);
    for (j = 0; j < simd_4xn_diag_size; j++)
    {
        nbat->simd_4xn_diagonal_j_minus_i[j] = j - 0.5;
    }

    snew_aligned(nbat->simd_2xnn_diagonal_j_minus_i, simd_width, NBNXN_MEM_ALIGN);
    for (j = 0; j < simd_width/2; j++)
    {
        /* The j-cluster size is half the SIMD width */
        nbat->simd_2xnn_diagonal_j_minus_i[j]              = j - 0.5;
        /* The next half of the SIMD width is for i + 1 */
        nbat->simd_2xnn_diagonal_j_minus_i[simd_width/2+j] = j - 1 - 0.5;
    }

    /* We use up to 32 bits for exclusion masking.
     * The same masks are used for the 4xN and 2x(N+N) kernels.
     * The masks are read either into epi32 SIMD registers or into
     * real SIMD registers (together with a cast).
     * In single precision this means the real and epi32 SIMD registers
     * are of equal size.
     * In double precision the epi32 registers can be smaller than
     * the real registers, so depending on the architecture, we might
     * need to use two, identical, 32-bit masks per real.
     */
    simd_excl_size = NBNXN_CPU_CLUSTER_I_SIZE*simd_width;
    snew_aligned(nbat->simd_exclusion_filter1, simd_excl_size,   NBNXN_MEM_ALIGN);
    snew_aligned(nbat->simd_exclusion_filter2, simd_excl_size*2, NBNXN_MEM_ALIGN);

    for (j = 0; j < simd_excl_size; j++)
    {
        /* Set the consecutive bits for masking pair exclusions */
        nbat->simd_exclusion_filter1[j]       = (1U << j);
        nbat->simd_exclusion_filter2[j*2 + 0] = (1U << j);
        nbat->simd_exclusion_filter2[j*2 + 1] = (1U << j);
    }

#if (defined GMX_SIMD_IBM_QPX)
    /* The QPX kernels shouldn't do the bit masking that is done on
     * x86, because the SIMD units lack bit-wise operations. Instead,
     * we generate a vector of all 2^4 possible ways an i atom
     * interacts with its 4 j atoms. Each array entry contains
     * simd_width signed ints that are read in a single SIMD
     * load. These ints must contain values that will be interpreted
     * as true and false when loaded in the SIMD floating-point
     * registers, ie. any positive or any negative value,
     * respectively. Each array entry encodes how this i atom will
     * interact with the 4 j atoms. Matching code exists in
     * set_ci_top_excls() to generate indices into this array. Those
     * indices are used in the kernels. */

    simd_excl_size = NBNXN_CPU_CLUSTER_I_SIZE*NBNXN_CPU_CLUSTER_I_SIZE;
    const int qpx_simd_width = GMX_SIMD_REAL_WIDTH;
    snew_aligned(simd_interaction_array, simd_excl_size * qpx_simd_width, NBNXN_MEM_ALIGN);
    for (j = 0; j < simd_excl_size; j++)
    {
        int index = j * qpx_simd_width;
        for (i = 0; i < qpx_simd_width; i++)
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
    int      i, j, nth;
    real     c6, c12, tol;
    char    *ptr;
    gmx_bool simple, bCombGeom, bCombLB, bSIMD;

    if (alloc == NULL)
    {
        nbat->alloc = nbnxn_alloc_aligned;
    }
    else
    {
        nbat->alloc = alloc;
    }
    if (free == NULL)
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
    if (ptr != NULL)
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
    for (i = 0; i < ntype; i++)
    {
        c6  = nbfp[(i*ntype+i)*2  ]/6.0;
        c12 = nbfp[(i*ntype+i)*2+1]/12.0;
        if (c6 > 0 && c12 > 0)
        {
            nbat->nbfp_comb[i*2  ] = pow(c12/c6, 1.0/6.0);
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

    for (i = 0; i < nbat->ntype; i++)
    {
        for (j = 0; j < nbat->ntype; j++)
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
                      gmx_within_tol(pow(c12/c6, 1.0/6.0), 0.5*(nbat->nbfp_comb[i*2]+nbat->nbfp_comb[j*2]), tol) &&
                      gmx_within_tol(0.25*c6*c6/c12, sqrt(nbat->nbfp_comb[i*2+1]*nbat->nbfp_comb[j*2+1]), tol)));
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
    nbat->type    = NULL;
    nbat->lj_comb = NULL;
    if (simple)
    {
        int pack_x;

        if (bSIMD)
        {
            pack_x = max(NBNXN_CPU_CLUSTER_I_SIZE,
                         nbnxn_kernel_to_cj_size(nb_kernel_type));
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
    nbat->q        = NULL;
    nbat->nenergrp = n_energygroups;
    if (!simple)
    {
        /* Energy groups not supported yet for super-sub lists */
        if (n_energygroups > 1 && fp != NULL)
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
    nbat->energrp = NULL;
    nbat->alloc((void **)&nbat->shift_vec, SHIFTS*sizeof(*nbat->shift_vec));
    nbat->xstride = (nbat->XFormat == nbatXYZQ ? STRIDE_XYZQ : DIM);
    nbat->fstride = (nbat->FFormat == nbatXYZQ ? STRIDE_XYZQ : DIM);
    nbat->x       = NULL;

#ifdef GMX_NBNXN_SIMD
    if (simple)
    {
        nbnxn_atomdata_init_simple_exclusion_masks(nbat);
    }
#endif

    /* Initialize the output data structures */
    nbat->nout    = nout;
    snew(nbat->out, nbat->nout);
    nbat->nalloc  = 0;
    for (i = 0; i < nbat->nout; i++)
    {
        nbnxn_atomdata_output_init(&nbat->out[i],
                                   nb_kernel_type,
                                   nbat->nenergrp, 1<<nbat->neg_2log,
                                   nbat->alloc);
    }
    nbat->buffer_flags.flag        = NULL;
    nbat->buffer_flags.flag_nalloc = 0;

    nth = gmx_omp_nthreads_get(emntNonbonded);

    ptr = getenv("GMX_USE_TREEREDUCE");
    if (ptr != NULL)
    {
        nbat->bUseTreeReduce = strtol(ptr, 0, 10);
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

static void copy_lj_to_nbat_lj_comb_x4(const real *ljparam_type,
                                       const int *type, int na,
                                       real *ljparam_at)
{
    int is, k, i;

    /* The LJ params follow the combination rule:
     * copy the params for the type array to the atom array.
     */
    for (is = 0; is < na; is += PACK_X4)
    {
        for (k = 0; k < PACK_X4; k++)
        {
            i = is + k;
            ljparam_at[is*2        +k] = ljparam_type[type[i]*2  ];
            ljparam_at[is*2+PACK_X4+k] = ljparam_type[type[i]*2+1];
        }
    }
}

static void copy_lj_to_nbat_lj_comb_x8(const real *ljparam_type,
                                       const int *type, int na,
                                       real *ljparam_at)
{
    int is, k, i;

    /* The LJ params follow the combination rule:
     * copy the params for the type array to the atom array.
     */
    for (is = 0; is < na; is += PACK_X8)
    {
        for (k = 0; k < PACK_X8; k++)
        {
            i = is + k;
            ljparam_at[is*2        +k] = ljparam_type[type[i]*2  ];
            ljparam_at[is*2+PACK_X8+k] = ljparam_type[type[i]*2+1];
        }
    }
}

/* Sets the atom type in nbnxn_atomdata_t */
static void nbnxn_atomdata_set_atomtypes(nbnxn_atomdata_t    *nbat,
                                         int                  ngrid,
                                         const nbnxn_search_t nbs,
                                         const int           *type)
{
    int                 g, i, ncz, ash;
    const nbnxn_grid_t *grid;

    for (g = 0; g < ngrid; g++)
    {
        grid = &nbs->grid[g];

        /* Loop over all columns and copy and fill */
        for (i = 0; i < grid->ncx*grid->ncy; i++)
        {
            ncz = grid->cxy_ind[i+1] - grid->cxy_ind[i];
            ash = (grid->cell0 + grid->cxy_ind[i])*grid->na_sc;

            copy_int_to_nbat_int(nbs->a+ash, grid->cxy_na[i], ncz*grid->na_sc,
                                 type, nbat->ntype-1, nbat->type+ash);
        }
    }
}

/* Sets the LJ combination rule parameters in nbnxn_atomdata_t */
static void nbnxn_atomdata_set_ljcombparams(nbnxn_atomdata_t    *nbat,
                                            int                  ngrid,
                                            const nbnxn_search_t nbs)
{
    int                 g, i, ncz, ash;
    const nbnxn_grid_t *grid;

    if (nbat->comb_rule != ljcrNONE)
    {
        for (g = 0; g < ngrid; g++)
        {
            grid = &nbs->grid[g];

            /* Loop over all columns and copy and fill */
            for (i = 0; i < grid->ncx*grid->ncy; i++)
            {
                ncz = grid->cxy_ind[i+1] - grid->cxy_ind[i];
                ash = (grid->cell0 + grid->cxy_ind[i])*grid->na_sc;

                if (nbat->XFormat == nbatX4)
                {
                    copy_lj_to_nbat_lj_comb_x4(nbat->nbfp_comb,
                                               nbat->type+ash, ncz*grid->na_sc,
                                               nbat->lj_comb+ash*2);
                }
                else if (nbat->XFormat == nbatX8)
                {
                    copy_lj_to_nbat_lj_comb_x8(nbat->nbfp_comb,
                                               nbat->type+ash, ncz*grid->na_sc,
                                               nbat->lj_comb+ash*2);
                }
            }
        }
    }
}

/* Sets the charges in nbnxn_atomdata_t *nbat */
static void nbnxn_atomdata_set_charges(nbnxn_atomdata_t    *nbat,
                                       int                  ngrid,
                                       const nbnxn_search_t nbs,
                                       const real          *charge)
{
    int                 g, cxy, ncz, ash, na, na_round, i, j;
    real               *q;
    const nbnxn_grid_t *grid;

    for (g = 0; g < ngrid; g++)
    {
        grid = &nbs->grid[g];

        /* Loop over all columns and copy and fill */
        for (cxy = 0; cxy < grid->ncx*grid->ncy; cxy++)
        {
            ash      = (grid->cell0 + grid->cxy_ind[cxy])*grid->na_sc;
            na       = grid->cxy_na[cxy];
            na_round = (grid->cxy_ind[cxy+1] - grid->cxy_ind[cxy])*grid->na_sc;

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
                                    int                  ngrid,
                                    const nbnxn_search_t nbs)
{
    real               *q;
    int                 stride_q, g, nsubc, c_offset, c, subc, i, ind;
    const nbnxn_grid_t *grid;

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

    for (g = 0; g < ngrid; g++)
    {
        grid = &nbs->grid[g];
        if (grid->bSimple)
        {
            nsubc = 1;
        }
        else
        {
            nsubc = GPU_NSUBCELL;
        }

        c_offset = grid->cell0*grid->na_sc;

        /* Loop over all columns and copy and fill */
        for (c = 0; c < grid->nc*nsubc; c++)
        {
            /* Does this cluster contain perturbed particles? */
            if (grid->fep[c] != 0)
            {
                for (i = 0; i < grid->na_c; i++)
                {
                    /* Is this a perturbed particle? */
                    if (grid->fep[c] & (1 << i))
                    {
                        ind = c_offset + c*grid->na_c + i;
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
    int i, j, sa, at;
    int comb;

    j = 0;
    for (i = 0; i < na; i += na_c)
    {
        /* Store na_c energy group numbers into one int */
        comb = 0;
        for (sa = 0; sa < na_c; sa++)
        {
            at = a[i+sa];
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
                                            int                  ngrid,
                                            const nbnxn_search_t nbs,
                                            const int           *atinfo)
{
    int                 g, i, ncz, ash;
    const nbnxn_grid_t *grid;

    if (nbat->nenergrp == 1)
    {
        return;
    }

    for (g = 0; g < ngrid; g++)
    {
        grid = &nbs->grid[g];

        /* Loop over all columns and copy and fill */
        for (i = 0; i < grid->ncx*grid->ncy; i++)
        {
            ncz = grid->cxy_ind[i+1] - grid->cxy_ind[i];
            ash = (grid->cell0 + grid->cxy_ind[i])*grid->na_sc;

            copy_egp_to_nbat_egps(nbs->a+ash, grid->cxy_na[i], ncz*grid->na_sc,
                                  nbat->na_c, nbat->neg_2log,
                                  atinfo, nbat->energrp+(ash>>grid->na_c_2log));
        }
    }
}

/* Sets all required atom parameter data in nbnxn_atomdata_t */
void nbnxn_atomdata_set(nbnxn_atomdata_t    *nbat,
                        int                  locality,
                        const nbnxn_search_t nbs,
                        const t_mdatoms     *mdatoms,
                        const int           *atinfo)
{
    int ngrid;

    if (locality == eatLocal)
    {
        ngrid = 1;
    }
    else
    {
        ngrid = nbs->ngrid;
    }

    nbnxn_atomdata_set_atomtypes(nbat, ngrid, nbs, mdatoms->typeA);

    nbnxn_atomdata_set_charges(nbat, ngrid, nbs, mdatoms->chargeA);

    if (nbs->bFEP)
    {
        nbnxn_atomdata_mask_fep(nbat, ngrid, nbs);
    }

    /* This must be done after masking types for FEP */
    nbnxn_atomdata_set_ljcombparams(nbat, ngrid, nbs);

    nbnxn_atomdata_set_energygroups(nbat, ngrid, nbs, atinfo);
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
        int g;

        for (g = g0; g < g1; g++)
        {
            const nbnxn_grid_t *grid;
            int                 cxy0, cxy1, cxy;

            grid = &nbs->grid[g];

            cxy0 = (grid->ncx*grid->ncy* th   +nth-1)/nth;
            cxy1 = (grid->ncx*grid->ncy*(th+1)+nth-1)/nth;

            for (cxy = cxy0; cxy < cxy1; cxy++)
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
                     * properly set before during ns.
                     */
                    na_fill = na;
                }
                copy_rvec_to_nbat_real(nbs->a+ash, na, na_fill, x,
                                       nbat->XFormat, nbat->x, ash,
                                       0, 0, 0);
            }
        }
    }
}

static void
nbnxn_atomdata_clear_reals(real * gmx_restrict dest,
                           int i0, int i1)
{
    int i;

    for (i = i0; i < i1; i++)
    {
        dest[i] = 0;
    }
}

static void
nbnxn_atomdata_reduce_reals(real * gmx_restrict dest,
                            gmx_bool bDestSet,
                            real ** gmx_restrict src,
                            int nsrc,
                            int i0, int i1)
{
    int i, s;

    if (bDestSet)
    {
        /* The destination buffer contains data, add to it */
        for (i = i0; i < i1; i++)
        {
            for (s = 0; s < nsrc; s++)
            {
                dest[i] += src[s][i];
            }
        }
    }
    else
    {
        /* The destination buffer is unitialized, set it first */
        for (i = i0; i < i1; i++)
        {
            dest[i] = src[0][i];
            for (s = 1; s < nsrc; s++)
            {
                dest[i] += src[s][i];
            }
        }
    }
}

static void
nbnxn_atomdata_reduce_reals_simd(real gmx_unused * gmx_restrict dest,
                                 gmx_bool gmx_unused bDestSet,
                                 real gmx_unused ** gmx_restrict src,
                                 int gmx_unused nsrc,
                                 int gmx_unused i0, int gmx_unused i1)
{
#ifdef GMX_NBNXN_SIMD
/* The SIMD width here is actually independent of that in the kernels,
 * but we use the same width for simplicity (usually optimal anyhow).
 */
    int             i, s;
    gmx_simd_real_t dest_SSE, src_SSE;

    if (bDestSet)
    {
        for (i = i0; i < i1; i += GMX_SIMD_REAL_WIDTH)
        {
            dest_SSE = gmx_simd_load_r(dest+i);
            for (s = 0; s < nsrc; s++)
            {
                src_SSE  = gmx_simd_load_r(src[s]+i);
                dest_SSE = gmx_simd_add_r(dest_SSE, src_SSE);
            }
            gmx_simd_store_r(dest+i, dest_SSE);
        }
    }
    else
    {
        for (i = i0; i < i1; i += GMX_SIMD_REAL_WIDTH)
        {
            dest_SSE = gmx_simd_load_r(src[0]+i);
            for (s = 1; s < nsrc; s++)
            {
                src_SSE  = gmx_simd_load_r(src[s]+i);
                dest_SSE = gmx_simd_add_r(dest_SSE, src_SSE);
            }
            gmx_simd_store_r(dest+i, dest_SSE);
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
    int         a, i, fa;
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

                for (a = a0; a < a1; a++)
                {
                    i = cell[a]*nbat->fstride;

                    f[a][XX] += fnb[i];
                    f[a][YY] += fnb[i+1];
                    f[a][ZZ] += fnb[i+2];
                }
            }
            else
            {
                for (a = a0; a < a1; a++)
                {
                    i = cell[a]*nbat->fstride;

                    for (fa = 0; fa < nfa; fa++)
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

                for (a = a0; a < a1; a++)
                {
                    i = X4_IND_A(cell[a]);

                    f[a][XX] += fnb[i+XX*PACK_X4];
                    f[a][YY] += fnb[i+YY*PACK_X4];
                    f[a][ZZ] += fnb[i+ZZ*PACK_X4];
                }
            }
            else
            {
                for (a = a0; a < a1; a++)
                {
                    i = X4_IND_A(cell[a]);

                    for (fa = 0; fa < nfa; fa++)
                    {
                        f[a][XX] += out[fa].f[i+XX*PACK_X4];
                        f[a][YY] += out[fa].f[i+YY*PACK_X4];
                        f[a][ZZ] += out[fa].f[i+ZZ*PACK_X4];
                    }
                }
            }
            break;
        case nbatX8:
            if (nfa == 1)
            {
                fnb = out[0].f;

                for (a = a0; a < a1; a++)
                {
                    i = X8_IND_A(cell[a]);

                    f[a][XX] += fnb[i+XX*PACK_X8];
                    f[a][YY] += fnb[i+YY*PACK_X8];
                    f[a][ZZ] += fnb[i+ZZ*PACK_X8];
                }
            }
            else
            {
                for (a = a0; a < a1; a++)
                {
                    i = X8_IND_A(cell[a]);

                    for (fa = 0; fa < nfa; fa++)
                    {
                        f[a][XX] += out[fa].f[i+XX*PACK_X8];
                        f[a][YY] += out[fa].f[i+YY*PACK_X8];
                        f[a][ZZ] += out[fa].f[i+ZZ*PACK_X8];
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

    int next_pow2 = 1<<(gmx_log2i(nth-1)+1);

    assert(nbat->nout == nth); /* tree-reduce currently only works for nout==nth */

    memset(nbat->syncStep, 0, sizeof(*(nbat->syncStep))*nth);

#pragma omp parallel num_threads(nth)
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
#else           /* TMPI_ATOMICS */
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
#ifdef GMX_NBNXN_SIMD
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
}


static void nbnxn_atomdata_add_nbat_f_to_f_stdreduce(const nbnxn_atomdata_t *nbat,
                                                     int                     nth)
{
    int th;
#pragma omp parallel for num_threads(nth) schedule(static)
    for (th = 0; th < nth; th++)
    {
        const nbnxn_buffer_flags_t *flags;
        int   b0, b1, b;
        int   i0, i1;
        int   nfptr;
        real *fptr[NBNXN_BUFFERFLAG_MAX_THREADS];
        int   out;

        flags = &nbat->buffer_flags;

        /* Calculate the cell-block range for our thread */
        b0 = (flags->nflag* th   )/nth;
        b1 = (flags->nflag*(th+1))/nth;

        for (b = b0; b < b1; b++)
        {
            i0 =  b   *NBNXN_BUFFERFLAG_SIZE*nbat->fstride;
            i1 = (b+1)*NBNXN_BUFFERFLAG_SIZE*nbat->fstride;

            nfptr = 0;
            for (out = 1; out < nbat->nout; out++)
            {
                if (bitmask_is_set(flags->flag[b], out))
                {
                    fptr[nfptr++] = nbat->out[out].f;
                }
            }
            if (nfptr > 0)
            {
#ifdef GMX_NBNXN_SIMD
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
}

/* Add the force array(s) from nbnxn_atomdata_t to f */
void nbnxn_atomdata_add_nbat_f_to_f(const nbnxn_search_t    nbs,
                                    int                     locality,
                                    const nbnxn_atomdata_t *nbat,
                                    rvec                   *f)
{
    int a0 = 0, na = 0;
    int nth, th;

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

    nth = gmx_omp_nthreads_get(emntNonbonded);

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
    for (th = 0; th < nth; th++)
    {
        nbnxn_atomdata_add_nbat_f_to_f_part(nbs, nbat,
                                            nbat->out,
                                            1,
                                            a0+((th+0)*na)/nth,
                                            a0+((th+1)*na)/nth,
                                            f);
    }

    nbs_cycle_stop(&nbs->cc[enbsCCreducef]);
}

/* Adds the shift forces from nbnxn_atomdata_t to fshift */
void nbnxn_atomdata_add_nbat_fshift_to_fshift(const nbnxn_atomdata_t *nbat,
                                              rvec                   *fshift)
{
    const nbnxn_atomdata_output_t *out;
    int  th;
    int  s;
    rvec sum;

    out = nbat->out;

    for (s = 0; s < SHIFTS; s++)
    {
        clear_rvec(sum);
        for (th = 0; th < nbat->nout; th++)
        {
            sum[XX] += out[th].fshift[s*DIM+XX];
            sum[YY] += out[th].fshift[s*DIM+YY];
            sum[ZZ] += out[th].fshift[s*DIM+ZZ];
        }
        rvec_inc(fshift[s], sum);
    }
}
