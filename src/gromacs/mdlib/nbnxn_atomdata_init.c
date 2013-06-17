/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <string.h>
#include "smalloc.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/simd/macros.h"
#include "vec.h"
#include "nbnxn_consts.h"
#include "nbnxn_internal.h"
#include "nbnxn_search.h"
#include "nbnxn_atomdata.h"

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

/* Determines the combination rule (or none) to be used, stores it,
 * and sets the LJ parameters required with the rule.
 */
static void set_combination_rule_data(nbnxn_atomdata_t *nbat)
{
    int  nt, i, j;
    real c6, c12;

    nt = nbat->ntype;

    switch (nbat->comb_rule)
    {
        case  ljcrGEOM:
            nbat->comb_rule = ljcrGEOM;

            for (i = 0; i < nt; i++)
            {
                /* Copy the diagonal from the nbfp matrix */
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
            break;
        default:
            gmx_incons("Unknown combination rule");
            break;
    }
}

void
nbnxn_atomdata_init_simple_exclusion_masks(nbnxn_atomdata_t *nbat)
{
    int       i, j;
    const int simd_width = GMX_SIMD_WIDTH_HERE;
    int       simd_excl_size;
    /* Set the diagonal cluster pair exclusion mask setup data.
     * In the kernel we check 0 < j - i to generate the masks.
     * Here we store j - i for generating the mask for the first i,
     * we substract 0.5 to avoid rounding issues.
     * In the kernel we can subtract 1 to generate the subsequent mask.
     */
    int       simd_4xn_diag_size;

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

#ifdef GMX_CPU_ACCELERATION_IBM_QPX
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

    real  simdFalse = -1, simdTrue = 1;
    real *simd_interaction_array;

    simd_excl_size = NBNXN_CPU_CLUSTER_I_SIZE*NBNXN_CPU_CLUSTER_I_SIZE;
    const int qpx_simd_width = 4;
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

/* Initializes an nbnxn_atomdata_t data structure */
void nbnxn_atomdata_init(FILE *fp,
                         nbnxn_atomdata_t *nbat,
                         int nb_kernel_type,
                         int ntype, const real *nbfp,
                         int n_energygroups,
                         int nout,
                         nbnxn_alloc_t *alloc,
                         nbnxn_free_t  *free)
{
    int      i, j;
    real     c6, c12, tol;
    char    *ptr;
    gmx_bool simple, bCombGeom, bCombLB;

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

    if (simple)
    {
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

        set_combination_rule_data(nbat);
    }
    else
    {
        nbat->comb_rule = ljcrNONE;

        nbat->free(nbat->nbfp_comb);
    }

    nbat->natoms  = 0;
    nbat->type    = NULL;
    nbat->lj_comb = NULL;
    if (simple)
    {
        int pack_x;

        switch (nb_kernel_type)
        {
            case nbnxnk4xN_SIMD_4xN:
            case nbnxnk4xN_SIMD_2xNN:
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
                break;
            default:
                nbat->XFormat = nbatXYZ;
                break;
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
}
