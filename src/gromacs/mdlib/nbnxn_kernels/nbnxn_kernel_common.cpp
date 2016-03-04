/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016, by the GROMACS development team, led by
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

#include "nbnxn_kernel_common.h"

#include <assert.h>

#include "gromacs/pbcutil/ishift.h"
#include "gromacs/simd/simd.h"

using namespace gmx; // TODO: Remove when this file is moved into gmx namespace

static void
clear_f_all(const nbnxn_atomdata_t *nbat, real *f)
{
#if GMX_SIMD_HAVE_REAL
    /* f is aligned and padded with  NBNXN_BUFFERFLAG_SIZE */
    assert(NBNXN_BUFFERFLAG_SIZE % GMX_SIMD_REAL_WIDTH == 0);
    SimdReal zero = setZero();
    for (int i = 0; i < nbat->natoms*nbat->fstride; i += GMX_SIMD_REAL_WIDTH)
    {
        store(f + i, zero);
    }
#else
    for (int i = 0; i < nbat->natoms*nbat->fstride; i++)
    {
        f[i] = 0;
    }
#endif
}

static void
clear_f_flagged(const nbnxn_atomdata_t *nbat, int output_index, real *f)
{
    const nbnxn_buffer_flags_t *flags;
    gmx_bitmask_t               our_flag;

    flags = &nbat->buffer_flags;

    bitmask_init_bit(&our_flag, output_index);

    for (int b = 0; b < flags->nflag; b++)
    {
        if (!bitmask_is_disjoint(flags->flag[b], our_flag))
        {
            int i0, i1;

            i0 = b*NBNXN_BUFFERFLAG_SIZE*nbat->fstride;
            i1 = i0 + NBNXN_BUFFERFLAG_SIZE*nbat->fstride;
#if GMX_SIMD_HAVE_REAL
            /* f is aligned and padded with  NBNXN_BUFFERFLAG_SIZE */
            assert(NBNXN_BUFFERFLAG_SIZE % GMX_SIMD_REAL_WIDTH == 0);
            SimdReal zero = setZero();
            for (int i = i0; i < i1; i += GMX_SIMD_REAL_WIDTH)
            {
                store(f + i, zero);
            }
#else
            for (int i = i0; i < i1; i++)
            {
                f[i] = 0;
            }
#endif
        }
    }
}

void
clear_f(const nbnxn_atomdata_t *nbat, int output_index, real *f)
{
    if (nbat->bUseBufferFlags)
    {
        clear_f_flagged(nbat, output_index, f);
    }
    else
    {
        clear_f_all(nbat, f);
    }
}

void
clear_fshift(real *fshift)
{
    int i;

    for (i = 0; i < SHIFTS*DIM; i++)
    {
        fshift[i] = 0;
    }
}

void
reduce_energies_over_lists(const nbnxn_atomdata_t     *nbat,
                           int                         nlist,
                           real                       *Vvdw,
                           real                       *Vc)
{
    int nb;
    int i, j, ind, indr;

    for (nb = 0; nb < nlist; nb++)
    {
        for (i = 0; i < nbat->nenergrp; i++)
        {
            /* Reduce the diagonal terms */
            ind        = i*nbat->nenergrp + i;
            Vvdw[ind] += nbat->out[nb].Vvdw[ind];
            Vc[ind]   += nbat->out[nb].Vc[ind];

            /* Reduce the off-diagonal terms */
            for (j = i+1; j < nbat->nenergrp; j++)
            {
                /* The output should contain only one off-diagonal part */
                ind        = i*nbat->nenergrp + j;
                indr       = j*nbat->nenergrp + i;
                Vvdw[ind] += nbat->out[nb].Vvdw[ind] + nbat->out[nb].Vvdw[indr];
                Vc[ind]   += nbat->out[nb].Vc[ind]   + nbat->out[nb].Vc[indr];
            }
        }
    }
}
