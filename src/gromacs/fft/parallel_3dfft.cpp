/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2005 David van der Spoel, Erik Lindahl, University of Groningen.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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

#include "parallel_3dfft.h"

#include <errno.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/fft/fft.h"
#include "gromacs/fft/fft5d.h"
#include "gromacs/math/gmxcomplex.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/smalloc.h"

struct gmx_parallel_3dfft  {
    fft5d_plan p1, p2;
};

int
gmx_parallel_3dfft_init   (gmx_parallel_3dfft_t     *    pfft_setup,
                           ivec                          ndata,
                           real     **                   real_data,
                           t_complex     **              complex_data,
                           MPI_Comm                      comm[2],
                           gmx_bool                      bReproducible,
                           int                           nthreads)
{
    int        rN      = ndata[2], M = ndata[1], K = ndata[0];
    int        flags   = FFT5D_REALCOMPLEX | FFT5D_ORDER_YZ; /* FFT5D_DEBUG */
    MPI_Comm   rcomm[] = {comm[1], comm[0]};
    int        Nb, Mb, Kb;                                   /* dimension for backtransform (in starting order) */
    t_complex *buf1, *buf2;                                  /*intermediate buffers - used internally.*/

    snew(*pfft_setup, 1);
    if (bReproducible)
    {
        flags |= FFT5D_NOMEASURE;
    }

    if (!(flags&FFT5D_ORDER_YZ))
    {
        Nb = M; Mb = K; Kb = rN;
    }
    else
    {
        Nb = K; Mb = rN; Kb = M;  /* currently always true because ORDER_YZ always set */
    }

    (*pfft_setup)->p1 = fft5d_plan_3d(rN, M, K, rcomm, flags, (t_complex**)real_data, complex_data, &buf1, &buf2, nthreads);

    (*pfft_setup)->p2 = fft5d_plan_3d(Nb, Mb, Kb, rcomm,
                                      (flags|FFT5D_BACKWARD|FFT5D_NOMALLOC)^FFT5D_ORDER_YZ, complex_data, (t_complex**)real_data, &buf1, &buf2, nthreads);

    return (*pfft_setup)->p1 != 0 && (*pfft_setup)->p2 != 0;
}


static int
fft5d_limits(fft5d_plan                p,
             ivec                      local_ndata,
             ivec                      local_offset,
             ivec                      local_size)
{
    local_offset[2] = 0;
    local_offset[1] = p->oM[0];  /*=p->coor[0]*p->MG/p->P[0]; */
    local_offset[0] = p->oK[0];  /*=p->coor[1]*p->KG/p->P[1]; */

    local_ndata[2] = p->rC[0];
    local_ndata[1] = p->pM[0];
    local_ndata[0] = p->pK[0];

    if ((!(p->flags&FFT5D_BACKWARD)) && (p->flags&FFT5D_REALCOMPLEX))
    {
        //C is length in multiples of complex local_size in multiples of real
        local_size[2] = p->C[0]*2;
    }
    else
    {
        local_size[2] = p->C[0];
    }
    local_size[1] = p->pM[0];
    local_size[0] = p->pK[0];
    return 0;
}

int
gmx_parallel_3dfft_real_limits(gmx_parallel_3dfft_t      pfft_setup,
                               ivec                      local_ndata,
                               ivec                      local_offset,
                               ivec                      local_size)
{
    return fft5d_limits(pfft_setup->p1, local_ndata, local_offset, local_size);
}

static void reorder_ivec_yzx(ivec v)
{
    int tmp;

    tmp   = v[0];
    v[XX] = v[2];
    v[ZZ] = v[1];
    v[YY] = tmp;
}

int
gmx_parallel_3dfft_complex_limits(gmx_parallel_3dfft_t      pfft_setup,
                                  ivec                      complex_order,
                                  ivec                      local_ndata,
                                  ivec                      local_offset,
                                  ivec                      local_size)
{
    int ret;

    /* For now everything is in-order, but prepare to save communication by avoiding transposes */
    complex_order[0] = 0;
    complex_order[1] = 1;
    complex_order[2] = 2;

    ret = fft5d_limits(pfft_setup->p2, local_ndata, local_offset, local_size);

    reorder_ivec_yzx(local_ndata);
    reorder_ivec_yzx(local_offset);
    reorder_ivec_yzx(local_size);

    return ret;
}


int
gmx_parallel_3dfft_execute(gmx_parallel_3dfft_t    pfft_setup,
                           enum gmx_fft_direction  dir,
                           int                     thread,
                           gmx_wallcycle_t         wcycle)
{
    if ((!(pfft_setup->p1->flags&FFT5D_REALCOMPLEX)) ^ (dir == GMX_FFT_FORWARD || dir == GMX_FFT_BACKWARD))
    {
        gmx_fatal(FARGS, "Invalid transform. Plan and execution don't match regarding reel/complex");
    }
    if (dir == GMX_FFT_FORWARD || dir == GMX_FFT_REAL_TO_COMPLEX)
    {
        fft5d_execute(pfft_setup->p1, thread, wcycle);
    }
    else
    {
        fft5d_execute(pfft_setup->p2, thread, wcycle);
    }
    return 0;
}

int
gmx_parallel_3dfft_destroy(gmx_parallel_3dfft_t    pfft_setup)
{
    if (pfft_setup)
    {
        fft5d_destroy(pfft_setup->p2);
        fft5d_destroy(pfft_setup->p1);
        sfree(pfft_setup);
    }
    return 0;
}
