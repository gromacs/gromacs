/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2016,2017,2018, by the GROMACS development team, led by
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
 *
 * \brief
 * CUDA kernel for GPU version of copy_rvec_to_nbat_real (nbnxn_atomdata.cpp).
 * Converts coordinate data from rvec to nb format.
 *
 * Tested on hardware and range of inputs covered by Gromacs test framework (15th Oct 2018)
 *  \author Alan Gray <alang@nvidia.com>
 *  \author Jon Vincent <jvincent@nvidia.com>
 */


__global__ void nbnxn_gpu_x_to_nbat_x_kernel(int         ncxy,
                                             float      *xnb,
                                             int         g,
                                             bool        FillLocal,
                                             const rvec *x,
                                             const int  *a,
                                             const int  *na_all,
                                             const int  *cxy_ind,
                                             int         cell0,
                                             int         na_sc)
{


    const float farAway = -1000000;

    /* map cell-level parallelism to y component of CUDA block index */
    int cxy = blockIdx.y;

    if (cxy < ncxy)
    {

        int na = na_all[cxy];
        int a0 = (cell0 + cxy_ind[cxy])*na_sc;
        int na_round;
        if (g == 0 && FillLocal)
        {
            na_round =
                (cxy_ind[cxy+1] - cxy_ind[cxy])*na_sc;
        }
        else
        {
            /* We fill only the real particle locations.
             * We assume the filling entries at the end have been
             * properly set before during pair-list generation.
             */
            na_round = na;
        }

        /* map parlallism within a cell to x component of CUDA block index linearized
         * with threads within a block */
        int i, j0;
        i = blockIdx.x*blockDim.x+threadIdx.x;

        j0 = a0*STRIDE_XYZQ;

        /* perform conversion of each element */
        if (i < na_round)
        {
            if (i < na)
            {

                xnb[j0+4*i]   = x[a[a0+i]][XX];
                xnb[j0+4*i+1] = x[a[a0+i]][YY];
                xnb[j0+4*i+2] = x[a[a0+i]][ZZ];
            }
            else
            {
                xnb[j0+4*i]   = farAway;
                xnb[j0+4*i+1] = farAway;
                xnb[j0+4*i+2] = farAway;
            }
        }
    }

}

/* CUDA kernel to add part of the force array(s) from nbnxn_atomdata_t to f */
__global__ void
nbnxn_gpu_add_nbat_f_to_f_kernel(const real* fnb, rvec* pme_f,
                                 rvec* f, const int* cell,
                                 const int a0, const int a1,
                                 const int stride, const bool addPmeF)
{

    /* map particle-level parallelism to 1D CUDA thread and block index */
    int a = a0+blockIdx.x*blockDim.x+threadIdx.x;

    /* perform addition for each particle*/
    if (a < a1)
    {

        int   i = cell[a]*stride;

        float f0, f1, f2;

        f0 = f[a][XX];
        f1 = f[a][YY];
        f2 = f[a][ZZ];

        f0 += fnb[i];
        f1 += fnb[i+1];
        f2 += fnb[i+2];

        if (addPmeF)
        {
            f0 += pme_f[a][XX];
            f1 += pme_f[a][YY];
            f2 += pme_f[a][ZZ];
        }

        f[a][XX] = f0;
        f[a][YY] = f1;
        f[a][ZZ] = f2;



    }

    return;
}
