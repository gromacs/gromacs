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
/*! \file
 * \brief functionality for performing update and constraints on the gpu
 * \author Jon Vincent <jvincent@nvidia.com> and Alan Gray <alang@nvidia.com>
 */

//
#include "gmxpre.h"
#include <assert.h>
#include <stdio.h>
#include <cmath>
#include <algorithm>
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/domdec_specatomcomm.h"
#include "gromacs/ewald/pme.cuh"
#include "gromacs/nbnxm/nbnxm_gpu.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/utility/gmxmpi.h"
//#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/nbnxm/cuda/gpuUpdateConstraintsCUDA.h"
#include "gromacs/nbnxm/cuda/gpuD2DCUDA.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/ewald/pme.cuh"
#include "gromacs/mdlib/lincs.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/gmxmpi.h"
//#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/pbcutil/pbc.h"


#include "omp.h"

#if defined(_MSVC)
#include <limits>
#endif

struct gpuUpdateConstraintsData{

    rvec  * invMassPerDim_d;
    real  * invmass_d;
    real  * rhs1_d;
    real  * rhs2_d;
    real  * sol_d;
    real  * blc_sol_d;
    rvec  * r_d;
    int   * bla_d;
    real  * blc_d;
    real  * bllen_d;
    real  * mlambda_d;
    int   * nlocat_d;
    real  * blcc_d;
    real  * blmf_d;
    int   * blnr_d;
    int   * blbnb_d;
    real  * vir_r_m_dr_d;
    t_pbc * pbc_d;
    rvec  * pbc_hbox_diag_d;
    matrix* pbc_box_d;
    rvec  * pbc_mhbox_diag_d;
    rvec *  pbc_fbox_diag_d;
    rvec *  pbc_tric_vec_d;
    ivec *  pbc_tric_shift_d;
    int   * settled_ow1_d;
    int   * settled_hw2_d;
    int   * settled_hw3_d;
    real  * settled_virfac_d;
    bool  * bErrorHasOccurred_d;
    real  * lambdaGroupArray;
    real  * lambdaGroupArray_d;


    cudaStream_t     stream;

    int              xsize;

    bool             bNS;
    bool             bNSNextStep;
    bool             bPressureStep;
    bool             copybackVelocity;

    gmx_nbnxn_gpu_t* gpu_nbv;

};

//TEMPORARY Global variable to store a copy of above struct within this module.
//TODO move this to a central location (e.g. gpu_nbv) and pass through fn args.
static gpuUpdateConstraintsData gpuUCDmod;




/*-------------------------------- CUDA kernels-------------------------------- */
/*------------------------------------------------------------------------------*/

// LINCS Constraints Kernels - see original CPU versions in lincs.cpp

__global__ void lincs_init_rhs_sol_kernel(int         ncons,
                                          t_pbc      *pbc,
                                          const rvec *x,
                                          const rvec *xp,
                                          rvec       *r,
                                          real       *rhs1,
                                          real       *sol,
                                          const int  *bla,
                                          const real *blc,
                                          const real *bllen,
                                          rvec      * pbc_hbox_diag,
                                          matrix    * pbc_box,
                                          rvec      * pbc_mhbox_diag,
                                          rvec     *  pbc_fbox_diag,
                                          rvec     *  pbc_tric_vec,
                                          ivec     *  pbc_tric_shift,
                                          bool        bpbc)
{

    int b = blockIdx.x*blockDim.x+threadIdx.x;

    if (b < ncons)
    {

        rvec dx;
        real mvb;

        if (bpbc)
        {
            pbc_dx_aiuc_gpu(pbc, x[bla[2*b]], x[bla[2*b+1]], dx, *pbc_hbox_diag,
                            *pbc_box, *pbc_mhbox_diag, *pbc_fbox_diag, pbc_tric_vec,
                            pbc_tric_shift);

            unitv_gpu(dx, r[b]);

            pbc_dx_aiuc_gpu(pbc, xp[bla[2*b]], xp[bla[2*b+1]], dx, *pbc_hbox_diag,
                            *pbc_box, *pbc_mhbox_diag, *pbc_fbox_diag, pbc_tric_vec,
                            pbc_tric_shift);

            mvb     = blc[b]*(iprod_gpu(r[b], dx) - bllen[b]);
            rhs1[b] = mvb;
            sol[b]  = mvb;
        }
        else
        {
            real tmp0, tmp1, tmp2, rlen, mvb;

            int  i       = bla[2*b];
            int  j       = bla[2*b+1];
            tmp0    = x[i][0] - x[j][0];
            tmp1    = x[i][1] - x[j][1];
            tmp2    = x[i][2] - x[j][2];
            rlen    = invsqrt(tmp0*tmp0 + tmp1*tmp1 + tmp2*tmp2);
            r[b][0] = rlen*tmp0;
            r[b][1] = rlen*tmp1;
            r[b][2] = rlen*tmp2;

            i       = bla[2*b];
            j       = bla[2*b+1];
            mvb     = blc[b]*(r[b][0]*(xp[i][0] - xp[j][0]) +
                              r[b][1]*(xp[i][1] - xp[j][1]) +
                              r[b][2]*(xp[i][2] - xp[j][2]) - bllen[b]);
            rhs1[b] = mvb;
            sol[b]  = mvb;

        }


    }

    return;

}


__global__ void lincs_construct_matrix_kernel(int         ncons,
                                              const real *blmf,
                                              const rvec *r,
                                              real       *blcc,
                                              const int  *blnr,
                                              const int  *blbnb)
{

    int b = blockIdx.x*blockDim.x+threadIdx.x;


    if (b < ncons)
    {

        for (int n = blnr[b]; n < blnr[b+1]; n++)
        {

            float* a = (float*) r[b];
            float* b = (float*) r[blbnb[n]];
            blcc[n] = blmf[n]*(a[XX]*b[XX]+a[YY]*b[YY]+a[ZZ]*b[ZZ]);

        }

    }

    return;

}

__global__ void lincs_matrix_expand_kernel(int         ncons,
                                           real       *rhs1,
                                           real       *rhs2,
                                           const real *blcc,
                                           real       *sol,
                                           const int  *blnr,
                                           const int  *blbnb)
{

    int b = blockIdx.x*blockDim.x+threadIdx.x;


    if (b < ncons)
    {

        real mvb;
        int  n;

        mvb = 0;

        for (n = blnr[b]; n < blnr[b+1]; n++)
        {

            mvb = mvb + blcc[n]*rhs1[blbnb[n]];

        }
        rhs2[b] = mvb;
        sol[b]  = sol[b] + mvb;
    }


    return;

}

__global__ void lincs_calc_mlambda_kernel(int         ncons,
                                          real       *mlambda,
                                          const real *blc,
                                          const real *sol)
{

    int b = blockIdx.x*blockDim.x+threadIdx.x;


    if (b < ncons)
    {

        mlambda[b] = blc[b]*sol[b];
    }

    return;

}


__global__ void lincs_update_atoms_kernel(int         ncons,
                                          const int  *bla,
                                          real        prefac,
                                          const real *fac,
                                          rvec       *r,
                                          const real *invmass,
                                          rvec       *x)
{
    int  i, j;
    real mvb, im1, im2, tmp0, tmp1, tmp2;

    int  b = blockIdx.x*blockDim.x+threadIdx.x;

    if (invmass != nullptr)
    {

        if (b < ncons)
        {
            {
                i        = bla[2*b];
                j        = bla[2*b+1];
                mvb      = prefac*fac[b];
                im1      = invmass[i];
                im2      = invmass[j];
                tmp0     = r[b][0]*mvb;
                tmp1     = r[b][1]*mvb;
                tmp2     = r[b][2]*mvb;

                atomicAdd(&(x[i][0]), -tmp0*im1);
                atomicAdd(&(x[i][1]), -tmp1*im1);
                atomicAdd(&(x[i][2]), -tmp2*im1);
                atomicAdd(&(x[j][0]), tmp0*im2);
                atomicAdd(&(x[j][1]), tmp1*im2);
                atomicAdd(&(x[j][2]), tmp2*im2);

            }
        }
    }
    else
    {

        if (b < ncons)
        {
            i        = bla[2*b];
            j        = bla[2*b+1];
            mvb      = prefac*fac[b];
            tmp0     = r[b][0]*mvb;
            tmp1     = r[b][1]*mvb;
            tmp2     = r[b][2]*mvb;

            atomicAdd(&(x[i][0]), -tmp0);
            atomicAdd(&(x[i][1]), -tmp1);
            atomicAdd(&(x[i][2]), -tmp2);
            atomicAdd(&(x[j][0]), tmp0);
            atomicAdd(&(x[j][1]), tmp1);
            atomicAdd(&(x[j][2]), tmp2);


        }
    }
}



__global__ void lincs_scale_mlambda_kernel(int ncons, real* mlambda, int* nlocat)
{

    int b = blockIdx.x*blockDim.x+threadIdx.x;

    if (b < ncons)
    {
        mlambda[b] *= 0.5*nlocat[b];
    }

}


__global__ void lincs_calc_dist_iter_kernel(int          ncons,
                                            const int   *bla,
                                            const rvec  *xp,
                                            const real  *bllen,
                                            const real  *blc,
                                            const t_pbc *pbc,
                                            real        *rhs,
                                            real        *sol,
                                            rvec       * pbc_hbox_diag,
                                            matrix     * pbc_box,
                                            rvec       * pbc_mhbox_diag,
                                            rvec      *  pbc_fbox_diag,
                                            rvec      *  pbc_tric_vec,
                                            ivec      *  pbc_tric_shift,
                                            bool         bpbc)
{

    int b = blockIdx.x*blockDim.x+threadIdx.x;

    if (b < ncons)
    {
        real len, len2, dlen2, mvb;
        rvec dx;

        len = bllen[b];
        if (bpbc)
        {
            pbc_dx_aiuc_gpu(pbc, xp[bla[2*b]], xp[bla[2*b+1]], dx,
                            *pbc_hbox_diag, *pbc_box, *pbc_mhbox_diag,
                            *pbc_fbox_diag, pbc_tric_vec, pbc_tric_shift);
        }
        else
        {
            rvec_sub_gpu(xp[bla[2*b]], xp[bla[2*b+1]], dx);
        }
        len2  = len*len;
        dlen2 = 2*len2 - norm2_gpu(dx);

        if (dlen2 > 0)
        {
            mvb = blc[b]*(len - dlen2*invsqrt(dlen2));
        }
        else
        {
            mvb = blc[b]*len;
        }
        rhs[b]  = mvb;
        sol[b]  = mvb;
    }
}


__global__ void lincs_update_blc_sol_mlambda_kernel(int         ncons,
                                                    real       *mlambda,
                                                    real       *blc_sol,
                                                    const real *blc,
                                                    const real *sol)
{

    int b = blockIdx.x*blockDim.x+threadIdx.x;

    if (b < ncons)
    {
        real mvb;

        mvb         = blc[b]*sol[b];
        blc_sol[b]  = mvb;
        mlambda[b] += mvb;
    }
}


__global__ void lincs_update_virial_kernel(int         ncons,
                                           real      * vir_r_m_dr,
                                           const real *mlambda,
                                           const real *bllen,
                                           const rvec *r)
{


    int i, j;

    int b = blockIdx.x*blockDim.x+threadIdx.x;

    /* Constraint virial */
    if (b < ncons)
    {
        real tmp0, tmp1;

        tmp0 = -bllen[b]*mlambda[b];
        for (i = 0; i < DIM; i++)
        {
            tmp1 = tmp0*r[b][i];
            for (j = 0; j < DIM; j++)
            {
                atomicAdd(&(vir_r_m_dr[i*DIM+j]), -tmp1*r[b][j]);
            }
        }
    }
    return;
}

// Update Kernel - see original CPU version in update.cpp

__global__ void updateMDLeapfrogSimple_kernel(int                 nrend,
                                              real                dt,
                                              real                dtPressureCouple,
                                              const rvec         *invMassPerDim,
                                              const real          lambdaGroup,
                                              const real         *lambdaGroupArray,
                                              const rvec         *x,
                                              rvec               *xprime,
                                              rvec               *v,
                                              const rvec         *f,
                                              NumTempScaleValues  numTempScaleValues)
{



    int a = blockIdx.x*blockDim.x+threadIdx.x;

    if (a < nrend)
    {

        for (int d = 0; d < DIM; d++)
        {
            real lambda;


            if (numTempScaleValues == NumTempScaleValues::single)
            {
                lambda = lambdaGroup;
            }

            if (numTempScaleValues == NumTempScaleValues::multiple)
            {
                lambda = lambdaGroupArray[a];
            }

            real vNew = lambda*v[a][d] + f[a][d]*invMassPerDim[a][d]*dt;

            v[a][d]      = vNew;
            xprime[a][d] = x[a][d] + vNew*dt;
        }
    }



}

// Settle Constraints Kernel - see original CPU version in settle.cpp
__global__ void settle_kernel(const int           nsettle,
                              const int         * settled_ow1,
                              const int         * settled_hw2,
                              const int         * settled_hw3,
                              const real        * settled_virfac,
                              const real          wh,
                              const real          rc,
                              const real          ra,
                              const real          rb,
                              const real          irc2,
                              const real          mO,
                              const real          mH,
                              const t_pbc       * pbc,
                              const real         *x,
                              real               *xprime,
                              real                invdt,
                              real * gmx_restrict v,
                              tensor              vir_r_m_dr,
                              bool               *bErrorHasOccurred,
                              bool                bCorrectVelocity,
                              bool                bCalcVirial,
                              rvec              * pbc_hbox_diag,
                              matrix            * pbc_box,
                              rvec              * pbc_mhbox_diag,
                              rvec             *  pbc_fbox_diag,
                              rvec             *  pbc_tric_vec,
                              ivec             *  pbc_tric_shift)
{
    /* ******************************************************************* */
    /*                                                                  ** */
    /*    Original code by Shuichi Miyamoto, last update Oct. 1, 1992   ** */
    /*                                                                  ** */
    /*    Algorithm changes by Berk Hess:                               ** */
    /*    2004-07-15 Convert COM to double precision to avoid drift     ** */
    /*    2006-10-16 Changed velocity update to use differences         ** */
    /*    2012-09-24 Use oxygen as reference instead of COM             ** */
    /*    2016-02    Complete rewrite of the code for SIMD              ** */
    /*                                                                  ** */
    /*    Reference for the SETTLE algorithm                            ** */
    /*           S. Miyamoto et al., J. Comp. Chem., 13, 952 (1992).    ** */
    /*                                                                  ** */
    /* ******************************************************************* */

    bool              bError      = bool(false);
    real              almost_zero = real(1e-12);

    int               i = blockIdx.x*blockDim.x+threadIdx.x;

    if (i < nsettle)
    {

        const int    *ow1 = settled_ow1 + i;
        const int    *hw2 = settled_hw2 + i;
        const int    *hw3 = settled_hw3 + i;

        real          x_ow1[DIM], x_hw2[DIM], x_hw3[DIM];
        x_ow1[XX] = x[3*(*ow1)];
        x_ow1[YY] = x[3*(*ow1)+1];
        x_ow1[ZZ] = x[3*(*ow1)+2];
        x_hw2[XX] = x[3*(*hw2)];
        x_hw2[YY] = x[3*(*hw2)+1];
        x_hw2[ZZ] = x[3*(*hw2)+2];
        x_hw3[XX] = x[3*(*hw3)];
        x_hw3[YY] = x[3*(*hw3)+1];
        x_hw3[ZZ] = x[3*(*hw3)+2];



        real xprime_ow1[DIM], xprime_hw2[DIM], xprime_hw3[DIM];
        xprime_ow1[XX] = xprime[3*(*ow1)];
        xprime_ow1[YY] = xprime[3*(*ow1)+1];
        xprime_ow1[ZZ] = xprime[3*(*ow1)+2];
        xprime_hw2[XX] = xprime[3*(*hw2)];
        xprime_hw2[YY] = xprime[3*(*hw2)+1];
        xprime_hw2[ZZ] = xprime[3*(*hw2)+2];
        xprime_hw3[XX] = xprime[3*(*hw3)];
        xprime_hw3[YY] = xprime[3*(*hw3)+1];
        xprime_hw3[ZZ] = xprime[3*(*hw3)+2];

        real dist21[DIM], dist31[DIM];
        real doh2[DIM], doh3[DIM];
        real sh_hw2[DIM], sh_hw3[DIM];

        pbc_dx_aiuc_gpu(pbc, x_hw2, x_ow1, dist21, *pbc_hbox_diag, *pbc_box,
                        *pbc_mhbox_diag, *pbc_fbox_diag, pbc_tric_vec, pbc_tric_shift);

        pbc_dx_aiuc_gpu(pbc, x_hw3, x_ow1, dist31, *pbc_hbox_diag, *pbc_box,
                        *pbc_mhbox_diag, *pbc_fbox_diag, pbc_tric_vec, pbc_tric_shift);


        /* Tedious way of doing pbc */
        pbc_dx_aiuc_gpu(pbc, xprime_hw2, xprime_ow1, doh2, *pbc_hbox_diag, *pbc_box,
                        *pbc_mhbox_diag, *pbc_fbox_diag, pbc_tric_vec, pbc_tric_shift);
        for (int d = 0; d < DIM; d++)
        {
            sh_hw2[d]     = xprime_hw2[d] - (xprime_ow1[d] + doh2[d]);
            xprime_hw2[d] = xprime_hw2[d] - sh_hw2[d];
        }
        pbc_dx_aiuc_gpu(pbc, xprime_hw3, xprime_ow1, doh3, *pbc_hbox_diag, *pbc_box,
                        *pbc_mhbox_diag, *pbc_fbox_diag, pbc_tric_vec, pbc_tric_shift);
        for (int d = 0; d < DIM; d++)
        {
            sh_hw3[d]     = xprime_hw3[d] - (xprime_ow1[d] + doh3[d]);
            xprime_hw3[d] = xprime_hw3[d] - sh_hw3[d];
        }

        /* Not calculating the center of mass using the oxygen position
         * and the O-H distances, as done below, will make SETTLE
         * the largest source of energy drift for simulations of water,
         * as then the oxygen coordinate is multiplied by 0.89 at every step,
         * which can then transfer a systematic rounding to the oxygen velocity.
         */

        real a1[DIM], com[DIM];
        for (int d = 0; d < DIM; d++)
        {
            a1[d]  = -(doh2[d] + doh3[d]) * wh;
            com[d] = xprime_ow1[d] - a1[d];
        }
        real b1[DIM];
        for (int d = 0; d < DIM; d++)
        {
            b1[d]  = xprime_hw2[d] - com[d];
        }
        real c1[DIM];
        for (int d = 0; d < DIM; d++)
        {
            c1[d]  = xprime_hw3[d] - com[d];
        }

        real xakszd = dist21[YY] * dist31[ZZ] - dist21[ZZ] * dist31[YY];
        real yakszd = dist21[ZZ] * dist31[XX] - dist21[XX] * dist31[ZZ];
        real zakszd = dist21[XX] * dist31[YY] - dist21[YY] * dist31[XX];
        real xaksxd = a1[YY] * zakszd - a1[ZZ] * yakszd;
        real yaksxd = a1[ZZ] * xakszd - a1[XX] * zakszd;
        real zaksxd = a1[XX] * yakszd - a1[YY] * xakszd;
        real xaksyd = yakszd * zaksxd - zakszd * yaksxd;
        real yaksyd = zakszd * xaksxd - xakszd * zaksxd;
        real zaksyd = xakszd * yaksxd - yakszd * xaksxd;

        real axlng = invsqrt(xaksxd * xaksxd + yaksxd * yaksxd + zaksxd * zaksxd);
        real aylng = invsqrt(xaksyd * xaksyd + yaksyd * yaksyd + zaksyd * zaksyd);
        real azlng = invsqrt(xakszd * xakszd + yakszd * yakszd + zakszd * zakszd);

        real trns1[DIM], trns2[DIM], trns3[DIM];

        trns1[XX] = xaksxd * axlng;
        trns2[XX] = yaksxd * axlng;
        trns3[XX] = zaksxd * axlng;
        trns1[YY] = xaksyd * aylng;
        trns2[YY] = yaksyd * aylng;
        trns3[YY] = zaksyd * aylng;
        trns1[ZZ] = xakszd * azlng;
        trns2[ZZ] = yakszd * azlng;
        trns3[ZZ] = zakszd * azlng;


        real b0d[2], c0d[2];

        for (int d = 0; d < 2; d++)
        {
            b0d[d] = trns1[d] * dist21[XX] + trns2[d] * dist21[YY] + trns3[d] * dist21[ZZ];
            c0d[d] = trns1[d] * dist31[XX] + trns2[d] * dist31[YY] + trns3[d] * dist31[ZZ];
        }

        real a1d_z, b1d[DIM], c1d[DIM];

        a1d_z = trns1[ZZ] * a1[XX] + trns2[ZZ] * a1[YY] + trns3[ZZ] * a1[ZZ];
        for (int d = 0; d < DIM; d++)
        {
            b1d[d] = trns1[d] * b1[XX] + trns2[d] * b1[YY] + trns3[d] * b1[ZZ];
            c1d[d] = trns1[d] * c1[XX] + trns2[d] * c1[YY] + trns3[d] * c1[ZZ];
        }

        real tmp, tmp2;

        real sinphi = a1d_z * invsqrt(ra*ra);
        tmp2     = 1.0 - sinphi * sinphi;

        /* If tmp2 gets close to or beyond zero we have severly distorted
         * water molecules and we should terminate the simulation.
         * Below we take the max with almost_zero to continue the loop.
         */
        bError   = bError || (tmp2 <= almost_zero);

        if (almost_zero > tmp2)
        {
            tmp2 = almost_zero;
        }

        tmp      = invsqrt(tmp2);
        real cosphi = tmp2*tmp;
        real sinpsi = (b1d[ZZ] - c1d[ZZ]) * irc2 * tmp;
        tmp2     = 1.0 - sinpsi * sinpsi;

        real cospsi = tmp2*invsqrt(tmp2);

        real a2d_y  =  ra * cosphi;
        real b2d_x  = -rc * cospsi;
        real t1     = -rb * cosphi;
        real t2     =  rc * sinpsi * sinphi;
        real b2d_y  =  t1 - t2;
        real c2d_y  =  t1 + t2;

        /*     --- Step3  al,be,ga            --- */
        real alpha  = b2d_x * (b0d[XX] - c0d[XX]) + b0d[YY] * b2d_y + c0d[YY] * c2d_y;
        real beta   = b2d_x * (c0d[YY] - b0d[YY]) + b0d[XX] * b2d_y + c0d[XX] * c2d_y;
        real gamma  = b0d[XX] * b1d[YY] - b1d[XX] * b0d[YY] + c0d[XX] * c1d[YY] - c1d[XX] * c0d[YY];
        real al2be2 = alpha * alpha + beta * beta;
        tmp2     = (al2be2 - gamma * gamma);
        real sinthe = (alpha * gamma - beta * tmp2*invsqrt(tmp2)) * invsqrt(al2be2*al2be2);

        /*  --- Step4  A3' --- */
        tmp2     = 1.0 - sinthe * sinthe;
        real costhe = tmp2*invsqrt(tmp2);

        real a3d[DIM], b3d[DIM], c3d[DIM];


        a3d[XX]  = -a2d_y * sinthe;
        a3d[YY]  = a2d_y * costhe;
        a3d[ZZ]  = a1d_z;
        b3d[XX]  = b2d_x * costhe - b2d_y * sinthe;
        b3d[YY]  = b2d_x * sinthe + b2d_y * costhe;
        b3d[ZZ]  = b1d[ZZ];
        c3d[XX]  = -b2d_x * costhe - c2d_y * sinthe;
        c3d[YY]  = -b2d_x * sinthe + c2d_y * costhe;
        c3d[ZZ]  = c1d[ZZ];

        /*    --- Step5  A3 --- */
        real a3[DIM], b3[DIM], c3[DIM];

        a3[XX] = trns1[XX]*a3d[XX] + trns1[YY]*a3d[YY] + trns1[ZZ]*a3d[ZZ];
        a3[YY] = trns2[XX]*a3d[XX] + trns2[YY]*a3d[YY] + trns2[ZZ]*a3d[ZZ];
        a3[ZZ] = trns3[XX]*a3d[XX] + trns3[YY]*a3d[YY] + trns3[ZZ]*a3d[ZZ];
        b3[XX] = trns1[XX]*b3d[XX] + trns1[YY]*b3d[YY] + trns1[ZZ]*b3d[ZZ];
        b3[YY] = trns2[XX]*b3d[XX] + trns2[YY]*b3d[YY] + trns2[ZZ]*b3d[ZZ];
        b3[ZZ] = trns3[XX]*b3d[XX] + trns3[YY]*b3d[YY] + trns3[ZZ]*b3d[ZZ];
        c3[XX] = trns1[XX]*c3d[XX] + trns1[YY]*c3d[YY] + trns1[ZZ]*c3d[ZZ];
        c3[YY] = trns2[XX]*c3d[XX] + trns2[YY]*c3d[YY] + trns2[ZZ]*c3d[ZZ];
        c3[ZZ] = trns3[XX]*c3d[XX] + trns3[YY]*c3d[YY] + trns3[ZZ]*c3d[ZZ];


        /* Compute and store the corrected new coordinate */
        for (int d = 0; d < DIM; d++)
        {
            xprime_ow1[d] = com[d] + a3[d];
        }
        for (int d = 0; d < DIM; d++)
        {
            xprime_hw2[d] = com[d] + b3[d] + sh_hw2[d];;
        }
        for (int d = 0; d < DIM; d++)
        {
            xprime_hw3[d] = com[d] + c3[d] + sh_hw3[d];
        }

        xprime[3*(*ow1)]   = xprime_ow1[XX];
        xprime[3*(*ow1)+1] = xprime_ow1[YY];
        xprime[3*(*ow1)+2] = xprime_ow1[ZZ];
        xprime[3*(*hw2)]   = xprime_hw2[XX];
        xprime[3*(*hw2)+1] = xprime_hw2[YY];
        xprime[3*(*hw2)+2] = xprime_hw2[ZZ];
        xprime[3*(*hw3)]   = xprime_hw3[XX];
        xprime[3*(*hw3)+1] = xprime_hw3[YY];
        xprime[3*(*hw3)+2] = xprime_hw3[ZZ];

        // cppcheck-suppress duplicateExpression
        if (bCorrectVelocity || bCalcVirial)
        {

            real da[DIM], db[DIM], dc[DIM];
            for (int d = 0; d < DIM; d++)
            {
                da[d] = a3[d] - a1[d];
            }
            for (int d = 0; d < DIM; d++)
            {
                db[d] = b3[d] - b1[d];
            }
            for (int d = 0; d < DIM; d++)
            {
                dc[d] = c3[d] - c1[d];
            }


            if (bCorrectVelocity)
            {
                real v_ow1[DIM], v_hw2[DIM], v_hw3[DIM];

                v_ow1[XX] = v[3*(*ow1)];
                v_ow1[YY] = v[3*(*ow1)+1];
                v_ow1[ZZ] = v[3*(*ow1)+2];
                v_hw2[XX] = v[3*(*hw2)];
                v_hw2[YY] = v[3*(*hw2)+1];
                v_hw2[ZZ] = v[3*(*hw2)+2];
                v_hw3[XX] = v[3*(*hw3)];
                v_hw3[YY] = v[3*(*hw3)+1];
                v_hw3[ZZ] = v[3*(*hw3)+2];

                /* Add the position correction divided by dt to the velocity */
                for (int d = 0; d < DIM; d++)
                {
                    v_ow1[d] = da[d]*invdt + v_ow1[d];
                }

                for (int d = 0; d < DIM; d++)
                {
                    v_hw2[d] = db[d]*invdt + v_hw2[d];

                }

                for (int d = 0; d < DIM; d++)
                {
                    v_hw3[d] = dc[d]*invdt + v_hw3[d];
                }

                v[3*(*ow1)]   = v_ow1[XX];
                v[3*(*ow1)+1] = v_ow1[YY];
                v[3*(*ow1)+2] = v_ow1[ZZ];
                v[3*(*hw2)]   = v_hw2[XX];
                v[3*(*hw2)+1] = v_hw2[YY];
                v[3*(*hw2)+2] = v_hw2[ZZ];
                v[3*(*hw3)]   = v_hw3[XX];
                v[3*(*hw3)+1] = v_hw3[YY];
                v[3*(*hw3)+2] = v_hw3[ZZ];

            }


            if (bCalcVirial)
            {
                /* Filter out the non-local settles */
                real filter = settled_virfac[i];
                real mOf    = filter*mO;
                real mHf    = filter*mH;

                real mdo[DIM], mdb[DIM], mdc[DIM];

                for (int d = 0; d < DIM; d++)
                {
                    mdb[d] = mHf*db[d];
                    mdc[d] = mHf*dc[d];
                    mdo[d] = mOf*da[d] + mdb[d] + mdc[d];
                }

                for (int d2 = 0; d2 < DIM; d2++)
                {
                    for (int d = 0; d < DIM; d++)
                    {

                        atomicAdd(&vir_r_m_dr[d2][d], -(x_ow1[d2]*mdo[d] +
                                                        dist21[d2]*mdb[d] +
                                                        dist31[d2]*mdc[d]));
                    }
                }
            }
        }
    }

    *bErrorHasOccurred = bError;

    return;
}



/*-------------------------------- End CUDA kernels-----------------------------*/

//#define USEGRAPH

static bool            graphCreated = false;
static bool            graphStep    = false;
#ifdef USEGRAPH
static cudaGraph_t     graph;
static cudaGraphExec_t instance;
#endif




/* externally visible function to perform update on GPU */
void
updateMDLeapfrogSimple_gpu(int                       nrend,
                           real                      dt,
                           real                      dtPressureCouple,
                           const rvec * gmx_restrict invMassPerDim,
                           const t_grp_tcstat      * tcstat,
                           const unsigned short    * cTC,
                           const rvec * gmx_restrict x,
                           rvec       * gmx_restrict xprime,
                           rvec       * gmx_restrict v,
                           const rvec * gmx_restrict f,
                           NumTempScaleValues        numTempScaleValues)
{



    gpuUpdateConstraintsData* gpuUCD = &gpuUCDmod;
    cudaStream_t              stream = gpuUCD->stream;

    real                      lambdaGroup = tcstat[0].lambda;

    const int                 threadsPerBlock = 128;

    rvec                    * x_d  = gpuUCD->gpu_nbv->xrvec;
    rvec                    * xp_d = gpuUCD->gpu_nbv->xprvec;
    rvec                    * f_d  = gpuUCD->gpu_nbv->frvec;
    rvec                    * v_d  = gpuUCD->gpu_nbv->vrvec;

    if (gpuUCD->bNS)
    {


        if (gpuUCD->lambdaGroupArray)
        {
            free(gpuUCD->lambdaGroupArray);
        }
        gpuUCD->lambdaGroupArray   = (real*) malloc(gpuUCD->xsize*sizeof(real));


        if (gpuUCD->lambdaGroupArray_d)
        {
            cudaFree(gpuUCD->lambdaGroupArray_d);
        }
        cudaMalloc(&gpuUCD->lambdaGroupArray_d, gpuUCD->xsize*sizeof(real));


        if (gpuUCD->invMassPerDim_d)
        {
            cudaFree(gpuUCD->invMassPerDim_d);
        }


        cudaMalloc(&gpuUCD->invMassPerDim_d, gpuUCD->xsize*sizeof(rvec));

        cudaMemcpy(gpuUCD->invMassPerDim_d, invMassPerDim, nrend*sizeof(rvec), cudaMemcpyHostToDevice);

        cudaCheckError();

        if (!gpuUCD->stream)
        {
            cudaStreamCreate(&gpuUCD->stream);
        }

        cudaMemcpy(v_d, v, nrend*sizeof(rvec), cudaMemcpyHostToDevice);

        cudaCheckError();

        if (!gpuUCD->stream)
        {
            cudaStreamCreate(&gpuUCD->stream);
        }

    }

    if (gpuUCD->bNS)
    {
        cudaMemcpy(x_d, x, gpuUCD->xsize*sizeof(rvec), cudaMemcpyHostToDevice);

        cudaCheckError();
        if (v != nullptr)
        {
            cudaMemcpy(v_d, v, gpuUCD->xsize*sizeof(rvec), cudaMemcpyHostToDevice);
        }
        cudaCheckError();
        cudaMemcpy(f_d, f, nrend*sizeof(rvec), cudaMemcpyHostToDevice);
    }

    if (numTempScaleValues == NumTempScaleValues::multiple)
    {
        for (int a = 0; a < nrend; a++)
        {
            gpuUCD->lambdaGroupArray[a] = tcstat[cTC[a]].lambda;
        }

        cudaMemcpy(gpuUCD->lambdaGroupArray_d, gpuUCD->lambdaGroupArray, nrend*sizeof(real), cudaMemcpyHostToDevice);

    }

    dim3 blocks((nrend+threadsPerBlock-1)/threadsPerBlock, 1, 1);
    dim3 threads(threadsPerBlock, 1, 1);

    updateMDLeapfrogSimple_kernel
    <<< blocks, threads, 0, stream>>>
    ( nrend,
      dt,
      dtPressureCouple,
      gpuUCD->invMassPerDim_d,
      lambdaGroup,
      gpuUCD->lambdaGroupArray_d,
      x_d,
      xp_d,
      v_d,
      f_d,
      numTempScaleValues);

    if (gpuUCD->bNS)
    {
        cudaMemcpy((void*) xprime, xp_d, gpuUCD->xsize*sizeof(rvec), cudaMemcpyDeviceToHost);
        if (v != nullptr)
        {
            cudaMemcpy((void*) v, v_d, gpuUCD->xsize*sizeof(rvec), cudaMemcpyDeviceToHost);
        }
    }
    cudaCheckError();

}




/* externally visible function to perform lincs on GPU */
void do_lincs_gpu(void            *lincsd_ptr,
                  t_pbc           *pbc,
                  rvec            *x,
                  rvec            *xp,
                  rvec           * v,
                  tensor           vir_r_m_dr,
                  real            *invmass,
                  const real       invdt,
                  gmx_bool         bCalcVir,
                  const t_commrec *cr,
                  matrix           box)
{


    Lincs   *lincsd = (Lincs*) lincsd_ptr;
    int     *bla, *blnr, *blbnb;
    rvec    *r;
    real    *blc, *blmf, *bllen;
    int     *nlocat;
    int      nc_alloc, ncc_alloc, nrec, ncons;

    bla        = lincsd->bla;
    r          = lincsd->tmpv;
    blnr       = lincsd->blnr;
    blbnb      = lincsd->blbnb;
    blc        = lincsd->blc;
    blmf       = lincsd->blmf;
    bllen      = lincsd->bllen;
    nlocat     = lincsd->nlocat;
    nc_alloc   = lincsd->nc_alloc;
    ncc_alloc  = lincsd->ncc_alloc;
    nrec       = lincsd->nOrder;
    ncons      = lincsd->nc_real;

    gpuUpdateConstraintsData* gpuUCD = &gpuUCDmod;
    cudaStream_t              stream = gpuUCD->stream;


    rvec                    * x_d  = gpuUCD->gpu_nbv->xrvec;
    rvec                    * xp_d = gpuUCD->gpu_nbv->xprvec;
    rvec                    * v_d  = gpuUCD->gpu_nbv->vrvec;


    bool move_x_setup_this_step = false;

    bool bpbc = (pbc != nullptr);

    graphStep = true;
    if (gpuUCD->bNSNextStep || gpuUCD->bNS || gpuUCD->copybackVelocity)
    {
        graphStep = false;
    }

    if (gpuUCD->bNS)
    {
        graphCreated = false;
    }

#ifdef USEGRAPH
    if (!graphCreated && graphStep)
    {
        cudaCheckError();
        cudaGraphCreate(&graph, 0); cudaCheckError();
        cudaStreamBeginCapture(stream); cudaCheckError();
    }
#endif

    if (gpuUCD->bNS)
    {

        if (gpuUCD->vir_r_m_dr_d)
        {
            cudaFree(gpuUCD->vir_r_m_dr_d);
        }
        if (gpuUCD->invmass_d)
        {
            cudaFree(gpuUCD->invmass_d);
        }
        if (gpuUCD->rhs1_d)
        {
            cudaFree(gpuUCD->rhs1_d);
        }
        if (gpuUCD->rhs2_d)
        {
            cudaFree(gpuUCD->rhs2_d);
        }
        if (gpuUCD->sol_d)
        {
            cudaFree(gpuUCD->sol_d);
        }
        if (gpuUCD->blc_sol_d)
        {
            cudaFree(gpuUCD->blc_sol_d);
        }
        if (gpuUCD->r_d)
        {
            cudaFree(gpuUCD->r_d);
        }
        if (gpuUCD->bla_d)
        {
            cudaFree(gpuUCD->bla_d);
        }
        if (gpuUCD->blc_d)
        {
            cudaFree(gpuUCD->blc_d);
        }
        if (gpuUCD->bllen_d)
        {
            cudaFree(gpuUCD->bllen_d);
        }
        if (gpuUCD->mlambda_d)
        {
            cudaFree(gpuUCD->mlambda_d);
        }
        if (gpuUCD->nlocat_d)
        {
            cudaFree(gpuUCD->nlocat_d);
        }
        if (gpuUCD->blcc_d)
        {
            cudaFree(gpuUCD->blcc_d);
        }
        if (gpuUCD->blnr_d)
        {
            cudaFree(gpuUCD->blnr_d);
        }
        if (gpuUCD->blmf_d)
        {
            cudaFree(gpuUCD->blmf_d);
        }
        if (gpuUCD->blbnb_d)
        {
            cudaFree(gpuUCD->blbnb_d);
        }
        if (gpuUCD->pbc_d)
        {
            cudaFree(gpuUCD->pbc_d);
        }
        if (gpuUCD->pbc_hbox_diag_d)
        {
            cudaFree(gpuUCD->pbc_hbox_diag_d);
        }
        if (gpuUCD->pbc_box_d)
        {
            cudaFree(gpuUCD->pbc_box_d);
        }
        if (gpuUCD->pbc_mhbox_diag_d)
        {
            cudaFree(gpuUCD->pbc_mhbox_diag_d);
        }
        if (gpuUCD->pbc_fbox_diag_d)
        {
            cudaFree(gpuUCD->pbc_fbox_diag_d);
        }
        if (gpuUCD->pbc_tric_vec_d)
        {
            cudaFree(gpuUCD->pbc_tric_vec_d);
        }
        if (gpuUCD->pbc_tric_shift_d)
        {
            cudaFree(gpuUCD->pbc_tric_shift_d);
        }


        cudaMalloc(&gpuUCD->invmass_d, gpuUCD->xsize*sizeof(real));
        cudaMalloc(&gpuUCD->rhs1_d, nc_alloc*sizeof(real));
        cudaMalloc(&gpuUCD->rhs2_d, nc_alloc*sizeof(real));
        cudaMalloc(&gpuUCD->sol_d, nc_alloc*sizeof(real));
        cudaMalloc(&gpuUCD->blc_sol_d, nc_alloc*sizeof(real));
        cudaMalloc(&gpuUCD->r_d, nc_alloc*sizeof(rvec));
        cudaMalloc(&gpuUCD->bla_d, 2*nc_alloc*sizeof(int));
        cudaMalloc(&gpuUCD->blc_d, nc_alloc*sizeof(real));
        cudaMalloc(&gpuUCD->bllen_d, nc_alloc*sizeof(real));
        cudaMalloc(&gpuUCD->mlambda_d, nc_alloc*sizeof(real));
        cudaMalloc(&gpuUCD->nlocat_d, nc_alloc*sizeof(int));
        cudaMalloc(&gpuUCD->blcc_d, ncc_alloc*sizeof(real));
        cudaMalloc(&gpuUCD->blmf_d, ncc_alloc*sizeof(real));
        cudaMalloc(&gpuUCD->blnr_d, (nc_alloc+1)*sizeof(int));
        cudaMalloc(&gpuUCD->blbnb_d, ncc_alloc*sizeof(int));
        cudaMalloc(&gpuUCD->vir_r_m_dr_d, DIM*DIM*sizeof(real));
        cudaMalloc(&gpuUCD->pbc_d, sizeof(t_pbc));
        cudaMalloc(&gpuUCD->pbc_hbox_diag_d, sizeof(rvec));
        cudaMalloc(&gpuUCD->pbc_box_d, sizeof(matrix));
        cudaMalloc(&gpuUCD->pbc_mhbox_diag_d, sizeof(rvec));
        cudaMalloc(&gpuUCD->pbc_fbox_diag_d, sizeof(rvec));
        cudaMalloc(&gpuUCD->pbc_tric_vec_d, MAX_NTRICVEC*sizeof(rvec));
        cudaMalloc(&gpuUCD->pbc_tric_shift_d, MAX_NTRICVEC*sizeof(ivec));

        if (bpbc)
        {
            cudaMemcpy(gpuUCD->pbc_hbox_diag_d, pbc->hbox_diag, sizeof(rvec), cudaMemcpyHostToDevice);
            cudaMemcpy(gpuUCD->pbc_box_d, pbc->box, sizeof(matrix), cudaMemcpyHostToDevice);
            cudaMemcpy(gpuUCD->pbc_mhbox_diag_d, pbc->mhbox_diag, sizeof(rvec), cudaMemcpyHostToDevice);
            cudaMemcpy(gpuUCD->pbc_fbox_diag_d, pbc->fbox_diag, sizeof(rvec), cudaMemcpyHostToDevice);
            cudaMemcpy(gpuUCD->pbc_tric_vec_d, pbc->tric_vec, MAX_NTRICVEC*sizeof(rvec), cudaMemcpyHostToDevice);
            cudaMemcpy(gpuUCD->pbc_tric_shift_d, pbc->tric_shift, MAX_NTRICVEC*sizeof(ivec), cudaMemcpyHostToDevice);
            cudaMemcpy(gpuUCD->pbc_d, pbc, sizeof(t_pbc), cudaMemcpyHostToDevice);
        }

        cudaMemcpy(gpuUCD->bla_d, bla, 2*nc_alloc*sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(gpuUCD->blc_d, blc, nc_alloc*sizeof(real), cudaMemcpyHostToDevice);
        cudaMemcpy(gpuUCD->bllen_d, bllen, nc_alloc*sizeof(real), cudaMemcpyHostToDevice);

        cudaMemcpy(gpuUCD->blnr_d, blnr, (nc_alloc+1)*sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(gpuUCD->blbnb_d, blbnb, ncc_alloc*sizeof(int), cudaMemcpyHostToDevice);

        cudaMemcpy(gpuUCD->blmf_d, blmf, ncc_alloc*sizeof(real), cudaMemcpyHostToDevice);

        cudaMemcpy(gpuUCD->r_d, r, nc_alloc*sizeof(rvec), cudaMemcpyHostToDevice);

        if (invmass != nullptr)
        {
            cudaMemcpy(gpuUCD->invmass_d, invmass, gpuUCD->xsize*sizeof(real), cudaMemcpyHostToDevice);
        }

        if (nlocat != nullptr)
        {
            cudaMemcpy(gpuUCD->nlocat_d, nlocat, nc_alloc*sizeof(int), cudaMemcpyHostToDevice);
        }

        cudaCheckError();


        move_x_setup_this_step = true;

    }




    cudaCheckError();

    if (gpuUCD->bNS)  //else these arrays already are present on GPU
    {
        cudaMemcpy(x_d, x, gpuUCD->xsize*sizeof(rvec), cudaMemcpyHostToDevice);
        cudaMemcpy(xp_d, xp, gpuUCD->xsize*sizeof(rvec), cudaMemcpyHostToDevice);
        if (v != nullptr)
        {
            cudaMemcpy(v_d, v, gpuUCD->xsize*sizeof(rvec), cudaMemcpyHostToDevice);
        }
    }


    const int    threadsPerBlock = 128;


    dim3 blocks((ncons+threadsPerBlock-1)/threadsPerBlock, 1, 1);
    dim3 threads(threadsPerBlock, 1, 1);

    cudaCheckError();


    //An update to GPU PBC data is required after every pressure coupling step
    if (gpuUCD->bPressureStep && bpbc)
    {
        cudaMemcpy(gpuUCD->pbc_hbox_diag_d, pbc->hbox_diag, sizeof(rvec), cudaMemcpyHostToDevice);
        cudaMemcpy(gpuUCD->pbc_box_d, pbc->box, sizeof(matrix), cudaMemcpyHostToDevice);
        cudaMemcpy(gpuUCD->pbc_mhbox_diag_d, pbc->mhbox_diag, sizeof(rvec), cudaMemcpyHostToDevice);
        cudaMemcpy(gpuUCD->pbc_fbox_diag_d, pbc->fbox_diag, sizeof(rvec), cudaMemcpyHostToDevice);
        cudaMemcpy(gpuUCD->pbc_tric_vec_d, pbc->tric_vec, MAX_NTRICVEC*sizeof(rvec), cudaMemcpyHostToDevice);
        cudaMemcpy(gpuUCD->pbc_tric_shift_d, pbc->tric_shift, MAX_NTRICVEC*sizeof(ivec), cudaMemcpyHostToDevice);
        cudaMemcpy(gpuUCD->pbc_d, pbc, sizeof(t_pbc), cudaMemcpyHostToDevice);
    }



    if (cr->dd && cr->dd->constraints && lincsd->bCommIter) //Perform device to device comms
    {
        printf("NVIDIA error: wrong input file in use\n"); fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, -1);

        gpuConstraintsD2D(x_d, cr, stream, box, move_x_setup_this_step);
        gpuConstraintsD2D(xp_d, cr, stream, box, false);

        cudaDeviceSynchronize();
        MPI_Barrier(cr->dd->mpi_comm_all);

    }

    cudaCheckError();



    if (!graphCreated || !graphStep)
    {

        lincs_init_rhs_sol_kernel
        <<< blocks, threads, 0, stream>>>
        (ncons, gpuUCD->pbc_d, x_d, xp_d, gpuUCD->r_d, gpuUCD->rhs1_d, gpuUCD->sol_d,
         gpuUCD->bla_d, gpuUCD->blc_d, gpuUCD->bllen_d, gpuUCD->pbc_hbox_diag_d,
         gpuUCD->pbc_box_d, gpuUCD->pbc_mhbox_diag_d, gpuUCD->pbc_fbox_diag_d,
         gpuUCD->pbc_tric_vec_d, gpuUCD->pbc_tric_shift_d, bpbc);


        lincs_construct_matrix_kernel
        <<< blocks, threads, 0, stream>>>
        (ncons, gpuUCD->blmf_d, gpuUCD->r_d, gpuUCD->blcc_d, gpuUCD->blnr_d, gpuUCD->blbnb_d);











        for (int rec = 0; rec < nrec; rec++)
        {

            lincs_matrix_expand_kernel
            <<< blocks, threads, 0, stream>>>
            (ncons, gpuUCD->rhs1_d, gpuUCD->rhs2_d, gpuUCD->blcc_d, gpuUCD->sol_d,
             gpuUCD->blnr_d, gpuUCD->blbnb_d);

            real *swap;

            swap           = gpuUCD->rhs1_d;
            gpuUCD->rhs1_d = gpuUCD->rhs2_d;
            gpuUCD->rhs2_d = swap;
        }




        lincs_calc_mlambda_kernel
        <<< blocks, threads, 0, stream>>>
        (ncons, gpuUCD->mlambda_d, gpuUCD->blc_d, gpuUCD->sol_d);

        lincs_update_atoms_kernel
        <<< blocks, threads, 0, stream>>>
        (ncons, gpuUCD->bla_d, 1.0, gpuUCD->mlambda_d, gpuUCD->r_d, gpuUCD->invmass_d, xp_d);

        for (int iter = 0; iter < lincsd->nIter; iter++)
        {

            if ((lincsd->bCommIter && DOMAINDECOMP(cr) && cr->dd &&
                 cr->dd->constraints && cr->dd->constraint_comm))
            {


                gpuConstraintsD2D(xp_d, cr, stream, box, false);

                cudaDeviceSynchronize();
                MPI_Barrier(cr->dd->mpi_comm_all);
            }


            lincs_calc_dist_iter_kernel
            <<< blocks, threads, 0, stream>>>
            (ncons, gpuUCD->bla_d, xp_d, gpuUCD->bllen_d, gpuUCD->blc_d, gpuUCD->pbc_d,
             gpuUCD->rhs1_d, gpuUCD->sol_d, gpuUCD->pbc_hbox_diag_d, gpuUCD->pbc_box_d,
             gpuUCD->pbc_mhbox_diag_d, gpuUCD->pbc_fbox_diag_d, gpuUCD->pbc_tric_vec_d,
             gpuUCD->pbc_tric_shift_d, bpbc);


            for (int rec = 0; rec < nrec; rec++)
            {


                lincs_matrix_expand_kernel
                <<< blocks, threads, 0, stream>>>
                (ncons, gpuUCD->rhs1_d, gpuUCD->rhs2_d, gpuUCD->blcc_d, gpuUCD->sol_d,
                 gpuUCD->blnr_d, gpuUCD->blbnb_d);

                real *swap;

                swap           = gpuUCD->rhs1_d;
                gpuUCD->rhs1_d = gpuUCD->rhs2_d;
                gpuUCD->rhs2_d = swap;
            }


            lincs_update_blc_sol_mlambda_kernel
            <<< blocks, threads, 0, stream>>>
            (ncons, gpuUCD->mlambda_d, gpuUCD->blc_sol_d, gpuUCD->blc_d, gpuUCD->sol_d );

            lincs_update_atoms_kernel
            <<< blocks, threads, 0, stream>>>
            (ncons, gpuUCD->bla_d, 1.0, gpuUCD->blc_sol_d, gpuUCD->r_d, gpuUCD->invmass_d, xp_d);
        }



        if (v != nullptr)
        {

            lincs_update_atoms_kernel
            <<< blocks, threads, 0, stream>>>
            (ncons, gpuUCD->bla_d, invdt, gpuUCD->mlambda_d, gpuUCD->r_d, gpuUCD->invmass_d, v_d);

        }

    }


    if (nlocat != nullptr && bCalcVir)
    {

        lincs_scale_mlambda_kernel
        <<< blocks, threads, 0, stream>>>
        (ncons, gpuUCD->mlambda_d, gpuUCD->nlocat_d);

    }




    if (bCalcVir)
    {
        cudaMemcpy(gpuUCD->vir_r_m_dr_d, vir_r_m_dr, DIM*DIM*sizeof(real), cudaMemcpyHostToDevice);
        lincs_update_virial_kernel
        <<< blocks, threads, 0, stream>>>
        (ncons, gpuUCD->vir_r_m_dr_d, gpuUCD->mlambda_d, gpuUCD->bllen_d, gpuUCD->r_d);

        cudaMemcpy(vir_r_m_dr, gpuUCD->vir_r_m_dr_d, DIM*DIM*sizeof(real), cudaMemcpyDeviceToHost);
    }

    cudaCheckError();

    if (gpuUCD->bNSNextStep || gpuUCD->bNS)
    {
        cudaMemcpy(x, x_d, gpuUCD->xsize*sizeof(rvec), cudaMemcpyDeviceToHost);
        cudaMemcpy(xp, xp_d, gpuUCD->xsize*sizeof(rvec), cudaMemcpyDeviceToHost);
    }


    if (gpuUCD->bNSNextStep || gpuUCD->bNS || gpuUCD->copybackVelocity)
    {
        if (v != nullptr)
        {
            cudaMemcpy(v, v_d, gpuUCD->xsize*sizeof(rvec), cudaMemcpyDeviceToHost);
        }
    }


    cudaCheckError();
    return;
}

/* externally visible function to perform settle on GPU */
void settle_gpu(void *settledptr,
                const t_pbc* pbc,
                const real *x, real *xprime,
                real invdt, real * gmx_restrict v,
                tensor vir_r_m_dr,
                bool *bErrorHasOccurred, bool bCorrectVelocity, bool bCalcVirial)

{

    settledata               *settled = (settledata*) settledptr;

    gpuUpdateConstraintsData* gpuUCD = &gpuUCDmod;
    cudaStream_t              stream = gpuUCD->stream;


    const int threadsPerBlock = 128;

    if (gpuUCD->bNS)
    {

        if (gpuUCD->bErrorHasOccurred_d)
        {
            cudaFree(gpuUCD->bErrorHasOccurred_d);
        }
        if (gpuUCD->settled_ow1_d)
        {
            cudaFree(gpuUCD->settled_ow1_d);
        }
        if (gpuUCD->settled_hw2_d)
        {
            cudaFree(gpuUCD->settled_hw2_d);
        }
        if (gpuUCD->settled_hw3_d)
        {
            cudaFree(gpuUCD->settled_hw3_d);
        }
        if (gpuUCD->settled_virfac_d)
        {
            cudaFree(gpuUCD->settled_virfac_d);
        }

        cudaMalloc(&gpuUCD->settled_ow1_d, settled->nsettle*sizeof(int));
        cudaMalloc(&gpuUCD->settled_hw2_d, settled->nsettle*sizeof(int));
        cudaMalloc(&gpuUCD->settled_hw3_d, settled->nsettle*sizeof(int));
        cudaMalloc(&gpuUCD->settled_virfac_d, settled->nsettle*sizeof(real));
        cudaMalloc(&gpuUCD->bErrorHasOccurred_d, sizeof(bool));


        if (gpuUCD->pbc_d)
        {
            cudaFree(gpuUCD->pbc_d);
        }
        if (gpuUCD->pbc_hbox_diag_d)
        {
            cudaFree(gpuUCD->pbc_hbox_diag_d);
        }
        if (gpuUCD->pbc_box_d)
        {
            cudaFree(gpuUCD->pbc_box_d);
        }
        if (gpuUCD->pbc_mhbox_diag_d)
        {
            cudaFree(gpuUCD->pbc_mhbox_diag_d);
        }
        if (gpuUCD->pbc_fbox_diag_d)
        {
            cudaFree(gpuUCD->pbc_fbox_diag_d);
        }
        if (gpuUCD->pbc_tric_vec_d)
        {
            cudaFree(gpuUCD->pbc_tric_vec_d);
        }
        if (gpuUCD->pbc_tric_shift_d)
        {
            cudaFree(gpuUCD->pbc_tric_shift_d);
        }

        cudaMalloc(&gpuUCD->pbc_d, sizeof(t_pbc));
        cudaMalloc(&gpuUCD->pbc_hbox_diag_d, sizeof(rvec));
        cudaMalloc(&gpuUCD->pbc_box_d, sizeof(matrix));
        cudaMalloc(&gpuUCD->pbc_mhbox_diag_d, sizeof(rvec));
        cudaMalloc(&gpuUCD->pbc_fbox_diag_d, sizeof(rvec));
        cudaMalloc(&gpuUCD->pbc_tric_vec_d, MAX_NTRICVEC*sizeof(rvec));
        cudaMalloc(&gpuUCD->pbc_tric_shift_d, MAX_NTRICVEC*sizeof(ivec));

        if (gpuUCD->vir_r_m_dr_d)
        {
            cudaFree(gpuUCD->vir_r_m_dr_d);
        }
        cudaMalloc(&gpuUCD->vir_r_m_dr_d, DIM*DIM*sizeof(real));

        cudaMemcpy(gpuUCD->settled_ow1_d, settled->ow1, settled->nsettle*sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(gpuUCD->settled_hw2_d, settled->hw2, settled->nsettle*sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(gpuUCD->settled_hw3_d, settled->hw3, settled->nsettle*sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(gpuUCD->settled_virfac_d, settled->virfac, settled->nsettle*sizeof(real), cudaMemcpyHostToDevice);
        cudaMemcpy(gpuUCD->pbc_hbox_diag_d, pbc->hbox_diag, sizeof(rvec), cudaMemcpyHostToDevice);
        cudaMemcpy(gpuUCD->pbc_box_d, pbc->box, sizeof(matrix), cudaMemcpyHostToDevice);
        cudaMemcpy(gpuUCD->pbc_mhbox_diag_d, pbc->mhbox_diag, sizeof(rvec), cudaMemcpyHostToDevice);
        cudaMemcpy(gpuUCD->pbc_fbox_diag_d, pbc->fbox_diag, sizeof(rvec), cudaMemcpyHostToDevice);
        cudaMemcpy(gpuUCD->pbc_tric_vec_d, pbc->tric_vec, MAX_NTRICVEC*sizeof(rvec), cudaMemcpyHostToDevice);
        cudaMemcpy(gpuUCD->pbc_tric_shift_d, pbc->tric_shift, MAX_NTRICVEC*sizeof(ivec), cudaMemcpyHostToDevice);
        cudaMemcpy(gpuUCD->pbc_d, pbc, sizeof(t_pbc), cudaMemcpyHostToDevice);

        cudaCheckError();
    }


    if (gpuUCD->bNS)
    {
        cudaMemcpy(gpuUCD->gpu_nbv->xrvec, x, gpuUCD->xsize*sizeof(rvec), cudaMemcpyHostToDevice);
        cudaMemcpy(gpuUCD->gpu_nbv->xprvec, xprime, gpuUCD->xsize*sizeof(rvec), cudaMemcpyHostToDevice);
        if (v != nullptr)
        {
            cudaMemcpy(gpuUCD->gpu_nbv->vrvec, v, gpuUCD->xsize*sizeof(rvec), cudaMemcpyHostToDevice);
        }
    }

    cudaCheckError();

    if (bCalcVirial)
    {
        cudaMemcpy(gpuUCD->vir_r_m_dr_d, vir_r_m_dr, DIM*DIM*sizeof(real), cudaMemcpyHostToDevice);
    }

    cudaCheckError();
    settleparam_t *p    = (settleparam_t*)  &settled->massw;

    dim3 blocks((settled->nsettle+threadsPerBlock-1)/threadsPerBlock, 1, 1);
    dim3 threads(threadsPerBlock, 1, 1);

    if (!graphCreated || !graphStep)
    {
        settle_kernel
        <<< blocks, threads, 0, stream>>>
        (
            settled->nsettle,
            gpuUCD->settled_ow1_d,
            gpuUCD->settled_hw2_d,
            gpuUCD->settled_hw3_d,
            gpuUCD->settled_virfac_d,
            p->wh,
            p->rc,
            p->ra,
            p->rb,
            p->irc2,
            p->mO,
            p->mH,
            gpuUCD->pbc_d,
            (real*) gpuUCD->gpu_nbv->xrvec,
            (real*) gpuUCD->gpu_nbv->xprvec,
            invdt,
            (real*) gpuUCD->gpu_nbv->vrvec,
            (real (*)[DIM])gpuUCD->vir_r_m_dr_d,
            gpuUCD->bErrorHasOccurred_d,
            bCorrectVelocity,
            bCalcVirial,
            gpuUCD->pbc_hbox_diag_d,
            gpuUCD->pbc_box_d,
            gpuUCD->pbc_mhbox_diag_d,
            gpuUCD->pbc_fbox_diag_d,
            gpuUCD->pbc_tric_vec_d, gpuUCD->
                pbc_tric_shift_d);



        cudaCheckError();

        cudaMemcpyAsync(bErrorHasOccurred, gpuUCD->bErrorHasOccurred_d, sizeof(bool), cudaMemcpyDeviceToHost, stream);
    }

    #ifdef USEGRAPH

    if (!graphCreated && graphStep)
    {
        cudaStreamEndCapture(stream, &graph); cudaCheckError();
        cudaGraphInstantiate(&instance, graph, NULL, NULL, 0); cudaCheckError();
        cudaGraphDestroy(graph); cudaCheckError();
        graphCreated = true;
    }

    if (graphCreated && graphStep)
    {
        cudaGraphLaunch(instance, stream);
    }
#endif

    cudaStreamSynchronize(stream);



    if (bCalcVirial)
    {
        cudaMemcpy(vir_r_m_dr, gpuUCD->vir_r_m_dr_d, DIM*DIM*sizeof(real), cudaMemcpyDeviceToHost);
    }


    if (gpuUCD->bNSNextStep || gpuUCD->bNS)
    {
        cudaMemcpy((void*) xprime, gpuUCD->gpu_nbv->xprvec, gpuUCD->xsize*sizeof(rvec), cudaMemcpyDeviceToHost);
    }


    if (gpuUCD->bNSNextStep || gpuUCD->bNS || gpuUCD->copybackVelocity)
    {

        if (v != nullptr)
        {
            cudaMemcpy((void*) v, gpuUCD->gpu_nbv->vrvec, gpuUCD->xsize*sizeof(rvec), cudaMemcpyDeviceToHost);
        }
    }

    return;
}

void gpuUpdateConstraintsSetSize(int size)
{

    gpuUpdateConstraintsData* gpuUCD = &gpuUCDmod;


    gpuUCD->xsize = size;
}

void gpuUpdateConstraintsSetTimestepInfo(bool bNS, bool bNSNextStep, bool copybackVelocity)
{

    gpuUpdateConstraintsData* gpuUCD = &gpuUCDmod;

    gpuUCD->bNS              = bNS;
    gpuUCD->bNSNextStep      = bNSNextStep;
    gpuUCD->copybackVelocity = copybackVelocity;

}

int gpuUpdateConstraintsCPUCopyRequired()
{

    gpuUpdateConstraintsData* gpuUCD = &gpuUCDmod;

    return (gpuUCD->bNS || gpuUCD->bNSNextStep);

}

void gpuUpdateConstraintsSetGpuNB(gmx_nbnxn_gpu_t* gpu_nbv)
{

    gpuUpdateConstraintsData* gpuUCD = &gpuUCDmod;

    gpuUCD->gpu_nbv = gpu_nbv;
}

void gpuUpdateConstraintsSetPressureCouplingStep(bool bPressureStep)
{

    gpuUpdateConstraintsData* gpuUCD = &gpuUCDmod;

    gpuUCD->bPressureStep = bPressureStep;
}

int gpuUpdateConstraintsGetSize()
{

    gpuUpdateConstraintsData* gpuUCD = &gpuUCDmod;

    return gpuUCD->xsize;
}

void gpuUpdateConstraintsCopyXPToXOnDevice()
{

    gpuUpdateConstraintsData* gpuUCD = &gpuUCDmod;

    cudaMemcpy(gpuUCD->gpu_nbv->xrvec, gpuUCD->gpu_nbv->xprvec,
               gpuUCD->xsize*sizeof(rvec), cudaMemcpyDeviceToDevice);

    return;
}

//TODO refactor this to a better place
#include "gromacs/ewald/pme.h"
#include "gromacs/ewald/pme.cuh"

void gpuUpdateConstraintsCopyXToPMEOnDevice(PmeGpu *pmeGpu)
{

    gpuUpdateConstraintsData* gpuUCD = &gpuUCDmod;

    cudaMemcpy((void*) pmeGpu->kernelParams->atoms.d_coordinates, gpuUCD->gpu_nbv->xrvec,
               pmeGpu->kernelParams->atoms.nAtoms*sizeof(rvec), cudaMemcpyDeviceToDevice);

    return;
}
