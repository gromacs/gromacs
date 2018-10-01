/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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
 * \brief Implements CUDA bonded functionality
 *
 * \author Jon Vincent <jvincent@nvidia.com>
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 * \author Berk Hess <hess@kth.se>
 * \author Szilárd Páll <pall.szilard@gmail.com>
 * \author Alan Gray <alang@nvidia.com>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \ingroup module_listed-forces
 */

#include "gmxpre.h"

#include <math_constants.h>

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gpu_vec.h"
#include "gromacs/listed-forces/listed-forces.h"
#include "gromacs/listed-forces/manage-threading.h"
#include "gromacs/math/units.h"
#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forcerec.h" // TODO try remove when integrating
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/pbcutil/gpu_pbc.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/gmxassert.h"

#if defined(_MSVC)
#include <limits>
#endif

// TODO remove after integrating
#define BO_MAX_RANKS 32
#include "gromacs/utility/gmxmpi.h"
//

//CUDA threads per block
#define TPB_BONDED 256

// packing some parameters for transfer.
// will be much more if we need free energy
// reals from we need for LJ14
enum {
    frFUDGEQQ, fr_icEPSFAC,
    fr_PAIRS_SCALE,
    PAIRS_NR
};

// ints we need for LJ14
enum {
    fr_PAIRS_STRIDE,
    md_NENERGRP,
    PAIRS_NI
};

/*-------------------------------- CUDA kernels-------------------------------- */
/*------------------------------------------------------------------------------*/


/*---------------- BONDED CUDA kernels--------------*/
__global__ void
reset_gpu_bonded_kernel(float *vtot, fvec *force,
                        const int size, float *energygrp_elec, float *energygrp_vdw,
                        const int nener, fvec *f_shift)
{
    int a = blockIdx.x * blockDim.x + threadIdx.x;

    if (a < F_NRE)
    {
        vtot[a] = 0.0f;
    }

    if (a < size)
    {
        force[a][XX] = 0.0f;
        force[a][YY] = 0.0f;
        force[a][ZZ] = 0.0f;
    }

    if (a < nener)
    {
        energygrp_vdw[a]  = 0.0f;
        energygrp_elec[a] = 0.0f;
    }

    if (a < SHIFTS)
    {
        f_shift[a][XX] = 0.0f;
        f_shift[a][YY] = 0.0f;
        f_shift[a][ZZ] = 0.0f;
    }
}


/* Harmonic */
__device__
static void harmonic_gpu(const float kA, const float xA, const float x, float *V, float *F)
{
    constexpr float half = 0.5f;
    float           dx, dx2;

    dx    = x-xA;
    dx2   = dx*dx;

    *F = -kA*dx;
    *V = half*kA*dx2;
}

template <bool calcEnerVir>
__global__
void bonds_gpu(float *vtot, const int nbonds, const int natoms,
               const t_iatom forceatoms[], const t_iparams forceparams[],
               const fvec x[], fvec force[], fvec fshift[],
               const PbcAiuc pbcAiuc)
{
    const int        i = blockIdx.x * blockDim.x + threadIdx.x;

    __shared__ float vtot_loc;
    __shared__ fvec  fshift_loc[SHIFTS];

    if (calcEnerVir)
    {
        if (threadIdx.x == 0)
        {
            vtot_loc = 0.0f;
        }
        if (threadIdx.x < SHIFTS)
        {
            fshift_loc[threadIdx.x][XX] = 0.0f;
            fshift_loc[threadIdx.x][YY] = 0.0f;
            fshift_loc[threadIdx.x][ZZ] = 0.0f;
        }
        __syncthreads();
    }

    if (i < nbonds)
    {
        int type = forceatoms[3*i];
        int ai   = forceatoms[3*i + 1];
        int aj   = forceatoms[3*i + 2];

        /* dx = xi - xj, corrected for periodic boundry conditions. */
        fvec  dx;
        int   ki = pbcDxAiuc<calcEnerVir>(pbcAiuc, x[ai], x[aj], dx);

        float dr2 = iprod_gpu(dx, dx);
        float dr  = sqrt(dr2);

        float vbond;
        float fbond;
        harmonic_gpu(forceparams[type].harmonic.krA,
                     forceparams[type].harmonic.rA,
                     dr, &vbond, &fbond);

        if (calcEnerVir)
        {
            atomicAdd(&vtot_loc, vbond);
        }

        if (dr2 != 0.0f)
        {
            fbond *= rsqrtf(dr2);

#pragma unroll
            for (int m = 0; m < DIM; m++)
            {
                float fij = fbond*dx[m];
                atomicAdd(&force[ai][m], fij);
                atomicAdd(&force[aj][m], -fij);
                if (calcEnerVir && ki != CENTRAL)
                {
                    atomicAdd(&fshift_loc[ki][m], fij);
                    atomicAdd(&fshift_loc[CENTRAL][m], -fij);
                }
            }
        }
    }

    if (calcEnerVir)
    {
        __syncthreads();
        if (threadIdx.x == 0)
        {
            atomicAdd(vtot, vtot_loc);
        }
        if (threadIdx.x < SHIFTS)
        {
            fvec_inc_atomic(fshift[threadIdx.x], fshift_loc[threadIdx.x]);
        }
    }
}

template <bool returnShift>
__device__
static float bond_angle_gpu(const fvec xi, const fvec xj, const fvec xk,
                            const PbcAiuc &pbcAiuc,
                            fvec r_ij, fvec r_kj, float *costh,
                            int *t1, int *t2)
/* Return value is the angle between the bonds i-j and j-k */
{
    *t1      = pbcDxAiuc<returnShift>(pbcAiuc, xi, xj, r_ij);
    *t2      = pbcDxAiuc<returnShift>(pbcAiuc, xk, xj, r_kj);

    *costh   = cos_angle_gpu(r_ij, r_kj);
    float th = acosf(*costh);

    return th;
}

template <bool calcEnerVir>
__global__
void angles_gpu(float *vtot, const int nbonds, const int natoms,
                const t_iatom forceatoms[], const t_iparams forceparams[],
                const fvec x[], fvec force[], fvec fshift[],
                const PbcAiuc pbcAiuc)
{
    __shared__ float vtot_loc;
    __shared__ fvec  fshift_loc[SHIFTS];

    const int        i = blockIdx.x*blockDim.x + threadIdx.x;

    if (calcEnerVir)
    {
        if (threadIdx.x == 0)
        {
            vtot_loc = 0.0f;
        }
        if (threadIdx.x < SHIFTS)
        {
            fshift_loc[threadIdx.x][XX] = 0.0f;
            fshift_loc[threadIdx.x][YY] = 0.0f;
            fshift_loc[threadIdx.x][ZZ] = 0.0f;
        }

        __syncthreads();
    }

    if (i < nbonds)
    {
        int   type = forceatoms[4*i];
        int   ai   = forceatoms[4*i + 1];
        int   aj   = forceatoms[4*i + 2];
        int   ak   = forceatoms[4*i + 3];

        fvec  r_ij;
        fvec  r_kj;
        float cos_theta;
        int   t1;
        int   t2;
        float theta =
            bond_angle_gpu<calcEnerVir>(x[ai], x[aj], x[ak], pbcAiuc,
                                        r_ij, r_kj, &cos_theta, &t1, &t2);

        float va;
        float dVdt;
        harmonic_gpu(forceparams[type].harmonic.krA,
                     forceparams[type].harmonic.rA*DEG2RAD,
                     theta, &va, &dVdt);

        if (calcEnerVir)
        {
            atomicAdd(&vtot_loc, va);
        }

        float cos_theta2 = cos_theta*cos_theta;
        if (cos_theta2 < 1.0f)
        {
            float st    = dVdt*rsqrtf(1.0f - cos_theta2);
            float sth   = st*cos_theta;
            float nrij2 = iprod_gpu(r_ij, r_ij);
            float nrkj2 = iprod_gpu(r_kj, r_kj);

            float nrij_1 = rsqrtf(nrij2);
            float nrkj_1 = rsqrtf(nrkj2);

            float cik = st*nrij_1*nrkj_1;
            float cii = sth*nrij_1*nrij_1;
            float ckk = sth*nrkj_1*nrkj_1;

            fvec  f_i;
            fvec  f_k;
            fvec  f_j;
            for (int m = 0; m < DIM; m++)
            {
                f_i[m]    = -(cik*r_kj[m] - cii*r_ij[m]);
                f_k[m]    = -(cik*r_ij[m] - ckk*r_kj[m]);
                f_j[m]    = -f_i[m] - f_k[m];
                atomicAdd(&force[ai][m], f_i[m]);
                atomicAdd(&force[aj][m], f_j[m]);
                atomicAdd(&force[ak][m], f_k[m]);
            }
            if (calcEnerVir)
            {
                fvec_inc_atomic(fshift_loc[t1], f_i);
                fvec_inc_atomic(fshift_loc[CENTRAL], f_j);
                fvec_inc_atomic(fshift_loc[t2], f_k);
            }
        }

    }

    if (calcEnerVir)
    {
        __syncthreads();

        if (threadIdx.x == 0)
        {
            atomicAdd(vtot, vtot_loc);
        }
        if (threadIdx.x < SHIFTS)
        {
            fvec_inc_atomic(fshift[threadIdx.x], fshift_loc[threadIdx.x]);
        }
    }
}

template <bool calcEnerVir>
__global__
void urey_bradley_gpu(float *vtot, const int nbonds, const int natoms,
                      const t_iatom forceatoms[], const t_iparams forceparams[],
                      const fvec x[], fvec force[], fvec fshift[],
                      const PbcAiuc pbcAiuc)
{
    __shared__ float vtot_loc;
    __shared__ fvec  fshift_loc[SHIFTS];

    const int        i = blockIdx.x*blockDim.x + threadIdx.x;

    if (calcEnerVir)
    {
        if (threadIdx.x == 0)
        {
            vtot_loc = 0.0f;
        }
        if (threadIdx.x < SHIFTS)
        {
            fshift_loc[threadIdx.x][XX] = 0.0f;
            fshift_loc[threadIdx.x][YY] = 0.0f;
            fshift_loc[threadIdx.x][ZZ] = 0.0f;
        }

        __syncthreads();
    }

    if (i < nbonds)
    {
        int   type  = forceatoms[4*i];
        int   ai    = forceatoms[4*i+1];
        int   aj    = forceatoms[4*i+2];
        int   ak    = forceatoms[4*i+3];

        float th0A  = forceparams[type].u_b.thetaA*DEG2RAD;
        float kthA  = forceparams[type].u_b.kthetaA;
        float r13A  = forceparams[type].u_b.r13A;
        float kUBA  = forceparams[type].u_b.kUBA;

        fvec  r_ij;
        fvec  r_kj;
        float cos_theta;
        int   t1;
        int   t2;
        float theta = bond_angle_gpu<calcEnerVir>(x[ai], x[aj], x[ak], pbcAiuc,
                                                  r_ij, r_kj, &cos_theta, &t1, &t2);

        float va;
        float dVdt;
        harmonic_gpu(kthA, th0A, theta, &va, &dVdt);

        if (calcEnerVir)
        {
            atomicAdd(&vtot_loc, va);
        }

        fvec  r_ik;
        int   ki = pbcDxAiuc<calcEnerVir>(pbcAiuc, x[ai], x[ak], r_ik);

        float dr2  = iprod_gpu(r_ik, r_ik);
        float dr   = dr2*rsqrtf(dr2);

        float vbond;
        float fbond;
        harmonic_gpu(kUBA, r13A, dr, &vbond, &fbond);

        float cos_theta2 = cos_theta*cos_theta;
        if (cos_theta2 < 1.0f)
        {
            float st    = dVdt*rsqrtf(1.0f - cos_theta2);
            float sth   = st*cos_theta;

            float nrkj2 = iprod_gpu(r_kj, r_kj);
            float nrij2 = iprod_gpu(r_ij, r_ij);

            float cik   = st*rsqrtf(nrkj2*nrij2);
            float cii   = sth/nrij2;
            float ckk   = sth/nrkj2;

            fvec  f_i;
            fvec  f_j;
            fvec  f_k;
            for (int m = 0; m < DIM; m++)
            {
                f_i[m]    = -(cik*r_kj[m]-cii*r_ij[m]);
                f_k[m]    = -(cik*r_ij[m]-ckk*r_kj[m]);
                f_j[m]    = -f_i[m]-f_k[m];
                atomicAdd(&force[ai][m], f_i[m]);
                atomicAdd(&force[aj][m], f_j[m]);
                atomicAdd(&force[ak][m], f_k[m]);
            }
            fvec_inc_atomic(fshift_loc[t1], f_i);
            fvec_inc_atomic(fshift_loc[CENTRAL], f_j);
            fvec_inc_atomic(fshift_loc[t2], f_k);
        }

        /* Time for the bond calculations */
        if (dr2 != 0.0f)
        {
            if (calcEnerVir)
            {
                atomicAdd(&vtot_loc, vbond);
            }

            fbond *= rsqrtf(dr2);

            for (int m = 0; m < DIM; m++)
            {
                float fik = fbond*r_ik[m];
                atomicAdd(&force[ai][m], fik);
                atomicAdd(&force[ak][m], -fik);

                if (calcEnerVir && ki != CENTRAL)
                {
                    atomicAdd(&fshift_loc[ki][m], fik);
                    atomicAdd(&fshift_loc[CENTRAL][m], -fik);
                }
            }
        }
    }

    if (calcEnerVir)
    {
        __syncthreads();

        if (threadIdx.x == 0)
        {
            atomicAdd(vtot, vtot_loc);
        }
        if (threadIdx.x < SHIFTS)
        {
            fvec_inc_atomic(fshift[threadIdx.x], fshift_loc[threadIdx.x]);
        }
    }
}

template <bool returnShift>
__device__
static float dih_angle_gpu(const fvec xi, const fvec xj, const fvec xk, const fvec xl,
                           const PbcAiuc &pbcAiuc,
                           fvec r_ij, fvec r_kj, fvec r_kl, fvec m, fvec n,
                           int *t1, int *t2, int *t3)
{
    *t1 = pbcDxAiuc<returnShift>(pbcAiuc, xi, xj, r_ij);
    *t2 = pbcDxAiuc<returnShift>(pbcAiuc, xk, xj, r_kj);
    *t3 = pbcDxAiuc<returnShift>(pbcAiuc, xk, xl, r_kl);

    cprod_gpu(r_ij, r_kj, m);
    cprod_gpu(r_kj, r_kl, n);
    float phi  = gmx_angle_gpu(m, n);
    float ipr  = iprod_gpu(r_ij, n);
    float sign = (ipr < 0.0f) ? -1.0f : 1.0f;
    phi        = sign*phi;

    return phi;
}


__device__
static void dopdihs_gpu(const float cpA, const float phiA, const int mult,
                        const float phi, float *V, float *F)
{
    float mdphi, sdphi;

    mdphi = mult*phi - phiA*DEG2RAD;
    sdphi = sinf(mdphi);
    *V    = cpA * (1.0f + cosf(mdphi));
    *F    = -cpA*mult*sdphi;
}

template <bool calcVir>
__device__
static void do_dih_fup_gpu(const int i, const int j, const int k, const int l, const int natoms,
                           const float ddphi, const fvec r_ij, const fvec r_kj, const fvec r_kl,
                           const fvec m, const fvec n, fvec force[], fvec fshift[],
                           const PbcAiuc &pbcAiuc,
                           const fvec x[], const int t1, const int t2, const int t3)
{
    float iprm  = iprod_gpu(m, m);
    float iprn  = iprod_gpu(n, n);
    float nrkj2 = iprod_gpu(r_kj, r_kj);
    float toler = nrkj2*GMX_REAL_EPS;
    if ((iprm > toler) && (iprn > toler))
    {
        float nrkj_1 = rsqrtf(nrkj2); // replacing std::invsqrt call
        float nrkj_2 = nrkj_1*nrkj_1;
        float nrkj   = nrkj2*nrkj_1;
        float a      = -ddphi*nrkj/iprm;
        fvec  f_i;
        svmul_gpu(a, m, f_i);
        float b      = ddphi*nrkj/iprn;
        fvec  f_l;
        svmul_gpu(b, n, f_l);
        float p      = iprod_gpu(r_ij, r_kj);
        p           *= nrkj_2;
        float q      = iprod_gpu(r_kl, r_kj);
        q           *= nrkj_2;
        fvec  uvec;
        svmul_gpu(p, f_i, uvec);
        fvec  vvec;
        svmul_gpu(q, f_l, vvec);
        fvec  svec;
        fvec_sub_gpu(uvec, vvec, svec);
        fvec  f_j;
        fvec_sub_gpu(f_i, svec, f_j);
        fvec  f_k;
        fvec_add_gpu(f_l, svec, f_k);
#pragma unroll
        for (int m = 0; (m < DIM); m++)
        {
            atomicAdd(&force[i][m], f_i[m]);
            atomicAdd(&force[j][m], -f_j[m]);
            atomicAdd(&force[k][m], -f_k[m]);
            atomicAdd(&force[l][m], f_l[m]);
        }

        if (calcVir)
        {
            fvec dx_jl;
            int  t3 = pbcDxAiuc<calcVir>(pbcAiuc, x[l], x[j], dx_jl);

            fvec_inc_atomic(fshift[t1], f_i);
            fvec_dec_atomic(fshift[CENTRAL], f_j);
            fvec_dec_atomic(fshift[t2], f_k);
            fvec_inc_atomic(fshift[t3], f_l);
        }
    }
}

template <bool calcEnerVir>
__global__
void  pdihs_gpu(float *vtot, const int nbonds, const int natoms,
                const t_iatom forceatoms[], const t_iparams forceparams[],
                const fvec x[], fvec f[], fvec fshift[],
                const PbcAiuc pbcAiuc)
{
    const int        i = blockIdx.x*blockDim.x + threadIdx.x;

    __shared__ float vtot_loc;
    __shared__ fvec  fshift_loc[SHIFTS];

    if (calcEnerVir)
    {
        if (threadIdx.x == 0)
        {
            vtot_loc = 0.0f;
        }
        if (threadIdx.x < SHIFTS)
        {
            fshift_loc[threadIdx.x][XX] = 0.0f;
            fshift_loc[threadIdx.x][YY] = 0.0f;
            fshift_loc[threadIdx.x][ZZ] = 0.0f;
        }
        __syncthreads();
    }

    if (i < nbonds)
    {
        int   type = forceatoms[5*i];
        int   ai   = forceatoms[5*i + 1];
        int   aj   = forceatoms[5*i + 2];
        int   ak   = forceatoms[5*i + 3];
        int   al   = forceatoms[5*i + 4];

        fvec  r_ij;
        fvec  r_kj;
        fvec  r_kl;
        fvec  m;
        fvec  n;
        int   t1;
        int   t2;
        int   t3;
        float phi  =
            dih_angle_gpu<calcEnerVir>(x[ai], x[aj], x[ak], x[al], pbcAiuc,
                                       r_ij, r_kj, r_kl, m, n, &t1, &t2, &t3);

        float vpd;
        float ddphi;
        dopdihs_gpu(forceparams[type].pdihs.cpA,
                    forceparams[type].pdihs.phiA,
                    forceparams[type].pdihs.mult,
                    phi, &vpd, &ddphi);

        atomicAdd(&vtot_loc, vpd);

        do_dih_fup_gpu<calcEnerVir>(ai, aj, ak, al, natoms,
                                    ddphi, r_ij, r_kj, r_kl, m, n,
                                    f, fshift_loc, pbcAiuc,
                                    x, t1, t2, t3);

    }

    if (calcEnerVir)
    {
        __syncthreads();

        if (threadIdx.x == 0)
        {
            atomicAdd(vtot, vtot_loc);
        }
        if (threadIdx.x < SHIFTS)
        {
            fvec_inc_atomic(fshift[threadIdx.x], fshift_loc[threadIdx.x]);
        }
    }
}

template <bool calcEnerVir>
__global__
void rbdihs_gpu(float *vtot, const int nbonds, const int natoms,
                const t_iatom forceatoms[], const t_iparams forceparams[],
                const fvec x[], fvec f[], fvec fshift[],
                const PbcAiuc pbcAiuc)
{
    constexpr float  c0 = 0.0f, c1 = 1.0f, c2 = 2.0f, c3 = 3.0f, c4 = 4.0f, c5 = 5.0f;

    __shared__ float vtot_loc;
    __shared__ fvec  fshift_loc[SHIFTS];

    const int        i = blockIdx.x*blockDim.x + threadIdx.x;

    if (calcEnerVir)
    {
        if (threadIdx.x == 0)
        {
            vtot_loc = 0.0f;
        }
        if (threadIdx.x < SHIFTS)
        {
            fshift_loc[threadIdx.x][XX] = 0.0f;
            fshift_loc[threadIdx.x][YY] = 0.0f;
            fshift_loc[threadIdx.x][ZZ] = 0.0f;
        }

        __syncthreads();
    }

    if (i < nbonds)
    {
        int   type = forceatoms[5*i];
        int   ai   = forceatoms[5*i+1];
        int   aj   = forceatoms[5*i+2];
        int   ak   = forceatoms[5*i+3];
        int   al   = forceatoms[5*i+4];

        fvec  r_ij;
        fvec  r_kj;
        fvec  r_kl;
        fvec  m;
        fvec  n;
        int   t1;
        int   t2;
        int   t3;
        float phi  =
            dih_angle_gpu<calcEnerVir>(x[ai], x[aj], x[ak], x[al], pbcAiuc,
                                       r_ij, r_kj, r_kl, m, n, &t1, &t2, &t3);

        /* Change to polymer convention */
        if (phi < c0)
        {
            phi += CUDART_PI_F;
        }
        else
        {
            phi -= CUDART_PI_F;

        }
        float cos_phi = cosf(phi);
        /* Beware of accuracy loss, cannot use 1-sqrt(cos^2) ! */
        float sin_phi = sinf(phi);

        float parm[NR_RBDIHS];
        for (int j = 0; j < NR_RBDIHS; j++)
        {
            parm[j]  = forceparams[type].rbdihs.rbcA[j];
        }
        /* Calculate cosine powers */
        /* Calculate the energy */
        /* Calculate the derivative */
        float v      = parm[0];
        float ddphi  = c0;
        float cosfac = c1;

        float rbp    = parm[1];
        ddphi       += rbp*cosfac;
        cosfac      *= cos_phi;
        rbp          = parm[2];
        ddphi       += c2*rbp*cosfac;
        cosfac      *= cos_phi;
        rbp          = parm[3];
        ddphi       += c3*rbp*cosfac;
        cosfac      *= cos_phi;
        rbp          = parm[4];
        ddphi       += c4*rbp*cosfac;
        cosfac      *= cos_phi;
        rbp          = parm[5];
        ddphi       += c5*rbp*cosfac;
        cosfac      *= cos_phi;

        ddphi = -ddphi*sin_phi;

        do_dih_fup_gpu<calcEnerVir>(ai, aj, ak, al, natoms,
                                    ddphi, r_ij, r_kj, r_kl, m, n,
                                    f, fshift_loc, pbcAiuc,
                                    x, t1, t2, t3);
        if (calcEnerVir)
        {
            atomicAdd(&vtot_loc, v);
        }
    }

    if (calcEnerVir)
    {
        __syncthreads();

        if (threadIdx.x == 0)
        {
            atomicAdd(vtot, vtot_loc);
        }
        if (threadIdx.x < SHIFTS)
        {
            fvec_inc_atomic(fshift[threadIdx.x], fshift_loc[threadIdx.x]);
        }
    }
}

__device__
static void make_dp_periodic_gpu(float *dp)
{
    /* dp cannot be outside (-pi,pi) */
    if (*dp >= CUDART_PI_F)
    {
        *dp -= 2*CUDART_PI_F;
    }
    else if (*dp < -CUDART_PI_F)
    {
        *dp += 2*CUDART_PI_F;
    }
}

template <bool calcEnerVir>
__global__
void  idihs_gpu(float *vtot, const int nbonds, const int natoms,
                const t_iatom forceatoms[], const t_iparams forceparams[],
                const fvec x[], fvec f[], fvec fshift[],
                const PbcAiuc pbcAiuc)
{
    const int        i = blockIdx.x*blockDim.x + threadIdx.x;

    __shared__ float vtot_loc;
    __shared__ fvec  fshift_loc[SHIFTS];

    if (calcEnerVir)
    {
        if (threadIdx.x == 0)
        {
            vtot_loc = 0.0f;
        }
        if (threadIdx.x < SHIFTS)
        {
            fshift_loc[threadIdx.x][XX] = 0.0f;
            fshift_loc[threadIdx.x][YY] = 0.0f;
            fshift_loc[threadIdx.x][ZZ] = 0.0f;
        }
        __syncthreads();
    }

    if (i < nbonds)
    {
        int   type = forceatoms[5*i];
        int   ai   = forceatoms[5*i + 1];
        int   aj   = forceatoms[5*i + 2];
        int   ak   = forceatoms[5*i + 3];
        int   al   = forceatoms[5*i + 4];

        fvec  r_ij;
        fvec  r_kj;
        fvec  r_kl;
        fvec  m;
        fvec  n;
        int   t1;
        int   t2;
        int   t3;
        float phi  =
            dih_angle_gpu<calcEnerVir>(x[ai], x[aj], x[ak], x[al], pbcAiuc,
                                       r_ij, r_kj, r_kl, m, n, &t1, &t2, &t3);

        /* phi can jump if phi0 is close to Pi/-Pi, which will cause huge
         * force changes if we just apply a normal harmonic.
         * Instead, we first calculate phi-phi0 and take it modulo (-Pi,Pi).
         * This means we will never have the periodicity problem, unless
         * the dihedral is Pi away from phiO, which is very unlikely due to
         * the potential.
         */
        float kA   = forceparams[type].harmonic.krA;
        float pA   = forceparams[type].harmonic.rA;

        float phi0 = pA*DEG2RAD;

        float dp   = phi - phi0;

        make_dp_periodic_gpu(&dp);

        float ddphi = -kA*dp;

        do_dih_fup_gpu<calcEnerVir>(ai, aj, ak, al, natoms,
                                    -ddphi, r_ij, r_kj, r_kl, m, n,
                                    f, fshift_loc, pbcAiuc,
                                    x, t1, t2, t3);
    }

    if (calcEnerVir)
    {
        __syncthreads();

        if (threadIdx.x == 0)
        {
            atomicAdd(vtot, vtot_loc);
        }
        if (threadIdx.x < SHIFTS)
        {
            fvec_inc_atomic(fshift[threadIdx.x], fshift_loc[threadIdx.x]);
        }
    }
}

template <bool calcEnerVir>
__global__
void pairs_gpu(int ftype, const int nbonds, const int natoms,
               const t_iatom iatoms[], const t_iparams iparams[],
               const fvec x[], fvec force[], fvec fshift[],
               const PbcAiuc pbcAiuc,
               const float md_chargeA[], const unsigned short md_cENER[],
               float energygrp_elec[], float energygrp_vdw[],
               const float scale_factor,
               int numEnergyGroups)
{
    const int       i = blockIdx.x*blockDim.x+threadIdx.x;

    __shared__ fvec fshift_loc[SHIFTS];

    if (calcEnerVir)
    {
        if (threadIdx.x < SHIFTS)
        {
            fshift_loc[threadIdx.x][XX] = 0.0f;
            fshift_loc[threadIdx.x][YY] = 0.0f;
            fshift_loc[threadIdx.x][ZZ] = 0.0f;
        }
        __syncthreads();
    }

    if (i <  nbonds)
    {
        int   itype = iatoms[3*i];
        int   ai    = iatoms[3*i + 1];
        int   aj    = iatoms[3*i + 2];

        float qq    = md_chargeA[ai]*md_chargeA[aj];
        float c6    = iparams[itype].lj14.c6A;
        float c12   = iparams[itype].lj14.c12A;

        /* Do we need to apply full periodic boundary conditions? */
        fvec  dr;
        int   fshift_index = pbcDxAiuc<calcEnerVir>(pbcAiuc, x[ai], x[aj], dr);

        float r2    = norm2_gpu(dr);
        float rinv  = rsqrtf(r2);
        float rinv2 = rinv*rinv;
        float rinv6 = rinv2*rinv2*rinv2;

        /* Calculate the Coulomb force * r */
        float velec = scale_factor*qq*rinv;

        /* Calculate the LJ force * r and add it to the Coulomb part */
        float fr    = (12.0f*c12*rinv6 - 6.0f*c6)*rinv6 + velec;

        float finvr = fr * rinv2;
        fvec  f;
        svmul_gpu(finvr, dr, f);

        /* Add the forces */
#pragma unroll
        for (int m = 0; m < DIM; m++)
        {
            atomicAdd(&force[ai][m], f[m]);
            atomicAdd(&force[aj][m], -f[m]);
        }

        if (calcEnerVir)
        {
            int gid = GID(md_cENER[ai], md_cENER[aj], numEnergyGroups);

            atomicAdd(&energygrp_elec[gid], velec);
            atomicAdd(&energygrp_vdw[gid], (c12*rinv6 - c6)*rinv6);

            if (fshift_index != CENTRAL)
            {
                fvec_inc_atomic(fshift_loc[fshift_index], f);
                fvec_dec_atomic(fshift_loc[CENTRAL], f);
            }
        }
    }

    if (calcEnerVir)
    {
        __syncthreads();
        if (threadIdx.x < SHIFTS)
        {
            fvec_inc_atomic(fshift[threadIdx.x], fshift_loc[threadIdx.x]);
        }
    }
}

/*-------------------------------- End CUDA kernels-----------------------------*/


static void setPbcAiuc(int           numPbcDim,
                       const matrix  box,
                       PbcAiuc      *pbcAiuc)
{
    if (numPbcDim > ZZ)
    {
        pbcAiuc->invBoxDiagZ = 1/box[ZZ][ZZ];
        pbcAiuc->boxZX       = box[ZZ][XX];
        pbcAiuc->boxZY       = box[ZZ][YY];
        pbcAiuc->boxZZ       = box[ZZ][ZZ];
    }
    else
    {
        pbcAiuc->invBoxDiagZ = 0;
        pbcAiuc->boxZX       = 0;
        pbcAiuc->boxZY       = 0;
        pbcAiuc->boxZZ       = 0;
    }
    if (numPbcDim > YY)
    {
        pbcAiuc->invBoxDiagY = 1/box[YY][YY];
        pbcAiuc->boxYX       = box[YY][XX];
        pbcAiuc->boxYY       = box[YY][YY];
    }
    else
    {
        pbcAiuc->invBoxDiagY = 0;
        pbcAiuc->boxYX       = 0;
        pbcAiuc->boxYY       = 0;
    }
    if (numPbcDim > XX)
    {
        pbcAiuc->invBoxDiagX = 1/box[XX][XX];
        pbcAiuc->boxXX       = box[XX][XX];
    }
    else
    {
        pbcAiuc->invBoxDiagX = 0;
        pbcAiuc->boxXX       = 0;
    }
}


/*------------------------------------------------------------------------------*/

/* static variables internal to module. Separate instances are needed
 * for each MPI rank (because, when using threadMPI, we need to avoid
 * interference and allow a global view of memory)
 */

/* variables to be integrated */

static float          *energygrp_elec_d[BO_MAX_RANKS];
static float          *energygrp_elec[BO_MAX_RANKS];
static float          *energygrp_vdw_d[BO_MAX_RANKS];
static float          *energygrp_vdw[BO_MAX_RANKS];

static cudaStream_t    streamBonded[BO_MAX_RANKS];

/*------------------------------------------------------------------------------*/
// bonded forces //

// initial setup (called once)
static void
init_gpu_bonded(GpuBondedLists    *gpuBondedLists,
                const int          rank,
                const t_idef      &idef,
                gmx_grppairener_t *grppener)
{
    cudaError_t stat;

    stat = cudaMalloc((void **)&gpuBondedLists->forceparamsDevice, sizeof(t_iparams)*idef.ntypes);
    CU_RET_ERR(stat, "cudaMalloc failed");

    stat = cudaMemcpyAsync(gpuBondedLists->forceparamsDevice, idef.iparams, idef.ntypes*sizeof(t_iparams), cudaMemcpyHostToDevice, streamBonded[rank]);
    CU_RET_ERR(stat, "cudaMemcpyAsync failed");

// device_iatom_alloc
// initiaized to 0 as static

// vtot
    gpuBondedLists->vtot.resize(F_NRE);
    stat = cudaMalloc(&gpuBondedLists->vtotDevice, F_NRE*sizeof(real));
    CU_RET_ERR(stat, "cudaMalloc failed");

// possibly these change size, but hopefully not
    stat = cudaMalloc(&energygrp_elec_d[rank], grppener->nener*sizeof(float));
    CU_RET_ERR(stat, "cudaMalloc failed");

    stat = cudaHostAlloc(&energygrp_elec[rank], grppener->nener*sizeof(float), cudaHostAllocPortable);
    CU_RET_ERR(stat, "cudaHostAlloc failed");

    stat = cudaMalloc(&energygrp_vdw_d[rank], grppener->nener*sizeof(float));
    CU_RET_ERR(stat, "cudaMalloc failed");


    stat = cudaHostAlloc(&energygrp_vdw[rank], grppener->nener*sizeof(float), cudaHostAllocPortable);
    CU_RET_ERR(stat, "cudaHostAlloc failed");

    stat = cudaMalloc(&gpuBondedLists->fshiftDevice, sizeof(fvec)*SHIFTS);
    CU_RET_ERR(stat, "cudaMalloc failed");

    gpuBondedLists->fshift.resize(SHIFTS);

    stat = cudaStreamCreate(&streamBonded[rank]);
    CU_RET_ERR(stat, "cudaStreamCreate failed");
}

//----//
// update after a neighbour list update step
void
update_gpu_bonded(GpuBondedLists *gpuBondedLists, const t_idef *idef,
                  const int xsize, const t_mdatoms *md, gmx_grppairener_t *grppener)
{
    GMX_RELEASE_ASSERT(gpuBondedLists, "Need a valid gpuBondedLists object");

    cudaError_t stat;
    int         rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Temporary hack to detect uninitialized data
    if (gpuBondedLists->forceparamsDevice == nullptr)
    {
        init_gpu_bonded(gpuBondedLists, rank, *idef, grppener);
    }

    // All x/f allocation below are temporary hacks which will be replaced with nbnxn buffers
    int  x_needed = (xsize > md->nr) ? xsize : md->nr;

    gpuBondedLists->f.resize(x_needed);

    int  n;
    int  nalloc = gpuBondedLists->nalloc;
    reallocateDeviceBuffer(&gpuBondedLists->xDevice, x_needed, &n, &nalloc, nullptr);
    nalloc = gpuBondedLists->nalloc;
    reallocateDeviceBuffer(&gpuBondedLists->fDevice, x_needed, &n, &nalloc, nullptr);
    nalloc = gpuBondedLists->nalloc;
    reallocateDeviceBuffer(&gpuBondedLists->qDevice, x_needed, &n, &nalloc, nullptr);
    nalloc = gpuBondedLists->nalloc;
    reallocateDeviceBuffer(&gpuBondedLists->cEnerDevice, x_needed, &n, &nalloc, nullptr);
    gpuBondedLists->nalloc = nalloc;

    if (gpuBondedLists->iLists[F_LJ14].size() > 0)
    {
        stat = cudaMemcpy(gpuBondedLists->cEnerDevice, md->cENER, md->nr*sizeof(unsigned short), cudaMemcpyHostToDevice);
        CU_RET_ERR(stat, "cudaMemcpy failed");
        stat = cudaMemcpy(gpuBondedLists->qDevice, md->chargeA, md->nr*sizeof(float), cudaMemcpyHostToDevice);
        CU_RET_ERR(stat, "cudaMemcpy failed");
    }

    // forceatoms
    for (int i = 0; i < F_NRE; i++)
    {
        const InteractionList &iList = gpuBondedLists->iLists[i];

        if (iList.size() > 0)
        {
            t_ilist &iListDevice = gpuBondedLists->iListsDevice[i];

            reallocateDeviceBuffer(&iListDevice.iatoms, iList.size(), &iListDevice.nr, &iListDevice.nalloc, nullptr);

            //copy  the data
            // Still uses idef ilists because of atom ordering
            stat = cudaMemcpy(iListDevice.iatoms, idef->il[i].iatoms, iList.size()*sizeof(t_iatom),
                              cudaMemcpyHostToDevice);
            CU_RET_ERR(stat, "cudaMemcpy failed");
        }
    }
}


/* -------------- */
// reset the internal variables
void
reset_gpu_bonded(GpuBondedLists *gpuBondedLists,
                 const int       size,
                 const int       nener)
{
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    dim3 blocks ( (size+TPB_BONDED-1)/TPB_BONDED, 1, 1);
    dim3 threads (TPB_BONDED, 1, 1);

    cudaError_t stat;

    fvec       *f_d_in    = gpuBondedLists->fDevice;
    fvec       *f_shift_d = gpuBondedLists->fshiftDevice;

    reset_gpu_bonded_kernel <<< blocks, threads, 0, streamBonded[rank]>>>
    (gpuBondedLists->vtotDevice, f_d_in, size, energygrp_elec_d[rank],
     energygrp_vdw_d[rank], nener, f_shift_d);

    stat = cudaGetLastError();
    CU_RET_ERR(stat, "reset bonded kernel failed");
}

/* -------------- */
void
do_bonded_gpu(t_forcerec *fr,
              int numEnergyGroups,
              const int flags, int natoms, fvec x[],
              const matrix box)
{
    if (!fr->gpuBondedLists->haveInteractions)
    {
        return;
    }

    int   nat1, nbonds;
    bool  bCalcEnerVir = (flags & (GMX_FORCE_VIRIAL | GMX_FORCE_ENERGY));

    int   rank  = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    fvec       *x_d_in = fr->gpuBondedLists->xDevice;
    fvec       *f_d_in = fr->gpuBondedLists->fDevice;

    cudaError_t stat;

    PbcAiuc     pbcAiuc;
    setPbcAiuc(fr->bMolPBC ? ePBC2npbcdim(fr->ePBC) : 0, box, &pbcAiuc);

    stat = cudaMemcpyAsync(x_d_in, x,
                           natoms*sizeof(fvec), cudaMemcpyHostToDevice, streamBonded[rank]);
    CU_RET_ERR(stat, "cudaMemcpy failed");

    t_iparams      *forceparams_d = fr->gpuBondedLists->forceparamsDevice;
    real           *charge_d      = fr->gpuBondedLists->qDevice;
    fvec           *f_shift_d     = fr->gpuBondedLists->fshiftDevice;
    unsigned short *md_cENER_d    = fr->gpuBondedLists->cEnerDevice;
    real           *vtot_d        = fr->gpuBondedLists->vtotDevice;

// launch kernels
// reordered for better overlap
    dim3 blocks, blocks_natoms;
    dim3 threads (TPB_BONDED, 1, 1);
    for (int ftype = 0; (ftype < F_NRE); ftype++)
    {
        const InteractionList &iList = fr->gpuBondedLists->iLists[ftype];

        if (iList.size() > 0)
        {
            nat1            = interaction_function[ftype].nratoms + 1;
            nbonds          = iList.size()/nat1;
            blocks.x        = (nbonds+TPB_BONDED-1)/TPB_BONDED;
            blocks_natoms.x = (natoms+TPB_BONDED-1)/TPB_BONDED;

            t_iatom *iatoms = fr->gpuBondedLists->iListsDevice[ftype].iatoms;

            if (ftype == F_PDIHS || ftype == F_PIDIHS)
            {
                if (bCalcEnerVir)
                {
                    pdihs_gpu<true> <<< blocks, threads, 0, streamBonded[rank]>>>
                    (&vtot_d[ftype], nbonds, natoms,
                     iatoms, forceparams_d,
                     x_d_in, f_d_in, f_shift_d,
                     pbcAiuc);
                }
                else
                {
                    pdihs_gpu<false> <<< blocks, threads, 0, streamBonded[rank]>>>
                    (nullptr, nbonds, natoms,
                     iatoms, forceparams_d,
                     x_d_in, f_d_in, nullptr,
                     pbcAiuc);
                }
            }
        }
    }

    for (int ftype = 0; (ftype < F_NRE); ftype++)
    {
        const InteractionList &iList = fr->gpuBondedLists->iLists[ftype];

        if (iList.size() > 0)
        {
            nat1            = interaction_function[ftype].nratoms + 1;
            nbonds          = iList.size()/nat1;
            blocks.x        = (nbonds+TPB_BONDED-1)/TPB_BONDED;
            blocks_natoms.x = (natoms+TPB_BONDED-1)/TPB_BONDED;
// in the main code they have a function pointer for this
// so we need something like that for final version
// note temp array here for output forces

            t_iatom *iatoms = fr->gpuBondedLists->iListsDevice[ftype].iatoms;

            if (ftype == F_BONDS)
            {
                if (bCalcEnerVir)
                {
                    bonds_gpu<true> <<< blocks, threads, 0, streamBonded[rank]>>>
                    (&vtot_d[ftype], nbonds, natoms,
                     iatoms, forceparams_d,
                     x_d_in, f_d_in, f_shift_d,
                     pbcAiuc);
                }
                else
                {
                    bonds_gpu<false> <<< blocks, threads, 0, streamBonded[rank]>>>
                    (nullptr, nbonds, natoms,
                     iatoms, forceparams_d,
                     x_d_in, f_d_in, nullptr,
                     pbcAiuc);
                }
            }

            if (ftype == F_ANGLES)
            {
                if (bCalcEnerVir)
                {
                    angles_gpu<true> <<< blocks, threads, 0, streamBonded[rank]>>>
                    (&vtot_d[ftype], nbonds, natoms,
                     iatoms, forceparams_d,
                     x_d_in, f_d_in, f_shift_d,
                     pbcAiuc);
                }
                else
                {
                    angles_gpu<false> <<< blocks, threads, 0, streamBonded[rank]>>>
                    (nullptr, nbonds, natoms,
                     iatoms, forceparams_d,
                     x_d_in, f_d_in, nullptr,
                     pbcAiuc);
                }
            }

            if (ftype == F_UREY_BRADLEY)
            {
                if (bCalcEnerVir)
                {
                    urey_bradley_gpu<true> <<< blocks, threads, 0, streamBonded[rank]>>>
                    (&vtot_d[ftype], nbonds, natoms,
                     iatoms, forceparams_d,
                     x_d_in, f_d_in, f_shift_d,
                     pbcAiuc);
                }
                else
                {
                    urey_bradley_gpu<false> <<< blocks, threads, 0, streamBonded[rank]>>>
                    (nullptr, nbonds, natoms,
                     iatoms, forceparams_d,
                     x_d_in, f_d_in, nullptr,
                     pbcAiuc);
                }
            }

            if (ftype == F_RBDIHS)
            {
                if (bCalcEnerVir)
                {
                    rbdihs_gpu<true> <<< blocks, threads, 0, streamBonded[rank]>>>
                    (&vtot_d[ftype], nbonds, natoms,
                     iatoms, forceparams_d,
                     x_d_in, f_d_in, f_shift_d,
                     pbcAiuc);
                }
                else
                {
                    rbdihs_gpu<false> <<< blocks, threads, 0, streamBonded[rank]>>>
                    (nullptr, nbonds, natoms,
                     iatoms, forceparams_d,
                     x_d_in, f_d_in, nullptr,
                     pbcAiuc);
                }
            }

            if (ftype == F_IDIHS)
            {
                if (bCalcEnerVir)
                {
                    idihs_gpu<true> <<< blocks, threads, 0, streamBonded[rank]>>>
                    (&vtot_d[ftype], nbonds, natoms,
                     iatoms, forceparams_d,
                     x_d_in, f_d_in, f_shift_d,
                     pbcAiuc);
                }
                else
                {
                    idihs_gpu<false> <<< blocks, threads, 0, streamBonded[rank]>>>
                    (nullptr, nbonds, natoms,
                     iatoms, forceparams_d,
                     x_d_in, f_d_in, nullptr,
                     pbcAiuc);
                }
            }

            if (ftype == F_LJ14)
            {
                if (bCalcEnerVir)
                {
                    pairs_gpu<true> <<< blocks, threads, 0, streamBonded[rank]>>>
                    (ftype, nbonds, natoms,
                     iatoms, forceparams_d,
                     x_d_in, f_d_in, f_shift_d,
                     pbcAiuc,
                     charge_d, md_cENER_d,
                     energygrp_elec_d[rank], energygrp_vdw_d[rank],
                     fr->ic->epsfac*fr->fudgeQQ,
                     numEnergyGroups);
                }
                else
                {
                    pairs_gpu<false> <<< blocks, threads, 0, streamBonded[rank]>>>
                    (ftype, nbonds, natoms,
                     iatoms, forceparams_d,
                     x_d_in, f_d_in, nullptr,
                     pbcAiuc,
                     charge_d, md_cENER_d,
                     energygrp_elec_d[rank], energygrp_vdw_d[rank],
                     fr->ic->epsfac*fr->fudgeQQ,
                     numEnergyGroups);
                }
            }
        }
    }
    stat = cudaGetLastError();
    CU_RET_ERR(stat, "kernels failed");
    return;
}

void
do_bonded_gpu_finalize(t_forcerec *fr, const int flags, const int natoms,
                       rvec *input_force, gmx_enerdata_t *enerd)
{
    if (!fr->gpuBondedLists->haveInteractions)
    {
        return;
    }

    bool               bCalcEnerVir = (flags & (GMX_FORCE_VIRIAL | GMX_FORCE_ENERGY));
    cudaError_t        stat;

    gmx_grppairener_t *grppener;
    grppener = &enerd->grpp;

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    fvec* f_d_in = fr->gpuBondedLists->fDevice;

//copy forces back
    rvec *f = as_rvec_array(fr->gpuBondedLists->f.data());
    stat = cudaMemcpyAsync(f, f_d_in,
                           natoms*sizeof(fvec), cudaMemcpyDeviceToHost, streamBonded[rank]);
    cudaStreamSynchronize(streamBonded[rank]);
    CU_RET_ERR(stat, "cudaMemcpy failed");

    for (int i = 0; i < natoms; i++)
    {
        rvec_inc(input_force[i], f[i]);
    }

    if (bCalcEnerVir)
    {
        // shift Forces
        rvec *fshift = as_rvec_array(fr->gpuBondedLists->fshift.data());
        stat = cudaMemcpyAsync(fshift, fr->gpuBondedLists->fshiftDevice,
                               SHIFTS*sizeof(fvec), cudaMemcpyDeviceToHost, streamBonded[rank]);
        cudaStreamSynchronize(streamBonded[rank]);
        CU_RET_ERR(stat, "cudaMemcpy failed");

        for (int i = 0; i < SHIFTS; i++)
        {
            rvec_inc(fr->fshift[i], fshift[i]);
        }

        // copy energies back
        real *vtot = fr->gpuBondedLists->vtot.data();
        stat = cudaMemcpyAsync(vtot, fr->gpuBondedLists->vtotDevice,
                               F_NRE*sizeof(float), cudaMemcpyDeviceToHost, streamBonded[rank]);
        cudaStreamSynchronize(streamBonded[rank]);
        CU_RET_ERR(stat, "cudaMemcpy failed");

        for (int i = 0; i < F_NRE; i++)
        {
            enerd->term[i] += vtot[i];
        }

        stat = cudaMemcpyAsync(energygrp_vdw[rank], energygrp_vdw_d[rank],
                               sizeof(float)*grppener->nener, cudaMemcpyDeviceToHost, streamBonded[rank]);
        CU_RET_ERR(stat, "cudaMemcpy failed");

        stat = cudaMemcpyAsync(energygrp_elec[rank], energygrp_elec_d[rank],
                               sizeof(float)*grppener->nener, cudaMemcpyDeviceToHost, streamBonded[rank]);
        CU_RET_ERR(stat, "cudaMemcpy failed");
        cudaStreamSynchronize(streamBonded[rank]);
        for (int i = 0; i < grppener->nener; i++)
        {
            grppener->ener[egCOUL14][i] += energygrp_elec[rank][i];
            grppener->ener[egLJ14][i]   +=  energygrp_vdw[rank][i];
        }
    }
    cudaStreamSynchronize(streamBonded[rank]); // needed if we have no copies for bufferops sync
}
/*------------------------------------------------------------------------------*/
