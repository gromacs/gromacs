/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
 * \ingroup module_listed_forces
 */

#include "gmxpre.h"

#include <math_constants.h>

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gpu_vec.cuh"
#include "gromacs/gpu_utils/gputraits.cuh"
#include "gromacs/listed_forces/gpubonded.h"
#include "gromacs/listed_forces/listed_forces.h"
#include "gromacs/math/units.h"
#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/pbcutil/gpu_pbc.cuh"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/gmxassert.h"

#include "gpubonded_impl.h"

#if defined(_MSVC)
#include <limits>
#endif

//CUDA threads per block
#define TPB_BONDED 256

/*-------------------------------- CUDA kernels-------------------------------- */
/*------------------------------------------------------------------------------*/


/*---------------- BONDED CUDA kernels--------------*/

/* Harmonic */
__device__ __forceinline__
static void harmonic_gpu(const float kA, const float xA, const float x, float *V, float *F)
{
    constexpr float half = 0.5f;
    float           dx, dx2;

    dx    = x-xA;
    dx2   = dx*dx;

    *F = -kA*dx;
    *V = half*kA*dx2;
}

template <bool calcVir, bool calcEner>
__global__
void bonds_gpu(float *vtot, const int nbonds,
               const t_iatom forceatoms[], const t_iparams forceparams[],
               const float4 xq[], fvec force[], fvec fshift[],
               const PbcAiuc pbcAiuc)
{
    const int        i = blockIdx.x * blockDim.x + threadIdx.x;

    __shared__ float vtot_loc;
    __shared__ fvec  fshift_loc[SHIFTS];

    if (calcVir || calcEner)
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
        int   ki = pbcDxAiuc<calcVir>(pbcAiuc, xq[ai], xq[aj], dx);

        float dr2 = iprod_gpu(dx, dx);
        float dr  = sqrt(dr2);

        float vbond;
        float fbond;
        harmonic_gpu(forceparams[type].harmonic.krA,
                     forceparams[type].harmonic.rA,
                     dr, &vbond, &fbond);

        if (calcEner)
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
                if (calcVir && ki != CENTRAL)
                {
                    atomicAdd(&fshift_loc[ki][m], fij);
                    atomicAdd(&fshift_loc[CENTRAL][m], -fij);
                }
            }
        }
    }

    if (calcVir || calcEner)
    {
        __syncthreads();
        if (calcEner && threadIdx.x == 0)
        {
            atomicAdd(vtot, vtot_loc);
        }
        if (calcVir && threadIdx.x < SHIFTS)
        {
            fvec_inc_atomic(fshift[threadIdx.x], fshift_loc[threadIdx.x]);
        }
    }
}

template <bool returnShift>
__device__ __forceinline__
static float bond_angle_gpu(const float4 xi, const float4 xj, const float4 xk,
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

template <bool calcVir, bool calcEner>
__global__
void angles_gpu(float *vtot, const int nbonds,
                const t_iatom forceatoms[], const t_iparams forceparams[],
                const float4 x[], fvec force[], fvec fshift[],
                const PbcAiuc pbcAiuc)
{
    __shared__ float vtot_loc;
    __shared__ fvec  fshift_loc[SHIFTS];

    const int        i = blockIdx.x*blockDim.x + threadIdx.x;

    if (calcVir || calcEner)
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
            bond_angle_gpu<calcVir>(x[ai], x[aj], x[ak], pbcAiuc,
                                    r_ij, r_kj, &cos_theta, &t1, &t2);

        float va;
        float dVdt;
        harmonic_gpu(forceparams[type].harmonic.krA,
                     forceparams[type].harmonic.rA*DEG2RAD,
                     theta, &va, &dVdt);

        if (calcEner)
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
            if (calcVir)
            {
                fvec_inc_atomic(fshift_loc[t1], f_i);
                fvec_inc_atomic(fshift_loc[CENTRAL], f_j);
                fvec_inc_atomic(fshift_loc[t2], f_k);
            }
        }

    }

    if (calcVir || calcEner)
    {
        __syncthreads();

        if (calcEner && threadIdx.x == 0)
        {
            atomicAdd(vtot, vtot_loc);
        }
        if (calcVir && threadIdx.x < SHIFTS)
        {
            fvec_inc_atomic(fshift[threadIdx.x], fshift_loc[threadIdx.x]);
        }
    }
}

template <bool calcVir, bool calcEner>
__global__
void urey_bradley_gpu(float *vtot, const int nbonds,
                      const t_iatom forceatoms[], const t_iparams forceparams[],
                      const float4 x[], fvec force[], fvec fshift[],
                      const PbcAiuc pbcAiuc)
{
    __shared__ float vtot_loc;
    __shared__ fvec  fshift_loc[SHIFTS];

    const int        i = blockIdx.x*blockDim.x + threadIdx.x;

    if (calcVir || calcEner)
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
        float theta = bond_angle_gpu<calcVir>(x[ai], x[aj], x[ak], pbcAiuc,
                                              r_ij, r_kj, &cos_theta, &t1, &t2);

        float va;
        float dVdt;
        harmonic_gpu(kthA, th0A, theta, &va, &dVdt);

        if (calcEner)
        {
            atomicAdd(&vtot_loc, va);
        }

        fvec  r_ik;
        int   ki = pbcDxAiuc<calcVir>(pbcAiuc, x[ai], x[ak], r_ik);

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
            if (calcEner)
            {
                atomicAdd(&vtot_loc, vbond);
            }

            fbond *= rsqrtf(dr2);

            for (int m = 0; m < DIM; m++)
            {
                float fik = fbond*r_ik[m];
                atomicAdd(&force[ai][m], fik);
                atomicAdd(&force[ak][m], -fik);

                if (calcVir && ki != CENTRAL)
                {
                    atomicAdd(&fshift_loc[ki][m], fik);
                    atomicAdd(&fshift_loc[CENTRAL][m], -fik);
                }
            }
        }
    }

    if (calcVir || calcEner)
    {
        __syncthreads();

        if (calcEner && threadIdx.x == 0)
        {
            atomicAdd(vtot, vtot_loc);
        }
        if (calcVir && threadIdx.x < SHIFTS)
        {
            fvec_inc_atomic(fshift[threadIdx.x], fshift_loc[threadIdx.x]);
        }
    }
}

template <bool returnShift, typename T>
__device__ __forceinline__
static float dih_angle_gpu(const T xi, const T xj, const T xk, const T xl,
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


__device__ __forceinline__
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
static void do_dih_fup_gpu(const int i, const int j, const int k, const int l,
                           const float ddphi, const fvec r_ij, const fvec r_kj, const fvec r_kl,
                           const fvec m, const fvec n, fvec force[], fvec fshift[],
                           const PbcAiuc &pbcAiuc,
                           const float4 x[], const int t1, const int t2, const int gmx_unused t3)
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

template <bool calcVir, bool calcEner>
__global__
void  pdihs_gpu(float *vtot, const int nbonds,
                const t_iatom forceatoms[], const t_iparams forceparams[],
                const float4 x[], fvec f[], fvec fshift[],
                const PbcAiuc pbcAiuc)
{
    const int        i = blockIdx.x*blockDim.x + threadIdx.x;

    __shared__ float vtot_loc;
    __shared__ fvec  fshift_loc[SHIFTS];

    if (calcVir || calcEner)
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
            dih_angle_gpu<calcVir>(x[ai], x[aj], x[ak], x[al], pbcAiuc,
                                   r_ij, r_kj, r_kl, m, n, &t1, &t2, &t3);

        float vpd;
        float ddphi;
        dopdihs_gpu(forceparams[type].pdihs.cpA,
                    forceparams[type].pdihs.phiA,
                    forceparams[type].pdihs.mult,
                    phi, &vpd, &ddphi);

        if (calcEner)
        {
            atomicAdd(&vtot_loc, vpd);
        }

        do_dih_fup_gpu<calcVir>(ai, aj, ak, al,
                                ddphi, r_ij, r_kj, r_kl, m, n,
                                f, fshift_loc, pbcAiuc,
                                x, t1, t2, t3);

    }

    if (calcVir || calcEner)
    {
        __syncthreads();

        if (calcEner && threadIdx.x == 0)
        {
            atomicAdd(vtot, vtot_loc);
        }
        if (calcVir && threadIdx.x < SHIFTS)
        {
            fvec_inc_atomic(fshift[threadIdx.x], fshift_loc[threadIdx.x]);
        }
    }
}

template <bool calcVir, bool calcEner>
__global__
void rbdihs_gpu(float *vtot, const int nbonds,
                const t_iatom forceatoms[], const t_iparams forceparams[],
                const float4 x[], fvec f[], fvec fshift[],
                const PbcAiuc pbcAiuc)
{
    constexpr float  c0 = 0.0f, c1 = 1.0f, c2 = 2.0f, c3 = 3.0f, c4 = 4.0f, c5 = 5.0f;

    __shared__ float vtot_loc;
    __shared__ fvec  fshift_loc[SHIFTS];

    const int        i = blockIdx.x*blockDim.x + threadIdx.x;

    if (calcVir || calcEner)
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
            dih_angle_gpu<calcVir>(x[ai], x[aj], x[ak], x[al], pbcAiuc,
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
        if (calcEner)
        {
            v       += cosfac*rbp;
        }
        rbp          = parm[2];
        ddphi       += c2*rbp*cosfac;
        cosfac      *= cos_phi;
        if (calcEner)
        {
            v       += cosfac*rbp;
        }
        rbp          = parm[3];
        ddphi       += c3*rbp*cosfac;
        cosfac      *= cos_phi;
        if (calcEner)
        {
            v       += cosfac*rbp;
        }
        rbp          = parm[4];
        ddphi       += c4*rbp*cosfac;
        cosfac      *= cos_phi;
        if (calcEner)
        {
            v       += cosfac*rbp;
        }
        rbp          = parm[5];
        ddphi       += c5*rbp*cosfac;
        cosfac      *= cos_phi;
        if (calcEner)
        {
            v       += cosfac*rbp;
        }

        ddphi = -ddphi*sin_phi;

        do_dih_fup_gpu<calcVir>(ai, aj, ak, al,
                                ddphi, r_ij, r_kj, r_kl, m, n,
                                f, fshift_loc, pbcAiuc,
                                x, t1, t2, t3);
        if (calcEner)
        {
            atomicAdd(&vtot_loc, v);
        }
    }

    if (calcVir || calcEner)
    {
        __syncthreads();

        if (calcEner && threadIdx.x == 0)
        {
            atomicAdd(vtot, vtot_loc);
        }
        if (calcVir && threadIdx.x < SHIFTS)
        {
            fvec_inc_atomic(fshift[threadIdx.x], fshift_loc[threadIdx.x]);
        }
    }
}

__device__ __forceinline__
static void make_dp_periodic_gpu(float *dp)
{
    /* dp cannot be outside (-pi,pi) */
    if (*dp >= CUDART_PI_F)
    {
        *dp -= 2.0f*CUDART_PI_F;
    }
    else if (*dp < -CUDART_PI_F)
    {
        *dp += 2.0f*CUDART_PI_F;
    }
}

template <bool calcVir, bool calcEner>
__global__
void  idihs_gpu(float *vtot, const int nbonds,
                const t_iatom forceatoms[], const t_iparams forceparams[],
                const float4 x[], fvec f[], fvec fshift[],
                const PbcAiuc pbcAiuc)
{
    const int        i = blockIdx.x*blockDim.x + threadIdx.x;

    __shared__ float vtot_loc;
    __shared__ fvec  fshift_loc[SHIFTS];

    if (calcVir || calcEner)
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
            dih_angle_gpu<calcVir>(x[ai], x[aj], x[ak], x[al], pbcAiuc,
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

        do_dih_fup_gpu<calcVir>(ai, aj, ak, al,
                                -ddphi, r_ij, r_kj, r_kl, m, n,
                                f, fshift_loc, pbcAiuc,
                                x, t1, t2, t3);

        if (calcEner)
        {
            atomicAdd(&vtot_loc, -0.5f*ddphi*dp);
        }
    }

    if (calcVir || calcEner)
    {
        __syncthreads();

        if (calcEner && threadIdx.x == 0)
        {
            atomicAdd(vtot, vtot_loc);
        }
        if (calcVir && threadIdx.x < SHIFTS)
        {
            fvec_inc_atomic(fshift[threadIdx.x], fshift_loc[threadIdx.x]);
        }
    }
}

template <bool calcVir, bool calcEner>
__global__
void pairs_gpu(const int nbonds,
               const t_iatom iatoms[], const t_iparams iparams[],
               const float4 xq[], fvec force[], fvec fshift[],
               const PbcAiuc pbcAiuc,
               const float scale_factor,
               float *vtotVdw, float *vtotElec)
{
    const int        i = blockIdx.x*blockDim.x+threadIdx.x;

    __shared__ float vtotVdw_loc;
    __shared__ float vtotElec_loc;
    __shared__ fvec  fshift_loc[SHIFTS];

    if (calcVir || calcEner)
    {
        if (threadIdx.x == 0)
        {
            vtotVdw_loc  = 0.0f;
            vtotElec_loc = 0.0f;
        }

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

        float qq    = xq[ai].w*xq[aj].w;
        float c6    = iparams[itype].lj14.c6A;
        float c12   = iparams[itype].lj14.c12A;

        /* Do we need to apply full periodic boundary conditions? */
        fvec  dr;
        int   fshift_index = pbcDxAiuc<calcVir>(pbcAiuc, xq[ai], xq[aj], dr);

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

        if (calcEner)
        {
            atomicAdd(&vtotVdw_loc, (c12*rinv6 - c6)*rinv6);
            atomicAdd(&vtotElec_loc, velec);
        }

        if (calcVir && fshift_index != CENTRAL)
        {
            fvec_inc_atomic(fshift_loc[fshift_index], f);
            fvec_dec_atomic(fshift_loc[CENTRAL], f);
        }
    }

    if (calcVir || calcEner)
    {
        __syncthreads();

        if (calcEner && threadIdx.x == 0)
        {
            atomicAdd(vtotVdw, vtotVdw_loc);
            atomicAdd(vtotElec, vtotElec_loc);
        }
        if (calcVir && threadIdx.x < SHIFTS)
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

namespace gmx
{

template <bool calcVir, bool calcEner>
void
GpuBonded::Impl::launchKernels(const t_forcerec *fr,
                               const matrix      box)
{
    GMX_ASSERT(haveInteractions_,
               "Cannot launch bonded GPU kernels unless bonded GPU work was scheduled");

    PbcAiuc       pbcAiuc;
    setPbcAiuc(fr->bMolPBC ? ePBC2npbcdim(fr->ePBC) : 0, box, &pbcAiuc);

    const t_iparams *forceparams_d = forceparamsDevice;
    float           *vtot_d        = vtotDevice;
    const float4    *xq_d          = xqDevice;
    fvec            *force_d       = forceDevice;
    fvec            *fshift_d      = fshiftDevice;

    for (int ftype : ftypesOnGpu)
    {
        const auto &iList = iLists[ftype];

        if (iList.size() > 0)
        {
            int                nat1   = interaction_function[ftype].nratoms + 1;
            int                nbonds = iList.size()/nat1;

            KernelLaunchConfig config;
            config.blockSize[0] = TPB_BONDED;
            config.blockSize[1] = 1;
            config.blockSize[2] = 1;
            config.gridSize[0]  = (nbonds + TPB_BONDED - 1)/TPB_BONDED;
            config.gridSize[1]  = 1;
            config.gridSize[2]  = 1;
            config.stream       = stream;

            const t_iatom *iatoms = iListsDevice[ftype].iatoms;

            if (ftype == F_PDIHS || ftype == F_PIDIHS)
            {
                auto       kernelPtr      = pdihs_gpu<calcVir, calcEner>;
                float     *ftypeEnergyPtr = vtot_d + ftype;
                const auto kernelArgs     = prepareGpuKernelArguments(kernelPtr, config,
                                                                      &ftypeEnergyPtr, &nbonds,
                                                                      &iatoms, &forceparams_d,
                                                                      &xq_d, &force_d, &fshift_d,
                                                                      &pbcAiuc);
                launchGpuKernel(kernelPtr, config, nullptr, "pdihs_gpu<calcVir, calcEner>", kernelArgs);
            }
        }
    }

    for (int ftype : ftypesOnGpu)
    {
        const auto &iList = iLists[ftype];

        if (iList.size() > 0)
        {
            int                nat1   = interaction_function[ftype].nratoms + 1;
            int                nbonds = iList.size()/nat1;

            const t_iatom     *iatoms = iListsDevice[ftype].iatoms;

            KernelLaunchConfig config;
            config.blockSize[0] = TPB_BONDED;
            config.blockSize[1] = 1;
            config.blockSize[2] = 1;
            config.gridSize[0]  = (nbonds + TPB_BONDED - 1)/TPB_BONDED;
            config.gridSize[1]  = 1;
            config.gridSize[2]  = 1;
            config.stream       = stream;

            float *ftypeEnergyPtr = vtot_d + ftype;
            // TODO consider using a map to assign the fn pointers to ftypes
            if (ftype == F_BONDS)
            {
                auto       kernelPtr  = bonds_gpu<calcVir, calcEner>;
                const auto kernelArgs = prepareGpuKernelArguments(kernelPtr, config,
                                                                  &ftypeEnergyPtr, &nbonds,
                                                                  &iatoms, &forceparams_d,
                                                                  &xq_d, &force_d, &fshift_d,
                                                                  &pbcAiuc);
                launchGpuKernel(kernelPtr, config, nullptr, "bonds_gpu<calcVir, calcEner>", kernelArgs);
            }

            if (ftype == F_ANGLES)
            {
                auto       kernelPtr  = angles_gpu<calcVir, calcEner>;
                const auto kernelArgs = prepareGpuKernelArguments(kernelPtr, config,
                                                                  &ftypeEnergyPtr, &nbonds,
                                                                  &iatoms, &forceparams_d,
                                                                  &xq_d, &force_d, &fshift_d,
                                                                  &pbcAiuc);
                launchGpuKernel(kernelPtr, config, nullptr, "angles_gpu<calcVir, calcEner>", kernelArgs);
            }

            if (ftype == F_UREY_BRADLEY)
            {
                auto       kernelPtr  = urey_bradley_gpu<calcVir, calcEner>;
                const auto kernelArgs = prepareGpuKernelArguments(kernelPtr, config,
                                                                  &ftypeEnergyPtr, &nbonds,
                                                                  &iatoms, &forceparams_d,
                                                                  &xq_d, &force_d, &fshift_d,
                                                                  &pbcAiuc);
                launchGpuKernel(kernelPtr, config, nullptr, "urey_bradley_gpu<calcVir, calcEner>", kernelArgs);
            }

            if (ftype == F_RBDIHS)
            {
                auto       kernelPtr  = rbdihs_gpu<calcVir, calcEner>;
                const auto kernelArgs = prepareGpuKernelArguments(kernelPtr, config,
                                                                  &ftypeEnergyPtr, &nbonds,
                                                                  &iatoms, &forceparams_d,
                                                                  &xq_d, &force_d, &fshift_d,
                                                                  &pbcAiuc);
                launchGpuKernel(kernelPtr, config, nullptr, "rbdihs_gpu<calcVir, calcEner>", kernelArgs);
            }

            if (ftype == F_IDIHS)
            {
                auto       kernelPtr  = idihs_gpu<calcVir, calcEner>;
                const auto kernelArgs = prepareGpuKernelArguments(kernelPtr, config,
                                                                  &ftypeEnergyPtr, &nbonds,
                                                                  &iatoms, &forceparams_d,
                                                                  &xq_d, &force_d, &fshift_d,
                                                                  &pbcAiuc);
                launchGpuKernel(kernelPtr, config, nullptr, "idihs_gpu<calcVir, calcEner>", kernelArgs);
            }

            if (ftype == F_LJ14)
            {
                auto       kernelPtr       = pairs_gpu<calcVir, calcEner>;
                float      scale_factor    = fr->ic->epsfac*fr->fudgeQQ;
                float     *lj14Energy      = vtot_d + F_LJ14;
                float     *coulomb14Energy = vtot_d + F_COUL14;
                const auto kernelArgs      = prepareGpuKernelArguments(kernelPtr, config,
                                                                       &nbonds,
                                                                       &iatoms, &forceparams_d,
                                                                       &xq_d, &force_d, &fshift_d,
                                                                       &pbcAiuc,
                                                                       &scale_factor,
                                                                       &lj14Energy, &coulomb14Energy);
                launchGpuKernel(kernelPtr, config, nullptr, "pairs_gpu<calcVir, calcEner>", kernelArgs);
            }
        }
    }
}

void
GpuBonded::launchKernels(const t_forcerec *fr,
                         int               forceFlags,
                         const matrix      box)
{
    if (forceFlags & GMX_FORCE_ENERGY)
    {
        // When we need the energy, we also need the virial
        impl_->launchKernels<true, true>
            (fr, box);
    }
    else if (forceFlags & GMX_FORCE_VIRIAL)
    {
        impl_->launchKernels<true, false>
            (fr, box);
    }
    else
    {
        impl_->launchKernels<false, false>
            (fr, box);
    }
}

} // namespace gmx
