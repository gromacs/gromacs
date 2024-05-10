/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
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

#include <cassert>

#include "gromacs/gpu_utils/cuda_arch_utils.cuh"
#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/gpu_utils/typecasts_cuda_hip.h"
#include "gromacs/gpu_utils/vectype_ops.cuh"
#include "gromacs/listed_forces/listed_forces_gpu.h"
#include "gromacs/math/units.h"
#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/pbcutil/pbc_aiuc_cuda.cuh"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/gmxassert.h"

#include "listed_forces_gpu_impl.h"

#if defined(_MSVC)
#    include <limits>
#endif


// \brief Staggered atomic force component accumulation to reduce clashes
//
// Reduce the number of atomic clashes by a theoretical max 3x by having consecutive threads
// accumulate different force components at the same time.
__device__ __forceinline__ void staggeredAtomicAddForce(float3* __restrict__ targetPtr, float3 f)
{
    int3 offset = make_int3(0, 1, 2);

    // Shift force components x, y, and z left by 2, 1, and 0, respectively
    // to end up with zxy, yzx, xyz on consecutive threads.
    f      = (threadIdx.x % 3 == 0) ? make_float3(f.y, f.z, f.x) : f;
    offset = (threadIdx.x % 3 == 0) ? make_int3(offset.y, offset.z, offset.x) : offset;
    f      = (threadIdx.x % 3 <= 1) ? make_float3(f.y, f.z, f.x) : f;
    offset = (threadIdx.x % 3 <= 1) ? make_int3(offset.y, offset.z, offset.x) : offset;

    atomicAdd(&targetPtr->x + offset.x, f.x);
    atomicAdd(&targetPtr->x + offset.y, f.y);
    atomicAdd(&targetPtr->x + offset.z, f.z);
}


/*-------------------------------- CUDA kernels-------------------------------- */
/*------------------------------------------------------------------------------*/

#define CUDA_DEG2RAD_F (CUDART_PI_F / 180.0F)

/*---------------- BONDED CUDA kernels--------------*/

/* Harmonic */
__device__ __forceinline__ static void
harmonic_gpu(const float kA, const float xA, const float x, float* V, float* F)
{
    constexpr float half = 0.5F;
    float           dx, dx2;

    dx  = x - xA;
    dx2 = dx * dx;

    *F = -kA * dx;
    *V = half * kA * dx2;
}

template<bool calcVir, bool calcEner>
__device__ __forceinline__ void bonds_gpu(const int i,
                                          float*    vtot_loc,

                                          const t_iatom   d_forceatoms[],
                                          const t_iparams d_forceparams[],
                                          const float4    gm_xq[],
                                          float3          gm_f[],
                                          float3          sm_fShiftLoc[],
                                          const PbcAiuc   pbcAiuc)
{
    const int3 bondData = *(reinterpret_cast<const int3*>(d_forceatoms + 3 * i));
    int        type     = bondData.x;
    int        ai       = bondData.y;
    int        aj       = bondData.z;

    /* dx = xi - xj, corrected for periodic boundary conditions. */
    float3 dx;
    int    ki = pbcDxAiuc<calcVir>(pbcAiuc, gm_xq[ai], gm_xq[aj], dx);

    float dr2 = norm2(dx);
    float dr  = sqrt(dr2);

    float vbond;
    float fbond;
    harmonic_gpu(d_forceparams[type].harmonic.krA, d_forceparams[type].harmonic.rA, dr, &vbond, &fbond);

    if (calcEner)
    {
        *vtot_loc += vbond;
    }

    if (dr2 != 0.0F)
    {
        fbond *= rsqrtf(dr2);

        float3 fij = fbond * dx;
        staggeredAtomicAddForce(&gm_f[ai], fij);
        staggeredAtomicAddForce(&gm_f[aj], -fij);
        if (calcVir && ki != gmx::c_centralShiftIndex)
        {
            staggeredAtomicAddForce(&sm_fShiftLoc[ki], fij);
            staggeredAtomicAddForce(&sm_fShiftLoc[gmx::c_centralShiftIndex], -fij);
        }
    }
}

template<bool returnShift>
__device__ __forceinline__ static float bond_angle_gpu(const float4   xi,
                                                       const float4   xj,
                                                       const float4   xk,
                                                       const PbcAiuc& pbcAiuc,
                                                       float3*        r_ij,
                                                       float3*        r_kj,
                                                       float*         costh,
                                                       int*           t1,
                                                       int*           t2)
/* Return value is the angle between the bonds i-j and j-k */
{
    *t1 = pbcDxAiuc<returnShift>(pbcAiuc, xi, xj, *r_ij);
    *t2 = pbcDxAiuc<returnShift>(pbcAiuc, xk, xj, *r_kj);

    *costh   = cos_angle(*r_ij, *r_kj);
    float th = acosf(*costh);

    return th;
}

template<bool calcVir, bool calcEner>
__device__ __forceinline__ void angles_gpu(const int i,
                                           float*    vtot_loc,

                                           const t_iatom   d_forceatoms[],
                                           const t_iparams d_forceparams[],
                                           const float4    gm_xq[],
                                           float3          gm_f[],
                                           float3          sm_fShiftLoc[],
                                           const PbcAiuc   pbcAiuc)
{
    const int4 angleData = *(reinterpret_cast<const int4*>(d_forceatoms + 4 * i));
    int        type      = angleData.x;
    int        ai        = angleData.y;
    int        aj        = angleData.z;
    int        ak        = angleData.w;

    float3 r_ij;
    float3 r_kj;
    float  cos_theta;
    int    t1;
    int    t2;
    float  theta = bond_angle_gpu<calcVir>(
            gm_xq[ai], gm_xq[aj], gm_xq[ak], pbcAiuc, &r_ij, &r_kj, &cos_theta, &t1, &t2);

    float va;
    float dVdt;
    harmonic_gpu(d_forceparams[type].harmonic.krA,
                 d_forceparams[type].harmonic.rA * CUDA_DEG2RAD_F,
                 theta,
                 &va,
                 &dVdt);

    if (calcEner)
    {
        *vtot_loc += va;
    }

    float cos_theta2 = cos_theta * cos_theta;
    if (cos_theta2 < 1.0F)
    {
        float st    = dVdt * rsqrtf(1.0F - cos_theta2);
        float sth   = st * cos_theta;
        float nrij2 = norm2(r_ij);
        float nrkj2 = norm2(r_kj);

        float nrij_1 = rsqrtf(nrij2);
        float nrkj_1 = rsqrtf(nrkj2);

        float cik = st * nrij_1 * nrkj_1;
        float cii = sth * nrij_1 * nrij_1;
        float ckk = sth * nrkj_1 * nrkj_1;

        float3 f_i = cii * r_ij - cik * r_kj;
        float3 f_k = ckk * r_kj - cik * r_ij;
        float3 f_j = -f_i - f_k;

        staggeredAtomicAddForce(&gm_f[ai], f_i);
        staggeredAtomicAddForce(&gm_f[aj], f_j);
        staggeredAtomicAddForce(&gm_f[ak], f_k);

        if (calcVir)
        {
            staggeredAtomicAddForce(&sm_fShiftLoc[t1], f_i);
            staggeredAtomicAddForce(&sm_fShiftLoc[gmx::c_centralShiftIndex], f_j);
            staggeredAtomicAddForce(&sm_fShiftLoc[t2], f_k);
        }
    }
}

template<bool calcVir, bool calcEner>
__device__ __forceinline__ void urey_bradley_gpu(const int i,
                                                 float*    vtot_loc,

                                                 const t_iatom   d_forceatoms[],
                                                 const t_iparams d_forceparams[],
                                                 const float4    gm_xq[],
                                                 float3          gm_f[],
                                                 float3          sm_fShiftLoc[],
                                                 const PbcAiuc   pbcAiuc)
{
    const int4 ubData = *(reinterpret_cast<const int4*>(d_forceatoms + 4 * i));
    int        type   = ubData.x;
    int        ai     = ubData.y;
    int        aj     = ubData.z;
    int        ak     = ubData.w;

    float th0A = d_forceparams[type].u_b.thetaA * CUDA_DEG2RAD_F;
    float kthA = d_forceparams[type].u_b.kthetaA;
    float r13A = d_forceparams[type].u_b.r13A;
    float kUBA = d_forceparams[type].u_b.kUBA;

    float3 r_ij;
    float3 r_kj;
    float  cos_theta;
    int    t1;
    int    t2;
    float  theta = bond_angle_gpu<calcVir>(
            gm_xq[ai], gm_xq[aj], gm_xq[ak], pbcAiuc, &r_ij, &r_kj, &cos_theta, &t1, &t2);

    float va;
    float dVdt;
    harmonic_gpu(kthA, th0A, theta, &va, &dVdt);

    if (calcEner)
    {
        *vtot_loc += va;
    }

    float3 r_ik;
    int    ki = pbcDxAiuc<calcVir>(pbcAiuc, gm_xq[ai], gm_xq[ak], r_ik);

    float dr2 = norm2(r_ik);
    float dr  = dr2 * rsqrtf(dr2);

    float vbond;
    float fbond;
    harmonic_gpu(kUBA, r13A, dr, &vbond, &fbond);

    float cos_theta2 = cos_theta * cos_theta;

    float3 f_i = make_float3(0.0F);
    float3 f_j = make_float3(0.0F);
    float3 f_k = make_float3(0.0F);

    if (cos_theta2 < 1.0F)
    {
        float st  = dVdt * rsqrtf(1.0F - cos_theta2);
        float sth = st * cos_theta;

        float nrkj2 = norm2(r_kj);
        float nrij2 = norm2(r_ij);

        float cik = st * rsqrtf(nrkj2 * nrij2);
        float cii = sth / nrij2;
        float ckk = sth / nrkj2;

        f_i = cii * r_ij - cik * r_kj;
        f_k = ckk * r_kj - cik * r_ij;
        f_j = -f_i - f_k;

        if (calcVir)
        {
            staggeredAtomicAddForce(&sm_fShiftLoc[t1], f_i);
            staggeredAtomicAddForce(&sm_fShiftLoc[gmx::c_centralShiftIndex], f_j);
            staggeredAtomicAddForce(&sm_fShiftLoc[t2], f_k);
        }
    }

    /* Time for the bond calculations */
    if (dr2 != 0.0F)
    {
        if (calcEner)
        {
            *vtot_loc += vbond;
        }

        fbond *= rsqrtf(dr2);

        float3 fik = fbond * r_ik;
        f_i += fik;
        f_k -= fik;

        if (calcVir && ki != gmx::c_centralShiftIndex)
        {
            staggeredAtomicAddForce(&sm_fShiftLoc[ki], fik);
            staggeredAtomicAddForce(&sm_fShiftLoc[gmx::c_centralShiftIndex], -fik);
        }
    }
    if ((cos_theta2 < 1.0F) || (dr2 != 0.0F))
    {
        staggeredAtomicAddForce(&gm_f[ai], f_i);
        staggeredAtomicAddForce(&gm_f[ak], f_k);
    }

    if (cos_theta2 < 1.0F)
    {
        staggeredAtomicAddForce(&gm_f[aj], f_j);
    }
}

template<bool returnShift, typename T>
__device__ __forceinline__ static float dih_angle_gpu(const T        xi,
                                                      const T        xj,
                                                      const T        xk,
                                                      const T        xl,
                                                      const PbcAiuc& pbcAiuc,
                                                      float3*        r_ij,
                                                      float3*        r_kj,
                                                      float3*        r_kl,
                                                      float3*        m,
                                                      float3*        n,
                                                      int*           t1,
                                                      int*           t2,
                                                      int*           t3)
{
    *t1 = pbcDxAiuc<returnShift>(pbcAiuc, xi, xj, *r_ij);
    *t2 = pbcDxAiuc<returnShift>(pbcAiuc, xk, xj, *r_kj);
    *t3 = pbcDxAiuc<returnShift>(pbcAiuc, xk, xl, *r_kl);

    *m         = cprod(*r_ij, *r_kj);
    *n         = cprod(*r_kj, *r_kl);
    float phi  = gmx_angle(*m, *n);
    float ipr  = iprod(*r_ij, *n);
    float sign = (ipr < 0.0F) ? -1.0F : 1.0F;
    phi        = sign * phi;

    return phi;
}
template<bool returnShift, typename T>
__device__ __forceinline__ static float dih_angle_gpu_sincos(const T        xi,
                                                             const T        xj,
                                                             const T        xk,
                                                             const T        xl,
                                                             const PbcAiuc& pbcAiuc,
                                                             float3*        r_ij,
                                                             float3*        r_kj,
                                                             float3*        r_kl,
                                                             float3*        m,
                                                             float3*        n,
                                                             int*           t1,
                                                             int*           t2,
                                                             float*         cosval)
{
    *t1 = pbcDxAiuc<returnShift>(pbcAiuc, xi, xj, *r_ij);
    *t2 = pbcDxAiuc<returnShift>(pbcAiuc, xk, xj, *r_kj);
    pbcDxAiuc<returnShift>(pbcAiuc, xk, xl, *r_kl);

    *m = cprod(*r_ij, *r_kj);
    *n = cprod(*r_kj, *r_kl);

    float3 w    = cprod(*m, *n);
    float  wLen = norm(w);
    float  s    = iprod(*m, *n);

    float mLenSq = norm2(*m);
    float nLenSq = norm2(*n);
    float mnInv  = rsqrtf(mLenSq * nLenSq);

    *cosval      = s * mnInv;
    float sinval = wLen * mnInv;

    float ipr  = iprod(*r_ij, *n);
    float sign = (ipr < 0.0F) ? -1.0F : 1.0F;
    return sign * sinval;
}


__device__ __forceinline__ static void
dopdihs_gpu(const float cpA, const float phiA, const int mult, const float phi, float* v, float* f)
{
    float mdphi, sdphi;

    mdphi = mult * phi - phiA * CUDA_DEG2RAD_F;
    sdphi = sinf(mdphi);
    *v    = cpA * (1.0F + cosf(mdphi));
    *f    = -cpA * mult * sdphi;
}

template<bool calcVir>
__device__ __forceinline__ static void do_dih_fup_gpu(const int      i,
                                                      const int      j,
                                                      const int      k,
                                                      const int      l,
                                                      const float    ddphi,
                                                      const float3   r_ij,
                                                      const float3   r_kj,
                                                      const float3   r_kl,
                                                      const float3   m,
                                                      const float3   n,
                                                      float3         gm_f[],
                                                      float3         sm_fShiftLoc[],
                                                      const PbcAiuc& pbcAiuc,
                                                      const float4   gm_xq[],
                                                      const int      t1,
                                                      const int      t2)
{
    float iprm  = norm2(m);
    float iprn  = norm2(n);
    float nrkj2 = norm2(r_kj);
    float toler = nrkj2 * GMX_REAL_EPS;
    if ((iprm > toler) && (iprn > toler))
    {
        float  nrkj_1 = rsqrtf(nrkj2); // replacing std::invsqrt call
        float  nrkj_2 = nrkj_1 * nrkj_1;
        float  nrkj   = nrkj2 * nrkj_1;
        float  a      = -ddphi * nrkj / iprm;
        float3 f_i    = a * m;
        float  b      = ddphi * nrkj / iprn;
        float3 f_l    = b * n;
        float  p      = iprod(r_ij, r_kj);
        p *= nrkj_2;
        float q = iprod(r_kl, r_kj);
        q *= nrkj_2;
        float3 uvec = p * f_i;
        float3 vvec = q * f_l;
        float3 svec = uvec - vvec;
        float3 f_j  = f_i - svec;
        float3 f_k  = f_l + svec;

        staggeredAtomicAddForce(&gm_f[i], f_i);
        staggeredAtomicAddForce(&gm_f[j], -f_j);
        staggeredAtomicAddForce(&gm_f[k], -f_k);
        staggeredAtomicAddForce(&gm_f[l], f_l);

        if (calcVir)
        {
            float3 dx_jl;
            int    t3 = pbcDxAiuc<calcVir>(pbcAiuc, gm_xq[l], gm_xq[j], dx_jl);

            staggeredAtomicAddForce(&sm_fShiftLoc[t1], f_i);
            staggeredAtomicAddForce(&sm_fShiftLoc[gmx::c_centralShiftIndex], -f_j);
            staggeredAtomicAddForce(&sm_fShiftLoc[t2], -f_k);
            staggeredAtomicAddForce(&sm_fShiftLoc[t3], f_l);
        }
    }
}

template<bool calcVir, bool calcEner>
__device__ __forceinline__ void pdihs_gpu(const int i,
                                          float*    vtot_loc,

                                          const t_iatom   d_forceatoms[],
                                          const t_iparams d_forceparams[],
                                          const float4    gm_xq[],
                                          float3          gm_f[],
                                          float3          sm_fShiftLoc[],
                                          const PbcAiuc   pbcAiuc)
{
    int type = d_forceatoms[5 * i];
    int ai   = d_forceatoms[5 * i + 1];
    int aj   = d_forceatoms[5 * i + 2];
    int ak   = d_forceatoms[5 * i + 3];
    int al   = d_forceatoms[5 * i + 4];

    float3 r_ij;
    float3 r_kj;
    float3 r_kl;
    float3 m;
    float3 n;
    int    t1;
    int    t2;
    int    t3;
    float  phi = dih_angle_gpu<calcVir>(
            gm_xq[ai], gm_xq[aj], gm_xq[ak], gm_xq[al], pbcAiuc, &r_ij, &r_kj, &r_kl, &m, &n, &t1, &t2, &t3);

    float vpd;
    float ddphi;
    dopdihs_gpu(d_forceparams[type].pdihs.cpA,
                d_forceparams[type].pdihs.phiA,
                d_forceparams[type].pdihs.mult,
                phi,
                &vpd,
                &ddphi);

    if (calcEner)
    {
        *vtot_loc += vpd;
    }

    do_dih_fup_gpu<calcVir>(
            ai, aj, ak, al, ddphi, r_ij, r_kj, r_kl, m, n, gm_f, sm_fShiftLoc, pbcAiuc, gm_xq, t1, t2);
}

template<bool calcVir, bool calcEner>
__device__ __forceinline__ void rbdihs_gpu(const int i,
                                           float*    vtot_loc,

                                           const t_iatom   d_forceatoms[],
                                           const t_iparams d_forceparams[],
                                           const float4    gm_xq[],
                                           float3          gm_f[],
                                           float3          sm_fShiftLoc[],
                                           const PbcAiuc   pbcAiuc)
{
    constexpr float c0 = 0.0F, c1 = 1.0F, c2 = 2.0F, c3 = 3.0F, c4 = 4.0F, c5 = 5.0F;

    int type = d_forceatoms[5 * i];
    int ai   = d_forceatoms[5 * i + 1];
    int aj   = d_forceatoms[5 * i + 2];
    int ak   = d_forceatoms[5 * i + 3];
    int al   = d_forceatoms[5 * i + 4];

    float3 r_ij;
    float3 r_kj;
    float3 r_kl;
    float3 m;
    float3 n;
    int    t1;
    int    t2;

    float cos_phi;
    // Changing the sign of sin and cos to convert to polymer convention
    float negative_sin_phi = dih_angle_gpu_sincos<calcVir>(
            gm_xq[ai], gm_xq[aj], gm_xq[ak], gm_xq[al], pbcAiuc, &r_ij, &r_kj, &r_kl, &m, &n, &t1, &t2, &cos_phi);
    cos_phi *= -1;


    float parm[NR_RBDIHS];
    for (int j = 0; j < NR_RBDIHS; j++)
    {
        parm[j] = d_forceparams[type].rbdihs.rbcA[j];
    }
    /* Calculate cosine powers */
    /* Calculate the energy */
    /* Calculate the derivative */
    float v      = parm[0];
    float ddphi  = c0;
    float cosfac = c1;

    float rbp = parm[1];
    ddphi += rbp * cosfac;
    cosfac *= cos_phi;
    if (calcEner)
    {
        v += cosfac * rbp;
    }
    rbp = parm[2];
    ddphi += c2 * rbp * cosfac;
    cosfac *= cos_phi;
    if (calcEner)
    {
        v += cosfac * rbp;
    }
    rbp = parm[3];
    ddphi += c3 * rbp * cosfac;
    cosfac *= cos_phi;
    if (calcEner)
    {
        v += cosfac * rbp;
    }
    rbp = parm[4];
    ddphi += c4 * rbp * cosfac;
    cosfac *= cos_phi;
    if (calcEner)
    {
        v += cosfac * rbp;
    }
    rbp = parm[5];
    ddphi += c5 * rbp * cosfac;
    cosfac *= cos_phi;
    if (calcEner)
    {
        v += cosfac * rbp;
    }

    ddphi = ddphi * negative_sin_phi;

    do_dih_fup_gpu<calcVir>(
            ai, aj, ak, al, ddphi, r_ij, r_kj, r_kl, m, n, gm_f, sm_fShiftLoc, pbcAiuc, gm_xq, t1, t2);
    if (calcEner)
    {
        *vtot_loc += v;
    }
}

__device__ __forceinline__ static void make_dp_periodic_gpu(float* dp)
{
    /* dp cannot be outside (-pi,pi) */
    if (*dp >= CUDART_PI_F)
    {
        *dp -= 2.0F * CUDART_PI_F;
    }
    else if (*dp < -CUDART_PI_F)
    {
        *dp += 2.0F * CUDART_PI_F;
    }
}

template<bool calcVir, bool calcEner>
__device__ __forceinline__ void idihs_gpu(const int i,
                                          float*    vtot_loc,

                                          const t_iatom   d_forceatoms[],
                                          const t_iparams d_forceparams[],
                                          const float4    gm_xq[],
                                          float3          gm_f[],
                                          float3          sm_fShiftLoc[],
                                          const PbcAiuc   pbcAiuc)
{
    int type = d_forceatoms[5 * i];
    int ai   = d_forceatoms[5 * i + 1];
    int aj   = d_forceatoms[5 * i + 2];
    int ak   = d_forceatoms[5 * i + 3];
    int al   = d_forceatoms[5 * i + 4];

    float3 r_ij;
    float3 r_kj;
    float3 r_kl;
    float3 m;
    float3 n;
    int    t1;
    int    t2;
    int    t3;
    float  phi = dih_angle_gpu<calcVir>(
            gm_xq[ai], gm_xq[aj], gm_xq[ak], gm_xq[al], pbcAiuc, &r_ij, &r_kj, &r_kl, &m, &n, &t1, &t2, &t3);

    /* phi can jump if phi0 is close to Pi/-Pi, which will cause huge
     * force changes if we just apply a normal harmonic.
     * Instead, we first calculate phi-phi0 and take it modulo (-Pi,Pi).
     * This means we will never have the periodicity problem, unless
     * the dihedral is Pi away from phiO, which is very unlikely due to
     * the potential.
     */
    float kA = d_forceparams[type].harmonic.krA;
    float pA = d_forceparams[type].harmonic.rA;

    float phi0 = pA * CUDA_DEG2RAD_F;

    float dp = phi - phi0;

    make_dp_periodic_gpu(&dp);

    float ddphi = -kA * dp;

    do_dih_fup_gpu<calcVir>(
            ai, aj, ak, al, -ddphi, r_ij, r_kj, r_kl, m, n, gm_f, sm_fShiftLoc, pbcAiuc, gm_xq, t1, t2);

    if (calcEner)
    {
        *vtot_loc += -0.5F * ddphi * dp;
    }
}

template<bool calcVir, bool calcEner>
__device__ __forceinline__ void pairs_gpu(const int i,

                                          const t_iatom   d_forceatoms[],
                                          const t_iparams iparams[],
                                          const float4    gm_xq[],
                                          float3          gm_f[],
                                          float3          sm_fShiftLoc[],
                                          const PbcAiuc   pbcAiuc,
                                          const float     scale_factor,
                                          float*          vtotVdw_loc,
                                          float*          vtotElec_loc)
{
    // TODO this should be made into a separate type, the GPU and CPU sizes should be compared
    const int3 pairData = *(reinterpret_cast<const int3*>(d_forceatoms + 3 * i));
    int        type     = pairData.x;
    int        ai       = pairData.y;
    int        aj       = pairData.z;

    float qq  = gm_xq[ai].w * gm_xq[aj].w;
    float c6  = iparams[type].lj14.c6A;
    float c12 = iparams[type].lj14.c12A;

    /* Do we need to apply full periodic boundary conditions? */
    float3 dr;
    int    fshift_index = pbcDxAiuc<calcVir>(pbcAiuc, gm_xq[ai], gm_xq[aj], dr);

    float r2    = norm2(dr);
    float rinv  = rsqrtf(r2);
    float rinv2 = rinv * rinv;
    float rinv6 = rinv2 * rinv2 * rinv2;

    /* Calculate the Coulomb force * r */
    float velec = scale_factor * qq * rinv;

    /* Calculate the LJ force * r and add it to the Coulomb part */
    float fr = (12.0F * c12 * rinv6 - 6.0F * c6) * rinv6 + velec;

    float  finvr = fr * rinv2;
    float3 f     = finvr * dr;

    /* Add the forces */
    staggeredAtomicAddForce(&gm_f[ai], f);
    staggeredAtomicAddForce(&gm_f[aj], -f);
    if (calcVir && fshift_index != gmx::c_centralShiftIndex)
    {
        staggeredAtomicAddForce(&sm_fShiftLoc[fshift_index], f);
        staggeredAtomicAddForce(&sm_fShiftLoc[gmx::c_centralShiftIndex], -f);
    }

    if (calcEner)
    {
        *vtotVdw_loc += (c12 * rinv6 - c6) * rinv6;
        *vtotElec_loc += velec;
    }
}

namespace gmx
{

template<bool calcVir, bool calcEner>
__global__ void bonded_kernel_gpu(BondedGpuKernelParameters kernelParams,
                                  BondedGpuKernelBuffers    kernelBuffers,
                                  float4*                   gm_xq,
                                  float3*                   gm_f,
                                  float3*                   gm_fShift)
{
    assert(blockDim.y == 1 && blockDim.z == 1);
    const int tid          = blockIdx.x * blockDim.x + threadIdx.x;
    float     vtot_loc     = 0.0F;
    float     vtotElec_loc = 0.0F; // Used only for F_LJ14

    extern __shared__ char sm_dynamicShmem[];
    char*                  sm_nextSlotPtr = sm_dynamicShmem;
    float3*                sm_fShiftLoc   = reinterpret_cast<float3*>(sm_nextSlotPtr);
    sm_nextSlotPtr += c_numShiftVectors * sizeof(float3);

    if (calcVir)
    {
        if (threadIdx.x < c_numShiftVectors)
        {
            sm_fShiftLoc[threadIdx.x] = make_float3(0.0F, 0.0F, 0.0F);
        }
        __syncthreads();
    }

    int  fType;
    bool threadComputedPotential = false;
#pragma unroll
    for (int j = 0; j < numFTypesOnGpu; j++)
    {
        if (tid >= kernelParams.fTypeRangeStart[j] && tid <= kernelParams.fTypeRangeEnd[j])
        {
            const int      numBonds = kernelParams.numFTypeBonds[j];
            int            fTypeTid = tid - kernelParams.fTypeRangeStart[j];
            const t_iatom* iatoms   = kernelBuffers.d_iatoms[j];
            fType                   = kernelParams.fTypesOnGpu[j];
            if (calcEner)
            {
                threadComputedPotential = true;
            }

            if (fTypeTid >= numBonds)
            {
                break;
            }


            switch (fType)
            {
                case F_BONDS:
                    bonds_gpu<calcVir, calcEner>(fTypeTid,
                                                 &vtot_loc,
                                                 iatoms,
                                                 kernelBuffers.d_forceParams,
                                                 gm_xq,
                                                 gm_f,
                                                 sm_fShiftLoc,
                                                 kernelParams.pbcAiuc);
                    break;
                case F_ANGLES:
                    angles_gpu<calcVir, calcEner>(fTypeTid,
                                                  &vtot_loc,
                                                  iatoms,
                                                  kernelBuffers.d_forceParams,
                                                  gm_xq,
                                                  gm_f,
                                                  sm_fShiftLoc,
                                                  kernelParams.pbcAiuc);
                    break;
                case F_UREY_BRADLEY:
                    urey_bradley_gpu<calcVir, calcEner>(fTypeTid,
                                                        &vtot_loc,
                                                        iatoms,
                                                        kernelBuffers.d_forceParams,
                                                        gm_xq,
                                                        gm_f,
                                                        sm_fShiftLoc,
                                                        kernelParams.pbcAiuc);
                    break;
                case F_PDIHS:
                case F_PIDIHS:
                    pdihs_gpu<calcVir, calcEner>(fTypeTid,
                                                 &vtot_loc,
                                                 iatoms,
                                                 kernelBuffers.d_forceParams,
                                                 gm_xq,
                                                 gm_f,
                                                 sm_fShiftLoc,
                                                 kernelParams.pbcAiuc);
                    break;
                case F_RBDIHS:
                    rbdihs_gpu<calcVir, calcEner>(fTypeTid,
                                                  &vtot_loc,
                                                  iatoms,
                                                  kernelBuffers.d_forceParams,
                                                  gm_xq,
                                                  gm_f,
                                                  sm_fShiftLoc,
                                                  kernelParams.pbcAiuc);
                    break;
                case F_IDIHS:
                    idihs_gpu<calcVir, calcEner>(fTypeTid,
                                                 &vtot_loc,
                                                 iatoms,
                                                 kernelBuffers.d_forceParams,
                                                 gm_xq,
                                                 gm_f,
                                                 sm_fShiftLoc,
                                                 kernelParams.pbcAiuc);
                    break;
                case F_LJ14:
                    pairs_gpu<calcVir, calcEner>(fTypeTid,
                                                 iatoms,
                                                 kernelBuffers.d_forceParams,
                                                 gm_xq,
                                                 gm_f,
                                                 sm_fShiftLoc,
                                                 kernelParams.pbcAiuc,
                                                 kernelParams.electrostaticsScaleFactor,
                                                 &vtot_loc,
                                                 &vtotElec_loc);
                    break;
            }
            break;
        }
    }

    if (threadComputedPotential)
    {
        float* vtot     = kernelBuffers.d_vTot + fType;
        float* vtotElec = kernelBuffers.d_vTot + F_COUL14;

        // Perform warp-local reduction
        vtot_loc += __shfl_down_sync(c_fullWarpMask, vtot_loc, 1);
        vtotElec_loc += __shfl_up_sync(c_fullWarpMask, vtotElec_loc, 1);
        if (threadIdx.x & 1)
        {
            vtot_loc = vtotElec_loc;
        }
#pragma unroll 4
        for (int i = 2; i < warpSize; i *= 2)
        {
            vtot_loc += __shfl_down_sync(c_fullWarpMask, vtot_loc, i);
        }

        // Write reduced warp results into global memory
        if (threadIdx.x % warpSize == 0)
        {
            atomicAdd(vtot, vtot_loc);
        }
        else if ((threadIdx.x % warpSize == 1) && (fType == F_LJ14))
        {
            atomicAdd(vtotElec, vtot_loc);
        }
    }
    /* Accumulate shift vectors from shared memory to global memory on the first c_numShiftVectors threads of the block. */
    if (calcVir)
    {
        __syncthreads();
        if (threadIdx.x < c_numShiftVectors)
        {
            staggeredAtomicAddForce(&gm_fShift[threadIdx.x], sm_fShiftLoc[threadIdx.x]);
        }
    }
}


/*-------------------------------- End CUDA kernels-----------------------------*/


template<bool calcVir, bool calcEner>
void ListedForcesGpu::Impl::launchKernel()
{
    GMX_ASSERT(haveInteractions_,
               "Cannot launch bonded GPU kernels unless bonded GPU work was scheduled");

    wallcycle_start_nocount(wcycle_, WallCycleCounter::LaunchGpuPp);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::LaunchGpuBonded);

    int fTypeRangeEnd = kernelParams_.fTypeRangeEnd[numFTypesOnGpu - 1];

    if (fTypeRangeEnd < 0)
    {
        return;
    }

    auto kernelPtr = bonded_kernel_gpu<calcVir, calcEner>;

    const auto kernelArgs = prepareGpuKernelArguments(
            kernelPtr, kernelLaunchConfig_, &kernelParams_, &kernelBuffers_, &d_xq_, &d_f_, &d_fShift_);

    launchGpuKernel(kernelPtr,
                    kernelLaunchConfig_,
                    deviceStream_,
                    nullptr,
                    "bonded_kernel_gpu<calcVir, calcEner>",
                    kernelArgs);

    wallcycle_sub_stop(wcycle_, WallCycleSubCounter::LaunchGpuBonded);
    wallcycle_stop(wcycle_, WallCycleCounter::LaunchGpuPp);
}

void ListedForcesGpu::launchKernel(const gmx::StepWorkload& stepWork)
{
    if (stepWork.computeEnergy)
    {
        // When we need the energy, we also need the virial
        impl_->launchKernel<true, true>();
    }
    else if (stepWork.computeVirial)
    {
        impl_->launchKernel<true, false>();
    }
    else
    {
        impl_->launchKernel<false, false>();
    }
}

} // namespace gmx
