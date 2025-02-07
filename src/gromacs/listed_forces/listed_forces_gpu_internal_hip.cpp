/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 * \brief Implements HIP bonded functionality
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

#include <cassert>

#include "hip/hip_math_constants.h"

#include "gromacs/gpu_utils/hip_sycl_kernel_utils.h"
#include "gromacs/gpu_utils/hiputils.h"
#include "gromacs/gpu_utils/typecasts_cuda_hip.h"
#include "gromacs/gpu_utils/vectype_ops_hip.h"
#include "gromacs/listed_forces/listed_forces_gpu.h"
#include "gromacs/math/units.h"
#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/pbcutil/pbc_aiuc_hip.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/gmxassert.h"

#include "listed_forces_gpu_impl.h"

#if defined(_MSVC)
#    include <limits>
#endif

/*-------------------------------- HIP kernels-------------------------------- */
/*------------------------------------------------------------------------------*/

#define HIP_DEG2RAD_F (HIP_PI / 180.0F)

/*---------------- BONDED HIP kernels--------------*/

/* There are some troubles optimizing the dynamic array
 * member access despite the fact that all the loops are unrolled.
 *
 * See https://developer.nvidia.com/blog/fast-dynamic-indexing-private-arrays-cuda/
 * for a details on why dynamic access is problematic.
 *
 * This wrapper avoid dynamic accesses into the array, replacing them
 * with a `switch` instead.
 */
template<typename T>
struct FTypeArray
{
    static_assert(gmx::numFTypesOnGpu == 8,
                  "Please update the member initializer list and the switch below");

    constexpr FTypeArray(const T in[gmx::numFTypesOnGpu]) :
        data{ in[0], in[1], in[2], in[3], in[4], in[5], in[6], in[7] }
    {
    }
    __device__ __forceinline__ constexpr T operator[](int index) const
    {
        // return values[index];

        switch (index)
        {
            case 0: return data[0];
            case 1: return data[1];
            case 2: return data[2];
            case 3: return data[3];
            case 4: return data[4];
            case 5: return data[5];
            case 6: return data[6];
            default: return data[7];
        }
    }
    T data[gmx::numFTypesOnGpu];
};


__device__ __forceinline__ float3 hipHeadSegmentedSum(float3 input, const bool flag)
{
    unsigned long long warpFlags = __ballot(flag);

    warpFlags >>= 1;
    unsigned int lane_id = (threadIdx.x & (warpSize - 1));

    warpFlags &= static_cast<unsigned long long>(-1)
                 ^ ((static_cast<unsigned long long>(1) << lane_id) - 1U);
    warpFlags |= static_cast<unsigned long long>(1) << (warpSize - 1U);
    unsigned int valid_items = __ffsll(warpFlags);

    float3 output = input;
#pragma unroll
    for (unsigned int offset = 1; offset < warpSize; offset *= 2)
    {
        float3 value = make_float3(__shfl_down(output.x, offset),
                                   __shfl_down(output.y, offset),
                                   __shfl_down(output.z, offset));
        if (lane_id + offset < valid_items)
        {
            output += value;
        }
    }
    return output;
}

__device__ __forceinline__ void storeForce(float3 gm_f[], int i, float3 f)
{
    // Combine forces of consecutive lanes that write forces for the same atom
    // but do it only if all lanes are active (the last block may have fewer lanes active or some
    // lanes may not pass tolerance conditions etc.)
    if (__popcll(__ballot(1)) == warpSize)
    {
        const int  prev_lane_i = __shfl_up(i, 1);
        const bool head        = (threadIdx.x & (warpSize - 1)) == 0 || i != prev_lane_i;

        f = hipHeadSegmentedSum(f, head);
        if (head)
        {
            // Reduce the number of conflicts that left after combining consecutive forces
            // by a factor of 3: different lanes write x, y and z in a different order
            int3 j = make_int3(0, 1, 2);
            f      = (threadIdx.x % 3 == 0) ? make_float3(f.y, f.z, f.x) : f;
            j      = (threadIdx.x % 3 == 0) ? make_int3(j.y, j.z, j.x) : j;
            f      = (threadIdx.x % 3 <= 1) ? make_float3(f.y, f.z, f.x) : f;
            j      = (threadIdx.x % 3 <= 1) ? make_int3(j.y, j.z, j.x) : j;

            atomicAdd(&gm_f[i].x + j.x, f.x);
            atomicAdd(&gm_f[i].x + j.y, f.y);
            atomicAdd(&gm_f[i].x + j.z, f.z);
        }
    }
    else
    {
        atomicAdd(&gm_f[i], f);
    }
}

/* Harmonic */

template<bool calcEner>
__device__ __forceinline__ static void
harmonic_gpu(const float kA, const float xA, const float x, float* V, float* F)
{
    constexpr float half = 0.5F;
    float           dx, dx2;

    dx  = x - xA;
    dx2 = dx * dx;

    *F = -kA * dx;
    if constexpr (calcEner)
    {
        *V = half * kA * dx2;
    }
}

template<bool calcVir, bool calcEner>
__device__ void bonds_gpu(const int       i,
                          float*          vtot_loc,
                          const int       numBonds,
                          const t_iatom   d_forceatoms[],
                          const t_iparams d_forceparams[],
                          const float4    gm_xq[],
                          float3          gm_f[],
                          float3          sm_fShiftLoc[],
                          const PbcAiuc   pbcAiuc)
{
    if (i < numBonds)
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
        harmonic_gpu<calcEner>(
                d_forceparams[type].harmonic.krA, d_forceparams[type].harmonic.rA, dr, &vbond, &fbond);

        if constexpr (calcEner)
        {
            *vtot_loc += vbond;
        }

        float3 fij = make_float3(0.0F);
        if (dr2 != 0.0F)
        {
            fbond *= __frsqrt_rn(dr2);

            fij = fbond * dx;
            if constexpr (calcVir)
            {
                if (ki != gmx::c_centralShiftIndex)
                {
                    atomicAdd(&sm_fShiftLoc[ki], fij);
                    atomicAdd(&sm_fShiftLoc[gmx::c_centralShiftIndex], -fij);
                }
            }
        }
        storeForce(gm_f, ai, fij);
        storeForce(gm_f, aj, -fij);
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
__device__ void angles_gpu(const int       i,
                           float*          vtot_loc,
                           const int       numBonds,
                           const t_iatom   d_forceatoms[],
                           const t_iparams d_forceparams[],
                           const float4    gm_xq[],
                           float3          gm_f[],
                           float3          sm_fShiftLoc[],
                           const PbcAiuc   pbcAiuc)
{
    if (i < numBonds)
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
        harmonic_gpu<calcEner>(d_forceparams[type].harmonic.krA,
                               d_forceparams[type].harmonic.rA * HIP_DEG2RAD_F,
                               theta,
                               &va,
                               &dVdt);

        if constexpr (calcEner)
        {
            *vtot_loc += va;
        }

        float cos_theta2 = cos_theta * cos_theta;
        if (cos_theta2 < 1.0F)
        {
            float st    = dVdt * __frsqrt_rn(1.0F - cos_theta2);
            float sth   = st * cos_theta;
            float nrij2 = norm2(r_ij);
            float nrkj2 = norm2(r_kj);

            float nrij_1 = __frsqrt_rn(nrij2);
            float nrkj_1 = __frsqrt_rn(nrkj2);

            float cik = st * nrij_1 * nrkj_1;
            float cii = sth * nrij_1 * nrij_1;
            float ckk = sth * nrkj_1 * nrkj_1;

            float3 f_i = cii * r_ij - cik * r_kj;
            float3 f_k = ckk * r_kj - cik * r_ij;
            float3 f_j = -f_i - f_k;

            storeForce(gm_f, ai, f_i);
            storeForce(gm_f, aj, f_j);
            storeForce(gm_f, ak, f_k);

            if constexpr (calcVir)
            {
                atomicAdd(&sm_fShiftLoc[t1], f_i);
                atomicAdd(&sm_fShiftLoc[gmx::c_centralShiftIndex], f_j);
                atomicAdd(&sm_fShiftLoc[t2], f_k);
            }
        }
    }
}

template<bool calcVir, bool calcEner>
__device__ void urey_bradley_gpu(const int       i,
                                 float*          vtot_loc,
                                 const int       numBonds,
                                 const t_iatom   d_forceatoms[],
                                 const t_iparams d_forceparams[],
                                 const float4    gm_xq[],
                                 float3          gm_f[],
                                 float3          sm_fShiftLoc[],
                                 const PbcAiuc   pbcAiuc)
{
    if (i < numBonds)
    {
        const int4 ubData = *(reinterpret_cast<const int4*>(d_forceatoms + 4 * i));
        int        type   = ubData.x;
        int        ai     = ubData.y;
        int        aj     = ubData.z;
        int        ak     = ubData.w;

        float th0A = d_forceparams[type].u_b.thetaA * HIP_DEG2RAD_F;
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
        harmonic_gpu<calcEner>(kthA, th0A, theta, &va, &dVdt);

        if constexpr (calcEner)
        {
            *vtot_loc += va;
        }

        float3 r_ik;
        int    ki = pbcDxAiuc<calcVir>(pbcAiuc, gm_xq[ai], gm_xq[ak], r_ik);

        float dr2 = norm2(r_ik);
        float dr  = dr2 * __frsqrt_rn(dr2);

        float vbond;
        float fbond;
        harmonic_gpu<calcEner>(kUBA, r13A, dr, &vbond, &fbond);

        float3 f_i = make_float3(0.0F);
        float3 f_j = make_float3(0.0F);
        float3 f_k = make_float3(0.0F);

        float cos_theta2 = cos_theta * cos_theta;
        if (cos_theta2 < 1.0F)
        {
            float st  = dVdt * __frsqrt_rn(1.0F - cos_theta2);
            float sth = st * cos_theta;

            float nrkj2 = norm2(r_kj);
            float nrij2 = norm2(r_ij);

            float cik = st * __frsqrt_rn(nrkj2 * nrij2);
            float cii = sth / nrij2;
            float ckk = sth / nrkj2;

            f_i = cii * r_ij - cik * r_kj;
            f_k = ckk * r_kj - cik * r_ij;
            f_j = -f_i - f_k;


            if constexpr (calcVir)
            {
                atomicAdd(&sm_fShiftLoc[t1], f_i);
                atomicAdd(&sm_fShiftLoc[gmx::c_centralShiftIndex], f_j);
                atomicAdd(&sm_fShiftLoc[t2], f_k);
            }
        }

        /* Time for the bond calculations */
        if (dr2 != 0.0F)
        {
            if constexpr (calcEner)
            {
                *vtot_loc += vbond;
            }

            fbond *= __frsqrt_rn(dr2);

            float3 fik = fbond * r_ik;
            f_i += fik;
            f_k -= fik;

            if constexpr (calcVir)
            {
                if (ki != gmx::c_centralShiftIndex)
                {
                    atomicAdd(&sm_fShiftLoc[ki], fik);
                    atomicAdd(&sm_fShiftLoc[gmx::c_centralShiftIndex], -fik);
                }
            }
        }

        storeForce(gm_f, ai, f_i);
        storeForce(gm_f, aj, f_j);
        storeForce(gm_f, ak, f_k);
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


__device__ __forceinline__ static void
dopdihs_gpu(const float cpA, const float phiA, const int mult, const float phi, float* v, float* f)
{
    float mdphi, sdphi;

    mdphi = mult * phi - phiA * HIP_DEG2RAD_F;
    sdphi = __sinf(mdphi);
    *v    = cpA * (1.0F + __cosf(mdphi));
    *f    = -cpA * mult * sdphi;
}

template<bool calcVir>
__device__ static void do_dih_fup_gpu(const int            i,
                                      const int            j,
                                      const int            k,
                                      const int            l,
                                      const float          ddphi,
                                      const float3         r_ij,
                                      const float3         r_kj,
                                      const float3         r_kl,
                                      const float3         m,
                                      const float3         n,
                                      float3               gm_f[],
                                      float3               sm_fShiftLoc[],
                                      const PbcAiuc&       pbcAiuc,
                                      const float4         gm_xq[],
                                      const int            t1,
                                      const int            t2,
                                      const int gmx_unused t3)
{
    float iprm  = norm2(m);
    float iprn  = norm2(n);
    float nrkj2 = norm2(r_kj);
    float toler = nrkj2 * GMX_REAL_EPS;
    if ((iprm > toler) && (iprn > toler))
    {
        float  nrkj_1 = __frsqrt_rn(nrkj2); // replacing std::invsqrt call
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

        storeForce(gm_f, i, f_i);
        storeForce(gm_f, j, -f_j);
        storeForce(gm_f, k, -f_k);
        storeForce(gm_f, l, f_l);

        if constexpr (calcVir)
        {
            float3 dx_jl;
            int    t3 = pbcDxAiuc<calcVir>(pbcAiuc, gm_xq[l], gm_xq[j], dx_jl);

            atomicAdd(&sm_fShiftLoc[t1], f_i);
            atomicAdd(&sm_fShiftLoc[gmx::c_centralShiftIndex], -f_j);
            atomicAdd(&sm_fShiftLoc[t2], -f_k);
            atomicAdd(&sm_fShiftLoc[t3], f_l);
        }
    }
}

template<bool calcVir, bool calcEner>
__device__ void pdihs_gpu(const int       i,
                          float*          vtot_loc,
                          const int       numBonds,
                          const t_iatom   d_forceatoms[],
                          const t_iparams d_forceparams[],
                          const float4    gm_xq[],
                          float3          gm_f[],
                          float3          sm_fShiftLoc[],
                          const PbcAiuc   pbcAiuc)
{
    if (i < numBonds)
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

        if constexpr (calcEner)
        {
            *vtot_loc += vpd;
        }

        do_dih_fup_gpu<calcVir>(
                ai, aj, ak, al, ddphi, r_ij, r_kj, r_kl, m, n, gm_f, sm_fShiftLoc, pbcAiuc, gm_xq, t1, t2, t3);
    }
}

template<bool calcVir, bool calcEner>
__device__ void rbdihs_gpu(const int       i,
                           float*          vtot_loc,
                           const int       numBonds,
                           const t_iatom   d_forceatoms[],
                           const t_iparams d_forceparams[],
                           const float4    gm_xq[],
                           float3          gm_f[],
                           float3          sm_fShiftLoc[],
                           const PbcAiuc   pbcAiuc)
{
    constexpr float c0 = 0.0F, c1 = 1.0F, c2 = 2.0F, c3 = 3.0F, c4 = 4.0F, c5 = 5.0F;

    if (i < numBonds)
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

        /* Change to polymer convention */
        if (phi < c0)
        {
            phi += HIP_PI;
        }
        else
        {
            phi -= HIP_PI;
        }
        float cos_phi = __cosf(phi);
        /* Beware of accuracy loss, cannot use 1-sqrt(cos^2) ! */
        float sin_phi = __sinf(phi);

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
        if constexpr (calcEner)
        {
            v += cosfac * rbp;
        }
        rbp = parm[2];
        ddphi += c2 * rbp * cosfac;
        cosfac *= cos_phi;
        if constexpr (calcEner)
        {
            v += cosfac * rbp;
        }
        rbp = parm[3];
        ddphi += c3 * rbp * cosfac;
        cosfac *= cos_phi;
        if constexpr (calcEner)
        {
            v += cosfac * rbp;
        }
        rbp = parm[4];
        ddphi += c4 * rbp * cosfac;
        cosfac *= cos_phi;
        if constexpr (calcEner)
        {
            v += cosfac * rbp;
        }
        rbp = parm[5];
        ddphi += c5 * rbp * cosfac;
        cosfac *= cos_phi;
        if constexpr (calcEner)
        {
            v += cosfac * rbp;
        }

        ddphi = -ddphi * sin_phi;

        do_dih_fup_gpu<calcVir>(
                ai, aj, ak, al, ddphi, r_ij, r_kj, r_kl, m, n, gm_f, sm_fShiftLoc, pbcAiuc, gm_xq, t1, t2, t3);
        if constexpr (calcEner)
        {
            *vtot_loc += v;
        }
    }
}

__device__ __forceinline__ static void make_dp_periodic_gpu(float* dp)
{
    /* dp cannot be outside (-pi,pi) */
    if (*dp >= HIP_PI)
    {
        *dp -= 2.0F * HIP_PI;
    }
    else if (*dp < -HIP_PI)
    {
        *dp += 2.0F * HIP_PI;
    }
}

template<bool calcVir, bool calcEner>
__device__ void idihs_gpu(const int       i,
                          float*          vtot_loc,
                          const int       numBonds,
                          const t_iatom   d_forceatoms[],
                          const t_iparams d_forceparams[],
                          const float4    gm_xq[],
                          float3          gm_f[],
                          float3          sm_fShiftLoc[],
                          const PbcAiuc   pbcAiuc)
{
    if (i < numBonds)
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

        float phi0 = pA * HIP_DEG2RAD_F;

        float dp = phi - phi0;

        make_dp_periodic_gpu(&dp);

        float ddphi = -kA * dp;

        do_dih_fup_gpu<calcVir>(
                ai, aj, ak, al, -ddphi, r_ij, r_kj, r_kl, m, n, gm_f, sm_fShiftLoc, pbcAiuc, gm_xq, t1, t2, t3);

        if constexpr (calcEner)
        {
            *vtot_loc += -0.5F * ddphi * dp;
        }
    }
}

template<bool calcVir, bool calcEner>
__device__ void pairs_gpu(const int       i,
                          const int       numBonds,
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
    if (i < numBonds)
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
        float rinv  = __frsqrt_rn(r2);
        float rinv2 = rinv * rinv;
        float rinv6 = rinv2 * rinv2 * rinv2;

        /* Calculate the Coulomb force * r */
        float velec = scale_factor * qq * rinv;

        /* Calculate the LJ force * r and add it to the Coulomb part */
        float fr = (12.0F * c12 * rinv6 - 6.0F * c6) * rinv6 + velec;

        float  finvr = fr * rinv2;
        float3 f     = finvr * dr;

        /* Add the forces */
        storeForce(gm_f, ai, f);
        storeForce(gm_f, aj, -f);
        if constexpr (calcVir)
        {
            if (fshift_index != gmx::c_centralShiftIndex)
            {
                atomicAdd(&sm_fShiftLoc[fshift_index], f);
                atomicAdd(&sm_fShiftLoc[gmx::c_centralShiftIndex], -f);
            }
        }

        if constexpr (calcEner)
        {
            *vtotVdw_loc += (c12 * rinv6 - c6) * rinv6;
            *vtotElec_loc += velec;
        }
    }
}

namespace gmx
{

template<bool calcVir, bool calcEner>
__launch_bounds__(c_threadsBondedPerBlock) __global__ void bonded_kernel_gpu(
        //! Periodic boundary data
        PbcAiuc pbcAiuc,
        //! Scale factor
        float electrostaticsScaleFactor,
        //! The bonded types on GPU
        const FTypeArray<int> fTypesOnGpu,
        //! The number of bonds for every function type
        const FTypeArray<int> numFTypeBonds,
        //! The start index in the range of each interaction type
        const FTypeArray<int> fTypeRangeStart,
        //! The end index in the range of each interaction type
        const FTypeArray<int> fTypeRangeEnd,
        //! Force parameters (on GPU)
        t_iparams* d_forceParams,
        //! Coordinates before the timestep (on GPU)
        const float4* gm_xq,
        //! Forces on atoms (on GPU)
        float3* gm_f,
        //! Force shifts on atoms (on GPU)
        float3* gm_fShift,
        //! Total Energy (on GPU)
        float* d_vTot,
        //! Interaction list atoms (on GPU)
        const FTypeArray<t_iatom*> d_iatoms)
{
    GMX_DEVICE_ASSERT(blockDim.y == 1 && blockDim.z == 1);
    const int tid          = blockIdx.x * blockDim.x + threadIdx.x;
    float     vtot_loc     = 0.0F;
    float     vtotVdw_loc  = 0.0F;
    float     vtotElec_loc = 0.0F;

    extern __shared__ float3 sm_dynamicShmem[];
    float3*                  sm_fShiftLoc = sm_dynamicShmem;

    if constexpr (calcVir)
    {
        if (threadIdx.x < c_numShiftVectors)
        {
            sm_fShiftLoc[threadIdx.x] = make_float3(0.0F, 0.0F, 0.0F);
        }
        __syncthreads();
    }

    int fType_shared_index = -1;
#pragma unroll
    for (int j = 0; j < numFTypesOnGpu; j++)
    {
        const int      numBonds = numFTypeBonds[j];
        const int      fTypeTid = tid - fTypeRangeStart[j];
        const t_iatom* iatoms   = d_iatoms[j];
        const int      fType    = fTypesOnGpu[j];
        const int      start    = fTypeRangeStart[j];
        const int      end      = fTypeRangeEnd[j];
        if (tid >= start && tid <= end)
        {
            fType_shared_index = j;

            switch (fType)
            {
                case F_BONDS:
                    bonds_gpu<calcVir, calcEner>(
                            fTypeTid, &vtot_loc, numBonds, iatoms, d_forceParams, gm_xq, gm_f, sm_fShiftLoc, pbcAiuc);
                    break;
                case F_ANGLES:
                    angles_gpu<calcVir, calcEner>(
                            fTypeTid, &vtot_loc, numBonds, iatoms, d_forceParams, gm_xq, gm_f, sm_fShiftLoc, pbcAiuc);
                    break;
                case F_UREY_BRADLEY:
                    urey_bradley_gpu<calcVir, calcEner>(
                            fTypeTid, &vtot_loc, numBonds, iatoms, d_forceParams, gm_xq, gm_f, sm_fShiftLoc, pbcAiuc);
                    break;
                case F_PDIHS:
                case F_PIDIHS:
                    pdihs_gpu<calcVir, calcEner>(
                            fTypeTid, &vtot_loc, numBonds, iatoms, d_forceParams, gm_xq, gm_f, sm_fShiftLoc, pbcAiuc);
                    break;
                case F_RBDIHS:
                    rbdihs_gpu<calcVir, calcEner>(
                            fTypeTid, &vtot_loc, numBonds, iatoms, d_forceParams, gm_xq, gm_f, sm_fShiftLoc, pbcAiuc);
                    break;
                case F_IDIHS:
                    idihs_gpu<calcVir, calcEner>(
                            fTypeTid, &vtot_loc, numBonds, iatoms, d_forceParams, gm_xq, gm_f, sm_fShiftLoc, pbcAiuc);
                    break;
                case F_LJ14:
                    pairs_gpu<calcVir, calcEner>(fTypeTid,
                                                 numBonds,
                                                 iatoms,
                                                 d_forceParams,
                                                 gm_xq,
                                                 gm_f,
                                                 sm_fShiftLoc,
                                                 pbcAiuc,
                                                 electrostaticsScaleFactor,
                                                 &vtotVdw_loc,
                                                 &vtotElec_loc);
                    break;
            }
            break;
        }
    }

    if constexpr (calcEner)
    {
#pragma unroll
        for (int j = 0; j < numFTypesOnGpu; j++)
        {
            int fType = fTypesOnGpu[j];
            if (__any(j == fType_shared_index))
            {
                float vtot_shuffle = j == fType_shared_index ? vtot_loc : 0.0f;
#pragma unroll
                for (unsigned int offset = (warpSize >> 1); offset > 0; offset >>= 1)
                {
                    vtot_shuffle += __shfl_down(vtot_shuffle, offset);
                }
                if ((threadIdx.x & (warpSize - 1)) == 0)
                {
                    atomicAdd((d_vTot + fType), vtot_shuffle);
                }
            }
        }

        float vtotVdw_shuffle  = vtotVdw_loc;
        float vtotElec_shuffle = vtotElec_loc;
#pragma unroll
        for (unsigned int offset = (warpSize >> 1); offset > 0; offset >>= 1)
        {
            vtotVdw_shuffle += __shfl_down(vtotVdw_shuffle, offset);
            vtotElec_shuffle += __shfl_down(vtotElec_shuffle, offset);
        }

        if ((threadIdx.x & (warpSize - 1)) == 0)
        { // One thread per warp accumulates partial sum into global sum
            atomicAdd(d_vTot + F_LJ14, vtotVdw_shuffle);
            atomicAdd(d_vTot + F_COUL14, vtotElec_shuffle);
        }
    }
    /* Accumulate shift vectors from shared memory to global memory on the first c_numShiftVectors threads of the block. */
    if constexpr (calcVir)
    {
        __syncthreads();
        if (threadIdx.x < c_numShiftVectors)
        {
            atomicAdd(&gm_fShift[threadIdx.x], sm_fShiftLoc[threadIdx.x]);
        }
    }
}


/*-------------------------------- End HIP kernels-----------------------------*/


template<bool calcVir, bool calcEner>
void ListedForcesGpu::Impl::launchKernel()
{
    GMX_ASSERT(haveInteractions_,
               "Cannot launch bonded GPU kernels unless bonded GPU work was scheduled");

    wallcycle_start_nocount(wcycle_, WallCycleCounter::LaunchGpuPp);
    wallcycle_sub_start(wcycle_, WallCycleSubCounter::LaunchGpuBonded);

    if (kernelParams_.fTypeRangeEnd[numFTypesOnGpu - 1] < 0)
    {
        return;
    }

    auto kernelPtr = bonded_kernel_gpu<calcVir, calcEner>;

    FTypeArray<int>      fTypesOnGpu(kernelParams_.fTypesOnGpu);
    FTypeArray<int>      numFTypeBonds(kernelParams_.numFTypeBonds);
    FTypeArray<int>      fTypeRangeStart(kernelParams_.fTypeRangeStart);
    FTypeArray<int>      fTypeRangeEnd(kernelParams_.fTypeRangeEnd);
    FTypeArray<t_iatom*> d_iatoms(kernelBuffers_.d_iatoms);
    const auto           kernelArgs = prepareGpuKernelArguments(kernelPtr,
                                                      kernelLaunchConfig_,
                                                      &kernelParams_.pbcAiuc,
                                                      &kernelParams_.electrostaticsScaleFactor,
                                                      &fTypesOnGpu,
                                                      &numFTypeBonds,
                                                      &fTypeRangeStart,
                                                      &fTypeRangeEnd,
                                                      &kernelBuffers_.d_forceParams,
                                                      &d_xq_,
                                                      &d_f_,
                                                      &d_fShift_,
                                                      &kernelBuffers_.d_vTot,
                                                      &d_iatoms);

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
