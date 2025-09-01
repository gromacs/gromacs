/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
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
 * \brief Implements generic GPU bonded functionality
 *
 * \author Andrey Alekseenko <al42and@gmail.com>
 * \author Jon Vincent <jvincent@nvidia.com>
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 * \author Berk Hess <hess@kth.se>
 * \author Szilárd Páll <pall.szilard@gmail.com>
 * \author Alan Gray <alang@nvidia.com>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \ingroup module_listed_forces
 */
#ifndef GMX_LISTED_FORCES_LISTED_FORCES_GPU_INTERNAL_SHARED_H
#define GMX_LISTED_FORCES_LISTED_FORCES_GPU_INTERNAL_SHARED_H

#include "config.h"

#if GMX_GPU_CUDA
#    include "gromacs/gpu_utils/cuda_kernel_utils.cuh"
#    include "gromacs/pbcutil/pbc_aiuc_cuda.cuh"
#elif GMX_GPU_SYCL
#    include "gromacs/gpu_utils/sycl_kernel_utils.h"
#    include "gromacs/gpu_utils/vectype_ops_sycl.h"
#    include "gromacs/pbcutil/pbc_aiuc_sycl.h"
#else
#    error Building GPU bonded for unsupported backend
#endif
#include "gromacs/gpu_utils/gputraits.h"
#include "gromacs/listed_forces/listed_forces_gpu.h"
#include "gromacs/math/units.h"

#include "listed_forces_gpu_impl.h"

#ifndef DOXYGEN

namespace gmx
{

static constexpr float c_deg2RadF = gmx::c_deg2Rad;
constexpr float        c_Pi       = M_PI;

/* Some SYCL targets have troubles optimizing the dynamic array
 * member access despite the fact that all the loops are unrolled.
 *
 * See https://developer.nvidia.com/blog/fast-dynamic-indexing-private-arrays-cuda/
 * for a details on why dynamic access is problematic.
 *
 * This seems to affect:
 * - hipSYCL 0.9.4 + Clang 14-15 for AMD (but not NVIDIA),
 * - IntelLLVM 2023-02 for NVIDIA and Arc (but to a much lesser extent PVC).
 *
 * This wrapper avoid dynamic accesses into the array, replacing them
 * with a `switch` instead.
 *
 * Based on the optimization by AMD/StreamHPC for their HIP port.
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
    constexpr T operator[](int idx) const
    {
        switch (idx)
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

/* Harmonic */
template<bool calcEner>
GMX_DEVICE_FUNC_ATTRIBUTE static inline void harmonic_gpu(const float             kA,
                                                          const float             xA,
                                                          const float             x,
                                                          DevicePrivatePtr<float> V,
                                                          DevicePrivatePtr<float> F)
{
    constexpr float half = 0.5F;
    float           dx   = x - xA;
    float           dx2  = dx * dx;

    *F = -kA * dx;
    if constexpr (calcEner)
    {
        *V = half * kA * dx2;
    }
}

template<bool calcVir, bool calcEner>
GMX_DEVICE_FUNC_ATTRIBUTE static inline void bonds_gpu(const int               i,
                                                       DevicePrivatePtr<float> vtot_loc,
                                                       const DeviceGlobalPtr<const t_iatom> gm_forceatoms,
                                                       const DeviceGlobalPtr<const t_iparams> gm_forceparams,
                                                       const DeviceGlobalPtr<const DeviceFloat4> gm_xq,
                                                       DeviceGlobalPtr<DeviceFloat3> gm_f,
                                                       DeviceLocalPtr<DeviceFloat3>  sm_fShiftLoc,
                                                       const PbcAiuc&                pbcAiuc,
                                                       const int                     localId)
{
    const int type = gm_forceatoms[3 * i];
    const int ai   = gm_forceatoms[3 * i + 1];
    const int aj   = gm_forceatoms[3 * i + 2];

    /* dx = xi - xj, corrected for periodic boundary conditions. */
    DeviceFloat3 dx;
    int          ki = pbcDxAiucGpu<calcVir>(pbcAiuc, gm_xq[ai], gm_xq[aj], dx);

    float dr2 = gmxDeviceNorm2(dx);
    float dr  = gmxDeviceSqrt(dr2);

    float vbond;
    float fbond;
    harmonic_gpu<calcEner>(
            gm_forceparams[type].harmonic.krA, gm_forceparams[type].harmonic.rA, dr, &vbond, &fbond);

    if constexpr (calcEner)
    {
        *vtot_loc += vbond;
    }

    if (dr2 != 0.0F)
    {
        fbond *= gmxDeviceRSqrt(dr2);

        DeviceFloat3 fij = fbond * dx;
        staggeredAtomicAddForce(&gm_f[ai], fij, localId);
        staggeredAtomicAddForce(&gm_f[aj], -fij, localId);
        if constexpr (calcVir)
        {
            if (ki != gmx::c_centralShiftIndex)
            {
                atomicFetchAddLocal(&sm_fShiftLoc[ki], fij);
                atomicFetchAddLocal(&sm_fShiftLoc[gmx::c_centralShiftIndex], -fij);
            }
        }
    }
}

template<bool returnShift>
GMX_DEVICE_FUNC_ATTRIBUTE static inline float bond_angle_gpu(const DeviceFloat4             xi,
                                                             const DeviceFloat4             xj,
                                                             const DeviceFloat4             xk,
                                                             const PbcAiuc&                 pbcAiuc,
                                                             DevicePrivatePtr<DeviceFloat3> r_ij,
                                                             DevicePrivatePtr<DeviceFloat3> r_kj,
                                                             DevicePrivatePtr<float>        costh,
                                                             DevicePrivatePtr<int>          t1,
                                                             DevicePrivatePtr<int>          t2)
{
    *t1 = pbcDxAiucGpu<returnShift>(pbcAiuc, xi, xj, *r_ij);
    *t2 = pbcDxAiucGpu<returnShift>(pbcAiuc, xk, xj, *r_kj);

    *costh = gmxDeviceCosAngle(*r_ij, *r_kj);
    // Return value is the angle between the bonds i-j and j-k
    return gmxDeviceAcos(*costh);
}

template<bool calcVir, bool calcEner>
GMX_DEVICE_FUNC_ATTRIBUTE static inline void angles_gpu(const int               i,
                                                        DevicePrivatePtr<float> vtot_loc,
                                                        const DeviceGlobalPtr<const t_iatom> gm_forceatoms,
                                                        const DeviceGlobalPtr<const t_iparams> gm_forceparams,
                                                        const DeviceGlobalPtr<const DeviceFloat4> gm_xq,
                                                        DeviceGlobalPtr<DeviceFloat3> gm_f,
                                                        DeviceLocalPtr<DeviceFloat3>  sm_fShiftLoc,
                                                        const PbcAiuc&                pbcAiuc,
                                                        const int                     localId)
{
    DeviceInt4 angleData = loadInt4(gm_forceatoms, i);
    const int  type      = angleData[0];
    const int  ai        = angleData[1];
    const int  aj        = angleData[2];
    const int  ak        = angleData[3];

    DeviceFloat3 r_ij;
    DeviceFloat3 r_kj;
    float        cos_theta;
    int          t1;
    int          t2;
    float        theta = bond_angle_gpu<calcVir>(
            gm_xq[ai], gm_xq[aj], gm_xq[ak], pbcAiuc, &r_ij, &r_kj, &cos_theta, &t1, &t2);

    float va;
    float dVdt;
    harmonic_gpu<calcEner>(
            gm_forceparams[type].harmonic.krA, gm_forceparams[type].harmonic.rA * c_deg2RadF, theta, &va, &dVdt);

    if constexpr (calcEner)
    {
        *vtot_loc += va;
    }

    float cos_theta2 = cos_theta * cos_theta;
    if (cos_theta2 < 1.0F)
    {
        float st    = dVdt * gmxDeviceRSqrt(1.0F - cos_theta2);
        float sth   = st * cos_theta;
        float nrij2 = gmxDeviceNorm2(r_ij);
        float nrkj2 = gmxDeviceNorm2(r_kj);

        float nrij_1 = gmxDeviceRSqrt(nrij2);
        float nrkj_1 = gmxDeviceRSqrt(nrkj2);

        float cik = st * nrij_1 * nrkj_1;
        float cii = sth * nrij_1 * nrij_1;
        float ckk = sth * nrkj_1 * nrkj_1;

        DeviceFloat3 f_i = cii * r_ij - cik * r_kj;
        DeviceFloat3 f_k = ckk * r_kj - cik * r_ij;
        DeviceFloat3 f_j = -f_i - f_k;

        staggeredAtomicAddForce(&gm_f[ai], f_i, localId);
        staggeredAtomicAddForce(&gm_f[aj], f_j, localId);
        staggeredAtomicAddForce(&gm_f[ak], f_k, localId);

        if constexpr (calcVir)
        {
            atomicFetchAddLocal(&sm_fShiftLoc[t1], f_i);
            atomicFetchAddLocal(&sm_fShiftLoc[gmx::c_centralShiftIndex], f_j);
            atomicFetchAddLocal(&sm_fShiftLoc[t2], f_k);
        }
    }
}

template<bool calcVir, bool calcEner>
GMX_DEVICE_FUNC_ATTRIBUTE static inline void
urey_bradley_gpu(const int                                 i,
                 DevicePrivatePtr<float>                   vtot_loc,
                 const DeviceGlobalPtr<const t_iatom>      gm_forceatoms,
                 const DeviceGlobalPtr<const t_iparams>    gm_forceparams,
                 const DeviceGlobalPtr<const DeviceFloat4> gm_xq,
                 DeviceGlobalPtr<DeviceFloat3>             gm_f,
                 DeviceLocalPtr<DeviceFloat3>              sm_fShiftLoc,
                 const PbcAiuc&                            pbcAiuc,
                 const int                                 localId)
{
    DeviceInt4 ubData = loadInt4(gm_forceatoms, i);
    const int  type   = ubData[0];
    const int  ai     = ubData[1];
    const int  aj     = ubData[2];
    const int  ak     = ubData[3];

    const float th0A = gm_forceparams[type].u_b.thetaA * c_deg2RadF;
    const float kthA = gm_forceparams[type].u_b.kthetaA;
    const float r13A = gm_forceparams[type].u_b.r13A;
    const float kUBA = gm_forceparams[type].u_b.kUBA;

    DeviceFloat3 r_ij;
    DeviceFloat3 r_kj;
    float        cos_theta;
    int          t1;
    int          t2;
    float        theta = bond_angle_gpu<calcVir>(
            gm_xq[ai], gm_xq[aj], gm_xq[ak], pbcAiuc, &r_ij, &r_kj, &cos_theta, &t1, &t2);

    float va;
    float dVdt;
    harmonic_gpu<calcEner>(kthA, th0A, theta, &va, &dVdt);

    if (calcEner)
    {
        *vtot_loc += va;
    }

    DeviceFloat3 r_ik;
    int          ki = pbcDxAiucGpu<calcVir>(pbcAiuc, gm_xq[ai], gm_xq[ak], r_ik);

    float dr2 = gmxDeviceNorm2(r_ik);
    float dr  = dr2 * gmxDeviceRSqrt(dr2);

    float vbond;
    float fbond;
    harmonic_gpu<calcEner>(kUBA, r13A, dr, &vbond, &fbond);

    float cos_theta2 = cos_theta * cos_theta;

    DeviceFloat3 f_i = { 0.0F, 0.0F, 0.0F };
    DeviceFloat3 f_j = { 0.0F, 0.0F, 0.0F };
    DeviceFloat3 f_k = { 0.0F, 0.0F, 0.0F };

    if (cos_theta2 < 1.0F)
    {
        float st  = dVdt * gmxDeviceRSqrt(1.0F - cos_theta2);
        float sth = st * cos_theta;

        float nrkj2 = gmxDeviceNorm2(r_kj);
        float nrij2 = gmxDeviceNorm2(r_ij);

        float cik = st * gmxDeviceRSqrt(nrkj2 * nrij2);
        float cii = sth / nrij2;
        float ckk = sth / nrkj2;

        f_i = cii * r_ij - cik * r_kj;
        f_k = ckk * r_kj - cik * r_ij;
        f_j = -f_i - f_k;

        if constexpr (calcVir)
        {
            atomicFetchAddLocal(&sm_fShiftLoc[t1], f_i);
            atomicFetchAddLocal(&sm_fShiftLoc[gmx::c_centralShiftIndex], f_j);
            atomicFetchAddLocal(&sm_fShiftLoc[t2], f_k);
        }
    }

    /* Time for the bond calculations */
    if (dr2 != 0.0F)
    {
        if constexpr (calcEner)
        {
            *vtot_loc += vbond;
        }

        fbond *= gmxDeviceRSqrt(dr2);
        DeviceFloat3 fik = fbond * r_ik;
        f_i += fik;
        f_k -= fik;


        if constexpr (calcVir)
        {
            if (ki != gmx::c_centralShiftIndex)
            {
                atomicFetchAddLocal(&sm_fShiftLoc[ki], fik);
                atomicFetchAddLocal(&sm_fShiftLoc[gmx::c_centralShiftIndex], -fik);
            }
        }
    }

    if ((cos_theta2 < 1.0F) || (dr2 != 0.0F))
    {
        staggeredAtomicAddForce(&gm_f[ai], f_i, localId);
        staggeredAtomicAddForce(&gm_f[ak], f_k, localId);
    }

    if (cos_theta2 < 1.0F)
    {
        staggeredAtomicAddForce(&gm_f[aj], f_j, localId);
    }
}

template<bool returnShift, typename T>
GMX_DEVICE_FUNC_ATTRIBUTE static inline float dih_angle_gpu(const T                        xi,
                                                            const T                        xj,
                                                            const T                        xk,
                                                            const T                        xl,
                                                            const PbcAiuc&                 pbcAiuc,
                                                            DevicePrivatePtr<DeviceFloat3> r_ij,
                                                            DevicePrivatePtr<DeviceFloat3> r_kj,
                                                            DevicePrivatePtr<DeviceFloat3> r_kl,
                                                            DevicePrivatePtr<DeviceFloat3> m,
                                                            DevicePrivatePtr<DeviceFloat3> n,
                                                            DevicePrivatePtr<int>          t1,
                                                            DevicePrivatePtr<int>          t2,
                                                            DevicePrivatePtr<int>          t3)
{
    *t1 = pbcDxAiucGpu<returnShift>(pbcAiuc, xi, xj, *r_ij);
    *t2 = pbcDxAiucGpu<returnShift>(pbcAiuc, xk, xj, *r_kj);
    *t3 = pbcDxAiucGpu<returnShift>(pbcAiuc, xk, xl, *r_kl);

    *m         = gmxDeviceCrossProd(*r_ij, *r_kj);
    *n         = gmxDeviceCrossProd(*r_kj, *r_kl);
    float phi  = gmxDeviceAngle(*m, *n);
    float ipr  = gmxDeviceInternalProd(*r_ij, *n);
    float sign = (ipr < 0.0F) ? -1.0F : 1.0F;
    phi        = sign * phi;

    return phi;
}

template<bool returnShift, typename T>
GMX_DEVICE_FUNC_ATTRIBUTE static inline float dih_angle_gpu_sincos(const T        xi,
                                                                   const T        xj,
                                                                   const T        xk,
                                                                   const T        xl,
                                                                   const PbcAiuc& pbcAiuc,
                                                                   DevicePrivatePtr<DeviceFloat3> r_ij,
                                                                   DevicePrivatePtr<DeviceFloat3> r_kj,
                                                                   DevicePrivatePtr<DeviceFloat3> r_kl,
                                                                   DevicePrivatePtr<DeviceFloat3> m,
                                                                   DevicePrivatePtr<DeviceFloat3> n,
                                                                   DevicePrivatePtr<int>   t1,
                                                                   DevicePrivatePtr<int>   t2,
                                                                   DevicePrivatePtr<float> cosval)
{
    *t1 = pbcDxAiucGpu<returnShift>(pbcAiuc, xi, xj, *r_ij);
    *t2 = pbcDxAiucGpu<returnShift>(pbcAiuc, xk, xj, *r_kj);
    pbcDxAiucGpu<returnShift>(pbcAiuc, xk, xl, *r_kl);

    *m = gmxDeviceCrossProd(*r_ij, *r_kj);
    *n = gmxDeviceCrossProd(*r_kj, *r_kl);

    DeviceFloat3 w = gmxDeviceCrossProd(*m, *n);

    float wlen = gmxDeviceSqrt(gmxDeviceNorm2(w));
    float s    = gmxDeviceInternalProd(*m, *n);

    float mLenSq = gmxDeviceNorm2(*m);
    float nLenSq = gmxDeviceNorm2(*n);
    float mnInv  = gmxDeviceRSqrt(mLenSq * nLenSq);

    *cosval      = s * mnInv;
    float sinval = wlen * mnInv;

    float ipr  = gmxDeviceInternalProd(*r_ij, *n);
    float sign = (ipr < 0.0F) ? -1.0F : 1.0F;

    return sign * sinval;
}

GMX_DEVICE_FUNC_ATTRIBUTE
static inline void dopdihs_gpu(const float             cpA,
                               const float             phiA,
                               const int               mult,
                               const float             phi,
                               DevicePrivatePtr<float> v,
                               DevicePrivatePtr<float> f)
{
    float mdphi = mult * phi - phiA * c_deg2RadF;
    float sdphi = gmxDeviceSin(mdphi);
    *v          = cpA * (1.0F + gmxDeviceCos(mdphi));
    *f          = -cpA * mult * sdphi;
}

template<bool calcVir>
GMX_DEVICE_FUNC_ATTRIBUTE static inline void do_dih_fup_gpu(const int                     i,
                                                            const int                     j,
                                                            const int                     k,
                                                            const int                     l,
                                                            const float                   ddphi,
                                                            const DeviceFloat3            r_ij,
                                                            const DeviceFloat3            r_kj,
                                                            const DeviceFloat3            r_kl,
                                                            const DeviceFloat3            m,
                                                            const DeviceFloat3            n,
                                                            DeviceGlobalPtr<DeviceFloat3> gm_f,
                                                            DeviceLocalPtr<DeviceFloat3> sm_fShiftLoc,
                                                            const PbcAiuc& pbcAiuc,
                                                            const DeviceGlobalPtr<const DeviceFloat4> gm_xq,
                                                            const int t1,
                                                            const int t2,
                                                            const int localId)
{
    float iprm  = gmxDeviceNorm2(m);
    float iprn  = gmxDeviceNorm2(n);
    float nrkj2 = gmxDeviceNorm2(r_kj);
    float toler = nrkj2 * GMX_REAL_EPS;
    if ((iprm > toler) && (iprn > toler))
    {
        float        nrkj_1 = gmxDeviceRSqrt(nrkj2);
        float        nrkj_2 = nrkj_1 * nrkj_1;
        float        nrkj   = nrkj2 * nrkj_1;
        float        a      = -ddphi * nrkj / iprm;
        DeviceFloat3 f_i    = a * m;
        float        b      = ddphi * nrkj / iprn;
        DeviceFloat3 f_l    = b * n;
        float        p      = gmxDeviceInternalProd(r_ij, r_kj);
        p *= nrkj_2;
        float q = gmxDeviceInternalProd(r_kl, r_kj);
        q *= nrkj_2;
        DeviceFloat3 uvec = p * f_i;
        DeviceFloat3 vvec = q * f_l;
        DeviceFloat3 svec = uvec - vvec;
        DeviceFloat3 f_j  = f_i - svec;
        DeviceFloat3 f_k  = f_l + svec;

        staggeredAtomicAddForce(&gm_f[i], f_i, localId);
        staggeredAtomicAddForce(&gm_f[j], -f_j, localId);
        staggeredAtomicAddForce(&gm_f[k], -f_k, localId);
        staggeredAtomicAddForce(&gm_f[l], f_l, localId);

        if constexpr (calcVir)
        {
            DeviceFloat3 dx_jl;
            int          t3 = pbcDxAiucGpu<calcVir>(pbcAiuc, gm_xq[l], gm_xq[j], dx_jl);

            atomicFetchAddLocal(&sm_fShiftLoc[t1], f_i);
            atomicFetchAddLocal(&sm_fShiftLoc[gmx::c_centralShiftIndex], -f_j);
            atomicFetchAddLocal(&sm_fShiftLoc[t2], -f_k);
            atomicFetchAddLocal(&sm_fShiftLoc[t3], f_l);
        }
    }
}

template<bool calcVir, bool calcEner>
GMX_DEVICE_FUNC_ATTRIBUTE static inline void pdihs_gpu(const int               i,
                                                       DevicePrivatePtr<float> vtot_loc,
                                                       const DeviceGlobalPtr<const t_iatom> gm_forceatoms,
                                                       const DeviceGlobalPtr<const t_iparams> gm_forceparams,
                                                       const DeviceGlobalPtr<const DeviceFloat4> gm_xq,
                                                       DeviceGlobalPtr<DeviceFloat3> gm_f,
                                                       DeviceLocalPtr<DeviceFloat3>  sm_fShiftLoc,
                                                       const PbcAiuc&                pbcAiuc,
                                                       const int                     localId)
{
    int type = gm_forceatoms[5 * i];
    int ai   = gm_forceatoms[5 * i + 1];
    int aj   = gm_forceatoms[5 * i + 2];
    int ak   = gm_forceatoms[5 * i + 3];
    int al   = gm_forceatoms[5 * i + 4];

    DeviceFloat3 r_ij;
    DeviceFloat3 r_kj;
    DeviceFloat3 r_kl;
    DeviceFloat3 m;
    DeviceFloat3 n;
    int          t1;
    int          t2;
    int          t3;
    float        phi = dih_angle_gpu<calcVir>(
            gm_xq[ai], gm_xq[aj], gm_xq[ak], gm_xq[al], pbcAiuc, &r_ij, &r_kj, &r_kl, &m, &n, &t1, &t2, &t3);

    float vpd;
    float ddphi;
    dopdihs_gpu(gm_forceparams[type].pdihs.cpA,
                gm_forceparams[type].pdihs.phiA,
                gm_forceparams[type].pdihs.mult,
                phi,
                &vpd,
                &ddphi);

    if constexpr (calcEner)
    {
        *vtot_loc += vpd;
    }

    do_dih_fup_gpu<calcVir>(
            ai, aj, ak, al, ddphi, r_ij, r_kj, r_kl, m, n, gm_f, sm_fShiftLoc, pbcAiuc, gm_xq, t1, t2, localId);
}

template<bool calcVir, bool calcEner>
GMX_DEVICE_FUNC_ATTRIBUTE static inline void rbdihs_gpu(const int               i,
                                                        DevicePrivatePtr<float> vtot_loc,
                                                        const DeviceGlobalPtr<const t_iatom> gm_forceatoms,
                                                        const DeviceGlobalPtr<const t_iparams> gm_forceparams,
                                                        const DeviceGlobalPtr<const DeviceFloat4> gm_xq,
                                                        DeviceGlobalPtr<DeviceFloat3> gm_f,
                                                        DeviceLocalPtr<DeviceFloat3>  sm_fShiftLoc,
                                                        const PbcAiuc&                pbcAiuc,
                                                        const int                     localId)
{
    constexpr float c0 = 0.0F, c1 = 1.0F, c2 = 2.0F, c3 = 3.0F, c4 = 4.0F, c5 = 5.0F;

    {
        int type = gm_forceatoms[5 * i];
        int ai   = gm_forceatoms[5 * i + 1];
        int aj   = gm_forceatoms[5 * i + 2];
        int ak   = gm_forceatoms[5 * i + 3];
        int al   = gm_forceatoms[5 * i + 4];

        DeviceFloat3 r_ij;
        DeviceFloat3 r_kj;
        DeviceFloat3 r_kl;
        DeviceFloat3 m;
        DeviceFloat3 n;
        int          t1;
        int          t2;

        // Changing the sign of sin and cos to convert to polymer convention
        float       cos_phi;
        const float negative_sin_phi = dih_angle_gpu_sincos<calcVir>(
                gm_xq[ai], gm_xq[aj], gm_xq[ak], gm_xq[al], pbcAiuc, &r_ij, &r_kj, &r_kl, &m, &n, &t1, &t2, &cos_phi);
        cos_phi *= -1;

        float parm[NR_RBDIHS];
#    pragma unroll NR_RBDIHS
        for (int j = 0; j < NR_RBDIHS; j++)
        {
            parm[j] = gm_forceparams[type].rbdihs.rbcA[j];
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

        ddphi = ddphi * negative_sin_phi;

        do_dih_fup_gpu<calcVir>(
                ai, aj, ak, al, ddphi, r_ij, r_kj, r_kl, m, n, gm_f, sm_fShiftLoc, pbcAiuc, gm_xq, t1, t2, localId);
        if constexpr (calcEner)
        {
            *vtot_loc += v;
        }
    }
}

//! Wrap angle from range [-3*pi; 3*pi) to [-pi; pi)
GMX_DEVICE_FUNC_ATTRIBUTE
static constexpr float wrapAngle(float a)
{
    constexpr float c_twoPi = 2.0F * c_Pi;
    if (a >= c_Pi)
    {
        return a - c_twoPi;
    }
    else if (a < -c_Pi)
    {
        return a + c_twoPi;
    }
    else
    {
        return a;
    }
}

template<bool calcVir, bool calcEner>
GMX_DEVICE_FUNC_ATTRIBUTE static inline void idihs_gpu(const int               i,
                                                       DevicePrivatePtr<float> vtot_loc,
                                                       const DeviceGlobalPtr<const t_iatom> gm_forceatoms,
                                                       const DeviceGlobalPtr<const t_iparams> gm_forceparams,
                                                       const DeviceGlobalPtr<const DeviceFloat4> gm_xq,
                                                       DeviceGlobalPtr<DeviceFloat3> gm_f,
                                                       DeviceLocalPtr<DeviceFloat3>  sm_fShiftLoc,
                                                       const PbcAiuc&                pbcAiuc,
                                                       const int                     localId)
{
    int type = gm_forceatoms[5 * i];
    int ai   = gm_forceatoms[5 * i + 1];
    int aj   = gm_forceatoms[5 * i + 2];
    int ak   = gm_forceatoms[5 * i + 3];
    int al   = gm_forceatoms[5 * i + 4];

    DeviceFloat3 r_ij;
    DeviceFloat3 r_kj;
    DeviceFloat3 r_kl;
    DeviceFloat3 m;
    DeviceFloat3 n;
    int          t1;
    int          t2;
    int          t3;
    float        phi = dih_angle_gpu<calcVir>(
            gm_xq[ai], gm_xq[aj], gm_xq[ak], gm_xq[al], pbcAiuc, &r_ij, &r_kj, &r_kl, &m, &n, &t1, &t2, &t3);

    /* phi can jump if phi0 is close to Pi/-Pi, which will cause huge
     * force changes if we just apply a normal harmonic.
     * Instead, we first calculate phi-phi0 and take it modulo (-Pi,Pi).
     * This means we will never have the periodicity problem, unless
     * the dihedral is Pi away from phiO, which is very unlikely due to
     * the potential.
     */
    float kA = gm_forceparams[type].harmonic.krA;
    float pA = gm_forceparams[type].harmonic.rA;

    float phi0 = pA * c_deg2RadF;

    float dp = wrapAngle(phi - phi0);

    float ddphi = -kA * dp;

    do_dih_fup_gpu<calcVir>(
            ai, aj, ak, al, -ddphi, r_ij, r_kj, r_kl, m, n, gm_f, sm_fShiftLoc, pbcAiuc, gm_xq, t1, t2, localId);

    if constexpr (calcEner)
    {
        *vtot_loc += -0.5F * ddphi * dp;
    }
}

template<bool calcVir, bool calcEner>
GMX_DEVICE_FUNC_ATTRIBUTE static inline void pairs_gpu(const int i,
                                                       const DeviceGlobalPtr<const t_iatom> gm_forceatoms,
                                                       const DeviceGlobalPtr<const t_iparams> gm_iparams,
                                                       const DeviceGlobalPtr<const DeviceFloat4> gm_xq,
                                                       DeviceGlobalPtr<DeviceFloat3> gm_f,
                                                       DeviceLocalPtr<DeviceFloat3>  sm_fShiftLoc,
                                                       const PbcAiuc&                pbcAiuc,
                                                       const float                   scale_factor,
                                                       DevicePrivatePtr<float>       vtot_loc,
                                                       DevicePrivatePtr<float>       vtotElec_loc,
                                                       const int                     localId)
{
    int type = gm_forceatoms[3 * i];
    int ai   = gm_forceatoms[3 * i + 1];
    int aj   = gm_forceatoms[3 * i + 2];

    float qq  = gm_xq[ai][3] * gm_xq[aj][3];
    float c6  = gm_iparams[type].lj14.c6A;
    float c12 = gm_iparams[type].lj14.c12A;

    /* Do we need to apply full periodic boundary conditions? */
    DeviceFloat3 dr;
    int          fshift_index = pbcDxAiucGpu<calcVir>(pbcAiuc, gm_xq[ai], gm_xq[aj], dr);

    float r2    = gmxDeviceNorm2(dr);
    float rinv  = gmxDeviceRSqrt(r2);
    float rinv2 = rinv * rinv;
    float rinv6 = rinv2 * rinv2 * rinv2;

    /* Calculate the Coulomb force * r */
    float velec = scale_factor * qq * rinv;

    /* Calculate the LJ force * r and add it to the Coulomb part */
    float fr = (12.0F * c12 * rinv6 - 6.0F * c6) * rinv6 + velec;

    float        finvr = fr * rinv2;
    DeviceFloat3 f     = finvr * dr;

    /* Add the forces */
    staggeredAtomicAddForce(&gm_f[ai], f, localId);
    staggeredAtomicAddForce(&gm_f[aj], -f, localId);
    if constexpr (calcVir)
    {
        if (fshift_index != gmx::c_centralShiftIndex)
        {
            atomicFetchAddLocal(&sm_fShiftLoc[fshift_index], f);
            atomicFetchAddLocal(&sm_fShiftLoc[gmx::c_centralShiftIndex], -f);
        }
    }

    // The elec and vdW contributions to the energy are separated only for the pairs
    // code, and combined later on.
    if constexpr (calcEner)
    {
        *vtot_loc += (c12 * rinv6 - c6) * rinv6;
        *vtotElec_loc += velec;
    }
}

} // namespace gmx

#endif

#endif
