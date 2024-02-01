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
 * \brief Implements SYCL bonded functionality
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

#include "gmxpre.h"

#include "gromacs/gpu_utils/devicebuffer_sycl.h"
#include "gromacs/gpu_utils/gmxsycl.h"
#include "gromacs/gpu_utils/sycl_kernel_utils.h"
#include "gromacs/listed_forces/listed_forces_gpu.h"
#include "gromacs/math/units.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/pbcutil/pbc_aiuc_sycl.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/gmxassert.h"

#include "listed_forces_gpu_impl.h"

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
        __builtin_assume(idx >= 0 && idx < gmx::numFTypesOnGpu);
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

// \brief Staggered atomic force component accummulation into global memory to reduce clashes
//
// Reduce the number of atomic clashes by a theoretical max 3x by having consecutive threads
// accumulate different force components at the same time.
static inline void staggeredAtomicAddForce(sycl::global_ptr<Float3> gm_f, Float3 f, const int localId)
{
    __builtin_assume(localId >= 0);

    using Int3  = sycl::int3;
    Int3 offset = { 0, 1, 2 };

    // Shift force components x (0), y (1), and z (2) left by 2, 1, and 0, respectively
    // to end up with zxy, yzx, xyz on consecutive threads.
    f      = (localId % 3 == 0) ? Float3(f[1], f[2], f[0]) : f;
    offset = (localId % 3 == 0) ? Int3(offset[1], offset[2], offset[0]) : offset;
    f      = (localId % 3 <= 1) ? Float3(f[1], f[2], f[0]) : f;
    offset = (localId % 3 <= 1) ? Int3(offset[1], offset[2], offset[0]) : offset;

    atomicFetchAdd(gm_f[0][offset[0]], f[0]);
    atomicFetchAdd(gm_f[0][offset[1]], f[1]);
    atomicFetchAdd(gm_f[0][offset[2]], f[2]);
}

/* Harmonic */
template<bool calcEner>
static inline void harmonic_gpu(const float              kA,
                                const float              xA,
                                const float              x,
                                sycl::private_ptr<float> V,
                                sycl::private_ptr<float> F)
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
static inline void bonds_gpu(const int                                  i,
                             sycl::private_ptr<float>                   vtot_loc,
                             const sycl::global_ptr<const t_iatom>      gm_forceatoms,
                             const sycl::global_ptr<const t_iparams>    gm_forceparams,
                             const sycl::global_ptr<const sycl::float4> gm_xq,
                             sycl::global_ptr<Float3>                   gm_f,
                             sycl::local_ptr<Float3>                    sm_fShiftLoc,
                             const PbcAiuc&                             pbcAiuc,
                             const int                                  localId)
{
    const int type = gm_forceatoms[3 * i];
    const int ai   = gm_forceatoms[3 * i + 1];
    const int aj   = gm_forceatoms[3 * i + 2];

    /* dx = xi - xj, corrected for periodic boundary conditions. */
    Float3 dx;
    int    ki = pbcDxAiucSycl<calcVir>(pbcAiuc, gm_xq[ai], gm_xq[aj], dx);

    float dr2 = dx.norm2();
    float dr  = sycl::sqrt(dr2);

    float vbond;
    float fbond;
    harmonic_gpu<calcEner>(
            gm_forceparams[type].harmonic.krA, gm_forceparams[type].harmonic.rA, dr, &vbond, &fbond);

    if (calcEner)
    {
        *vtot_loc += vbond;
    }

    if (dr2 != 0.0F)
    {
        fbond *= sycl::rsqrt(dr2);

        Float3 fij = fbond * dx;
        staggeredAtomicAddForce(&gm_f[ai], fij, localId);
        staggeredAtomicAddForce(&gm_f[aj], -fij, localId);
        if (calcVir && ki != gmx::c_centralShiftIndex)
        {
            atomicFetchAddLocal(sm_fShiftLoc[ki], fij);
            atomicFetchAddLocal(sm_fShiftLoc[gmx::c_centralShiftIndex], -fij);
        }
    }
}

/*! \brief Cosine of an angle between two vectors.
 *
 * Computes cosine using the following formula:
 *
 *                  ax*bx + ay*by + az*bz
 * cos-vec (a,b) =  ---------------------
 *                      ||a|| * ||b||
 *
 * This function also makes sure that the cosine does not leave the [-1, 1]
 * interval, which can happen due to numerical errors.
 *
 * \param[in] a  First vector.
 * \param[in] b  Second vector.
 * \returns Cosine between a and b.
 */
static inline float cos_angle(const Float3& a, const Float3& b)
{
    float cosval;

    float ipa  = a.norm2();
    float ipb  = b.norm2();
    float ip   = a.dot(b);
    float ipab = ipa * ipb;
    if (ipab > 0.0F)
    {
        cosval = ip * sycl::rsqrt(ipab);
    }
    else
    {
        cosval = 1.0F;
    }
    if (cosval > 1.0F)
    {
        return 1.0F;
    }
    if (cosval < -1.0F)
    {
        return -1.0F;
    }

    return cosval;
}

/*! \brief Compute the angle between two vectors.
 *
 * Uses atan( |axb| / a.b ) formula.
 *
 * \param[in] a  First vector.
 * \param[in] b  Second vector.
 * \returns Angle between vectors in radians.
 */
static inline float angle(const Float3& a, const Float3& b)
{
    Float3 w = a.cross(b);

    float wlen = sycl::sqrt(w.norm2());
    float s    = a.dot(b);

    return sycl::atan2(wlen, s);
}

template<bool returnShift>
static float bond_angle_gpu(const sycl::float4        xi,
                            const sycl::float4        xj,
                            const sycl::float4        xk,
                            const PbcAiuc&            pbcAiuc,
                            sycl::private_ptr<Float3> r_ij,
                            sycl::private_ptr<Float3> r_kj,
                            sycl::private_ptr<float>  costh,
                            sycl::private_ptr<int>    t1,
                            sycl::private_ptr<int>    t2)
{
    *t1 = pbcDxAiucSycl<returnShift>(pbcAiuc, xi, xj, *r_ij);
    *t2 = pbcDxAiucSycl<returnShift>(pbcAiuc, xk, xj, *r_kj);

    *costh = cos_angle(*r_ij, *r_kj);
    // Return value is the angle between the bonds i-j and j-k
    return sycl::acos(*costh);
}

template<bool calcVir, bool calcEner>
static void angles_gpu(const int                                  i,
                       sycl::private_ptr<float>                   vtot_loc,
                       const sycl::global_ptr<const t_iatom>      gm_forceatoms,
                       const sycl::global_ptr<const t_iparams>    gm_forceparams,
                       const sycl::global_ptr<const sycl::float4> gm_xq,
                       sycl::global_ptr<Float3>                   gm_f,
                       sycl::local_ptr<Float3>                    sm_fShiftLoc,
                       const PbcAiuc&                             pbcAiuc,
                       const int                                  localId)
{
    sycl::int4 angleData;
    angleData.load(i, gm_forceatoms);
    const int type = angleData.x();
    const int ai   = angleData.y();
    const int aj   = angleData.z();
    const int ak   = angleData.w();

    Float3 r_ij;
    Float3 r_kj;
    float  cos_theta;
    int    t1;
    int    t2;
    float  theta = bond_angle_gpu<calcVir>(
            gm_xq[ai], gm_xq[aj], gm_xq[ak], pbcAiuc, &r_ij, &r_kj, &cos_theta, &t1, &t2);

    float va;
    float dVdt;
    harmonic_gpu<calcEner>(
            gm_forceparams[type].harmonic.krA, gm_forceparams[type].harmonic.rA * c_deg2RadF, theta, &va, &dVdt);

    if (calcEner)
    {
        *vtot_loc += va;
    }

    float cos_theta2 = cos_theta * cos_theta;
    if (cos_theta2 < 1.0F)
    {
        float st    = dVdt * sycl::rsqrt(1.0F - cos_theta2);
        float sth   = st * cos_theta;
        float nrij2 = r_ij.norm2();
        float nrkj2 = r_kj.norm2();

        float nrij_1 = sycl::rsqrt(nrij2);
        float nrkj_1 = sycl::rsqrt(nrkj2);

        float cik = st * nrij_1 * nrkj_1;
        float cii = sth * nrij_1 * nrij_1;
        float ckk = sth * nrkj_1 * nrkj_1;

        Float3 f_i = cii * r_ij - cik * r_kj;
        Float3 f_k = ckk * r_kj - cik * r_ij;
        Float3 f_j = -f_i - f_k;

        staggeredAtomicAddForce(&gm_f[ai], f_i, localId);
        staggeredAtomicAddForce(&gm_f[aj], f_j, localId);
        staggeredAtomicAddForce(&gm_f[ak], f_k, localId);

        if (calcVir)
        {
            atomicFetchAddLocal(sm_fShiftLoc[t1], f_i);
            atomicFetchAddLocal(sm_fShiftLoc[gmx::c_centralShiftIndex], f_j);
            atomicFetchAddLocal(sm_fShiftLoc[t2], f_k);
        }
    }
}

template<bool calcVir, bool calcEner>
static void urey_bradley_gpu(const int                                  i,
                             sycl::private_ptr<float>                   vtot_loc,
                             const sycl::global_ptr<const t_iatom>      gm_forceatoms,
                             const sycl::global_ptr<const t_iparams>    gm_forceparams,
                             const sycl::global_ptr<const sycl::float4> gm_xq,
                             sycl::global_ptr<Float3>                   gm_f,
                             sycl::local_ptr<Float3>                    sm_fShiftLoc,
                             const PbcAiuc&                             pbcAiuc,
                             const int                                  localId)
{
    sycl::int4 ubData;
    ubData.load(i, gm_forceatoms);
    const int type = ubData.x();
    const int ai   = ubData.y();
    const int aj   = ubData.z();
    const int ak   = ubData.w();

    const float th0A = gm_forceparams[type].u_b.thetaA * c_deg2RadF;
    const float kthA = gm_forceparams[type].u_b.kthetaA;
    const float r13A = gm_forceparams[type].u_b.r13A;
    const float kUBA = gm_forceparams[type].u_b.kUBA;

    Float3 r_ij;
    Float3 r_kj;
    float  cos_theta;
    int    t1;
    int    t2;
    float  theta = bond_angle_gpu<calcVir>(
            gm_xq[ai], gm_xq[aj], gm_xq[ak], pbcAiuc, &r_ij, &r_kj, &cos_theta, &t1, &t2);

    float va;
    float dVdt;
    harmonic_gpu<calcEner>(kthA, th0A, theta, &va, &dVdt);

    if (calcEner)
    {
        *vtot_loc += va;
    }

    Float3 r_ik;
    int    ki = pbcDxAiucSycl<calcVir>(pbcAiuc, gm_xq[ai], gm_xq[ak], r_ik);

    float dr2 = r_ik.norm2();
    float dr  = dr2 * sycl::rsqrt(dr2);

    float vbond;
    float fbond;
    harmonic_gpu<calcEner>(kUBA, r13A, dr, &vbond, &fbond);

    float cos_theta2 = cos_theta * cos_theta;

    Float3 f_i = { 0.0F, 0.0F, 0.0F };
    Float3 f_j = { 0.0F, 0.0F, 0.0F };
    Float3 f_k = { 0.0F, 0.0F, 0.0F };

    if (cos_theta2 < 1.0F)
    {
        float st  = dVdt * sycl::rsqrt(1.0F - cos_theta2);
        float sth = st * cos_theta;

        float nrkj2 = r_kj.norm2();
        float nrij2 = r_ij.norm2();

        float cik = st * sycl::rsqrt(nrkj2 * nrij2);
        float cii = sth / nrij2;
        float ckk = sth / nrkj2;

        f_i = cii * r_ij - cik * r_kj;
        f_k = ckk * r_kj - cik * r_ij;
        f_j = -f_i - f_k;

        if (calcVir)
        {
            atomicFetchAddLocal(sm_fShiftLoc[t1], f_i);
            atomicFetchAddLocal(sm_fShiftLoc[gmx::c_centralShiftIndex], f_j);
            atomicFetchAddLocal(sm_fShiftLoc[t2], f_k);
        }
    }

    /* Time for the bond calculations */
    if (dr2 != 0.0F)
    {
        if (calcEner)
        {
            *vtot_loc += vbond;
        }

        fbond *= sycl::rsqrt(dr2);
        Float3 fik = fbond * r_ik;
        f_i += fik;
        f_k -= fik;


        if (calcVir && ki != gmx::c_centralShiftIndex)
        {
            atomicFetchAddLocal(sm_fShiftLoc[ki], fik);
            atomicFetchAddLocal(sm_fShiftLoc[gmx::c_centralShiftIndex], -fik);
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
static float dih_angle_gpu(const T                   xi,
                           const T                   xj,
                           const T                   xk,
                           const T                   xl,
                           const PbcAiuc&            pbcAiuc,
                           sycl::private_ptr<Float3> r_ij,
                           sycl::private_ptr<Float3> r_kj,
                           sycl::private_ptr<Float3> r_kl,
                           sycl::private_ptr<Float3> m,
                           sycl::private_ptr<Float3> n,
                           sycl::private_ptr<int>    t1,
                           sycl::private_ptr<int>    t2,
                           sycl::private_ptr<int>    t3)
{
    *t1 = pbcDxAiucSycl<returnShift>(pbcAiuc, xi, xj, *r_ij);
    *t2 = pbcDxAiucSycl<returnShift>(pbcAiuc, xk, xj, *r_kj);
    *t3 = pbcDxAiucSycl<returnShift>(pbcAiuc, xk, xl, *r_kl);

    *m         = r_ij->cross(*r_kj);
    *n         = r_kj->cross(*r_kl);
    float phi  = angle(*m, *n);
    float ipr  = r_ij->dot(*n);
    float sign = (ipr < 0.0F) ? -1.0F : 1.0F;
    phi        = sign * phi;

    return phi;
}

template<bool returnShift, typename T>
static float dih_angle_gpu_sincos(const T                   xi,
                                  const T                   xj,
                                  const T                   xk,
                                  const T                   xl,
                                  const PbcAiuc&            pbcAiuc,
                                  sycl::private_ptr<Float3> r_ij,
                                  sycl::private_ptr<Float3> r_kj,
                                  sycl::private_ptr<Float3> r_kl,
                                  sycl::private_ptr<Float3> m,
                                  sycl::private_ptr<Float3> n,
                                  sycl::private_ptr<int>    t1,
                                  sycl::private_ptr<int>    t2,
                                  sycl::private_ptr<float>  cosval)
{
    *t1 = pbcDxAiucSycl<returnShift>(pbcAiuc, xi, xj, *r_ij);
    *t2 = pbcDxAiucSycl<returnShift>(pbcAiuc, xk, xj, *r_kj);
    pbcDxAiucSycl<returnShift>(pbcAiuc, xk, xl, *r_kl);

    *m = r_ij->cross(*r_kj);
    *n = r_kj->cross(*r_kl);

    Float3 w = m->cross(*n);

    float wlen = sycl::sqrt(w.norm2());
    float s    = m->dot(*n);

    float mLenSq = m->norm2();
    float nLenSq = n->norm2();
    float mnInv  = sycl::rsqrt(mLenSq * nLenSq);

    *cosval      = s * mnInv;
    float sinval = wlen * mnInv;

    float ipr  = r_ij->dot(*n);
    float sign = (ipr < 0.0F) ? -1.0F : 1.0F;

    return sign * sinval;
}


static void dopdihs_gpu(const float              cpA,
                        const float              phiA,
                        const int                mult,
                        const float              phi,
                        sycl::private_ptr<float> v,
                        sycl::private_ptr<float> f)
{
    float mdphi = mult * phi - phiA * c_deg2RadF;
    float sdphi = sycl::sin(mdphi);
    *v          = cpA * (1.0F + sycl::cos(mdphi));
    *f          = -cpA * mult * sdphi;
}

template<bool calcVir>
static void do_dih_fup_gpu(const int                            i,
                           const int                            j,
                           const int                            k,
                           const int                            l,
                           const float                          ddphi,
                           const Float3                         r_ij,
                           const Float3                         r_kj,
                           const Float3                         r_kl,
                           const Float3                         m,
                           const Float3                         n,
                           sycl::global_ptr<Float3>             gm_f,
                           sycl::local_ptr<Float3>              sm_fShiftLoc,
                           const PbcAiuc&                       pbcAiuc,
                           const sycl::global_ptr<const Float4> gm_xq,
                           const int                            t1,
                           const int                            t2,
                           const int                            localId)
{
    float iprm  = m.norm2();
    float iprn  = n.norm2();
    float nrkj2 = r_kj.norm2();
    float toler = nrkj2 * GMX_REAL_EPS;
    if ((iprm > toler) && (iprn > toler))
    {
        float  nrkj_1 = sycl::rsqrt(nrkj2);
        float  nrkj_2 = nrkj_1 * nrkj_1;
        float  nrkj   = nrkj2 * nrkj_1;
        float  a      = -ddphi * nrkj / iprm;
        Float3 f_i    = a * m;
        float  b      = ddphi * nrkj / iprn;
        Float3 f_l    = b * n;
        float  p      = r_ij.dot(r_kj);
        p *= nrkj_2;
        float q = r_kl.dot(r_kj);
        q *= nrkj_2;
        Float3 uvec = p * f_i;
        Float3 vvec = q * f_l;
        Float3 svec = uvec - vvec;
        Float3 f_j  = f_i - svec;
        Float3 f_k  = f_l + svec;

        staggeredAtomicAddForce(&gm_f[i], f_i, localId);
        staggeredAtomicAddForce(&gm_f[j], -f_j, localId);
        staggeredAtomicAddForce(&gm_f[k], -f_k, localId);
        staggeredAtomicAddForce(&gm_f[l], f_l, localId);

        if constexpr (calcVir)
        {
            Float3 dx_jl;
            int    t3 = pbcDxAiucSycl<calcVir>(pbcAiuc, gm_xq[l], gm_xq[j], dx_jl);

            atomicFetchAddLocal(sm_fShiftLoc[t1], f_i);
            atomicFetchAddLocal(sm_fShiftLoc[gmx::c_centralShiftIndex], -f_j);
            atomicFetchAddLocal(sm_fShiftLoc[t2], -f_k);
            atomicFetchAddLocal(sm_fShiftLoc[t3], f_l);
        }
    }
}

template<bool calcVir, bool calcEner>
static void pdihs_gpu(const int                                  i,
                      sycl::private_ptr<float>                   vtot_loc,
                      const sycl::global_ptr<const t_iatom>      gm_forceatoms,
                      const sycl::global_ptr<const t_iparams>    gm_forceparams,
                      const sycl::global_ptr<const sycl::float4> gm_xq,
                      sycl::global_ptr<Float3>                   gm_f,
                      sycl::local_ptr<Float3>                    sm_fShiftLoc,
                      const PbcAiuc&                             pbcAiuc,
                      const int                                  localId)
{
    int type = gm_forceatoms[5 * i];
    int ai   = gm_forceatoms[5 * i + 1];
    int aj   = gm_forceatoms[5 * i + 2];
    int ak   = gm_forceatoms[5 * i + 3];
    int al   = gm_forceatoms[5 * i + 4];

    Float3 r_ij;
    Float3 r_kj;
    Float3 r_kl;
    Float3 m;
    Float3 n;
    int    t1;
    int    t2;
    int    t3;
    float  phi = dih_angle_gpu<calcVir>(
            gm_xq[ai], gm_xq[aj], gm_xq[ak], gm_xq[al], pbcAiuc, &r_ij, &r_kj, &r_kl, &m, &n, &t1, &t2, &t3);

    float vpd;
    float ddphi;
    dopdihs_gpu(gm_forceparams[type].pdihs.cpA,
                gm_forceparams[type].pdihs.phiA,
                gm_forceparams[type].pdihs.mult,
                phi,
                &vpd,
                &ddphi);

    if (calcEner)
    {
        *vtot_loc += vpd;
    }

    do_dih_fup_gpu<calcVir>(
            ai, aj, ak, al, ddphi, r_ij, r_kj, r_kl, m, n, gm_f, sm_fShiftLoc, pbcAiuc, gm_xq, t1, t2, localId);
}

template<bool calcVir, bool calcEner>
static void rbdihs_gpu(const int                                  i,
                       sycl::private_ptr<float>                   vtot_loc,
                       const sycl::global_ptr<const t_iatom>      gm_forceatoms,
                       const sycl::global_ptr<const t_iparams>    gm_forceparams,
                       const sycl::global_ptr<const sycl::float4> gm_xq,
                       sycl::global_ptr<Float3>                   gm_f,
                       sycl::local_ptr<Float3>                    sm_fShiftLoc,
                       const PbcAiuc&                             pbcAiuc,
                       const int                                  localId)
{
    constexpr float c0 = 0.0F, c1 = 1.0F, c2 = 2.0F, c3 = 3.0F, c4 = 4.0F, c5 = 5.0F;

    {
        int type = gm_forceatoms[5 * i];
        int ai   = gm_forceatoms[5 * i + 1];
        int aj   = gm_forceatoms[5 * i + 2];
        int ak   = gm_forceatoms[5 * i + 3];
        int al   = gm_forceatoms[5 * i + 4];

        Float3 r_ij;
        Float3 r_kj;
        Float3 r_kl;
        Float3 m;
        Float3 n;
        int    t1;
        int    t2;

        // Changing the sign of sin and cos to convert to polymer convention
        float       cos_phi;
        const float negative_sin_phi = dih_angle_gpu_sincos<calcVir>(
                gm_xq[ai], gm_xq[aj], gm_xq[ak], gm_xq[al], pbcAiuc, &r_ij, &r_kj, &r_kl, &m, &n, &t1, &t2, &cos_phi);
        cos_phi *= -1;

        float parm[NR_RBDIHS];
#pragma unroll NR_RBDIHS
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
                ai, aj, ak, al, ddphi, r_ij, r_kj, r_kl, m, n, gm_f, sm_fShiftLoc, pbcAiuc, gm_xq, t1, t2, localId);
        if (calcEner)
        {
            *vtot_loc += v;
        }
    }
}

//! Wrap angle from range [-3*pi; 3*pi) to [-pi; pi)
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
static void idihs_gpu(const int                                  i,
                      sycl::private_ptr<float>                   vtot_loc,
                      const sycl::global_ptr<const t_iatom>      gm_forceatoms,
                      const sycl::global_ptr<const t_iparams>    gm_forceparams,
                      const sycl::global_ptr<const sycl::float4> gm_xq,
                      sycl::global_ptr<Float3>                   gm_f,
                      sycl::local_ptr<Float3>                    sm_fShiftLoc,
                      const PbcAiuc&                             pbcAiuc,
                      const int                                  localId)
{
    int type = gm_forceatoms[5 * i];
    int ai   = gm_forceatoms[5 * i + 1];
    int aj   = gm_forceatoms[5 * i + 2];
    int ak   = gm_forceatoms[5 * i + 3];
    int al   = gm_forceatoms[5 * i + 4];

    Float3 r_ij;
    Float3 r_kj;
    Float3 r_kl;
    Float3 m;
    Float3 n;
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
    float kA = gm_forceparams[type].harmonic.krA;
    float pA = gm_forceparams[type].harmonic.rA;

    float phi0 = pA * c_deg2RadF;

    float dp = wrapAngle(phi - phi0);

    float ddphi = -kA * dp;

    do_dih_fup_gpu<calcVir>(
            ai, aj, ak, al, -ddphi, r_ij, r_kj, r_kl, m, n, gm_f, sm_fShiftLoc, pbcAiuc, gm_xq, t1, t2, localId);

    if (calcEner)
    {
        *vtot_loc += -0.5F * ddphi * dp;
    }
}

template<bool calcVir, bool calcEner>
static void pairs_gpu(const int                                  i,
                      const sycl::global_ptr<const t_iatom>      gm_forceatoms,
                      const sycl::global_ptr<const t_iparams>    gm_iparams,
                      const sycl::global_ptr<const sycl::float4> gm_xq,
                      sycl::global_ptr<Float3>                   gm_f,
                      sycl::local_ptr<Float3>                    sm_fShiftLoc,
                      const PbcAiuc&                             pbcAiuc,
                      const float                                scale_factor,
                      sycl::private_ptr<float>                   vtotVdw_loc,
                      sycl::private_ptr<float>                   vtotElec_loc,
                      const int                                  localId)
{
    int type = gm_forceatoms[3 * i];
    int ai   = gm_forceatoms[3 * i + 1];
    int aj   = gm_forceatoms[3 * i + 2];

    float qq  = gm_xq[ai].w() * gm_xq[aj].w();
    float c6  = gm_iparams[type].lj14.c6A;
    float c12 = gm_iparams[type].lj14.c12A;

    /* Do we need to apply full periodic boundary conditions? */
    Float3 dr;
    int    fshift_index = pbcDxAiucSycl<calcVir>(pbcAiuc, gm_xq[ai], gm_xq[aj], dr);

    float r2    = dr.norm2();
    float rinv  = sycl::rsqrt(r2);
    float rinv2 = rinv * rinv;
    float rinv6 = rinv2 * rinv2 * rinv2;

    /* Calculate the Coulomb force * r */
    float velec = scale_factor * qq * rinv;

    /* Calculate the LJ force * r and add it to the Coulomb part */
    float fr = (12.0F * c12 * rinv6 - 6.0F * c6) * rinv6 + velec;

    float  finvr = fr * rinv2;
    Float3 f     = finvr * dr;

    /* Add the forces */
    staggeredAtomicAddForce(&gm_f[ai], f, localId);
    staggeredAtomicAddForce(&gm_f[aj], -f, localId);
    if (calcVir && fshift_index != gmx::c_centralShiftIndex)
    {
        atomicFetchAddLocal(sm_fShiftLoc[fshift_index], f);
        atomicFetchAddLocal(sm_fShiftLoc[gmx::c_centralShiftIndex], -f);
    }

    if (calcEner)
    {
        *vtotVdw_loc += (c12 * rinv6 - c6) * rinv6;
        *vtotElec_loc += velec;
    }
}

template<bool calcVir, bool calcEner>
class BondedKernel;

namespace gmx
{

using sycl::access::fence_space;
using mode = sycl::access_mode;

template<bool calcVir, bool calcEner>
auto bondedKernel(sycl::handler&                   cgh,
                  const BondedGpuKernelParameters& kernelParams,
                  const DeviceBuffer<t_iatom>      gm_iatoms_[numFTypesOnGpu],
                  float* __restrict__ gm_vTot,
                  const t_iparams* __restrict__ gm_forceParams_,
                  const sycl::float4* __restrict__ gm_xq_,
                  Float3* __restrict__ gm_f_,
                  Float3* __restrict__ gm_fShift_)
{
    sycl::global_ptr<const t_iatom> gm_iatomsTemp[numFTypesOnGpu];
    for (int i = 0; i < numFTypesOnGpu; i++)
    {
        gm_iatomsTemp[i] = gm_iatoms_[i].get_pointer();
    }
    const FTypeArray<sycl::global_ptr<const t_iatom>> gm_iatoms(gm_iatomsTemp);

    const FTypeArray<int> fTypeRangeStart(kernelParams.fTypeRangeStart);
    const FTypeArray<int> fTypeRangeEnd(kernelParams.fTypeRangeEnd);
    const FTypeArray<int> numFTypeBonds(kernelParams.numFTypeBonds);

    const auto electrostaticsScaleFactor = kernelParams.electrostaticsScaleFactor;

    sycl::local_accessor<Float3, 1> sm_fShiftLoc{ sycl::range<1>(c_numShiftVectors), cgh };

    const PbcAiuc pbcAiuc = kernelParams.pbcAiuc;

    return [=](sycl::nd_item<1> itemIdx) {
        sycl::global_ptr<const t_iparams> gm_forceParams = gm_forceParams_;
        sycl::global_ptr<const Float4>    gm_xq          = gm_xq_;
        sycl::global_ptr<Float3>          gm_f           = gm_f_;
        sycl::global_ptr<Float3>          gm_fShift      = gm_fShift_;

        const int tid          = itemIdx.get_global_linear_id();
        const int localId      = itemIdx.get_local_linear_id();
        float     vtot_loc     = 0.0F;
        float     vtotElec_loc = 0.0F;

        if constexpr (calcVir)
        {
            if (localId < c_numShiftVectors)
            {
                sm_fShiftLoc[localId] = { 0.0F, 0.0F, 0.0F };
            }
            itemIdx.barrier(fence_space::local_space);
        }

        int  fType;
        bool threadComputedPotential = false;
#pragma unroll
        for (int j = 0; j < numFTypesOnGpu; j++)
        {
            if (tid >= fTypeRangeStart[j] && tid <= fTypeRangeEnd[j])
            {
                const int                             numBonds = numFTypeBonds[j];
                const int                             fTypeTid = tid - fTypeRangeStart[j];
                const sycl::global_ptr<const t_iatom> iatoms   = gm_iatoms[j];
                fType                                          = fTypesOnGpu[j];

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
                        bonds_gpu<calcVir, calcEner>(
                                fTypeTid, &vtot_loc, iatoms, gm_forceParams, gm_xq, gm_f, sm_fShiftLoc, pbcAiuc, localId);
                        break;
                    case F_ANGLES:
                        angles_gpu<calcVir, calcEner>(
                                fTypeTid, &vtot_loc, iatoms, gm_forceParams, gm_xq, gm_f, sm_fShiftLoc, pbcAiuc, localId);
                        break;
                    case F_UREY_BRADLEY:
                        urey_bradley_gpu<calcVir, calcEner>(
                                fTypeTid, &vtot_loc, iatoms, gm_forceParams, gm_xq, gm_f, sm_fShiftLoc, pbcAiuc, localId);
                        break;
                    case F_PDIHS:
                    case F_PIDIHS:
                        pdihs_gpu<calcVir, calcEner>(
                                fTypeTid, &vtot_loc, iatoms, gm_forceParams, gm_xq, gm_f, sm_fShiftLoc, pbcAiuc, localId);
                        break;
                    case F_RBDIHS:
                        rbdihs_gpu<calcVir, calcEner>(
                                fTypeTid, &vtot_loc, iatoms, gm_forceParams, gm_xq, gm_f, sm_fShiftLoc, pbcAiuc, localId);
                        break;
                    case F_IDIHS:
                        idihs_gpu<calcVir, calcEner>(
                                fTypeTid, &vtot_loc, iatoms, gm_forceParams, gm_xq, gm_f, sm_fShiftLoc, pbcAiuc, localId);
                        break;
                    case F_LJ14:
                        pairs_gpu<calcVir, calcEner>(fTypeTid,
                                                     iatoms,
                                                     gm_forceParams,
                                                     gm_xq,
                                                     gm_f,
                                                     sm_fShiftLoc,
                                                     pbcAiuc,
                                                     electrostaticsScaleFactor,
                                                     &vtot_loc,
                                                     &vtotElec_loc,
                                                     localId);
                        break;
                }
                break;
            }
        }

        if (calcEner && threadComputedPotential)
        {
            subGroupBarrier(itemIdx); // Should not be needed, but https://github.com/AdaptiveCpp/AdaptiveCpp/issues/823
            sycl::sub_group sg = itemIdx.get_sub_group();
            vtot_loc           = sycl::reduce_over_group(sg, vtot_loc, sycl::plus<float>());
            vtotElec_loc       = sycl::reduce_over_group(sg, vtotElec_loc, sycl::plus<float>());
            if (sg.leader())
            {
                atomicFetchAdd(gm_vTot[fType], vtot_loc);
                if (fType == F_LJ14)
                {
                    atomicFetchAdd(gm_vTot[F_COUL14], vtotElec_loc);
                }
            }
        }
        /* Accumulate shift vectors from shared memory to global memory on the first c_numShiftVectors threads of the block. */
        if constexpr (calcVir)
        {
            itemIdx.barrier(fence_space::local_space);
            if (localId < c_numShiftVectors)
            {
                const Float3 tmp = sm_fShiftLoc[localId];
                atomicFetchAdd(gm_fShift[localId], tmp);
            }
        }
    };
}


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

    using kernelNameType = BondedKernel<calcVir, calcEner>;

    const sycl::nd_range<1> rangeAll(kernelLaunchConfig_.blockSize[0] * kernelLaunchConfig_.gridSize[0],
                                     kernelLaunchConfig_.blockSize[0]);
    sycl::queue             q = deviceStream_.stream();

    q.submit(GMX_SYCL_DISCARD_EVENT[&](sycl::handler & cgh) {
        auto kernel = bondedKernel<calcVir, calcEner>(cgh,
                                                      kernelParams_,
                                                      kernelBuffers_.d_iatoms,
                                                      kernelBuffers_.d_vTot.get_pointer(),
                                                      kernelBuffers_.d_forceParams.get_pointer(),
                                                      d_xq_.get_pointer(),
                                                      d_f_.get_pointer(),
                                                      d_fShift_.get_pointer());
        cgh.parallel_for<kernelNameType>(rangeAll, kernel);
    });

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
