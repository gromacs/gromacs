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
/*! \file
 * \brief functionality for performing update and constraints on the gpu
 * \author Jon Vincent <jvincent@nvidia.com> and Alan Gray <alang@nvidia.com>
 */

#include <stdio.h>
#include <math_constants.h>

#include "gromacs/ewald/pme.cuh"
#include "gromacs/mdlib/nbnxn_gpu.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/gpu_utils/cudautils.cuh"

#include "gpuBufferOpsCUDA.h"
#include "gpuBondedCUDA.h"

#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/math/units.h"
#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/mdlib/force_flags.h"

// this needs to be in an include file probably
// duplicated from  pbcutil/pbc.cpp
enum {
    epbcdxRECTANGULAR = 1, epbcdxTRICLINIC,
    epbcdx2D_RECT,       epbcdx2D_TRIC,
    epbcdx1D_RECT,       epbcdx1D_TRIC,
    epbcdxSCREW_RECT,    epbcdxSCREW_TRIC,
    epbcdxNOPBC,         epbcdxUNSUPPORTED
};

#include "gromacs/mdlib/nbnxn_cuda/gpu_vec.h"

#if defined(_MSVC)
#include <limits>
#endif

static gmx_bool
ftype_is_bonded_potential(int ftype)
{
    return (((interaction_function[ftype].flags & IF_BOND) != 0u) &&
            !(ftype == F_CONNBONDS || ftype == F_POSRES || ftype == F_FBPOSRES));
}

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
                        int size, float *energygrp_elec, float *energygrp_vdw,
                        int nener, fvec *f_shift)
{
    int a = blockIdx.x * blockDim.x + threadIdx.x;

    if (a < F_NRE)
    {
        vtot[a] = 0.0;
    }

    if (a < size)
    {
        force[a][XX] = 0.0;
        force[a][YY] = 0.0;
        force[a][ZZ] = 0.0;
    }

    if (a < nener)
    {
        energygrp_vdw[a]  = 0.0;
        energygrp_elec[a] = 0.0;
    }

    if (a < SHIFTS)
    {
        f_shift[a][XX] = 0.0;
        f_shift[a][YY] = 0.0;
        f_shift[a][ZZ] = 0.0;
    }
}

__device__
static int pbc_fvec_sub_gpu(const t_pbc *pbc, const fvec xi, const fvec xj, fvec dx,
                            const fvec *pbc_hbox_diag, const matrix *pbc_box,
                            const fvec *pbc_mhbox_diag, const fvec *pbc_fbox_diag,
                            const fvec *pbc_tric_vec, const ivec* pbc_tric_shift)
{
    /* Need to be careful here, gpu pointer needs to be null if cpu pointer is null */
    if (pbc)
    {
        return pbc_dx_aiuc_gpu(pbc, xi, xj, dx, *pbc_hbox_diag, *pbc_box, *pbc_mhbox_diag,
                               *pbc_fbox_diag, pbc_tric_vec, pbc_tric_shift);
    }
    else
    {
        fvec_sub_gpu(xi, xj, dx);
        return CENTRAL;
    }
}

/* Harmonic */
__device__
static void harmonic_gpu(float kA, float xA, float x, float *V, float *F)
{
    constexpr float half = 0.5f;
    float       dx, dx2;

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
               const t_pbc *pbc, const fvec *pbc_hbox_diag, const matrix *pbc_box,
               const fvec *pbc_mhbox_diag, const fvec  *pbc_fbox_diag,
               const fvec *pbc_tric_vec, const ivec *pbc_tric_shift)
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
        fvec dx;
        int  ki =
            pbc_fvec_sub_gpu(pbc, x[ai], x[aj], dx,
                             pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                             pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift);

        float dr2 = iprod_gpu(dx, dx);
        float dr  = sqrt(dr2);

        float vbond;
        float fbond;
        harmonic_gpu(forceparams[type].harmonic.krA,
                     forceparams[type].harmonic.rA,
                     dr, &vbond, &fbond);
        if (dr2 == 0.0f)
        {
            return;
        }

        if (calcEnerVir)
        {
            atomicAdd(&vtot_loc, vbond);
        }

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
static float bond_angle_gpu(const fvec xi, const fvec xj, const fvec xk, const t_pbc *pbc,
                            const fvec *pbc_hbox_diag, const matrix *pbc_box,
                            const fvec *pbc_mhbox_diag, const fvec  *pbc_fbox_diag,
                            const fvec*  pbc_tric_vec,  const ivec*  pbc_tric_shift,
                            fvec r_ij, fvec r_kj, float *costh,
                            int *t1, int *t2)
/* Return value is the angle between the bonds i-j and j-k */
{
    float th;

    *t1 = pbc_fvec_sub_gpu(pbc, xi, xj, r_ij,
                           pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                           pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift);
    *t2 = pbc_fvec_sub_gpu(pbc, xk, xj, r_kj,
                           pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                           pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift);

    *costh = cos_angle_gpu(r_ij, r_kj);
    th     = acosf(*costh);

    return th;
}

template <bool calcEnerVir>
__global__
void angles_gpu(float *vtot, int nbonds, int natoms,
                const t_iatom forceatoms[], const t_iparams forceparams[],
                const fvec x[], fvec force[], fvec fshift[],
                const t_pbc *pbc, const fvec *pbc_hbox_diag, const matrix *pbc_box,
                const fvec *pbc_mhbox_diag, const fvec  *pbc_fbox_diag,
                const fvec*  pbc_tric_vec,  const ivec*  pbc_tric_shift)
{
    __shared__ float vtot_loc;
    __shared__ fvec  fshift_loc[SHIFTS];

    const int        i = blockIdx.x*blockDim.x+threadIdx.x;

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
        int type = forceatoms[4*i];
        int ai   = forceatoms[4*i + 1];
        int aj   = forceatoms[4*i + 2];
        int ak   = forceatoms[4*i + 3];

        fvec  r_ij;
        fvec  r_kj;
        float cos_theta;
        int   t1;
        int   t2;
        float theta =
            bond_angle_gpu(x[ai], x[aj], x[ak], pbc,
                           pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                           pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift,
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

            fvec f_i;
            fvec f_k;
            fvec f_j;
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

__global__
void urey_bradley_gpu(float *vtot, int nbonds, int natoms,
                      const t_iatom forceatoms[], const t_iparams forceparams[],
                      const fvec x[], fvec force[], fvec fshift[],
                      const t_pbc *pbc, const fvec *pbc_hbox_diag, const matrix *pbc_box,
                      const fvec *pbc_mhbox_diag, const fvec  *pbc_fbox_diag,
                      const fvec*  pbc_tric_vec,  const ivec*  pbc_tric_shift)
{
    int              i, m, ai, aj, ak, t1, t2, type, ki;
    fvec             r_ij, r_kj, r_ik;
    float            cos_theta, cos_theta2, theta;
    float            dVdt, va, dr, dr2, vbond, fbond, fik;
    float            kthA, th0A, kUBA, r13A;

    __shared__ float vtot_loc;
    __shared__ fvec  fshift_loc[SHIFTS];


    i = blockIdx.x*blockDim.x+threadIdx.x;
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
    if (i < nbonds)
    {
        type  = forceatoms[4*i];
        ai    = forceatoms[4*i+1];
        aj    = forceatoms[4*i+2];
        ak    = forceatoms[4*i+3];
        th0A  = forceparams[type].u_b.thetaA*DEG2RAD;
        kthA  = forceparams[type].u_b.kthetaA;
        r13A  = forceparams[type].u_b.r13A;
        kUBA  = forceparams[type].u_b.kUBA;

        theta  = bond_angle_gpu(x[ai], x[aj], x[ak], pbc,
                                pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                                pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift,
                                r_ij, r_kj, &cos_theta, &t1, &t2);

        harmonic_gpu(kthA, th0A, theta, &va, &dVdt);
        atomicAdd(&vtot_loc,va);
        
        ki = pbc_fvec_sub_gpu(pbc, x[ai], x[ak], r_ik,
                              pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                              pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift);
        dr2  = iprod_gpu(r_ik, r_ik);
        dr   = dr2*rsqrtf(dr2);

        harmonic_gpu(kUBA, r13A, dr, &vbond, &fbond);

        cos_theta2 = cos_theta*cos_theta;
        if (cos_theta2 < 1.0f)
        {
            float st, sth;
            float cik, cii, ckk;
            float nrkj2, nrij2;
            fvec  f_i, f_j, f_k;

            st  = dVdt*rsqrtf(1.0f - cos_theta2);
            sth = st*cos_theta;

            nrkj2 = iprod_gpu(r_kj, r_kj);
            nrij2 = iprod_gpu(r_ij, r_ij);

            cik = st*rsqrtf(nrkj2*nrij2);
            cii = sth/nrij2;
            ckk = sth/nrkj2;

            for (m = 0; (m < DIM); m++)
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
            atomicAdd(&vtot_loc, vbond);
            fbond *= rsqrtf(dr2);

            for (m = 0; (m < DIM); m++)
            {
                fik = fbond*r_ik[m];
                atomicAdd(&force[ai][m], fik);
                atomicAdd(&force[ak][m], -fik);
                if (ki != CENTRAL)
                {
                    atomicAdd(&fshift_loc[ki][m], fik);
                    atomicAdd(&fshift_loc[CENTRAL][m], -fik);
                }
            }
        }
    }
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

__global__
void urey_bradley_gpu_noener(int nbonds, int natoms,
                             const t_iatom forceatoms[], const t_iparams forceparams[],
                             const fvec x[], fvec force[],
                             const t_pbc *pbc, const fvec *pbc_hbox_diag, const matrix *pbc_box,
                             const fvec *pbc_mhbox_diag, const fvec  *pbc_fbox_diag,
                             const fvec*  pbc_tric_vec,  const ivec*  pbc_tric_shift)
{
    int   i, m, ai, aj, ak, t1, t2, type;
    fvec  r_ij, r_kj, r_ik;
    float cos_theta, cos_theta2, theta;
    float dVdt, va, dr, dr2, vbond, fbond, fik;
    float kthA, th0A, kUBA, r13A;

    i = blockIdx.x*blockDim.x+threadIdx.x;
    if (i < nbonds)
    {
        type  = forceatoms[4*i];
        ai    = forceatoms[4*i+1];
        aj    = forceatoms[4*i+2];
        ak    = forceatoms[4*i+3];
        th0A  = forceparams[type].u_b.thetaA*DEG2RAD;
        kthA  = forceparams[type].u_b.kthetaA;
        r13A  = forceparams[type].u_b.r13A;
        kUBA  = forceparams[type].u_b.kUBA;

        theta  = bond_angle_gpu(x[ai], x[aj], x[ak], pbc,
                                pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                                pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift,
                                r_ij, r_kj, &cos_theta, &t1, &t2);

        harmonic_gpu(kthA, th0A, theta, &va, &dVdt);

        pbc_fvec_sub_gpu(pbc, x[ai], x[ak], r_ik,
                         pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                         pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift);
        dr2  = iprod_gpu(r_ik, r_ik);
        dr   = dr2*rsqrtf(dr2);

        harmonic_gpu(kUBA, r13A, dr, &vbond, &fbond);

        cos_theta2 = cos_theta*cos_theta;
        if (cos_theta2 < 1.0f)
        {
            float st, sth;
            float cik, cii, ckk;
            float nrkj2, nrij2;
            fvec  f_i, f_j, f_k;

            st  = dVdt*rsqrtf(1.0f - cos_theta2);
            sth = st*cos_theta;

            nrkj2 = iprod_gpu(r_kj, r_kj);
            nrij2 = iprod_gpu(r_ij, r_ij);

            cik = st*rsqrtf(nrkj2*nrij2);
            cii = sth/nrij2;
            ckk = sth/nrkj2;

            for (m = 0; (m < DIM); m++)
            {
                f_i[m]    = -(cik*r_kj[m]-cii*r_ij[m]);
                f_k[m]    = -(cik*r_ij[m]-ckk*r_kj[m]);
                f_j[m]    = -f_i[m]-f_k[m];
                atomicAdd(&force[ai][m], f_i[m]);
                atomicAdd(&force[aj][m], f_j[m]);
                atomicAdd(&force[ak][m], f_k[m]);
            }
        }

        /* Time for the bond calculations */
        if (dr2 != 0.0f)
        {
            fbond *= rsqrtf(dr2);

            for (m = 0; (m < DIM); m++)
            {
                fik = fbond*r_ik[m];
                atomicAdd(&force[ai][m], fik);
                atomicAdd(&force[ak][m], -fik);
            }
        }
    }
}

__device__
static float dih_angle_gpu(const fvec xi, const fvec xj, const fvec xk, const fvec xl,
                           const t_pbc *pbc,  const fvec *pbc_hbox_diag, const matrix *pbc_box,
                           const fvec *pbc_mhbox_diag, const fvec  *pbc_fbox_diag,
                           const fvec*  pbc_tric_vec,  const ivec*  pbc_tric_shift,
                           fvec r_ij, fvec r_kj, fvec r_kl, fvec m, fvec n,
                           int *t1, int *t2, int *t3)
{
    *t1 = pbc_fvec_sub_gpu(pbc, xi, xj, r_ij,
                           pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                           pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift);
    *t2 = pbc_fvec_sub_gpu(pbc, xk, xj, r_kj,
                           pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                           pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift);
    *t3 = pbc_fvec_sub_gpu(pbc, xk, xl, r_kl,
                           pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                           pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift);

    cprod_gpu(r_ij, r_kj, m);
    cprod_gpu(r_kj, r_kl, n);
    float phi  = gmx_angle_gpu(m, n);
    float ipr  = iprod_gpu(r_ij, n);
    float sign = (ipr < 0.0f) ? -1.0f : 1.0f;
    phi        = sign*phi;

    return phi;
}


__device__
static void dopdihs_gpu(float cpA, float phiA, int mult, float phi, float *V, float *F)
{
    float mdphi, sdphi;

    mdphi = mult*phi - phiA*DEG2RAD;
    sdphi = sinf(mdphi);
    *V    = cpA * (1.0f + cosf(mdphi));
    *F    = -cpA*mult*sdphi;
}

template <bool calcVir>
__device__
static void do_dih_fup_gpu(int i, int j, int k, int l, int natoms,
                           float ddphi, fvec r_ij, fvec r_kj, fvec r_kl,
                           fvec m, fvec n, fvec force[], fvec fshift[],
                           const t_pbc *pbc,  const fvec *pbc_hbox_diag, const matrix *pbc_box,
                           const fvec *pbc_mhbox_diag, const fvec  *pbc_fbox_diag,
                           const fvec*  pbc_tric_vec,  const ivec*  pbc_tric_shift,
                           const fvec x[], int t1, int t2, int t3)
{
    float iprm  = iprod_gpu(m, m);
    float iprn  = iprod_gpu(n, n);
    float nrkj2 = iprod_gpu(r_kj, r_kj);
    float toler = nrkj2*GMX_REAL_EPS;
    if ((iprm > toler) && (iprn > toler))
    {
        float nrkj_1 = 1.0f/sqrt(nrkj2); // replacing std::invsqrt call
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
            int  t3 =
                pbc_fvec_sub_gpu(pbc, x[l], x[j], dx_jl,
                                 pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                                 pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift);

            fvec_inc_atomic(fshift[t1], f_i);
            fvec_dec_atomic(fshift[CENTRAL], f_j);
            fvec_dec_atomic(fshift[t2], f_k);
            fvec_inc_atomic(fshift[t3], f_l);
        }
    }
}

template <bool calcEnerVir>
__global__
void  pdihs_gpu(float *vtot, int nbonds, int natoms,
                const t_iatom forceatoms[], const t_iparams forceparams[],
                const fvec x[], fvec f[], fvec fshift[],
                const t_pbc *pbc, const fvec *pbc_hbox_diag, const matrix *pbc_box,
                const fvec *pbc_mhbox_diag, const fvec  *pbc_fbox_diag,
                const fvec*  pbc_tric_vec,  const ivec*  pbc_tric_shift)
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
            dih_angle_gpu(x[ai], x[aj], x[ak], x[al], pbc,
                          pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                          pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift,
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
                                    f, fshift_loc, pbc,
                                    pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                                    pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift,
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
void rbdihs_gpu(float *vtot, int nbonds, int natoms,
                const t_iatom forceatoms[], const t_iparams forceparams[],
                const fvec x[], fvec f[], fvec fshift[],
                const t_pbc *pbc,  const fvec *pbc_hbox_diag, const matrix *pbc_box,
                const fvec *pbc_mhbox_diag, const fvec  *pbc_fbox_diag,
                const fvec*  pbc_tric_vec,  const ivec*  pbc_tric_shift)
{
    constexpr float  c0 = 0.0f, c1 = 1.0f, c2 = 2.0f, c3 = 3.0f, c4 = 4.0f, c5 = 5.0f;

    __shared__ float vtot_loc;
    __shared__ fvec  fshift_loc[SHIFTS];

    const int i = blockIdx.x*blockDim.x + threadIdx.x;

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
            dih_angle_gpu(x[ai], x[aj], x[ak], x[al], pbc,
                          pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                          pbc_fbox_diag, pbc_tric_vec, pbc_tric_shift,
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
                                    f, fshift_loc, pbc,
                                    pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                                    pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift,
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
    return;
}

template <bool calcEnerVir>
__global__
void  idihs_gpu(float *vtot, int nbonds, int natoms,
                const t_iatom forceatoms[], const t_iparams forceparams[],
                const fvec x[], fvec f[], fvec fshift[],
                const t_pbc *pbc,  const fvec *pbc_hbox_diag, const matrix *pbc_box,
                const fvec *pbc_mhbox_diag, const fvec  *pbc_fbox_diag,
                const fvec*  pbc_tric_vec,  const ivec*  pbc_tric_shift)
{
    const int i = blockIdx.x*blockDim.x + threadIdx.x;

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
            dih_angle_gpu(x[ai], x[aj], x[ak], x[al], pbc,
                          pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                          pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift,
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
                                    f, fshift_loc, pbc,
                                    pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                                    pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift,
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

__device__
static float evaluate_single(float r2, float tabscale, float *vftab, float tableStride,
                             float qq, float c6, float c12, float *velec, float *vvdw)
{
    float      rinv, r, rtab, eps, eps2, Y, F, Geps, Heps2, Fp, VVe, FFe, VVd, FFd, VVr, FFr, fscal;
    int        ntab;

    /* Do the tabulated interactions - first table lookup */
    rinv             = rsqrtf(r2);
    r                = r2*rinv;
    rtab             = r*tabscale;
    ntab             = static_cast<int>(rtab);
    eps              = rtab-ntab;
    eps2             = eps*eps;
    ntab             = tableStride*ntab;
    /* Electrostatics */
    Y                = vftab[ntab];
    F                = vftab[ntab+1];
    Geps             = eps*vftab[ntab+2];
    Heps2            = eps2*vftab[ntab+3];
    Fp               = F+Geps+Heps2;
    VVe              = Y+eps*Fp;
    FFe              = Fp+Geps+2.0f*Heps2;
    /* Dispersion */
    Y                = vftab[ntab+4];
    F                = vftab[ntab+5];
    Geps             = eps*vftab[ntab+6];
    Heps2            = eps2*vftab[ntab+7];
    Fp               = F+Geps+Heps2;
    VVd              = Y+eps*Fp;
    FFd              = Fp+Geps+2.0f*Heps2;
    /* Repulsion */
    Y                = vftab[ntab+8];
    F                = vftab[ntab+9];
    Geps             = eps*vftab[ntab+10];
    Heps2            = eps2*vftab[ntab+11];
    Fp               = F+Geps+Heps2;
    VVr              = Y+eps*Fp;
    FFr              = Fp+Geps+2.0f*Heps2;

    *velec           = qq*VVe;
    *vvdw            = c6*VVd+c12*VVr;
    fscal            = -(qq*FFe+c6*FFd+c12*FFr)*tabscale*rinv;

    return fscal;
}

__global__
void pairs_gpu(int ftype, int nbonds, int natoms,
               const t_iatom iatoms[], const t_iparams iparams[],
               const fvec x[], fvec force[], fvec fshift[],
               const t_pbc *pbc,  const fvec *pbc_hbox_diag, const matrix *pbc_box,
               const fvec *pbc_mhbox_diag, const fvec  *pbc_fbox_diag,
               const fvec*  pbc_tric_vec,  const ivec*  pbc_tric_shift,
               float md_chargeA[], unsigned short md_cENER[],
               float energygrp_elec[], float energygrp_vdw[],
               float fr_pairsTable_data[],
               float packed_R[], int packed_I[])
{
    float             qq, c6, c12;
    fvec              dx;
    int               i, itype, ai, aj, gid;
    int               fshift_index;
    float             r2;
    float             fscal, velec, vvdw;

    const float       epsfac = packed_R[fr_icEPSFAC];

    __shared__ fvec  fshift_loc[SHIFTS];
    
    i = blockIdx.x*blockDim.x+threadIdx.x;
    if (threadIdx.x < SHIFTS)
    {
        fshift_loc[threadIdx.x][XX] = 0.0f;
        fshift_loc[threadIdx.x][YY] = 0.0f;
        fshift_loc[threadIdx.x][ZZ] = 0.0f;
    }

    __syncthreads();

    if (i <  nbonds)
    {
        itype = iatoms[3*i];
        ai    = iatoms[3*i+1];
        aj    = iatoms[3*i+2];
        gid   = GID(md_cENER[ai], md_cENER[aj], packed_I[md_NENERGRP]);

        qq               = md_chargeA[ai]*md_chargeA[aj]*epsfac*packed_R[frFUDGEQQ];
        c6               = iparams[itype].lj14.c6A;
        c12              = iparams[itype].lj14.c12A;

        /* To save flops in the optimized kernels, c6/c12 have 6.0/12.0 derivative prefactors
         * included in the general nfbp array now. This means the tables are scaled down by the
         * same factor, so when we use the original c6/c12 parameters from iparams[] they must
         * be scaled up.
         */
        c6  *= 6.0f;
        c12 *= 12.0f;

        /* Do we need to apply full periodic boundary conditions? */
        fshift_index = pbc_fvec_sub_gpu(pbc, x[ai], x[aj], dx,
                                        pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                                        pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift);

        r2 = norm2_gpu(dx);

        /* Evaluate tabulated interaction without free energy */
        fscal = evaluate_single(r2, packed_R[fr_PAIRS_SCALE],
                                fr_pairsTable_data, packed_I[fr_PAIRS_STRIDE],
                                qq, c6, c12, &velec, &vvdw);

        atomicAdd(&energygrp_elec[gid], velec);
        atomicAdd(&energygrp_vdw[gid], vvdw);

        svmul_gpu(fscal, dx, dx);

        /* Add the forces */
        for (int m = 0; m < DIM; m++)
        {
            atomicAdd(&force[ai][m], dx[m]);
            atomicAdd(&force[aj][m], -dx[m]);
        }

        if (fshift_index != CENTRAL)
        {
            fvec_inc_atomic(fshift_loc[fshift_index], dx);
            fvec_dec_atomic(fshift_loc[CENTRAL], dx);
        }
    }
    __syncthreads();
    if (threadIdx.x < SHIFTS)
    {
        fvec_inc_atomic(fshift[threadIdx.x], fshift_loc[threadIdx.x]);
    }
}

__global__
void pairs_gpu_noener(int ftype, int nbonds, int natoms,
                      const t_iatom iatoms[], const t_iparams iparams[],
                      const fvec x[], fvec force[], fvec fshift[],
                      const t_pbc *pbc,  const fvec *pbc_hbox_diag, const matrix *pbc_box,
                      const fvec *pbc_mhbox_diag, const fvec  *pbc_fbox_diag,
                      const fvec*  pbc_tric_vec,  const ivec*  pbc_tric_shift,
                      float md_chargeA[], unsigned short md_cENER[],
                      float energygrp_elec[], float energygrp_vdw[],
                      float fr_pairsTable_data[],
                      float packed_R[], int packed_I[])
{
    float             qq, c6, c12;
    fvec              dx;
    int               i, itype, ai, aj;
    float             r2;
    float             fscal, velec, vvdw;

    const float       epsfac = packed_R[fr_icEPSFAC];

    i = blockIdx.x*blockDim.x+threadIdx.x;

    if (i <  nbonds)
    {
        itype = iatoms[3*i];
        ai    = iatoms[3*i+1];
        aj    = iatoms[3*i+2];

        qq               = md_chargeA[ai]*md_chargeA[aj]*epsfac*packed_R[frFUDGEQQ];
        c6               = iparams[itype].lj14.c6A;
        c12              = iparams[itype].lj14.c12A;

        /* To save flops in the optimized kernels, c6/c12 have 6.0/12.0 derivative prefactors
         * included in the general nfbp array now. This means the tables are scaled down by the
         * same factor, so when we use the original c6/c12 parameters from iparams[] they must
         * be scaled up.
         */
        c6  *= 6.0f;
        c12 *= 12.0f;

        /* Do we need to apply full periodic boundary conditions? */
        pbc_fvec_sub_gpu(pbc, x[ai], x[aj], dx,
                         pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                         pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift);

        r2 = norm2_gpu(dx);

        /* Evaluate tabulated interaction without free energy */
        fscal = evaluate_single(r2, packed_R[fr_PAIRS_SCALE],
                                fr_pairsTable_data, packed_I[fr_PAIRS_STRIDE],
                                qq, c6, c12, &velec, &vvdw);
        svmul_gpu(fscal, dx, dx);

        /* Add the forces */
        for (int m = 0; m < DIM; m++)
        {
            atomicAdd(&force[ai][m], dx[m]);
            atomicAdd(&force[aj][m], -dx[m]);
        }

    }
}

// ----- //
__global__
void reset_next(int natoms, int force_next[] )
{
    int  i = blockIdx.x*blockDim.x+threadIdx.x;
    if (i < natoms)
    {
        force_next[i] = 0;
    }
}

/*-------------------------------- End CUDA kernels-----------------------------*/



/*------------------------------------------------------------------------------*/

/* static variables internal to module. Separate instances are needed
 * for each MPI rank (because, when using threadMPI, we need to avoid
 * interference and allow a global view of memory)
 */

/* variables to be integrated */
#ifdef STANDALONE
static fvec *x_bonded_d[BO_MAX_RANKS];
static fvec *f_bonded_d[BO_MAX_RANKS];
#endif

static fvec *f_bonded[BO_MAX_RANKS];
static fvec *f_shift_d[BO_MAX_RANKS];
static fvec *f_shift[BO_MAX_RANKS];

static float vtot[BO_MAX_RANKS][F_NRE];

static float *vtot_d[BO_MAX_RANKS];
static t_iparams *forceparams_d[BO_MAX_RANKS];

// TODO : we have duplicate data for constraints
static t_pbc  *pbc_d[BO_MAX_RANKS];
static fvec   *pbc_hbox_diag_d[BO_MAX_RANKS];
static matrix *pbc_box_d[BO_MAX_RANKS];
static fvec   *pbc_mhbox_diag_d[BO_MAX_RANKS];
static fvec   *pbc_fbox_diag_d[BO_MAX_RANKS];
static fvec   *pbc_tric_vec_d[BO_MAX_RANKS];
static ivec   *pbc_tric_shift_d[BO_MAX_RANKS];

static t_iatom *iatoms_d[BO_MAX_RANKS][F_NRE];
static int      device_iatom_alloc[BO_MAX_RANKS][F_NRE];

int            xbonded_size_max[BO_MAX_RANKS];
static bool    bonded_init_done[BO_MAX_RANKS];
static int     ntypes_alloc[BO_MAX_RANKS];

static float  *pairs_Rpacked[BO_MAX_RANKS]; // packed array of relavant values
static float  *pairs_Rpacked_d[BO_MAX_RANKS];
static int    *pairs_Ipacked[BO_MAX_RANKS];
static int    *pairs_Ipacked_d[BO_MAX_RANKS];

static float  *md_chargeA_d[BO_MAX_RANKS];

static float  *energygrp_elec_d[BO_MAX_RANKS];
static float  *energygrp_elec[BO_MAX_RANKS];
static float  *energygrp_vdw_d[BO_MAX_RANKS];
static float  *energygrp_vdw[BO_MAX_RANKS];

static float  *fr_pairsTable_data_d[BO_MAX_RANKS];

static unsigned short *md_cENER_d[BO_MAX_RANKS];

static int    *force_next_d[BO_MAX_RANKS][F_NRE];
static int    *force_max_d[BO_MAX_RANKS];
static int     force_max[BO_MAX_RANKS];

static cudaStream_t    streamBonded[BO_MAX_RANKS];

static bool    copy_required_this_step[BO_MAX_RANKS];

/*------------------------------------------------------------------------------*/
// bonded forces //

// initial setup (called once)
static void
init_gpu_bonded( t_pbc *pbc, int rank, gmx_grppairener_t *grppener, const t_forcerec *fr)
{
    cudaError_t stat;

// device_iatom_alloc
// initiaized to 0 as static

    pbc_d[rank] = NULL;
// null status is used as a marker for no pbc
    if (pbc)
    {
        stat = cudaMalloc(&pbc_d[rank], sizeof(t_pbc));
        CU_RET_ERR(stat, "cudaMalloc failed");
        stat = cudaMalloc(&pbc_hbox_diag_d[rank], sizeof(fvec));
        CU_RET_ERR(stat, "cudaMalloc failed");
        stat = cudaMalloc(&pbc_box_d[rank], sizeof(matrix));
        CU_RET_ERR(stat, "cudaMalloc failed");
        stat = cudaMalloc(&pbc_mhbox_diag_d[rank], sizeof(fvec));
        CU_RET_ERR(stat, "cudaMalloc failed");
        stat = cudaMalloc(&pbc_fbox_diag_d[rank], sizeof(fvec));
        CU_RET_ERR(stat, "cudaMalloc failed");
        stat = cudaMalloc(&pbc_tric_vec_d[rank], MAX_NTRICVEC*sizeof(fvec));
        CU_RET_ERR(stat, "cudaMalloc failed");
        stat = cudaMalloc(&pbc_tric_shift_d[rank], MAX_NTRICVEC*sizeof(ivec));
        CU_RET_ERR(stat, "cudaMalloc failed");
    }
// vtot
    stat = cudaMalloc(&vtot_d[rank],F_NRE*sizeof(real));
    CU_RET_ERR(stat, "cudaMalloc failed");

// ener arrays
    stat = cudaMalloc(&pairs_Rpacked_d[rank], PAIRS_NR*sizeof(float));
    CU_RET_ERR(stat, "cudaMalloc failed");

// packed arrays
    stat = cudaMalloc(&pairs_Rpacked_d[rank], PAIRS_NR*sizeof(float));
    CU_RET_ERR(stat, "cudaMalloc failed");
    pairs_Rpacked[rank] = (float *) malloc(sizeof(float)*PAIRS_NR);

    stat = cudaMalloc(&pairs_Ipacked_d[rank], PAIRS_NI*sizeof(int));
    CU_RET_ERR(stat, "cudaMalloc failed");
    pairs_Ipacked[rank] = (int *) malloc(sizeof(int)*PAIRS_NI);

// possibly these change size, but hopefully not
    stat = cudaMalloc(&energygrp_elec_d[rank], grppener->nener*sizeof(float));
    CU_RET_ERR(stat, "cudaMalloc failed");

    stat = cudaHostAlloc(&energygrp_elec[rank], grppener->nener*sizeof(float), cudaHostAllocPortable);
    CU_RET_ERR(stat, "cudaHostAlloc failed");

    stat = cudaMalloc(&energygrp_vdw_d[rank], grppener->nener*sizeof(float));
    CU_RET_ERR(stat, "cudaMalloc failed");


    stat = cudaHostAlloc(&energygrp_vdw[rank], grppener->nener*sizeof(float), cudaHostAllocPortable);
    CU_RET_ERR(stat, "cudaHostAlloc failed");


    stat = cudaMalloc(&fr_pairsTable_data_d[rank], (fr->pairsTable->n + 1)*fr->pairsTable->stride*sizeof(float));
    CU_RET_ERR(stat, "cudaMalloc failed");

    stat = cudaMemcpy(fr_pairsTable_data_d[rank], fr->pairsTable->data,
                      (fr->pairsTable->n + 1)*fr->pairsTable->stride*sizeof(float), cudaMemcpyHostToDevice);
    CU_RET_ERR(stat, "cudaMemcpy failed");

    stat = cudaMalloc(&f_shift_d[rank], sizeof(fvec)*SHIFTS);
    CU_RET_ERR(stat, "cudaMalloc failed");

    f_shift[rank] = (fvec *) malloc(sizeof(fvec)*SHIFTS);

    stat = cudaMalloc(&force_max_d[rank], sizeof(int));
    CU_RET_ERR(stat, "cudaMalloc failed");

    stat = cudaStreamCreate(&streamBonded[rank]);
    CU_RET_ERR(stat, "cudaStreamCreate failed");
}

//----//
// update after a neighbour list update step
void
update_gpu_bonded(const t_idef *idef, const t_forcerec *fr, matrix box,
                  int xsize, const t_mdatoms *md, gmx_grppairener_t *grppener)
{
    cudaError_t stat;
    int         rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const float alloc_factor = 1.2;
    int         xsize_alloc;

    copy_required_this_step[rank] = true;

    t_pbc *pbc_null = NULL;
    t_pbc  pbc;
    set_pbc(&pbc, fr->ePBC, box);

    if (fr->bMolPBC)
    {
        pbc_null = &pbc;
    }
    else
    {
        pbc_null = NULL;
    }                            // there is probaly a better way to do this

    if (!bonded_init_done[rank]) //static variables init to 0 i.e. false
    {
        init_gpu_bonded(pbc_null, rank, grppener, fr);
        bonded_init_done[rank] = true;
    }

    bool x_resized = false;
    xsize_alloc = 0;
    int  x_needed = (xsize > md->nr) ? xsize : md->nr;
    // x + F_bonded_host (temporary for testing), others needed in final version
    if (x_needed > xbonded_size_max[rank])
    {
        x_resized   = true;
        xsize_alloc = alloc_factor*x_needed;
        if (md_chargeA_d[rank])
        {
#ifdef STANDALONE
            stat = cudaFree(x_bonded_d[rank]);
            CU_RET_ERR(stat, "cudaFree failed");
            stat = cudaFree(f_bonded_d[rank]);
            CU_RET_ERR(stat, "cudaFree failed");
#endif
            free(f_bonded[rank]);
            stat = cudaFree(md_cENER_d[rank]);
            CU_RET_ERR(stat, "cudaFree failed");
            stat = cudaFree(md_chargeA_d[rank]);
            CU_RET_ERR(stat, "cudaFree failed");
        }
#ifdef STANDALONE
        stat = cudaMalloc(&x_bonded_d[rank], sizeof(fvec)*xsize_alloc);
        CU_RET_ERR(stat, "cudaMalloc failed");
        stat = cudaMalloc(&f_bonded_d[rank], sizeof(fvec)*xsize_alloc);
        CU_RET_ERR(stat, "cudaMalloc failed");
#endif
        f_bonded[rank] = (fvec *) malloc(sizeof(fvec)*xsize_alloc);
        stat           = cudaMalloc(&md_cENER_d[rank], sizeof(unsigned short)*xsize_alloc);
        CU_RET_ERR(stat, "cudaMalloc failed");
        stat = cudaMalloc(&md_chargeA_d[rank], sizeof(float)*xsize_alloc);
        CU_RET_ERR(stat, "cudaMalloc failed");

        xbonded_size_max[rank] = xsize_alloc;
    }

    stat = cudaMemcpy(md_cENER_d[rank], md->cENER, md->nr*sizeof(unsigned short), cudaMemcpyHostToDevice);
    CU_RET_ERR(stat, "cudaMemcpy failed");
    stat = cudaMemcpy(md_chargeA_d[rank], md->chargeA, md->nr*sizeof(float), cudaMemcpyHostToDevice);
    CU_RET_ERR(stat, "cudaMemcpy failed");

    //forceparams
    int              ntypes = idef->ntypes;
    const t_iparams *iparams;
    iparams = idef->iparams;

    if (ntypes > ntypes_alloc[rank])
    {
        if (forceparams_d[rank])
        {
            stat = cudaFree(forceparams_d[rank]);
            CU_RET_ERR(stat, "cudaFree failed");
        }
        stat = cudaMalloc(&forceparams_d[rank], sizeof(t_iparams)*ntypes);
        CU_RET_ERR(stat, "cudaMalloc failed");
        ntypes_alloc[rank] = ntypes;
    }

    stat = cudaMemcpyAsync(forceparams_d[rank], iparams, ntypes*sizeof(t_iparams), cudaMemcpyHostToDevice, streamBonded[rank]);

    // forceatoms
    for (int i = 0; i < F_NRE; i++)
    {
        if (idef->il[i].nr > 0 && ftype_is_bonded_potential(i))
        {
            if (x_resized)
            {
                if (force_next_d[rank][i])
                {
                    stat = cudaFree(force_next_d[rank][i]);
                    CU_RET_ERR(stat, "cudaFree failed");
                }
                stat = cudaMalloc(&force_next_d[rank][i], sizeof(int)*xsize_alloc);
                CU_RET_ERR(stat, "cudaMalloc failed");
            }

            //floatlocate if we do not have enough memory
            if (idef->il[i].nalloc > device_iatom_alloc[rank][i])
            {
                if (iatoms_d[rank][i])
                {
                    stat = cudaFree(iatoms_d[rank][i]);
                    CU_RET_ERR(stat, "cudaFree failed");
                }
                stat = cudaMalloc(&iatoms_d[rank][i], idef->il[i].nalloc*sizeof(t_iatom));
                CU_RET_ERR(stat, "cudaMalloc failed");
                device_iatom_alloc[rank][i] = idef->il[i].nalloc;
            }
            //copy  the data
            stat = cudaMemcpy(iatoms_d[rank][i], idef->il[i].iatoms, idef->il[i].nr*sizeof(t_iatom),
                              cudaMemcpyHostToDevice);
            CU_RET_ERR(stat, "cudaMemcpy failed");
        }
    }
    //pbc
    //this can be null
    if (pbc_null)
    {
        stat = cudaMemcpy(pbc_d[rank], pbc_null, sizeof(t_pbc), cudaMemcpyHostToDevice);
        CU_RET_ERR(stat, "cudaMemcpy failed");
        stat = cudaMemcpy(pbc_hbox_diag_d[rank], pbc_null->hbox_diag, sizeof(fvec), cudaMemcpyHostToDevice);
        CU_RET_ERR(stat, "cudaMemcpy failed");
        stat = cudaMemcpy(pbc_box_d[rank], pbc_null->box, sizeof(matrix), cudaMemcpyHostToDevice);
        CU_RET_ERR(stat, "cudaMemcpy failed");
        stat = cudaMemcpy(pbc_mhbox_diag_d[rank], pbc_null->mhbox_diag, sizeof(fvec), cudaMemcpyHostToDevice);
        CU_RET_ERR(stat, "cudaMemcpy failed");
        stat = cudaMemcpy(pbc_fbox_diag_d[rank], pbc_null->fbox_diag, sizeof(fvec), cudaMemcpyHostToDevice);
        CU_RET_ERR(stat, "cudaMemcpy failed");
        stat = cudaMemcpy(pbc_tric_vec_d[rank], pbc_null->tric_vec, MAX_NTRICVEC*sizeof(fvec), cudaMemcpyHostToDevice);
        CU_RET_ERR(stat, "cudaMemcpy failed");
        stat = cudaMemcpy(pbc_tric_shift_d[rank], pbc_null->tric_shift, MAX_NTRICVEC*sizeof(ivec), cudaMemcpyHostToDevice);
        CU_RET_ERR(stat, "cudaMemcpy failed");
    }

    // Packed arrays for LJ14

    pairs_Rpacked[rank][frFUDGEQQ]      = fr->fudgeQQ;
    pairs_Rpacked[rank][fr_icEPSFAC]    = fr->ic->epsfac;
    pairs_Rpacked[rank][fr_PAIRS_SCALE] = fr->pairsTable->scale;

    stat = cudaMemcpy(pairs_Rpacked_d[rank], pairs_Rpacked[rank],
                      sizeof(float)*PAIRS_NR, cudaMemcpyHostToDevice);
    CU_RET_ERR(stat, "cudaMemcpy failed");

    pairs_Ipacked[rank][fr_PAIRS_STRIDE] = fr->pairsTable->stride;
    pairs_Ipacked[rank][md_NENERGRP]     = md->nenergrp;

    stat = cudaMemcpy(pairs_Ipacked_d[rank], pairs_Ipacked[rank],
                      sizeof(int)*PAIRS_NI, cudaMemcpyHostToDevice);
    CU_RET_ERR(stat, "cudaMemcpy failed");

    dim3 blocks, blocks_reset;
    dim3 threads (TPB_BONDED, 1, 1);
    blocks_reset.x = (xbonded_size_max[rank]+TPB_BONDED-1)/TPB_BONDED;

    stat = cudaMemcpyAsync(&force_max[rank], force_max_d[rank],
                           sizeof(int), cudaMemcpyDeviceToHost, streamBonded[rank]);
    cudaStreamSynchronize(streamBonded[rank]);
    stat = cudaGetLastError();
    CU_RET_ERR(stat, "Async error");

    cudaStreamSynchronize(streamBonded[rank]);
}


/* -------------- */
// reset the internal variables
void
reset_gpu_bonded(const int size, const int nener)
{
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    dim3 blocks ( (size+TPB_BONDED-1)/TPB_BONDED, 1, 1);
    dim3 threads (TPB_BONDED, 1, 1);

    cudaError_t stat;

#ifdef STANDALONE
    fvec* f_d_in = f_bonded_d[rank];
#else
    fvec* f_d_in = gpuBufferOpsGetFPtr();
#endif

    reset_gpu_bonded_kernel <<< blocks, threads, 0, streamBonded[rank]>>>
    (vtot_d[rank], f_d_in, size, energygrp_elec_d[rank],
     energygrp_vdw_d[rank], nener, f_shift_d[rank]);

    stat = cudaGetLastError();
    CU_RET_ERR(stat, "reset bonded kernel failed");
}

/* -------------- */
void
do_bonded_gpu(t_forcerec *fr, const t_inputrec *ir, const t_idef *idef,
              int flags,  const t_graph *graph, int natoms, fvec x[])
{
    int nat1, nbonds;
    bool  bCalcEnerVir = (flags & (GMX_FORCE_VIRIAL | GMX_FORCE_ENERGY));

    bool abort = false;
    int  rank  = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#ifdef STANDALONE
    fvec* x_d_in = x_bonded_d[rank];
    fvec* f_d_in = f_bonded_d[rank];
#else
    fvec* x_d_in = gpuBufferOpsGetXPtr();
    fvec* f_d_in = gpuBufferOpsGetFPtr();
#endif

    cudaError_t stat;

// sanity check

    if (fr->bQMMM)
    {
        abort = true;
    }
    if (ir->nwall)
    {
        abort = true;
    }
    if (ir->implicit_solvent)
    {
        abort = true;
    }
    if ((fr->cutoff_scheme == ecutsGROUP) && (flags & GMX_FORCE_NONBONDED))
    {
        abort = true;
    }
    if (graph)
    {
        abort = true;
    }
    if (flags & GMX_FORCE_DHDL)
    {
        abort = true;
    }
//     if (fr->efep != efepNO) abort=true;

    if (abort)
    {
        printf("calculation not supported on GPU\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    bool do_copy_this_step = false;

#ifdef STANDALONE
    copy_required_this_step[rank] = true;
#endif

    if (copy_required_this_step[rank])
    {
        do_copy_this_step = true;
    }

    if (do_copy_this_step)
    {
        stat = cudaMemcpyAsync(x_d_in, x,
                               natoms*sizeof(fvec), cudaMemcpyHostToDevice, streamBonded[rank]);
        CU_RET_ERR(stat, "cudaMemcpy failed");
    }

// launch kernels
// reordered for better overlap
    dim3 blocks, blocks_natoms;
    dim3 threads (TPB_BONDED, 1, 1);
    for (int ftype = 0; (ftype < F_NRE); ftype++)
    {
        if (idef->il[ftype].nr > 0 && ftype_is_bonded_potential(ftype))
        {
            nat1            = interaction_function[ftype].nratoms + 1;
            nbonds          = idef->il[ftype].nr/nat1;
            blocks.x        = (nbonds+TPB_BONDED-1)/TPB_BONDED;
            blocks_natoms.x = (natoms+TPB_BONDED-1)/TPB_BONDED;
            if (ftype == F_PDIHS || ftype == F_PIDIHS)
            {
                if (bCalcEnerVir)
                {
                    pdihs_gpu<true> <<< blocks, threads, 0, streamBonded[rank]>>>
                    (&vtot_d[rank][ftype], nbonds, natoms,
                    iatoms_d[rank][ftype], forceparams_d[rank],
                    x_d_in, f_d_in, f_shift_d[rank],
                    pbc_d[rank], pbc_hbox_diag_d[rank],
                    pbc_box_d[rank], pbc_mhbox_diag_d[rank],
                    pbc_fbox_diag_d[rank], pbc_tric_vec_d[rank],
                    pbc_tric_shift_d[rank]);
                }
                else
                {
                    pdihs_gpu<false> <<< blocks, threads, 0, streamBonded[rank]>>>
                    (nullptr, nbonds, natoms,
                    iatoms_d[rank][ftype], forceparams_d[rank],
                    x_d_in, f_d_in, nullptr,
                    pbc_d[rank], pbc_hbox_diag_d[rank],
                    pbc_box_d[rank], pbc_mhbox_diag_d[rank],
                    pbc_fbox_diag_d[rank], pbc_tric_vec_d[rank],
                    pbc_tric_shift_d[rank]);
                }
            }
        }
    }

    for (int ftype = 0; (ftype < F_NRE); ftype++)
    {
        if (idef->il[ftype].nr > 0 && ftype_is_bonded_potential(ftype))
        {
            nat1            = interaction_function[ftype].nratoms + 1;
            nbonds          = idef->il[ftype].nr/nat1;
            blocks.x        = (nbonds+TPB_BONDED-1)/TPB_BONDED;
            blocks_natoms.x = (natoms+TPB_BONDED-1)/TPB_BONDED;
// in the main code they have a function pointer for this
// so we need something like that for final version
// note temp array here for output forces

            if (ftype == F_BONDS)
            {
                if (bCalcEnerVir)
                {
                    bonds_gpu<true> <<< blocks, threads, 0, streamBonded[rank]>>>
                    (&vtot_d[rank][ftype], nbonds, natoms,
                    iatoms_d[rank][ftype], forceparams_d[rank],
                    x_d_in, f_d_in, f_shift_d[rank],
                    pbc_d[rank], pbc_hbox_diag_d[rank],
                    pbc_box_d[rank], pbc_mhbox_diag_d[rank],
                    pbc_fbox_diag_d[rank], pbc_tric_vec_d[rank],
                    pbc_tric_shift_d[rank]);
                }
                else
                {
                    bonds_gpu<false> <<< blocks, threads, 0, streamBonded[rank]>>>
                    (nullptr, nbonds, natoms,
                    iatoms_d[rank][ftype], forceparams_d[rank],
                    x_d_in, f_d_in, nullptr,
                    pbc_d[rank], pbc_hbox_diag_d[rank],
                    pbc_box_d[rank], pbc_mhbox_diag_d[rank],
                    pbc_fbox_diag_d[rank], pbc_tric_vec_d[rank],
                    pbc_tric_shift_d[rank]);
                }
            }

            if (ftype == F_ANGLES)
            {
                if (bCalcEnerVir)
                {
                    angles_gpu<true> <<< blocks, threads, 0, streamBonded[rank]>>>
                    (&vtot_d[rank][ftype], nbonds, natoms,
                    iatoms_d[rank][ftype], forceparams_d[rank],
                    x_d_in, f_d_in, f_shift_d[rank],
                    pbc_d[rank], pbc_hbox_diag_d[rank],
                    pbc_box_d[rank], pbc_mhbox_diag_d[rank],
                    pbc_fbox_diag_d[rank], pbc_tric_vec_d[rank],
                    pbc_tric_shift_d[rank]);
                }
                else
                {
                    angles_gpu<false> <<< blocks, threads, 0, streamBonded[rank]>>>
                    (nullptr, nbonds, natoms,
                    iatoms_d[rank][ftype], forceparams_d[rank],
                    x_d_in, f_d_in, nullptr,
                    pbc_d[rank], pbc_hbox_diag_d[rank],
                    pbc_box_d[rank], pbc_mhbox_diag_d[rank],
                    pbc_fbox_diag_d[rank], pbc_tric_vec_d[rank],
                    pbc_tric_shift_d[rank]);
                }
            }

            if (ftype == F_UREY_BRADLEY)
            {
                if (bCalcEnerVir)
                {
                    urey_bradley_gpu <<< blocks, threads, 0, streamBonded[rank]>>>
                    (&vtot_d[rank][ftype], nbonds, natoms,
                    iatoms_d[rank][ftype], forceparams_d[rank],
                    x_d_in, f_d_in, f_shift_d[rank],
                    pbc_d[rank], pbc_hbox_diag_d[rank],
                    pbc_box_d[rank], pbc_mhbox_diag_d[rank],
                    pbc_fbox_diag_d[rank], pbc_tric_vec_d[rank],
                    pbc_tric_shift_d[rank]);
                }
                else
                {
                    urey_bradley_gpu_noener <<< blocks, threads, 0, streamBonded[rank]>>>
                    (nbonds, natoms,
                    iatoms_d[rank][ftype], forceparams_d[rank],
                    x_d_in, f_d_in,
                    pbc_d[rank], pbc_hbox_diag_d[rank],
                    pbc_box_d[rank], pbc_mhbox_diag_d[rank],
                    pbc_fbox_diag_d[rank], pbc_tric_vec_d[rank],
                    pbc_tric_shift_d[rank]);
                }
            }

            if (ftype == F_RBDIHS)
            {
                if (bCalcEnerVir)
                {
                    rbdihs_gpu<true> <<< blocks, threads, 0, streamBonded[rank]>>>
                    (&vtot_d[rank][ftype], nbonds, natoms,
                    iatoms_d[rank][ftype], forceparams_d[rank],
                    x_d_in, f_d_in, f_shift_d[rank],
                    pbc_d[rank], pbc_hbox_diag_d[rank],
                    pbc_box_d[rank], pbc_mhbox_diag_d[rank],
                    pbc_fbox_diag_d[rank], pbc_tric_vec_d[rank],
                    pbc_tric_shift_d[rank]);
                }
                else
                {
                    rbdihs_gpu<false> <<< blocks, threads, 0, streamBonded[rank]>>>
                    (nullptr, nbonds, natoms,
                    iatoms_d[rank][ftype], forceparams_d[rank],
                    x_d_in, f_d_in, nullptr,
                    pbc_d[rank], pbc_hbox_diag_d[rank],
                    pbc_box_d[rank], pbc_mhbox_diag_d[rank],
                    pbc_fbox_diag_d[rank], pbc_tric_vec_d[rank],
                    pbc_tric_shift_d[rank]);
                }
            }

            if (ftype == F_IDIHS)
            {
                if (bCalcEnerVir)
                {
                    idihs_gpu<true> <<< blocks, threads, 0, streamBonded[rank]>>>
                    (&vtot_d[rank][ftype], nbonds, natoms,
                    iatoms_d[rank][ftype], forceparams_d[rank],
                    x_d_in, f_d_in, f_shift_d[rank],
                    pbc_d[rank], pbc_hbox_diag_d[rank],
                    pbc_box_d[rank], pbc_mhbox_diag_d[rank],
                    pbc_fbox_diag_d[rank], pbc_tric_vec_d[rank],
                    pbc_tric_shift_d[rank]);
                }
                else
                {
                    idihs_gpu<false> <<< blocks, threads, 0, streamBonded[rank]>>>
                    (nullptr, nbonds, natoms,
                    iatoms_d[rank][ftype], forceparams_d[rank],
                     x_d_in, f_d_in, nullptr,
                    pbc_d[rank], pbc_hbox_diag_d[rank],
                    pbc_box_d[rank], pbc_mhbox_diag_d[rank],
                    pbc_fbox_diag_d[rank], pbc_tric_vec_d[rank],
                    pbc_tric_shift_d[rank]);
                }
            }

            if (ftype == F_LJ14)
            {
                if (bCalcEnerVir)
                {
                    pairs_gpu <<< blocks, threads, 0, streamBonded[rank]>>>
                    (ftype, nbonds, natoms,
                    iatoms_d[rank][ftype], forceparams_d[rank],
                    x_d_in, f_d_in, f_shift_d[rank],
                    pbc_d[rank], pbc_hbox_diag_d[rank],
                    pbc_box_d[rank], pbc_mhbox_diag_d[rank],
                    pbc_fbox_diag_d[rank], pbc_tric_vec_d[rank],
                    pbc_tric_shift_d[rank],
                    md_chargeA_d[rank], md_cENER_d[rank],
                    energygrp_elec_d[rank], energygrp_vdw_d[rank],
                    fr_pairsTable_data_d[rank],
                    pairs_Rpacked_d[rank], pairs_Ipacked_d[rank]);
                }
                else
                {
                    pairs_gpu_noener <<< blocks, threads, 0, streamBonded[rank]>>>
                    (ftype, nbonds, natoms,
                    iatoms_d[rank][ftype], forceparams_d[rank],
                    x_d_in, f_d_in, f_shift_d[rank],
                    pbc_d[rank], pbc_hbox_diag_d[rank],
                    pbc_box_d[rank], pbc_mhbox_diag_d[rank],
                    pbc_fbox_diag_d[rank], pbc_tric_vec_d[rank],
                    pbc_tric_shift_d[rank],
                    md_chargeA_d[rank], md_cENER_d[rank],
                    energygrp_elec_d[rank], energygrp_vdw_d[rank],
                    fr_pairsTable_data_d[rank],
                    pairs_Rpacked_d[rank], pairs_Ipacked_d[rank]);
                }
            }
        }
    }
    stat = cudaGetLastError();
    CU_RET_ERR(stat, "kernels failed");
    return;
}

void
do_bonded_gpu_finalize(t_forcerec *fr, int flags, int natoms,
                       rvec *input_force, gmx_enerdata_t *enerd)
{
    int                i;
    bool               bCalcEnerVir = (flags & (GMX_FORCE_VIRIAL | GMX_FORCE_ENERGY));
    cudaError_t        stat;

    gmx_grppairener_t *grppener;
    grppener = &enerd->grpp;

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#ifdef STANDALONE
    fvec* f_d_in = f_bonded_d[rank];
#else
    fvec* f_d_in = gpuBufferOpsGetFPtr();
#endif

    bool do_copy_this_step = false;

    if (copy_required_this_step[rank])
    {
        do_copy_this_step             = true;
        copy_required_this_step[rank] = false;
    }

//copy forces back
    if (do_copy_this_step)
    {
        stat = cudaMemcpyAsync(f_bonded[rank], f_d_in,
                               natoms*sizeof(fvec), cudaMemcpyDeviceToHost, streamBonded[rank]);
        cudaStreamSynchronize(streamBonded[rank]);
        CU_RET_ERR(stat, "cudaMemcpy failed");

        for (i = 0; i < natoms; i++)
        {
            rvec_inc(input_force[i], f_bonded[rank][i]);
        }
    }

    if (bCalcEnerVir)
    {
        // shift Forces
        stat = cudaMemcpyAsync(f_shift[rank], f_shift_d[rank],
                               SHIFTS*sizeof(fvec), cudaMemcpyDeviceToHost, streamBonded[rank]);
        cudaStreamSynchronize(streamBonded[rank]);
        CU_RET_ERR(stat, "cudaMemcpy failed");

        for (i = 0; i < SHIFTS; i++)
        {
            rvec_inc(fr->fshift[i], f_shift[rank][i]);
        }

        // copy energies back
        stat = cudaMemcpyAsync(vtot[rank], vtot_d[rank],
                               F_NRE*sizeof(float), cudaMemcpyDeviceToHost, streamBonded[rank]);
        cudaStreamSynchronize(streamBonded[rank]);
        CU_RET_ERR(stat, "cudaMemcpy failed");

        for (i = 0; i < F_NRE; i++)
        {
            enerd->term[i] += vtot[rank][i];
        }

        stat = cudaMemcpyAsync(energygrp_vdw[rank], energygrp_vdw_d[rank],
                               sizeof(float)*grppener->nener, cudaMemcpyDeviceToHost, streamBonded[rank]);
        CU_RET_ERR(stat, "cudaMemcpy failed");

        stat = cudaMemcpyAsync(energygrp_elec[rank], energygrp_elec_d[rank],
                               sizeof(float)*grppener->nener, cudaMemcpyDeviceToHost, streamBonded[rank]);
        CU_RET_ERR(stat, "cudaMemcpy failed");
        cudaStreamSynchronize(streamBonded[rank]);
        for (i = 0; i < grppener->nener; i++)
        {
            grppener->ener[egCOUL14][i] += energygrp_elec[rank][i];
            grppener->ener[egLJ14][i]   +=  energygrp_vdw[rank][i];
        }
    }
    cudaStreamSynchronize(streamBonded[rank]); // needed if we have no copies for bufferops sync
}
/*------------------------------------------------------------------------------*/
