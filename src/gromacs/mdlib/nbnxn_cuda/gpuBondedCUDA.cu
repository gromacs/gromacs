/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
#include "gromacs/utility/real.h"
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
    return
        ((interaction_function[ftype].flags & IF_BOND) != 0u) &&
        !(ftype == F_CONNBONDS || ftype == F_POSRES || ftype == F_FBPOSRES);
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
reset_gpu_bonded_kernel(real *vtot, real *dvd_lambda, rvec *force , int size, 
          real *energygrp_elec, real *energygrp_vdw, int nener, rvec *f_shift )
{
  int a=blockIdx.x*blockDim.x+threadIdx.x;
 
  if(a<F_NRE){
    vtot[a]=0.0;
  }

  if(a<egNR){
   dvd_lambda[a]=0.0;
  }

  if(a<size){
    force[a][XX]=0.0;
    force[a][YY]=0.0;
    force[a][ZZ]=0.0;
  }

  if(a<nener){
   energygrp_vdw[a]=0.0;
   energygrp_elec[a]=0.0;

  }

  if(a<SHIFTS)
  {
    f_shift[a][XX]=0.0;
    f_shift[a][YY]=0.0;
    f_shift[a][ZZ]=0.0;
  }
}

__device__
int pbc_rvec_sub_gpu(const t_pbc *pbc, const rvec xi, const rvec xj, rvec dx,
                     const rvec *pbc_hbox_diag, const matrix *pbc_box,
                     const rvec *pbc_mhbox_diag, const rvec  *pbc_fbox_diag,
                    const rvec*  pbc_tric_vec,  const ivec*  pbc_tric_shift)
{
    if (pbc) //need to be careful here, gpu pointer needs to be
            // null if cpu pointer is null
    {
        return pbc_dx_aiuc_gpu(pbc, xi, xj, dx, *pbc_hbox_diag, *pbc_box, *pbc_mhbox_diag,
            *pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift );
    }
    else
    { 
        rvec_sub_gpu(xi, xj, dx);
        return CENTRAL;
   }
}

/* Harmonic */

__device__ real harmonic_gpu(real kA, real kB, real xA, real xB, real x, real lambda,
                     real *V, real *F)
{
    const real half = 0.5;
    real       L1, kk, x0, dx, dx2;
    real       v, f, dvdlambda;

    L1    = 1.0-lambda;
    kk    = L1*kA+lambda*kB;
    x0    = L1*xA+lambda*xB;

    dx    = x-x0;
    dx2   = dx*dx;

    f          = -kk*dx;
    v          = half*kk*dx2;
    dvdlambda  = half*(kB-kA)*dx2 + (xA-xB)*kk*dx;

    *F    = f;
    *V    = v;

    return dvdlambda;

    /* That was 19 flops */
}

__global__ 
void bonds_gpu(real *vtot, int nbonds, int natoms,
               const t_iatom forceatoms[], const t_iparams forceparams[],
               const rvec x[], rvec f_distributed[], rvec fshift[], 
               const t_pbc *pbc,  const rvec *pbc_hbox_diag, const matrix *pbc_box, 
               const rvec *pbc_mhbox_diag, const rvec  *pbc_fbox_diag,  
               const rvec*  pbc_tric_vec,  const ivec*  pbc_tric_shift,
               real lambda, real *dvdlambda,
               int force_mapping[] ) 
{

  int i=blockIdx.x*blockDim.x+threadIdx.x;

  int  m, ki, ai, aj, type;
  int li,lj;
  float dr, dr2, fbond=0, vbond=0, fij;
  rvec dx;
  __shared__ real vtot_loc;
  __shared__ rvec fshift_loc[SHIFTS];

  if(threadIdx.x == 0) vtot_loc=0.0;
  if(threadIdx.x<SHIFTS) {
    fshift_loc[threadIdx.x][XX]=0.0;
    fshift_loc[threadIdx.x][YY]=0.0;
    fshift_loc[threadIdx.x][ZZ]=0.0;
  }
   __syncthreads();

  // probably makes sense to do this with more calculations
  // per thread
  if(i< nbonds) {
     type = forceatoms[3*i];
     ai   = forceatoms[3*i+1];
     aj   = forceatoms[3*i+2];
     li   = force_mapping[2*i];
     lj   = force_mapping[2*i+1];

// dx = xi - xj, corrected for periodic boundry conditions.
     ki   = pbc_rvec_sub_gpu(pbc, x[ai], x[aj], dx, 
            pbc_hbox_diag, pbc_box, pbc_mhbox_diag, 
            pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift);
// inner product
     dr2  = iprod_gpu(dx, dx);
     dr   = sqrt(dr2);

     harmonic_gpu(forceparams[type].harmonic.krA,
                              forceparams[type].harmonic.krB,
                              forceparams[type].harmonic.rA,
                              forceparams[type].harmonic.rB,
                              dr, lambda, &vbond, &fbond); 

    if (dr2 == 0.0 ) return; //will be continue if we had a loop

     atomicAdd(&vtot_loc,vbond);

     fbond *= rsqrtf(dr2);

#pragma unroll
    for (m = 0; (m < DIM); m++)  
        {
            fij                 = fbond*dx[m];
            f_distributed[ai+natoms*li][m] = fij;
            f_distributed[aj+natoms*lj][m] = -1.0f*fij;
            if(ki != CENTRAL) {
              atomicAdd(&fshift_loc[ki][m],fij);
              atomicAdd(&fshift_loc[CENTRAL][m],-1.0f*fij);
            }
        }
  }
  __syncthreads();
  if(threadIdx.x == 0)  atomicAdd(vtot,vtot_loc);
  if(threadIdx.x <SHIFTS) rvec_inc_atomic(fshift[threadIdx.x],fshift_loc[threadIdx.x]);
  return;
}

__global__
void bonds_gpu_noener(real *vtot, int nbonds, int natoms,
               const t_iatom forceatoms[], const t_iparams forceparams[],
               const rvec x[], rvec f_distributed[], rvec fshift[],
               const t_pbc *pbc, const rvec *pbc_hbox_diag, const matrix *pbc_box,
               const rvec *pbc_mhbox_diag, const rvec  *pbc_fbox_diag,
               const rvec*  pbc_tric_vec,  const ivec*  pbc_tric_shift,
               real lambda, real *dvdlambda,
               int force_mapping[] )
{

  int i=blockIdx.x*blockDim.x+threadIdx.x;

  int  m, ai, aj, type;
  int li,lj;
  float dr, dr2, fbond=0, vbond=0, fij;
  rvec dx;

  if(i< nbonds) {
     type = forceatoms[3*i];
     ai   = forceatoms[3*i+1];
     aj   = forceatoms[3*i+2];
     li   = force_mapping[2*i];
     lj   = force_mapping[2*i+1];
   
// dx = xi - xj, corrected for periodic boundry conditions.
     pbc_rvec_sub_gpu(pbc, x[ai], x[aj], dx,
            pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
            pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift);
  
// inner product
     dr2  = iprod_gpu(dx, dx);
     dr   = sqrt(dr2);

     harmonic_gpu(forceparams[type].harmonic.krA,
                              forceparams[type].harmonic.krB,
                              forceparams[type].harmonic.rA,
                              forceparams[type].harmonic.rB,
                              dr, lambda, &vbond, &fbond);

    if (dr2 == 0.0 ) return; //will be continue if we had a loop

     fbond *= rsqrtf(dr2);

#pragma unroll
    for (m = 0; (m < DIM); m++)
        {
            fij                 = fbond*dx[m];
            f_distributed[ai+natoms*li][m] = fij;
            f_distributed[aj+natoms*lj][m] = -1.0f*fij;
        }
  }
  return;
}

__device__
real bond_angle_gpu(const rvec xi, const rvec xj, const rvec xk, const t_pbc *pbc,
                    const rvec *pbc_hbox_diag, const matrix *pbc_box,
                    const rvec *pbc_mhbox_diag, const rvec  *pbc_fbox_diag,
                    const rvec*  pbc_tric_vec,  const ivec*  pbc_tric_shift,
                    rvec r_ij, rvec r_kj, real *costh,
                    int *t1, int *t2)
/* Return value is the angle between the bonds i-j and j-k */
{
    /* 41 FLOPS */
    real th;

    *t1 = pbc_rvec_sub_gpu(pbc, xi, xj, r_ij,
                     pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                     pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift); 
    *t2 = pbc_rvec_sub_gpu(pbc, xk, xj, r_kj,
                     pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                     pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift); 

    *costh = cos_angle_gpu(r_ij, r_kj);        /* 25                */
    th     = acosf(*costh);            /* 10                */
    /* 41 TOTAL */
    return th;
}

__global__
void angles_gpu(real *vtot, int nbonds, int natoms,
            const t_iatom forceatoms[], const t_iparams forceparams[],
            const rvec x[], rvec f_distributed[], rvec fshift[],
            const t_pbc *pbc,  const rvec *pbc_hbox_diag, const matrix *pbc_box,
            const rvec *pbc_mhbox_diag, const rvec  *pbc_fbox_diag,
            const rvec*  pbc_tric_vec,  const ivec*  pbc_tric_shift,
            real lambda, real *dvdlambda,
            int force_mapping[]) //,
{
   __shared__ real vtot_loc;
   __shared__ rvec fshift_loc[SHIFTS];
    int  i, ai, aj, ak, t1, t2, type;
    int li,lj,lk;
    rvec r_ij, r_kj;
    real cos_theta, cos_theta2, theta, dVdt, va;

    i=blockIdx.x*blockDim.x+threadIdx.x;
    if(threadIdx.x == 0) vtot_loc=0;
    if(threadIdx.x<SHIFTS) {
      fshift_loc[threadIdx.x][XX]=0.0;
      fshift_loc[threadIdx.x][YY]=0.0;
      fshift_loc[threadIdx.x][ZZ]=0.0;
    }

    __syncthreads();
    if (i < nbonds )
    {
        type = forceatoms[4*i];
        ai   = forceatoms[4*i+1];
        aj   = forceatoms[4*i+2];
        ak   = forceatoms[4*i+3];
        li   = force_mapping[3*i];
        lj   = force_mapping[3*i+1];
        lk   = force_mapping[3*i+2];

        theta  = bond_angle_gpu(x[ai], x[aj], x[ak], pbc,
                       pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                       pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift,
                       r_ij, r_kj, &cos_theta, &t1, &t2);

        harmonic_gpu(forceparams[type].harmonic.krA,
                               forceparams[type].harmonic.krB,
                               forceparams[type].harmonic.rA*DEG2RAD,
                               forceparams[type].harmonic.rB*DEG2RAD,
                               theta, lambda, &va, &dVdt);  
/*  21  */
        atomicAdd(&vtot_loc,va);

        cos_theta2 = cos_theta*cos_theta;
        if (cos_theta2 < 1)
        {
            int  m;
            real st, sth;
            real cik, cii, ckk;
            real nrkj2, nrij2;
            real nrkj_1, nrij_1;
            rvec f_i, f_j, f_k;

            st  = dVdt*rsqrtf(1 - cos_theta2); /*  12             */
            sth = st*cos_theta;                      /*   1             */

            nrij2 = iprod_gpu(r_ij, r_ij);      /*   5              */
            nrkj2 = iprod_gpu(r_kj, r_kj);      /*   5              */

            nrij_1 = rsqrtf(nrij2);   /*  10              */
            nrkj_1 = rsqrtf(nrkj2);   /*  10              */

            cik = st*nrij_1*nrkj_1;         /*   2              */
            cii = sth*nrij_1*nrij_1;        /*   2              */
            ckk = sth*nrkj_1*nrkj_1;        /*   2              */

            for (m = 0; m < DIM; m++)
            {           /*  39          */
                f_i[m]    = -(cik*r_kj[m] - cii*r_ij[m]);
                f_k[m]    = -(cik*r_ij[m] - ckk*r_kj[m]);
                f_j[m]    = -f_i[m] - f_k[m];
                f_distributed[ai+natoms*li][m] = f_i[m];
                f_distributed[aj+natoms*lj][m] = f_j[m];
                f_distributed[ak+natoms*lk][m] = f_k[m];
            }
            rvec_inc_atomic(fshift_loc[t1], f_i);
            rvec_inc_atomic(fshift_loc[CENTRAL], f_j);
            rvec_inc_atomic(fshift_loc[t2], f_k);
        }                                           /* 161 TOTAL        */
    }
    __syncthreads();
    if(threadIdx.x == 0)  atomicAdd(vtot,vtot_loc);
    if(threadIdx.x <SHIFTS)  rvec_inc_atomic(fshift[threadIdx.x],fshift_loc[threadIdx.x]);
    return;
}

__global__
void angles_gpu_noener(real *vtot, int nbonds, int natoms,
            const t_iatom forceatoms[], const t_iparams forceparams[],
            const rvec x[], rvec f_distributed[], rvec fshift[],
            const t_pbc *pbc, const rvec *pbc_hbox_diag, const matrix *pbc_box,
            const rvec *pbc_mhbox_diag, const rvec  *pbc_fbox_diag,
            const rvec*  pbc_tric_vec,  const ivec*  pbc_tric_shift,
            real lambda, real *dvdlambda,
            int force_mapping[]) //,
{
    int  i, ai, aj, ak, t1, t2, type;
    int li,lj,lk;
    rvec r_ij, r_kj;
    real cos_theta, cos_theta2, theta, dVdt, va;

    i=blockIdx.x*blockDim.x+threadIdx.x;

    if (i < nbonds )
    {
        type = forceatoms[4*i];
        ai   = forceatoms[4*i+1];
        aj   = forceatoms[4*i+2];
        ak   = forceatoms[4*i+3];
        li   = force_mapping[3*i];
        lj   = force_mapping[3*i+1];
        lk   = force_mapping[3*i+2];

        theta  = bond_angle_gpu(x[ai], x[aj], x[ak], pbc,
                     pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                     pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift,
                     r_ij, r_kj, &cos_theta, &t1, &t2);

        harmonic_gpu(forceparams[type].harmonic.krA,
                               forceparams[type].harmonic.krB,
                               forceparams[type].harmonic.rA*DEG2RAD,
                               forceparams[type].harmonic.rB*DEG2RAD,
                               theta, lambda, &va, &dVdt);
/*  21  */

        cos_theta2 = cos_theta*cos_theta;
        if (cos_theta2 < 1)
        {
            int  m;
            real st, sth;
            real cik, cii, ckk;
            real nrkj2, nrij2;
            real nrkj_1, nrij_1;
            rvec f_i, f_j, f_k;

            st  = dVdt*rsqrtf(1 - cos_theta2); /*  12             */
            sth = st*cos_theta;                      /*   1             */
            nrij2 = iprod_gpu(r_ij, r_ij);      /*   5              */
            nrkj2 = iprod_gpu(r_kj, r_kj);      /*   5              */

            nrij_1 = rsqrtf(nrij2);   /*  10              */
            nrkj_1 = rsqrtf(nrkj2);   /*  10              */

            cik = st*nrij_1*nrkj_1;         /*   2              */
            cii = sth*nrij_1*nrij_1;        /*   2              */
            ckk = sth*nrkj_1*nrkj_1;        /*   2              */

            for (m = 0; m < DIM; m++)
            {           /*  39          */
                f_i[m]    = -(cik*r_kj[m] - cii*r_ij[m]);
                f_k[m]    = -(cik*r_ij[m] - ckk*r_kj[m]);
                f_j[m]    = -f_i[m] - f_k[m];
                f_distributed[ai+natoms*li][m] = f_i[m];
                f_distributed[aj+natoms*lj][m] = f_j[m];
                f_distributed[ak+natoms*lk][m] = f_k[m];
            }
        }                                           /* 161 TOTAL        */
    }
    return;
}



__global__
void urey_bradley_gpu(real *vtot, int nbonds,int natoms,
                  const t_iatom forceatoms[], const t_iparams forceparams[],
                  const rvec x[], rvec f_distributed[], rvec fshift[],
                  const t_pbc *pbc, const rvec *pbc_hbox_diag, const matrix *pbc_box,
                  const rvec *pbc_mhbox_diag, const rvec  *pbc_fbox_diag,
                  const rvec*  pbc_tric_vec,  const ivec*  pbc_tric_shift,
                  real lambda, real *dvdlambda,
                  int force_mapping[])
{
    int  i, m, ai, aj, ak, t1, t2, type, ki;
    int li,lj,lk;
    rvec r_ij, r_kj, r_ik;
    real cos_theta, cos_theta2, theta;
    real dVdt, va, dr, dr2, vbond, fbond, fik;
    real kthA, th0A, kUBA, r13A, kthB, th0B, kUBB, r13B;
    real dvdlambda_loc;

    __shared__ real vtot_loc;
    __shared__ rvec fshift_loc[SHIFTS];
 
    i=blockIdx.x*blockDim.x+threadIdx.x;
    if(threadIdx.x == 0) vtot_loc=0;
    if( threadIdx.x <SHIFTS) {
      fshift_loc[threadIdx.x][XX]=0.0;
      fshift_loc[threadIdx.x][YY]=0.0;
      fshift_loc[threadIdx.x][ZZ]=0.0;
    }
    __syncthreads();
    if (i < nbonds )
    {
        type  = forceatoms[4*i];
        ai    = forceatoms[4*i+1];
        aj    = forceatoms[4*i+2];
        ak    = forceatoms[4*i+3];
        li   = force_mapping[3*i];
        lj   = force_mapping[3*i+1];
        lk   = force_mapping[3*i+2];
        th0A  = forceparams[type].u_b.thetaA*DEG2RAD;
        kthA  = forceparams[type].u_b.kthetaA;
        r13A  = forceparams[type].u_b.r13A;
        kUBA  = forceparams[type].u_b.kUBA;
        th0B  = forceparams[type].u_b.thetaB*DEG2RAD;
        kthB  = forceparams[type].u_b.kthetaB;
        r13B  = forceparams[type].u_b.r13B;
        kUBB  = forceparams[type].u_b.kUBB;

        theta  = bond_angle_gpu(x[ai], x[aj], x[ak], pbc,
                    pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                    pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift,
                    r_ij, r_kj, &cos_theta, &t1, &t2);                

        dvdlambda_loc  = harmonic_gpu(kthA, kthB, th0A, th0B, theta, lambda, &va, &dVdt); /*  21  */
        atomicAdd(&vtot_loc,va);

        ki   = pbc_rvec_sub_gpu(pbc, x[ai], x[ak], r_ik,
                     pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                     pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift); 

        dr2  = iprod_gpu(r_ik, r_ik);                                                   /*   5              */
        dr   = dr2*rsqrtf(dr2);                                               /*  10              */

        dvdlambda_loc += harmonic_gpu(kUBA, kUBB, r13A, r13B, dr, lambda, &vbond, &fbond); /*  19  */

        cos_theta2 = cos_theta*cos_theta;                                        /*   1              */
        if (cos_theta2 < 1)
        {
            real st, sth;
            real cik, cii, ckk;
            real nrkj2, nrij2;
            rvec f_i, f_j, f_k;

            st  = dVdt*rsqrtf(1 - cos_theta2); /*  12             */
            sth = st*cos_theta;                      /*   1             */

            nrkj2 = iprod_gpu(r_kj, r_kj);  /*   5          */
            nrij2 = iprod_gpu(r_ij, r_ij);

            cik = st*rsqrtf(nrkj2*nrij2); /*  12          */
            cii = sth/nrij2;                    /*  10          */
            ckk = sth/nrkj2;                    /*  10          */

            for (m = 0; (m < DIM); m++)         /*  39          */
            {
                f_i[m]    = -(cik*r_kj[m]-cii*r_ij[m]);
                f_k[m]    = -(cik*r_ij[m]-ckk*r_kj[m]);
                f_j[m]    = -f_i[m]-f_k[m];
                f_distributed[ai+natoms*li][m] = f_i[m];
                f_distributed[aj+natoms*lj][m] = f_j[m];
                f_distributed[ak+natoms*lk][m] = f_k[m];
            }
            rvec_inc_atomic(fshift_loc[t1], f_i);
            rvec_inc_atomic(fshift_loc[CENTRAL], f_j);
            rvec_inc_atomic(fshift_loc[t2], f_k);
        } else {
           for (m = 0; (m < DIM); m++)  
           {
             f_distributed[ai+natoms*li][m] = 0.0;
             f_distributed[aj+natoms*lj][m] = 0.0;
             f_distributed[ak+natoms*lk][m] = 0.0;
           }
        }

                             /* 161 TOTAL    */
        /* Time for the bond calculations */
        if (dr2 != 0.0)
        {
         
        
          atomicAdd(&vtot_loc,vbond); /* 1*/
          fbond *= rsqrtf(dr2); /*   6              */

          for (m = 0; (m < DIM); m++)     /*  15          */
          {
            fik                 = fbond*r_ik[m];
              f_distributed[ai+natoms*li][m] += fik;
              f_distributed[ak+natoms*lk][m] -= fik;
              if(ki != CENTRAL) {
                atomicAdd(&fshift_loc[ki][m],fik);
                atomicAdd(&fshift_loc[CENTRAL][m],-1.0f*fik);
             }
          }
        }
    }
    __syncthreads();
    if(threadIdx.x == 0) atomicAdd(vtot,vtot_loc); 
    if( threadIdx.x <SHIFTS)  rvec_inc_atomic(fshift[threadIdx.x],fshift_loc[threadIdx.x]);
    return ;
}

__global__
void urey_bradley_gpu_noener(real *vtot, int nbonds,int natoms,
                  const t_iatom forceatoms[], const t_iparams forceparams[],
                  const rvec x[], rvec f_distributed[], rvec fshift[],
                  const t_pbc *pbc, const rvec *pbc_hbox_diag, const matrix *pbc_box,
                  const rvec *pbc_mhbox_diag, const rvec  *pbc_fbox_diag,
                  const rvec*  pbc_tric_vec,  const ivec*  pbc_tric_shift,
                  real lambda, real *dvdlambda,
                  int force_mapping[])
{
    int  i, m, ai, aj, ak, t1, t2, type;
    int li,lj,lk;
    rvec r_ij, r_kj, r_ik;
    real cos_theta, cos_theta2, theta;
    real dVdt, va, dr, dr2, vbond, fbond, fik;
    real kthA, th0A, kUBA, r13A, kthB, th0B, kUBB, r13B;
    real dvdlambda_loc;

    i=blockIdx.x*blockDim.x+threadIdx.x;
    if (i < nbonds )
    {
        type  = forceatoms[4*i];
        ai    = forceatoms[4*i+1];
        aj    = forceatoms[4*i+2];
        ak    = forceatoms[4*i+3];
        li   = force_mapping[3*i];
        lj   = force_mapping[3*i+1];
        lk   = force_mapping[3*i+2];
        th0A  = forceparams[type].u_b.thetaA*DEG2RAD;
        kthA  = forceparams[type].u_b.kthetaA;
        r13A  = forceparams[type].u_b.r13A;
        kUBA  = forceparams[type].u_b.kUBA;
        th0B  = forceparams[type].u_b.thetaB*DEG2RAD;
        kthB  = forceparams[type].u_b.kthetaB;
        r13B  = forceparams[type].u_b.r13B;
        kUBB  = forceparams[type].u_b.kUBB;

        theta  = bond_angle_gpu(x[ai], x[aj], x[ak], pbc,
                       pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                       pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift,
                       r_ij, r_kj, &cos_theta, &t1, &t2);

        dvdlambda_loc  = harmonic_gpu(kthA, kthB, th0A, th0B, theta, lambda, &va, &dVdt); /*  21  */

        pbc_rvec_sub_gpu(pbc, x[ai], x[ak], r_ik,
                     pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                     pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift);                               /*   3      */
        dr2  = iprod_gpu(r_ik, r_ik);                                                   /*   5              */
        dr   = dr2*rsqrtf(dr2);                                               /*  10              */

        dvdlambda_loc += harmonic_gpu(kUBA, kUBB, r13A, r13B, dr, lambda, &vbond, &fbond); /*  19  */

        cos_theta2 = cos_theta*cos_theta;                                        /*   1              */
        if (cos_theta2 < 1)
        {
            real st, sth;
            real cik, cii, ckk;
            real nrkj2, nrij2;
            rvec f_i, f_j, f_k;

            st  = dVdt*rsqrtf(1 - cos_theta2); /*  12             */
            sth = st*cos_theta;                      /*   1             */

            nrkj2 = iprod_gpu(r_kj, r_kj);  /*   5          */
            nrij2 = iprod_gpu(r_ij, r_ij);

            cik = st*rsqrtf(nrkj2*nrij2); /*  12          */
            cii = sth/nrij2;                    /*  10          */
            ckk = sth/nrkj2;                    /*  10          */

            for (m = 0; (m < DIM); m++)         /*  39          */
            {
                f_i[m]    = -(cik*r_kj[m]-cii*r_ij[m]);
                f_k[m]    = -(cik*r_ij[m]-ckk*r_kj[m]);
                f_j[m]    = -f_i[m]-f_k[m];
                f_distributed[ai+natoms*li][m] = f_i[m];
                f_distributed[aj+natoms*lj][m] = f_j[m];
                f_distributed[ak+natoms*lk][m] = f_k[m];
            }
        } else {
           for (m = 0; (m < DIM); m++)
           {
             f_distributed[ai+natoms*li][m] = 0.0;
             f_distributed[aj+natoms*lj][m] = 0.0;
             f_distributed[ak+natoms*lk][m] = 0.0;
           }
        }

                             /* 161 TOTAL    */
        /* Time for the bond calculations */
        if (dr2 != 0.0)
        {
          fbond *= rsqrtf(dr2); /*   6              */

          for (m = 0; (m < DIM); m++)     /*  15          */
          {
            fik                 = fbond*r_ik[m];
              f_distributed[ai+natoms*li][m] += fik;
              f_distributed[ak+natoms*lk][m] -= fik;
          }
        }
    }
    return ;
}

__device__
real dih_angle_gpu(const rvec xi, const rvec xj, const rvec xk, const rvec xl,
               const t_pbc *pbc,  const rvec *pbc_hbox_diag, const matrix *pbc_box,
               const rvec *pbc_mhbox_diag, const rvec  *pbc_fbox_diag,
               const rvec*  pbc_tric_vec,  const ivec*  pbc_tric_shift,
               rvec r_ij, rvec r_kj, rvec r_kl, rvec m, rvec n,
               int *t1, int *t2, int *t3)
{
    *t1 = pbc_rvec_sub_gpu(pbc, xi, xj, r_ij,
                     pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                     pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift); /*  3        */
    *t2 = pbc_rvec_sub_gpu(pbc, xk, xj, r_kj,
                     pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                     pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift); /*  3                */
    *t3 = pbc_rvec_sub_gpu(pbc, xk, xl, r_kl,
                     pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                     pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift); /*  3                */

    cprod_gpu(r_ij, r_kj, m);                  /*  9        */
    cprod_gpu(r_kj, r_kl, n);                  /*  9                */
    real phi  = gmx_angle_gpu(m, n);           /* 49 (assuming 25 for atan2) */
    real ipr  = iprod_gpu(r_ij, n);            /*  5        */
    real sign = (ipr < 0.0) ? -1.0 : 1.0;
    phi       = sign*phi;                  /*  1                */
    /* 82 TOTAL */
    return phi;
}


__device__
real dopdihs_gpu(real cpA, real cpB, real phiA, real phiB, int mult,
                    real phi, real lambda, real *V, real *F)
{
    real v, dvdlambda, mdphi, v1, sdphi, ddphi;
    real L1   = 1.0 - lambda;
    real ph0  = (L1*phiA + lambda*phiB)*DEG2RAD;
    real dph0 = (phiB - phiA)*DEG2RAD;
    real cp   = L1*cpA + lambda*cpB;

    mdphi =  mult*phi - ph0;
    sdphi = sinf(mdphi); //requires float
    ddphi = -cp*mult*sdphi;
    v1    = 1.0 + cosf(mdphi); //requires float
    v     = cp*v1;

    dvdlambda  = (cpB - cpA)*v1 + cp*dph0*sdphi;

    *V = v;
    *F = ddphi;

    return dvdlambda;

    /* That was 40 flops */
}

__device__
void do_dih_fup_gpu(int i, int j, int k, int l, 
                    int li, int lj, int lk, int ll,int natoms,
                real ddphi,
                rvec r_ij, rvec r_kj, rvec r_kl,
                rvec m, rvec n, rvec f_distributed[], rvec fshift[],
                const t_pbc *pbc,  const rvec *pbc_hbox_diag, const matrix *pbc_box,
                const rvec *pbc_mhbox_diag, const rvec  *pbc_fbox_diag,
                const rvec*  pbc_tric_vec,  const ivec*  pbc_tric_shift,
                const rvec x[], int t1, int t2, int t3)
{
    /* 143 FLOPS */
    rvec f_i, f_j, f_k, f_l;
    rvec uvec, vvec, svec, dx_jl;
    real iprm, iprn, nrkj, nrkj2, nrkj_1, nrkj_2;
    real a, b, p, q, toler;

    iprm  = iprod_gpu(m, m);       /*  5    */
    iprn  = iprod_gpu(n, n);       /*  5    */
    nrkj2 = iprod_gpu(r_kj, r_kj); /*  5    */
    toler = nrkj2*GMX_REAL_EPS;
    if ((iprm > toler) && (iprn > toler))
    {
        nrkj_1 = 1.0f/sqrt(nrkj2);    /* 10     */  // replacing std::invsqrt call
        nrkj_2 = nrkj_1*nrkj_1;       /*  1     */
        nrkj   = nrkj2*nrkj_1;        /*  1     */
        a      = -ddphi*nrkj/iprm;    /* 11     */
        svmul_gpu(a, m, f_i);             /*  3     */
        b     = ddphi*nrkj/iprn;      /* 11     */
        svmul_gpu(b, n, f_l);             /*  3  */
        p     = iprod_gpu(r_ij, r_kj);    /*  5     */
        p    *= nrkj_2;               /*  1     */
        q     = iprod_gpu(r_kl, r_kj);    /*  5     */
        q    *= nrkj_2;               /*  1     */
        svmul_gpu(p, f_i, uvec);          /*  3     */
        svmul_gpu(q, f_l, vvec);          /*  3     */
        rvec_sub_gpu(uvec, vvec, svec);   /*  3     */
        rvec_sub_gpu(f_i, svec, f_j);     /*  3     */
        rvec_add_gpu(f_l, svec, f_k);     /*  3     */
#pragma unroll
        for (int m = 0; (m < DIM); m++)
        {
           f_distributed[i+natoms*li][m] = f_i[m];
           f_distributed[j+natoms*lj][m] = -1.0*f_j[m];
           f_distributed[k+natoms*lk][m] = -1.0*f_k[m];
           f_distributed[l+natoms*ll][m] = f_l[m];
        }
          
        if (pbc)
        {
            t3 = pbc_rvec_sub_gpu(pbc, x[l], x[j], dx_jl,
                     pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                     pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift);
        }
        else
        {
            t3 = CENTRAL;
        }

        rvec_inc_atomic(fshift[t1], f_i);
        rvec_dec_atomic(fshift[CENTRAL], f_j);
        rvec_dec_atomic(fshift[t2], f_k);
        rvec_inc_atomic(fshift[t3], f_l);
    }
    /* 112 TOTAL    */
}

__device__
void do_dih_fup_gpu_noener(int i, int j, int k, int l,
                    int li, int lj, int lk, int ll,int natoms,
                real ddphi,
                rvec r_ij, rvec r_kj, rvec r_kl,
                rvec m, rvec n, rvec f_distributed[], rvec fshift[],
                const t_pbc *pbc,  const rvec *pbc_hbox_diag, const matrix *pbc_box,
                const rvec *pbc_mhbox_diag, const rvec  *pbc_fbox_diag,
                const rvec*  pbc_tric_vec,  const ivec*  pbc_tric_shift,
                const rvec x[], int t1, int t2, int gmx_unused t3)
{
    /* 143 FLOPS */
    rvec f_i, f_j, f_k, f_l;
    rvec uvec, vvec, svec, dx_jl;
    real iprm, iprn, nrkj, nrkj2, nrkj_1, nrkj_2;
    real a, b, p, q, toler;

    iprm  = iprod_gpu(m, m);       /*  5    */
    iprn  = iprod_gpu(n, n);       /*  5    */
    nrkj2 = iprod_gpu(r_kj, r_kj); /*  5    */
    toler = nrkj2*GMX_REAL_EPS;
    if ((iprm > toler) && (iprn > toler))
    {
        nrkj_1 = 1.0f/sqrt(nrkj2);    /* 10     */  // replacing std::invsqrt call
        nrkj_2 = nrkj_1*nrkj_1;       /*  1     */
        nrkj   = nrkj2*nrkj_1;        /*  1     */
        a      = -ddphi*nrkj/iprm;    /* 11     */
        svmul_gpu(a, m, f_i);             /*  3     */
        b     = ddphi*nrkj/iprn;      /* 11     */
        svmul_gpu(b, n, f_l);             /*  3  */
        p     = iprod_gpu(r_ij, r_kj);    /*  5     */
        p    *= nrkj_2;               /*  1     */
        q     = iprod_gpu(r_kl, r_kj);    /*  5     */
        q    *= nrkj_2;               /*  1     */
        svmul_gpu(p, f_i, uvec);          /*  3     */
        svmul_gpu(q, f_l, vvec);          /*  3     */
        rvec_sub_gpu(uvec, vvec, svec);   /*  3     */
        rvec_sub_gpu(f_i, svec, f_j);     /*  3     */
        rvec_add_gpu(f_l, svec, f_k);     /*  3     */
#pragma unroll
        for (int m = 0; (m < DIM); m++)
        {
           f_distributed[i+natoms*li][m] = f_i[m];
           f_distributed[j+natoms*lj][m] = -1.0*f_j[m];
           f_distributed[k+natoms*lk][m] = -1.0*f_k[m];
           f_distributed[l+natoms*ll][m] = f_l[m];
        }
        if (pbc)
        {
            pbc_rvec_sub_gpu(pbc, x[l], x[j], dx_jl,
                     pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                     pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift);
        }

    }
    /* 112 TOTAL    */
}

__global__
void  pdihs_gpu(real *vtot, int nbonds, int natoms,
           const t_iatom forceatoms[], const t_iparams forceparams[],
           const rvec x[], rvec f[], rvec fshift[],
           const t_pbc *pbc,  const rvec *pbc_hbox_diag, const matrix *pbc_box,
           const rvec *pbc_mhbox_diag, const rvec  *pbc_fbox_diag,
           const rvec*  pbc_tric_vec,  const ivec*  pbc_tric_shift,
           real lambda, real *dvdlambda,
           int force_mapping[] ) 
{
    int  i, type, ai, aj, ak, al;
    int li,lj,lk,ll;
    int  t1, t2, t3;
    rvec r_ij, r_kj, r_kl, m, n;
    real phi, ddphi, vpd;

    i=blockIdx.x*blockDim.x+threadIdx.x;
    __shared__ real vtot_loc;
    __shared__ rvec fshift_loc[SHIFTS];
    if(threadIdx.x == 0) vtot_loc=0;
    if(threadIdx.x < SHIFTS) {
      fshift_loc[threadIdx.x][XX]=0.0;
      fshift_loc[threadIdx.x][YY]=0.0;
      fshift_loc[threadIdx.x][ZZ]=0.0;
    } 
    __syncthreads();
     
    if(i<nbonds)
    {
        type = forceatoms[5*i];
        ai   = forceatoms[5*i+1];
        aj   = forceatoms[5*i+2];
        ak   = forceatoms[5*i+3];
        al   = forceatoms[5*i+4];
        li   = force_mapping[4*i];
        lj   = force_mapping[4*i+1];
        lk   = force_mapping[4*i+2];
        ll   = force_mapping[4*i+3];

        phi = dih_angle_gpu(x[ai], x[aj], x[ak], x[al], pbc, 
                     pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                     pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift,
                     r_ij, r_kj, r_kl, m, n,
                     &t1, &t2, &t3);  /*  84      */
        dopdihs_gpu(forceparams[type].pdihs.cpA,
                              forceparams[type].pdihs.cpB,
                              forceparams[type].pdihs.phiA,
                              forceparams[type].pdihs.phiB,
                              forceparams[type].pdihs.mult,
                              phi, lambda, &vpd, &ddphi);

        atomicAdd(&vtot_loc,vpd);
        do_dih_fup_gpu(ai, aj, ak, al, li, lj,lk,ll,natoms,
                   ddphi, r_ij, r_kj, r_kl, m, n,
                   f, fshift_loc, pbc, 
                   pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                   pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift, 
                   x, t1, t2, t3); /* 112            */

    } /* 223 TOTAL  */
    __syncthreads();
    if(threadIdx.x == 0) atomicAdd(vtot,vtot_loc);
    if(threadIdx.x < SHIFTS) rvec_inc_atomic(fshift[threadIdx.x],fshift_loc[threadIdx.x]);
    return;
}

__global__
void  pdihs_gpu_noener(real *vtot, int nbonds, int natoms,
           const t_iatom forceatoms[], const t_iparams forceparams[],
           const rvec x[], rvec f[], rvec fshift[],
           const t_pbc *pbc, const rvec *pbc_hbox_diag, const matrix *pbc_box,
           const rvec *pbc_mhbox_diag, const rvec  *pbc_fbox_diag,
           const rvec*  pbc_tric_vec,  const ivec*  pbc_tric_shift,
           real lambda, real *dvdlambda,
           int force_mapping[] ) 
{
    int  i, type, ai, aj, ak, al;
    int li,lj,lk,ll;
    int  t1, t2, t3;
    rvec r_ij, r_kj, r_kl, m, n;
    real phi, ddphi, vpd;

    i=blockIdx.x*blockDim.x+threadIdx.x;

    if(i<nbonds)
    {
        type = forceatoms[5*i];
        ai   = forceatoms[5*i+1];
        aj   = forceatoms[5*i+2];
        ak   = forceatoms[5*i+3];
        al   = forceatoms[5*i+4];
        li   = force_mapping[4*i];
        lj   = force_mapping[4*i+1];
        lk   = force_mapping[4*i+2];
        ll   = force_mapping[4*i+3];

        phi = dih_angle_gpu(x[ai], x[aj], x[ak], x[al], pbc, 
                     pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                     pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift,
                     r_ij, r_kj, r_kl, m, n,
                     &t1, &t2, &t3);  /*  84      */
        dopdihs_gpu(forceparams[type].pdihs.cpA,
                              forceparams[type].pdihs.cpB,
                              forceparams[type].pdihs.phiA,
                              forceparams[type].pdihs.phiB,
                              forceparams[type].pdihs.mult,
                              phi, lambda, &vpd, &ddphi);

        do_dih_fup_gpu_noener(ai, aj, ak, al, li, lj,lk,ll,natoms,
                   ddphi, r_ij, r_kj, r_kl, m, n,
                   f, fshift, pbc, 
                   pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                   pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift,
                   x, t1, t2, t3); /* 112            */

    } /* 223 TOTAL  */
    return;
}

__global__
void rbdihs_gpu(real *vtot, int nbonds, int natoms,
            const t_iatom forceatoms[], const t_iparams forceparams[],
            const rvec x[], rvec f[], rvec fshift[],
            const t_pbc *pbc,  const rvec *pbc_hbox_diag, const matrix *pbc_box,
            const rvec *pbc_mhbox_diag, const rvec  *pbc_fbox_diag,
            const rvec*  pbc_tric_vec,  const ivec*  pbc_tric_shift,
            real lambda, real *dvdlambda,
            int force_mapping[])
{
    const real c0 = 0.0, c1 = 1.0, c2 = 2.0, c3 = 3.0, c4 = 4.0, c5 = 5.0;
    int        type, ai, aj, ak, al, i, j;
    int li,lj,lk,ll;
    int        t1, t2, t3;
    rvec       r_ij, r_kj, r_kl, m, n;
    real       parmA[NR_RBDIHS];
    real       parmB[NR_RBDIHS];
    real       parm[NR_RBDIHS];
    real       cos_phi, phi, rbp, rbpBA;
    real       v, ddphi, sin_phi;
    real       cosfac;
    real       L1        = 1.0-lambda;
    real       dvdl_term = 0;

    __shared__ real vtot_loc;
    __shared__ rvec fshift_loc[SHIFTS];

    i=blockIdx.x*blockDim.x+threadIdx.x;
    
     if(threadIdx.x == 0) vtot_loc=0;
     if(threadIdx.x < SHIFTS) {
       fshift_loc[threadIdx.x][XX]=0.0;
       fshift_loc[threadIdx.x][YY]=0.0;
       fshift_loc[threadIdx.x][ZZ]=0.0;
     }
    __syncthreads();

    if (i < nbonds )
    {
        type = forceatoms[5*i];
        ai   = forceatoms[5*i+1];
        aj   = forceatoms[5*i+2];
        ak   = forceatoms[5*i+3];
        al   = forceatoms[5*i+4];
        li   = force_mapping[4*i];
        lj   = force_mapping[4*i+1];
        lk   = force_mapping[4*i+2];
        ll   = force_mapping[4*i+3];

        phi = dih_angle_gpu(x[ai], x[aj], x[ak], x[al], pbc, 
                     pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                     pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift,
                     r_ij, r_kj, r_kl, m, n,
                     &t1, &t2, &t3);  /*  84         */

        /* Change to polymer convention */
        if (phi < c0)
        {
            phi += M_PI;
        }
        else
        {
            phi -= M_PI;    /*   1              */

        }
        cos_phi = cosf(phi);
        /* Beware of accuracy loss, cannot use 1-sqrt(cos^2) ! */
        sin_phi = sinf(phi);

        for (j = 0; (j < NR_RBDIHS); j++)
        {
            parmA[j] = forceparams[type].rbdihs.rbcA[j];
            parmB[j] = forceparams[type].rbdihs.rbcB[j];
            parm[j]  = L1*parmA[j]+lambda*parmB[j];
        }
        /* Calculate cosine powers */
        /* Calculate the energy */
        /* Calculate the derivative */

        v            = parm[0];
        dvdl_term   += (parmB[0]-parmA[0]);
        ddphi        = c0;
        cosfac       = c1;

        rbp          = parm[1];
        rbpBA        = parmB[1]-parmA[1];
        ddphi       += rbp*cosfac;
        cosfac      *= cos_phi;
        v           += cosfac*rbp;
        dvdl_term   += cosfac*rbpBA;
        rbp          = parm[2];
        rbpBA        = parmB[2]-parmA[2];
        ddphi       += c2*rbp*cosfac;
        cosfac      *= cos_phi;
        v           += cosfac*rbp;
        dvdl_term   += cosfac*rbpBA;
        rbp          = parm[3];
        rbpBA        = parmB[3]-parmA[3];
        ddphi       += c3*rbp*cosfac;
        cosfac      *= cos_phi;
        v           += cosfac*rbp;
        dvdl_term   += cosfac*rbpBA;
        rbp          = parm[4];
        rbpBA        = parmB[4]-parmA[4];
        ddphi       += c4*rbp*cosfac;
        cosfac      *= cos_phi;
        v           += cosfac*rbp;
        dvdl_term   += cosfac*rbpBA;
        rbp          = parm[5];
        rbpBA        = parmB[5]-parmA[5];
        ddphi       += c5*rbp*cosfac;
        cosfac      *= cos_phi;
        v           += cosfac*rbp;
        dvdl_term   += cosfac*rbpBA;

        ddphi = -ddphi*sin_phi;         /*  11          */

        do_dih_fup_gpu(ai, aj, ak, al, 
                   li, lj ,lk ,ll, natoms,
                   ddphi, r_ij, r_kj, r_kl, m, n,
                   f, fshift_loc, pbc, 
                   pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                   pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift, 
                   x, t1, t2, t3); /* 112            */
         atomicAdd(&vtot_loc,v);
    }
    __syncthreads();
    if(threadIdx.x == 0) atomicAdd(vtot,vtot_loc);
    if(threadIdx.x < SHIFTS) rvec_inc_atomic(fshift[threadIdx.x],fshift_loc[threadIdx.x]);
    return;
}

__global__
void rbdihs_gpu_noener(real *vtot, int nbonds, int natoms,
            const t_iatom forceatoms[], const t_iparams forceparams[],
            const rvec x[], rvec f[], rvec fshift[],
            const t_pbc *pbc,  const rvec *pbc_hbox_diag, const matrix *pbc_box,
            const rvec *pbc_mhbox_diag, const rvec  *pbc_fbox_diag,
            const rvec*  pbc_tric_vec,  const ivec*  pbc_tric_shift,
            real lambda, real *dvdlambda,
            int force_mapping[])
{
    const real c0 = 0.0, c1 = 1.0, c2 = 2.0, c3 = 3.0, c4 = 4.0, c5 = 5.0;
    int        type, ai, aj, ak, al, i, j;
    int li,lj,lk,ll;
    int        t1, t2, t3;
    rvec       r_ij, r_kj, r_kl, m, n;
    real       parmA[NR_RBDIHS];
    real       parmB[NR_RBDIHS];
    real       parm[NR_RBDIHS];
    real       cos_phi, phi, rbp, rbpBA;
    real       v, ddphi, sin_phi;
    real       cosfac;
    real       L1        = 1.0-lambda;
    real       dvdl_term = 0;

    i=blockIdx.x*blockDim.x+threadIdx.x;

    if (i < nbonds )
    {
        type = forceatoms[5*i];
        ai   = forceatoms[5*i+1];
        aj   = forceatoms[5*i+2];
        ak   = forceatoms[5*i+3];
        al   = forceatoms[5*i+4];
        li   = force_mapping[4*i];
        lj   = force_mapping[4*i+1];
        lk   = force_mapping[4*i+2];
        ll   = force_mapping[4*i+3];

        phi = dih_angle_gpu(x[ai], x[aj], x[ak], x[al], pbc, 
                     pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                     pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift,
                     r_ij, r_kj, r_kl, m, n,
                     &t1, &t2, &t3);  /*  84         */

        /* Change to polymer convention */
        if (phi < c0)
        {
            phi += M_PI;
        }
        else
        {
            phi -= M_PI;    /*   1              */

        }
        cos_phi = cosf(phi);
        /* Beware of accuracy loss, cannot use 1-sqrt(cos^2) ! */
        sin_phi = sinf(phi);

        for (j = 0; (j < NR_RBDIHS); j++)
        {
            parmA[j] = forceparams[type].rbdihs.rbcA[j];
            parmB[j] = forceparams[type].rbdihs.rbcB[j];
            parm[j]  = L1*parmA[j]+lambda*parmB[j];
        }
        /* Calculate cosine powers */
        /* Calculate the energy */
        /* Calculate the derivative */

        v            = parm[0];
        dvdl_term   += (parmB[0]-parmA[0]);
        ddphi        = c0;
        cosfac       = c1;

        rbp          = parm[1];
        rbpBA        = parmB[1]-parmA[1];
        ddphi       += rbp*cosfac;
        cosfac      *= cos_phi;
        v           += cosfac*rbp;
        dvdl_term   += cosfac*rbpBA;
        rbp          = parm[2];
        rbpBA        = parmB[2]-parmA[2];
        ddphi       += c2*rbp*cosfac;
        cosfac      *= cos_phi;
        v           += cosfac*rbp;
        dvdl_term   += cosfac*rbpBA;
        rbp          = parm[3];
        rbpBA        = parmB[3]-parmA[3];
        ddphi       += c3*rbp*cosfac;
        cosfac      *= cos_phi;
        v           += cosfac*rbp;
        dvdl_term   += cosfac*rbpBA;
        rbp          = parm[4];
        rbpBA        = parmB[4]-parmA[4];
        ddphi       += c4*rbp*cosfac;
        cosfac      *= cos_phi;
        v           += cosfac*rbp;
        dvdl_term   += cosfac*rbpBA;
        rbp          = parm[5];
        rbpBA        = parmB[5]-parmA[5];
        ddphi       += c5*rbp*cosfac;
        cosfac      *= cos_phi;
        v           += cosfac*rbp;
        dvdl_term   += cosfac*rbpBA;

        ddphi = -ddphi*sin_phi;         /*  11          */

        do_dih_fup_gpu_noener(ai, aj, ak, al,
                   li, lj ,lk ,ll, natoms,
                   ddphi, r_ij, r_kj, r_kl, m, n,
                   f, fshift, pbc,
                   pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                   pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift,
                   x, t1, t2, t3); /* 112            */
    }
    return;
}


__device__
void make_dp_periodic_gpu(real *dp)  /* 1 flop? */
{
    /* dp cannot be outside (-pi,pi) */
    if (*dp >= M_PI)
    {
        *dp -= 2*M_PI;
    }
    else if (*dp < -M_PI)
    {
        *dp += 2*M_PI;
    }
    return;
}

__global__
void  idihs_gpu(real *vtot, int nbonds, int natoms,
           const t_iatom forceatoms[], const t_iparams forceparams[],
           const rvec x[], rvec f[], rvec fshift[],
           const t_pbc *pbc,  const rvec *pbc_hbox_diag, const matrix *pbc_box,
           const rvec *pbc_mhbox_diag, const rvec  *pbc_fbox_diag,
           const rvec*  pbc_tric_vec,  const ivec*  pbc_tric_shift,
           real lambda, real *dvdlambda,
           int force_mapping[] )
{
    int  i, type, ai, aj, ak, al;
    int li,lj,lk,ll;
    int  t1, t2, t3;
    real phi, phi0, dphi0, ddphi;
    rvec r_ij, r_kj, r_kl, m, n;
    real L1, kk, dp, dp2, kA, kB, pA, pB, dvdl_term;

    L1        = 1.0-lambda;
    dvdl_term = 0;
    i=blockIdx.x*blockDim.x+threadIdx.x;
    __shared__ real vtot_loc;
    __shared__ rvec fshift_loc[SHIFTS];
    if(threadIdx.x == 0) vtot_loc=0;
    if(threadIdx.x<SHIFTS) {
      fshift_loc[threadIdx.x][XX]=0.0;
      fshift_loc[threadIdx.x][YY]=0.0;
      fshift_loc[threadIdx.x][ZZ]=0.0;
    }
    __syncthreads();

    if ( i < nbonds )
    {
        type = forceatoms[5*i];
        ai   = forceatoms[5*i+1];
        aj   = forceatoms[5*i+2];
        ak   = forceatoms[5*i+3];
        al   = forceatoms[5*i+4];
        li   = force_mapping[4*i];
        lj   = force_mapping[4*i+1];
        lk   = force_mapping[4*i+2];
        ll   = force_mapping[4*i+3];
   
        phi = dih_angle_gpu(x[ai], x[aj], x[ak], x[al], pbc, 
                     pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                     pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift,
                     r_ij, r_kj, r_kl, m, n,
                     &t1, &t2, &t3);  /*  84         */

        /* phi can jump if phi0 is close to Pi/-Pi, which will cause huge
         * force changes if we just apply a normal harmonic.
         * Instead, we first calculate phi-phi0 and take it modulo (-Pi,Pi).
         * This means we will never have the periodicity problem, unless
         * the dihedral is Pi away from phiO, which is very unlikely due to
         * the potential.
         */
        kA = forceparams[type].harmonic.krA;
        kB = forceparams[type].harmonic.krB;
        pA = forceparams[type].harmonic.rA;
        pB = forceparams[type].harmonic.rB;

        kk    = L1*kA + lambda*kB;
        phi0  = (L1*pA + lambda*pB)*DEG2RAD;
        dphi0 = (pB - pA)*DEG2RAD;

        dp = phi-phi0;

        make_dp_periodic_gpu(&dp);

        dp2 = dp*dp;

        atomicAdd(&vtot_loc,0.5*kk*dp2);
        ddphi = -kk*dp;

        dvdl_term += 0.5*(kB - kA)*dp2 - kk*dphi0*dp;

        do_dih_fup_gpu(ai, aj, ak, al, 
                   li,lj,lk,ll,natoms,
                   -ddphi, r_ij, r_kj, r_kl, m, n,
                   f, fshift_loc, pbc, 
                   pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                   pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift, 
                   x, t1, t2, t3); /* 112            */
        /* 218 TOTAL    */
    }
    __syncthreads();
    if(threadIdx.x == 0) atomicAdd(vtot,vtot_loc);
    if(threadIdx.x<SHIFTS) rvec_inc_atomic(fshift[threadIdx.x],fshift_loc[threadIdx.x]);
    return;
}

__global__
void  idihs_gpu_noener(real *vtot, int nbonds, int natoms,
           const t_iatom forceatoms[], const t_iparams forceparams[],
           const rvec x[], rvec f[], rvec fshift[],
           const t_pbc *pbc,  const rvec *pbc_hbox_diag, const matrix *pbc_box,
           const rvec *pbc_mhbox_diag, const rvec  *pbc_fbox_diag,
           const rvec*  pbc_tric_vec,  const ivec*  pbc_tric_shift,
           real lambda, real *dvdlambda,
           int force_mapping[] )
{
    int  i, type, ai, aj, ak, al;
    int li,lj,lk,ll;
    int  t1, t2, t3;
    real phi, phi0, dphi0, ddphi;
    rvec r_ij, r_kj, r_kl, m, n;
    real L1, kk, dp, dp2, kA, kB, pA, pB, dvdl_term;

    L1        = 1.0-lambda;
    dvdl_term = 0;
    i=blockIdx.x*blockDim.x+threadIdx.x;
    if ( i < nbonds )
    {
        type = forceatoms[5*i];
        ai   = forceatoms[5*i+1];
        aj   = forceatoms[5*i+2];
        ak   = forceatoms[5*i+3];
        al   = forceatoms[5*i+4];
        li   = force_mapping[4*i];
        lj   = force_mapping[4*i+1];
        lk   = force_mapping[4*i+2];
        ll   = force_mapping[4*i+3];

        phi = dih_angle_gpu(x[ai], x[aj], x[ak], x[al], pbc, 
                     pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                     pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift,
                     r_ij, r_kj, r_kl, m, n,
                     &t1, &t2, &t3);  /*  84         */

        /* phi can jump if phi0 is close to Pi/-Pi, which will cause huge
         * force changes if we just apply a normal harmonic.
         * Instead, we first calculate phi-phi0 and take it modulo (-Pi,Pi).
         * This means we will never have the periodicity problem, unless
         * the dihedral is Pi away from phiO, which is very unlikely due to
         * the potential.
         */
        kA = forceparams[type].harmonic.krA;
        kB = forceparams[type].harmonic.krB;
        pA = forceparams[type].harmonic.rA;
        pB = forceparams[type].harmonic.rB;

        kk    = L1*kA + lambda*kB;
        phi0  = (L1*pA + lambda*pB)*DEG2RAD;
        dphi0 = (pB - pA)*DEG2RAD;

        dp = phi-phi0;

        make_dp_periodic_gpu(&dp);
        dp2 = dp*dp;

        ddphi = -kk*dp;

        dvdl_term += 0.5*(kB - kA)*dp2 - kk*dphi0*dp;

        do_dih_fup_gpu_noener(ai, aj, ak, al,
                   li,lj,lk,ll,natoms,
                   -ddphi, r_ij, r_kj, r_kl, m, n,
                   f, fshift, pbc,
                   pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                   pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift,
                   x, t1, t2, t3); /* 112            */
        /* 218 TOTAL    */
    }
    return;
}

__device__ real
evaluate_single(real r2, real tabscale, real *vftab, real tableStride,
                real qq, real c6, real c12, real *velec, real *vvdw)
{
    real       rinv, r, rtab, eps, eps2, Y, F, Geps, Heps2, Fp, VVe, FFe, VVd, FFd, VVr, FFr, fscal;
    int        ntab;

    /* Do the tabulated interactions - first table lookup */
    rinv             = rsqrtf(r2); // float only 
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
    FFe              = Fp+Geps+2.0*Heps2;
    /* Dispersion */
    Y                = vftab[ntab+4];
    F                = vftab[ntab+5];
    Geps             = eps*vftab[ntab+6];
    Heps2            = eps2*vftab[ntab+7];
    Fp               = F+Geps+Heps2;
    VVd              = Y+eps*Fp;
    FFd              = Fp+Geps+2.0*Heps2;
    /* Repulsion */
    Y                = vftab[ntab+8];
    F                = vftab[ntab+9];
    Geps             = eps*vftab[ntab+10];
    Heps2            = eps2*vftab[ntab+11];
    Fp               = F+Geps+Heps2;
    VVr              = Y+eps*Fp;
    FFr              = Fp+Geps+2.0*Heps2;

    *velec           = qq*VVe;
    *vvdw            = c6*VVd+c12*VVr;
    fscal            = -(qq*FFe+c6*FFd+c12*FFr)*tabscale*rinv;

    return fscal;
}



// Just now for F_LJ14 only
__global__
void pairs_gpu(int ftype, int nbonds,int natoms,
                 const t_iatom iatoms[], const t_iparams iparams[],
                 const rvec x[], rvec f[], rvec fshift[],
                 const t_pbc *pbc,  const rvec *pbc_hbox_diag, const matrix *pbc_box,
                 const rvec *pbc_mhbox_diag, const rvec  *pbc_fbox_diag,
                 const rvec*  pbc_tric_vec,  const ivec*  pbc_tric_shift,
                 real md_chargeA[], unsigned short md_cENER[],
                 real energygrp_elec[], real energygrp_vdw[],
                 real fr_pairsTable_data[],
                 real packed_R[],int packed_I[] ,int force_mapping[]
                 )
{
  __shared__ rvec fshift_loc[SHIFTS];
  real             qq, c6, c12;
  rvec             dx;
  int i,itype,ai,aj,gid;
  int li,lj;
  int              fshift_index;
  real             r2;
  real             fscal, velec, vvdw;

  const real epsfac = packed_R[fr_icEPSFAC];

  i=blockIdx.x*blockDim.x+threadIdx.x;  
  if(threadIdx.x<SHIFTS) {
      fshift_loc[threadIdx.x][XX]=0.0;
      fshift_loc[threadIdx.x][YY]=0.0;
      fshift_loc[threadIdx.x][ZZ]=0.0;
  }

    __syncthreads();
  
  if (i <  nbonds )
  {
        itype = iatoms[3*i];
        ai    = iatoms[3*i+1];
        aj    = iatoms[3*i+2];
        li   = force_mapping[2*i];
        lj   = force_mapping[2*i+1];
        gid   = GID(md_cENER[ai], md_cENER[aj], packed_I[md_NENERGRP]);

        /* Get parameters */
            // F_LJ14 Only
                qq               = md_chargeA[ai]*md_chargeA[aj]*epsfac*packed_R[frFUDGEQQ];
                c6               = iparams[itype].lj14.c6A;
                c12              = iparams[itype].lj14.c12A;

        /* To save flops in the optimized kernels, c6/c12 have 6.0/12.0 derivative prefactors
         * included in the general nfbp array now. This means the tables are scaled down by the
         * same factor, so when we use the original c6/c12 parameters from iparams[] they must
         * be scaled up.
         */
        c6  *= 6.0;
        c12 *= 12.0;

        /* Do we need to apply full periodic boundary conditions? */
         fshift_index =  pbc_rvec_sub_gpu(pbc, x[ai], x[aj], dx,
                     pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                     pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift);

        r2           = norm2_gpu(dx);
// TODO : warn if r2 >= fr->pairsTable->r*fr->pairsTable->r

          // Evaluate tabulated interaction without free energy //
            fscal            = evaluate_single(r2, packed_R[fr_PAIRS_SCALE], 
                               fr_pairsTable_data, packed_I[fr_PAIRS_STRIDE],
                                               qq, c6, c12, &velec, &vvdw);
        atomicAdd(&energygrp_elec[gid],velec);
        atomicAdd(&energygrp_vdw[gid],vvdw); 
     
        svmul_gpu(fscal, dx, dx);

        /* Add the forces */
          for (int m = 0; m < DIM; m++)
          {
             f[ai+natoms*li][m] = dx[m];
             f[aj+natoms*lj][m] = -1.0*dx[m];
          }

        if (fshift_index != CENTRAL)
        {
            rvec_inc_atomic(fshift_loc[fshift_index], dx);
            rvec_dec_atomic(fshift_loc[CENTRAL], dx);
        } 
    }
    __syncthreads();
    if(threadIdx.x<SHIFTS) rvec_inc_atomic(fshift[threadIdx.x],fshift_loc[threadIdx.x]);
    return;
}

__global__
void pairs_gpu_noener(int ftype, int nbonds,int natoms,
                 const t_iatom iatoms[], const t_iparams iparams[],
                 const rvec x[], rvec f[], rvec fshift[],
                 const t_pbc *pbc,  const rvec *pbc_hbox_diag, const matrix *pbc_box,
                 const rvec *pbc_mhbox_diag, const rvec  *pbc_fbox_diag,
                 const rvec*  pbc_tric_vec,  const ivec*  pbc_tric_shift,
                 real md_chargeA[], unsigned short md_cENER[],
                 real energygrp_elec[], real energygrp_vdw[],
                 real fr_pairsTable_data[],
                 real packed_R[],int packed_I[] ,int force_mapping[]
                 )
{
  real             qq, c6, c12;
  rvec             dx;
  int              i,itype,ai,aj;
  int              li,lj;
  real             r2;
  real             fscal, velec, vvdw;

  const real epsfac = packed_R[fr_icEPSFAC];

  i=blockIdx.x*blockDim.x+threadIdx.x;

  if (i <  nbonds )
  {
        itype = iatoms[3*i];
        ai    = iatoms[3*i+1];
        aj    = iatoms[3*i+2];
        li   = force_mapping[2*i];
        lj   = force_mapping[2*i+1];

            // F_LJ14 ONLY
                qq               = md_chargeA[ai]*md_chargeA[aj]*epsfac*packed_R[frFUDGEQQ];
                c6               = iparams[itype].lj14.c6A;
                c12              = iparams[itype].lj14.c12A;

        /* To save flops in the optimized kernels, c6/c12 have 6.0/12.0 derivative prefactors
         * included in the general nfbp array now. This means the tables are scaled down by the
         * same factor, so when we use the original c6/c12 parameters from iparams[] they must
         * be scaled up.
         */
        c6  *= 6.0;
        c12 *= 12.0;

        /* Do we need to apply full periodic boundary conditions? */
         pbc_rvec_sub_gpu(pbc, x[ai], x[aj], dx,
                     pbc_hbox_diag, pbc_box, pbc_mhbox_diag,
                     pbc_fbox_diag,  pbc_tric_vec, pbc_tric_shift);

        r2           = norm2_gpu(dx);
// TODO : warn if r2 >= fr->pairsTable->r*fr->pairsTable->r

          // Evaluate tabulated interaction without free energy //
        fscal            = evaluate_single(r2, packed_R[fr_PAIRS_SCALE],
                               fr_pairsTable_data, packed_I[fr_PAIRS_STRIDE],
                                               qq, c6, c12, &velec, &vvdw);
        svmul_gpu(fscal, dx, dx);

        /* Add the forces */
        for (int m = 0; m < DIM; m++)
        {
           f[ai+natoms*li][m] = dx[m];
           f[aj+natoms*lj][m] = -1.0*dx[m];
        }

    }
    return;
}

// ----- //
__global__ 
void reset_next(int natoms, int force_next[] )
{
  int  i=blockIdx.x*blockDim.x+threadIdx.x;
  if(i < natoms)
  {
    force_next[i]=0;
  }
}

__global__
void calc_force_mapping ( int nbonds, int atoms_per_bond, 
                         const t_iatom iatoms[], int force_mapping[],
                         int force_next[], int *force_max )
{
  int i,j,a,p;
  i=blockIdx.x*blockDim.x+threadIdx.x;
  if (i < nbonds ) 
  {
    for(j=0;j<atoms_per_bond;j++)
    {
      a = iatoms[(atoms_per_bond+1)*i+j+1];
      p = atomicAdd(&force_next[a],1);
      force_mapping[atoms_per_bond*i+j]=p;
      atomicMax(force_max,p);
    } 
    assert(*force_max<BO_MAX_FORCE);
  } 
}

__global__ 
void consolidate_forces (int natoms, rvec force[], rvec force_distributed[], int force_next[])
{
  int i,j;
  rvec force_loc = {0.0f, 0.0f, 0.0f };

  i=blockIdx.x*blockDim.x+threadIdx.x;
  if (i <  natoms )
  {
     for(j=0;j<force_next[i];j++) // could do the max as well, should compare.
     {
       force_loc[XX]+=force_distributed[j*natoms+i][XX];
       force_loc[YY]+=force_distributed[j*natoms+i][YY];
       force_loc[ZZ]+=force_distributed[j*natoms+i][ZZ];
     }
  rvec_inc_gpu(force[i],force_loc);
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
static rvec* x_bonded_d[BO_MAX_RANKS];
static rvec* f_bonded_d[BO_MAX_RANKS]; 
#endif

static rvec* f_bonded[BO_MAX_RANKS];
static rvec* f_shift_d[BO_MAX_RANKS]; 
static rvec* f_shift[BO_MAX_RANKS];

static real vtot[BO_MAX_RANKS][F_NRE];

static real *vtot_d[BO_MAX_RANKS]; 
static t_iparams *forceparams_d[BO_MAX_RANKS];
static real *dvdlambda_d[BO_MAX_RANKS];

// TODO : we have duplicate data for constraints
static t_pbc *pbc_d[BO_MAX_RANKS];
static rvec* pbc_hbox_diag_d[BO_MAX_RANKS];
static matrix* pbc_box_d[BO_MAX_RANKS];
static rvec* pbc_mhbox_diag_d[BO_MAX_RANKS];
static rvec*  pbc_fbox_diag_d[BO_MAX_RANKS];
static rvec*  pbc_tric_vec_d[BO_MAX_RANKS];
static ivec*  pbc_tric_shift_d[BO_MAX_RANKS];

static t_iatom *iatoms_d[BO_MAX_RANKS][F_NRE];
static int device_iatom_alloc[BO_MAX_RANKS][F_NRE]; 

int xbonded_size_max[BO_MAX_RANKS]; 
static bool bonded_init_done[BO_MAX_RANKS]; 
static int ntypes_alloc[BO_MAX_RANKS];

static real *pairs_Rpacked[BO_MAX_RANKS]; // packed array of relavant values
static real *pairs_Rpacked_d[BO_MAX_RANKS];
static int  *pairs_Ipacked[BO_MAX_RANKS];
static int  *pairs_Ipacked_d[BO_MAX_RANKS]; 

static real *md_chargeA_d[BO_MAX_RANKS];

static real *energygrp_elec_d[BO_MAX_RANKS];
static real *energygrp_elec[BO_MAX_RANKS];
static real *energygrp_vdw_d[BO_MAX_RANKS];
static real *energygrp_vdw[BO_MAX_RANKS];

static real *fr_pairsTable_data_d[BO_MAX_RANKS];

static unsigned short *md_cENER_d[BO_MAX_RANKS];

static rvec *force_distributed_d[BO_MAX_RANKS];
static int force_distributed_allocated[BO_MAX_RANKS];
static int *force_next_d[BO_MAX_RANKS][F_NRE];
static int *force_mapping_d[BO_MAX_RANKS][F_NRE];
static int force_mapping_allocated[BO_MAX_RANKS][F_NRE];
static int *force_max_d[BO_MAX_RANKS];
static int force_max[BO_MAX_RANKS];

static cudaStream_t streamBonded[BO_MAX_RANKS];

static bool copy_required_this_step[BO_MAX_RANKS];

/*------------------------------------------------------------------------------*/
// bonded forces //

// initial setup (called once)
void
init_gpu_bonded( t_pbc *pbc,int rank, gmx_grppairener_t *grppener, const t_forcerec *fr)
{  
  cudaError_t stat;

// device_iatom_alloc 
// initiaized to 0 as static

// pbc
   pbc_d[rank] = NULL;
// null status is used as a marker for no pbc  
   if(pbc) {
     stat = cudaMalloc(&pbc_d[rank],sizeof(t_pbc));
     CU_RET_ERR(stat, "cudaMalloc failed");
     stat = cudaMalloc(&pbc_hbox_diag_d[rank],sizeof(rvec));
     CU_RET_ERR(stat, "cudaMalloc failed");
     stat = cudaMalloc(&pbc_box_d[rank],sizeof(matrix));
     CU_RET_ERR(stat, "cudaMalloc failed");
     stat = cudaMalloc(&pbc_mhbox_diag_d[rank],sizeof(rvec));
     CU_RET_ERR(stat, "cudaMalloc failed");
     stat = cudaMalloc(&pbc_fbox_diag_d[rank],sizeof(rvec));
     CU_RET_ERR(stat, "cudaMalloc failed");
     stat = cudaMalloc(&pbc_tric_vec_d[rank],MAX_NTRICVEC*sizeof(rvec));
     CU_RET_ERR(stat, "cudaMalloc failed");
     stat = cudaMalloc(&pbc_tric_shift_d[rank],MAX_NTRICVEC*sizeof(ivec));
     CU_RET_ERR(stat, "cudaMalloc failed");
   }
// malloc dvdlambda
   stat = cudaMalloc(&dvdlambda_d[rank],efptNR*sizeof(real));
   CU_RET_ERR(stat, "cudaMalloc failed");
// vtot
   stat = cudaMalloc(&vtot_d[rank],F_NRE*sizeof(real));
   CU_RET_ERR(stat, "cudaMalloc failed");

// ener arrays
  stat = cudaMalloc(&pairs_Rpacked_d[rank],PAIRS_NR*sizeof(real));
   CU_RET_ERR(stat, "cudaMalloc failed");

// packed arrays
   stat = cudaMalloc(&pairs_Rpacked_d[rank],PAIRS_NR*sizeof(real));
   CU_RET_ERR(stat, "cudaMalloc failed");   
   pairs_Rpacked[rank] = (real *) malloc(sizeof(real)*PAIRS_NR);

   stat = cudaMalloc(&pairs_Ipacked_d[rank],PAIRS_NI*sizeof(int));
   CU_RET_ERR(stat, "cudaMalloc failed");
   pairs_Ipacked[rank] = (int *) malloc(sizeof(int)*PAIRS_NI);  
 
   stat = cudaMalloc(&energygrp_elec_d[rank],grppener->nener*sizeof(real));
   CU_RET_ERR(stat, "cudaMalloc failed");

   stat = cudaHostAlloc(&energygrp_elec[rank],grppener->nener*sizeof(real),cudaHostAllocPortable);
   CU_RET_ERR(stat, "cudaHostAlloc failed");

   stat = cudaMalloc(&energygrp_vdw_d[rank],grppener->nener*sizeof(real));
   CU_RET_ERR(stat, "cudaMalloc failed");
  

   stat = cudaHostAlloc(&energygrp_vdw[rank],grppener->nener*sizeof(real),cudaHostAllocPortable);
   CU_RET_ERR(stat, "cudaHostAlloc failed");   

  stat = cudaMalloc(&fr_pairsTable_data_d[rank],(fr->pairsTable->n + 1)*fr->pairsTable->stride*sizeof(real));
  CU_RET_ERR(stat, "cudaMalloc failed");

  stat = cudaMemcpy(fr_pairsTable_data_d[rank],fr->pairsTable->data,
          (fr->pairsTable->n + 1)*fr->pairsTable->stride*sizeof(real),cudaMemcpyHostToDevice);
  CU_RET_ERR(stat, "cudaMemcpy failed");  

  stat = cudaMalloc(&f_shift_d[rank],sizeof(rvec)*SHIFTS);
  CU_RET_ERR(stat, "cudaMalloc failed");

  f_shift[rank] = (rvec *) malloc(sizeof(rvec)*SHIFTS);

  stat = cudaMalloc(&force_max_d[rank],sizeof(int));
  CU_RET_ERR(stat, "cudaMalloc failed");

  stat = cudaStreamCreate(&streamBonded[rank]);
  CU_RET_ERR(stat, "cudaStreamCreate failed");
}

//----//
// update after a neighbour list update step
void 
update_gpu_bonded( const t_idef *idef, const t_forcerec *fr, matrix box,
                  int xsize, const t_mdatoms *md, const real *lambda,
                  gmx_grppairener_t *grppener )
{
  
  cudaError_t stat;
  int nat,nbonds;
  int rank=0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  const real alloc_factor=1.2;
  int xsize_alloc;

  copy_required_this_step[rank]=true;

  t_pbc *pbc_null = NULL;
  t_pbc pbc;
  set_pbc(&pbc, fr->ePBC, box);

  if (fr->bMolPBC)
  {
    pbc_null=&pbc;
  }
  else
  {
    pbc_null = NULL;
  } 

  if(!bonded_init_done[rank]) // currently static so inits to 0
  {
    init_gpu_bonded(pbc_null, rank, grppener, fr);
    bonded_init_done[rank]=true;
  }

bool x_resized=false;
xsize_alloc=0;
int x_needed = (xsize > md->nr) ? xsize : md->nr ;
  if(x_needed > xbonded_size_max[rank]) 
  {
    x_resized=true;
    xsize_alloc=alloc_factor*x_needed;
    if(md_chargeA_d[rank]) 
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
    stat = cudaMalloc(&x_bonded_d[rank],sizeof(rvec)*xsize_alloc);
    CU_RET_ERR(stat, "cudaMalloc failed");
    stat = cudaMalloc(&f_bonded_d[rank],sizeof(rvec)*xsize_alloc);
    CU_RET_ERR(stat, "cudaMalloc failed");
#endif

    f_bonded[rank] = (rvec *) malloc(sizeof(rvec)*xsize_alloc);
    stat = cudaMalloc(&md_cENER_d[rank],sizeof(unsigned short)*xsize_alloc); 
    CU_RET_ERR(stat, "cudaMalloc failed");
    stat = cudaMalloc(&md_chargeA_d[rank],sizeof(real)*xsize_alloc); 
    CU_RET_ERR(stat, "cudaMalloc failed");

    xbonded_size_max[rank]=xsize_alloc;
  }

  stat = cudaMemcpy(md_cENER_d[rank],md->cENER,md->nr*sizeof(unsigned short),cudaMemcpyHostToDevice);
  CU_RET_ERR(stat, "cudaMemcpy failed");
  stat = cudaMemcpy(md_chargeA_d[rank],md->chargeA,md->nr*sizeof(real),cudaMemcpyHostToDevice);
  CU_RET_ERR(stat, "cudaMemcpy failed");

//forceparams
  int ntypes = idef->ntypes;
  const t_iparams *iparams;
  iparams = idef->iparams;

  if(ntypes > ntypes_alloc[rank]){
     if( forceparams_d[rank]) {
        stat = cudaFree(forceparams_d[rank]);
        CU_RET_ERR(stat, "cudaFree failed");
     }
     stat = cudaMalloc(&forceparams_d[rank],sizeof(t_iparams)*ntypes);
     CU_RET_ERR(stat, "cudaMalloc failed");
     ntypes_alloc[rank]=ntypes;
  }
  
  stat = cudaMemcpyAsync(forceparams_d[rank],iparams,ntypes*sizeof(t_iparams),cudaMemcpyHostToDevice,streamBonded[rank]);

// forceatoms
  for(int i=0;i<F_NRE;i++)
  {
    if(idef->il[i].nr > 0 && ftype_is_bonded_potential(i) ) 
    {
      if(x_resized)
      {
        if(force_next_d[rank][i])
        {
          stat = cudaFree(force_next_d[rank][i]);
          CU_RET_ERR(stat, "cudaFree failed");
        }
        stat = cudaMalloc(&force_next_d[rank][i],sizeof(int)*xsize_alloc);
        CU_RET_ERR(stat, "cudaMalloc failed");
      }     
     
      //reallocate if we do not have enough memory
      if(idef->il[i].nalloc > device_iatom_alloc[rank][i])
      { 
        if(iatoms_d[rank][i]) {
          stat = cudaFree(iatoms_d[rank][i]);
          CU_RET_ERR(stat, "cudaFree failed");
        }
        stat = cudaMalloc(&iatoms_d[rank][i],idef->il[i].nalloc*sizeof(t_iatom));
        CU_RET_ERR(stat, "cudaMalloc failed");
        device_iatom_alloc[rank][i] = idef->il[i].nalloc;
      }
      //copy  the data
      stat = cudaMemcpy(iatoms_d[rank][i],idef->il[i].iatoms,idef->il[i].nr*sizeof(t_iatom)
         ,cudaMemcpyHostToDevice);
      CU_RET_ERR(stat, "cudaMemcpy failed");
    }
  }
//pbc
//this can be null
  if(pbc_null) 
  {
    stat = cudaMemcpy(pbc_d[rank],pbc_null,sizeof(t_pbc),cudaMemcpyHostToDevice);
    CU_RET_ERR(stat, "cudaMemcpy failed"); 
    stat = cudaMemcpy(pbc_hbox_diag_d[rank],pbc_null->hbox_diag,sizeof(rvec),cudaMemcpyHostToDevice);
    CU_RET_ERR(stat, "cudaMemcpy failed");
    stat = cudaMemcpy(pbc_box_d[rank],pbc_null->box,sizeof(matrix),cudaMemcpyHostToDevice);
    CU_RET_ERR(stat, "cudaMemcpy failed");
    stat = cudaMemcpy(pbc_mhbox_diag_d[rank],pbc_null->mhbox_diag,sizeof(rvec),cudaMemcpyHostToDevice);
    CU_RET_ERR(stat, "cudaMemcpy failed");
    stat = cudaMemcpy(pbc_fbox_diag_d[rank],pbc_null->fbox_diag,sizeof(rvec),cudaMemcpyHostToDevice);
    CU_RET_ERR(stat, "cudaMemcpy failed");
    stat = cudaMemcpy(pbc_tric_vec_d[rank],pbc_null->tric_vec,MAX_NTRICVEC*sizeof(rvec),cudaMemcpyHostToDevice);
    CU_RET_ERR(stat, "cudaMemcpy failed");
    stat = cudaMemcpy(pbc_tric_shift_d[rank],pbc_null->tric_shift,MAX_NTRICVEC*sizeof(ivec),cudaMemcpyHostToDevice);
    CU_RET_ERR(stat, "cudaMemcpy failed");
  } 

// Packed arrays for LJ14

   pairs_Rpacked[rank][frFUDGEQQ]=fr->fudgeQQ;
   pairs_Rpacked[rank][fr_icEPSFAC]=fr->ic->epsfac;
   pairs_Rpacked[rank][fr_PAIRS_SCALE]=fr->pairsTable->scale;

   stat = cudaMemcpy(pairs_Rpacked_d[rank],pairs_Rpacked[rank],
                     sizeof(real)*PAIRS_NR,cudaMemcpyHostToDevice);
   CU_RET_ERR(stat, "cudaMemcpy failed");

   pairs_Ipacked[rank][fr_PAIRS_STRIDE]=fr->pairsTable->stride;
   pairs_Ipacked[rank][md_NENERGRP]=md->nenergrp;

   stat = cudaMemcpy(pairs_Ipacked_d[rank],pairs_Ipacked[rank],
                     sizeof(int)*PAIRS_NI,cudaMemcpyHostToDevice);
   CU_RET_ERR(stat, "cudaMemcpy failed");

  dim3 blocks,blocks_reset;
  dim3 threads (TPB_BONDED,1,1);
  blocks_reset.x= (xbonded_size_max[rank]+TPB_BONDED-1)/TPB_BONDED;

// distributed arrays to remove atomics
  for (int ftype = 0; (ftype < F_NRE) ; ftype++)
  {
    if(idef->il[ftype].nr > 0 && ftype_is_bonded_potential(ftype))
    {
      nat      = interaction_function[ftype].nratoms;
      nbonds    = idef->il[ftype].nr/(nat+1);
// resize if needed
      if( idef->il[ftype].nr > force_mapping_allocated[rank][ftype] )
      {
        if ( force_mapping_allocated[rank][ftype] != 0) {
          stat = cudaFree(force_mapping_d[rank][ftype]);
          CU_RET_ERR(stat, "cudaFree failed");
        }
        stat = cudaMalloc(&force_mapping_d[rank][ftype],sizeof(int)*nbonds*nat);
        CU_RET_ERR(stat, "cudaMalloc failed");
        force_mapping_allocated[rank][ftype] = nbonds*nat;
      }
      blocks.x = (nbonds+TPB_BONDED-1)/TPB_BONDED;
      reset_next<<<blocks_reset,threads,0,streamBonded[rank]>>>(xbonded_size_max[rank],force_next_d[rank][ftype]);
      calc_force_mapping<<<blocks,threads,0,streamBonded[rank]>>> ( nbonds , nat, iatoms_d[rank][ftype],
         force_mapping_d[rank][ftype] , force_next_d[rank][ftype], force_max_d[rank]);       
    }
  }
  stat = cudaMemcpyAsync(&force_max[rank],force_max_d[rank],
           sizeof(int),cudaMemcpyDeviceToHost,streamBonded[rank]);
  cudaStreamSynchronize(streamBonded[rank]);
  stat = cudaGetLastError();
  CU_RET_ERR(stat, "Async error");
 
// reallocate force_distributed_d if needed   
  if(xsize*(force_max[rank]+1)>force_distributed_allocated[rank]) {
    if(force_distributed_d[rank]) {
       stat = cudaFree(force_distributed_d[rank]);
       CU_RET_ERR(stat, "cudaFree failed");
    }
    stat = cudaMalloc(&force_distributed_d[rank],sizeof(rvec)*xsize*(force_max[rank]+1));
    CU_RET_ERR(stat, "cudaMalloc failed");
    force_distributed_allocated[rank]=xsize*(force_max[rank]+1);
  }
  cudaStreamSynchronize(streamBonded[rank]); 
}


/* -------------- */
// reset the internal variables
void 
reset_gpu_bonded(const int size, const int nener )
{
  int rank=0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  dim3 blocks ( (size+TPB_BONDED-1)/TPB_BONDED,1,1);
  dim3 threads (TPB_BONDED,1,1);

  cudaError_t stat;

#ifdef STANDALONE
  rvec* f_d_in = f_bonded_d[rank];
#else
  rvec* f_d_in = gpuBufferOpsGetFPtr();
#endif

  reset_gpu_bonded_kernel <<<blocks,threads,0,streamBonded[rank]>>>
    (vtot_d[rank], dvdlambda_d[rank], f_d_in, size, energygrp_elec_d[rank],
     energygrp_vdw_d[rank],nener,f_shift_d[rank] );

  stat = cudaGetLastError();
  CU_RET_ERR(stat, "reset bonded kernel failed");
}

/* -------------- */
void 
do_bonded_gpu(t_forcerec *fr, const t_inputrec *ir, const t_idef *idef, 
              int flags,  const t_graph *graph, int natoms, rvec x[], 
              real *lambda, const t_mdatoms *md, 
              rvec *input_force, t_lambda *fepvals, gmx_enerdata_t *enerd)
{
  int nat1,nbonds,efptFTYPE;

  bool  bCalcEnerVir = (flags & (GMX_FORCE_VIRIAL | GMX_FORCE_ENERGY));

  bool abort=false;
  int rank=0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

#ifdef STANDALONE
  rvec* x_d_in = x_bonded_d[rank];
  rvec* f_d_in = f_bonded_d[rank];
#else
  rvec* x_d_in = gpuBufferOpsGetXPtr();
  rvec* f_d_in = gpuBufferOpsGetFPtr();
#endif  

  cudaError_t stat;

// sanity check

  if (fr->bQMMM) abort=true;
  if (ir->nwall) abort=true;
  if (ir->implicit_solvent) abort=true;
  if ((fr->cutoff_scheme == ecutsGROUP) && (flags & GMX_FORCE_NONBONDED)) abort=true;
  if (graph)  abort=true; 
  if (fepvals->n_lambda > 0  && (flags & GMX_FORCE_DHDL) ) abort=true;
  if (flags & GMX_FORCE_DHDL) abort=true;
  if (fr->efep != efepNO) abort=true;

  if(abort) 
  {
    printf("calculation not supported on GPU\n");
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  bool do_copy_this_step=false;

#ifdef STANDALONE
  copy_required_this_step[rank]=true;  
#endif

  if(copy_required_this_step[rank]){
      do_copy_this_step=true;
  }

  if(do_copy_this_step){
    stat = cudaMemcpyAsync(x_d_in,x,
               natoms*sizeof(rvec),cudaMemcpyHostToDevice,streamBonded[rank]);
    CU_RET_ERR(stat, "cudaMemcpy failed");
  }

// launch kernels
// reordered for better overlap  
  dim3 blocks,blocks_natoms;
  dim3 threads (TPB_BONDED,1,1);
  for (int ftype = 0; (ftype < F_NRE) ; ftype++)
  {
    if(idef->il[ftype].nr > 0 && ftype_is_bonded_potential(ftype))
    {
      nat1      = interaction_function[ftype].nratoms + 1;
      nbonds    = idef->il[ftype].nr/nat1;
      blocks.x = (nbonds+TPB_BONDED-1)/TPB_BONDED;
      blocks_natoms.x = (natoms+TPB_BONDED-1)/TPB_BONDED;
      if (IS_RESTRAINT_TYPE(ftype))
      {
        efptFTYPE = efptRESTRAINT;
      }
      else
      {
        efptFTYPE = efptBONDED;
      }
       if(ftype == F_PDIHS || ftype == F_PIDIHS ) {
        if(bCalcEnerVir) {
         pdihs_gpu <<<blocks,threads,0,streamBonded[rank]>>>
           (&vtot_d[rank][ftype],nbonds,natoms,
           iatoms_d[rank][ftype],forceparams_d[rank],
           x_d_in, force_distributed_d[rank], f_shift_d[rank],
           pbc_d[rank], pbc_hbox_diag_d[rank],
           pbc_box_d[rank],pbc_mhbox_diag_d[rank],
           pbc_fbox_diag_d[rank],pbc_tric_vec_d[rank],      
           pbc_tric_shift_d[rank],
           lambda[efptFTYPE], &dvdlambda_d[rank][efptFTYPE],
           force_mapping_d[rank][ftype]);
        } else {
           pdihs_gpu_noener <<<blocks,threads,0,streamBonded[rank]>>>
           (&vtot_d[rank][ftype],nbonds,natoms,
           iatoms_d[rank][ftype],forceparams_d[rank]  ,
           x_d_in, force_distributed_d[rank], f_shift_d[rank],
           pbc_d[rank],pbc_hbox_diag_d[rank],
           pbc_box_d[rank],pbc_mhbox_diag_d[rank],
           pbc_fbox_diag_d[rank],pbc_tric_vec_d[rank],
           pbc_tric_shift_d[rank], 
           lambda[efptFTYPE], &dvdlambda_d[rank][efptFTYPE],
           force_mapping_d[rank][ftype]);
        }
        consolidate_forces<<<blocks_natoms,threads,0,streamBonded[rank]>>>
        (natoms, f_d_in, force_distributed_d[rank],force_next_d[rank][ftype]);
      }

    }
  } 
  
  for (int ftype = 0; (ftype < F_NRE) ; ftype++)
  {
    if(idef->il[ftype].nr > 0 && ftype_is_bonded_potential(ftype)) 
    {
      nat1      = interaction_function[ftype].nratoms + 1;
      nbonds    = idef->il[ftype].nr/nat1;
      blocks.x = (nbonds+TPB_BONDED-1)/TPB_BONDED;
      blocks_natoms.x = (natoms+TPB_BONDED-1)/TPB_BONDED;
      if (IS_RESTRAINT_TYPE(ftype))
      {
        efptFTYPE = efptRESTRAINT;
      }
      else
      {
        efptFTYPE = efptBONDED;
      }
// in the main code they have a function pointer for this
// so we need something like that for final version    

      if(ftype == F_BONDS) {
        if(bCalcEnerVir) {
         bonds_gpu <<<blocks,threads,0,streamBonded[rank]>>> 
          (&vtot_d[rank][ftype],nbonds,natoms,
          iatoms_d[rank][ftype],forceparams_d[rank]  ,
          x_d_in, force_distributed_d[rank], f_shift_d[rank], 
          pbc_d[rank], pbc_hbox_diag_d[rank],
          pbc_box_d[rank],pbc_mhbox_diag_d[rank],
          pbc_fbox_diag_d[rank],pbc_tric_vec_d[rank],
          pbc_tric_shift_d[rank],
          lambda[efptFTYPE], &dvdlambda_d[rank][efptFTYPE],
          force_mapping_d[rank][ftype]);
        } else {
          bonds_gpu_noener <<<blocks,threads,0,streamBonded[rank]>>>
          (&vtot_d[rank][ftype],nbonds,natoms,
          iatoms_d[rank][ftype],forceparams_d[rank]  ,
          x_d_in, force_distributed_d[rank], f_shift_d[rank],
          pbc_d[rank], pbc_hbox_diag_d[rank],
           pbc_box_d[rank],pbc_mhbox_diag_d[rank],
           pbc_fbox_diag_d[rank],pbc_tric_vec_d[rank],
           pbc_tric_shift_d[rank],
          lambda[efptFTYPE], &dvdlambda_d[rank][efptFTYPE],
          force_mapping_d[rank][ftype]);

        }

        consolidate_forces<<<blocks_natoms,threads,0,streamBonded[rank]>>>
        (natoms, f_d_in, force_distributed_d[rank],force_next_d[rank][ftype]);
      }

      if(ftype == F_ANGLES) {
       if(bCalcEnerVir) {
        angles_gpu <<<blocks,threads,0,streamBonded[rank]>>>
          (&vtot_d[rank][ftype],nbonds,natoms,
          iatoms_d[rank][ftype],forceparams_d[rank]  ,
          x_d_in, force_distributed_d[rank], f_shift_d[rank],
          pbc_d[rank], pbc_hbox_diag_d[rank],
           pbc_box_d[rank],pbc_mhbox_diag_d[rank],
          pbc_fbox_diag_d[rank],pbc_tric_vec_d[rank],
          pbc_tric_shift_d[rank],
          lambda[efptFTYPE], &dvdlambda_d[rank][efptFTYPE],
          force_mapping_d[rank][ftype]);
        } else {
          angles_gpu_noener <<<blocks,threads,0,streamBonded[rank]>>>
          (&vtot_d[rank][ftype],nbonds,natoms,
          iatoms_d[rank][ftype],forceparams_d[rank]  ,
          x_d_in, force_distributed_d[rank], f_shift_d[rank],
          pbc_d[rank], pbc_hbox_diag_d[rank],
          pbc_box_d[rank],pbc_mhbox_diag_d[rank],
          pbc_fbox_diag_d[rank],pbc_tric_vec_d[rank],
          pbc_tric_shift_d[rank],
          lambda[efptFTYPE], &dvdlambda_d[rank][efptFTYPE],
          force_mapping_d[rank][ftype]);

        }
 
        consolidate_forces<<<blocks_natoms,threads,0,streamBonded[rank]>>>
        (natoms, f_d_in, force_distributed_d[rank],force_next_d[rank][ftype]);
      }

      if(ftype == F_UREY_BRADLEY) {
        if(bCalcEnerVir) {
         urey_bradley_gpu <<<blocks,threads,0,streamBonded[rank]>>>
          (&vtot_d[rank][ftype],nbonds,natoms,
          iatoms_d[rank][ftype],forceparams_d[rank]  ,
          x_d_in, force_distributed_d[rank], f_shift_d[rank],
          pbc_d[rank], pbc_hbox_diag_d[rank],
          pbc_box_d[rank],pbc_mhbox_diag_d[rank],
          pbc_fbox_diag_d[rank],pbc_tric_vec_d[rank],
          pbc_tric_shift_d[rank], 
          lambda[efptFTYPE], &dvdlambda_d[rank][efptFTYPE],
          force_mapping_d[rank][ftype]);
        } else {
         urey_bradley_gpu_noener <<<blocks,threads,0,streamBonded[rank]>>>
          (&vtot_d[rank][ftype],nbonds,natoms,
          iatoms_d[rank][ftype],forceparams_d[rank]  ,
          x_d_in, force_distributed_d[rank], f_shift_d[rank],
          pbc_d[rank], pbc_hbox_diag_d[rank],
          pbc_box_d[rank],pbc_mhbox_diag_d[rank],
          pbc_fbox_diag_d[rank],pbc_tric_vec_d[rank],
          pbc_tric_shift_d[rank],       
          lambda[efptFTYPE], &dvdlambda_d[rank][efptFTYPE],
          force_mapping_d[rank][ftype]);
       }

        consolidate_forces<<<blocks_natoms,threads,0,streamBonded[rank]>>>
        (natoms, f_d_in, force_distributed_d[rank],force_next_d[rank][ftype]);
      }

       if(ftype == F_RBDIHS  ) {
         if(bCalcEnerVir) {
         rbdihs_gpu <<<blocks,threads,0,streamBonded[rank]>>>
           (&vtot_d[rank][ftype],nbonds,natoms,
           iatoms_d[rank][ftype],forceparams_d[rank]  ,
           x_d_in, force_distributed_d[rank], f_shift_d[rank],
           pbc_d[rank], pbc_hbox_diag_d[rank],
           pbc_box_d[rank],pbc_mhbox_diag_d[rank],
           pbc_fbox_diag_d[rank],pbc_tric_vec_d[rank],
           pbc_tric_shift_d[rank],      
           lambda[efptFTYPE], &dvdlambda_d[rank][efptFTYPE],
           force_mapping_d[rank][ftype]);
         } else {
           rbdihs_gpu_noener <<<blocks,threads,0,streamBonded[rank]>>>
           (&vtot_d[rank][ftype],nbonds,natoms,
           iatoms_d[rank][ftype],forceparams_d[rank]  ,
           x_d_in, force_distributed_d[rank], f_shift_d[rank],
           pbc_d[rank], pbc_hbox_diag_d[rank],
           pbc_box_d[rank],pbc_mhbox_diag_d[rank],
           pbc_fbox_diag_d[rank],pbc_tric_vec_d[rank],
           pbc_tric_shift_d[rank],
           lambda[efptFTYPE], &dvdlambda_d[rank][efptFTYPE],
           force_mapping_d[rank][ftype]);
         }

        consolidate_forces<<<blocks_natoms,threads,0,streamBonded[rank]>>>
        (natoms, f_d_in, force_distributed_d[rank],force_next_d[rank][ftype]);
      }

      if(ftype == F_IDIHS  ) {
         if(bCalcEnerVir) {
         idihs_gpu <<<blocks,threads,0,streamBonded[rank]>>>
           (&vtot_d[rank][ftype],nbonds,natoms,
           iatoms_d[rank][ftype],forceparams_d[rank]  ,
           x_d_in, force_distributed_d[rank], f_shift_d[rank],
           pbc_d[rank], pbc_hbox_diag_d[rank],
           pbc_box_d[rank],pbc_mhbox_diag_d[rank],
           pbc_fbox_diag_d[rank],pbc_tric_vec_d[rank],
           pbc_tric_shift_d[rank],
           lambda[efptFTYPE], &dvdlambda_d[rank][efptFTYPE],
           force_mapping_d[rank][ftype]);
         } else {
           idihs_gpu_noener <<<blocks,threads,0,streamBonded[rank]>>>
           (&vtot_d[rank][ftype],nbonds,natoms,
           iatoms_d[rank][ftype],forceparams_d[rank]  ,
           x_d_in, force_distributed_d[rank], f_shift_d[rank],
           pbc_d[rank], pbc_hbox_diag_d[rank],
           pbc_box_d[rank],pbc_mhbox_diag_d[rank],
           pbc_fbox_diag_d[rank],pbc_tric_vec_d[rank],
           pbc_tric_shift_d[rank], 
           lambda[efptFTYPE], &dvdlambda_d[rank][efptFTYPE],
           force_mapping_d[rank][ftype]);
         }
        consolidate_forces<<<blocks_natoms,threads,0,streamBonded[rank]>>>
        (natoms, f_d_in, force_distributed_d[rank],force_next_d[rank][ftype]);

      }  
     
      if(ftype == F_LJ14) {
         if(bCalcEnerVir) {
           pairs_gpu <<<blocks,threads,0,streamBonded[rank]>>>
            (ftype,nbonds,natoms,
            iatoms_d[rank][ftype],forceparams_d[rank]  ,
            x_d_in, force_distributed_d[rank], f_shift_d[rank],
            pbc_d[rank], pbc_hbox_diag_d[rank],
            pbc_box_d[rank],pbc_mhbox_diag_d[rank],
            pbc_fbox_diag_d[rank],pbc_tric_vec_d[rank],
            pbc_tric_shift_d[rank], 
            md_chargeA_d[rank],md_cENER_d[rank],
            energygrp_elec_d[rank], energygrp_vdw_d[rank],
            fr_pairsTable_data_d[rank],
            pairs_Rpacked_d[rank],pairs_Ipacked_d[rank],
            force_mapping_d[rank][ftype]
            );
          } else {
            pairs_gpu_noener <<<blocks,threads,0,streamBonded[rank]>>>
            (ftype,nbonds,natoms,
            iatoms_d[rank][ftype],forceparams_d[rank]  ,
            x_d_in, force_distributed_d[rank], f_shift_d[rank],
            pbc_d[rank], pbc_hbox_diag_d[rank],
            pbc_box_d[rank],pbc_mhbox_diag_d[rank],
            pbc_fbox_diag_d[rank],pbc_tric_vec_d[rank],
            pbc_tric_shift_d[rank], 
            md_chargeA_d[rank],md_cENER_d[rank],
            energygrp_elec_d[rank], energygrp_vdw_d[rank],
            fr_pairsTable_data_d[rank],
            pairs_Rpacked_d[rank],pairs_Ipacked_d[rank],
            force_mapping_d[rank][ftype]
            );

          }
          consolidate_forces<<<blocks_natoms,threads,0,streamBonded[rank]>>>
           (natoms, f_d_in, force_distributed_d[rank],force_next_d[rank][ftype]);
      } 
    }
  }  

  stat = cudaGetLastError();
  CU_RET_ERR(stat, "kernels failed");
  return;
}

void
do_bonded_gpu_finalize(t_forcerec *fr,
              int flags, int natoms, rvec *input_force, gmx_enerdata_t *enerd)
{
   int i;
   bool bCalcEnerVir = (flags & (GMX_FORCE_VIRIAL | GMX_FORCE_ENERGY));
   cudaError_t stat;
   
   gmx_grppairener_t *grppener;
   grppener=&enerd->grpp ;

   int rank=0;
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);

#ifdef STANDALONE
  rvec* f_d_in = f_bonded_d[rank];
#else
  rvec* f_d_in = gpuBufferOpsGetFPtr();
#endif

  bool do_copy_this_step=false;

  if(copy_required_this_step[rank]){
      do_copy_this_step=true;
      copy_required_this_step[rank]=false;
  }

//copy forces back
  if(do_copy_this_step){
    stat = cudaMemcpyAsync(f_bonded[rank],f_d_in,
               natoms*sizeof(rvec),cudaMemcpyDeviceToHost,streamBonded[rank]);
    cudaStreamSynchronize(streamBonded[rank]);
    CU_RET_ERR(stat, "cudaMemcpy failed");
  
    for( i=0 ; i< natoms ; i++)
    {
       rvec_inc(input_force[i],f_bonded[rank][i]);
    }
  }
 
  if(bCalcEnerVir) {
  // shift Forces
    stat = cudaMemcpyAsync(f_shift[rank],f_shift_d[rank],
             SHIFTS*sizeof(rvec),cudaMemcpyDeviceToHost,streamBonded[rank]);
    cudaStreamSynchronize(streamBonded[rank]);
  CU_RET_ERR(stat, "cudaMemcpy failed");

  for (i = 0; i < SHIFTS; i++)
  {
    rvec_inc(fr->fshift[i], f_shift[rank][i]);
  }

  // copy energies back
  stat = cudaMemcpyAsync(vtot[rank],vtot_d[rank],
             F_NRE*sizeof(real),cudaMemcpyDeviceToHost,streamBonded[rank]);
  cudaStreamSynchronize(streamBonded[rank]);
  CU_RET_ERR(stat, "cudaMemcpy failed");

   for (i = 0; i < F_NRE; i++)
   {
     enerd->term[i] += vtot[rank][i];
   }

   stat = cudaMemcpyAsync(energygrp_vdw[rank],energygrp_vdw_d[rank],
        sizeof(real)*grppener->nener,cudaMemcpyDeviceToHost,streamBonded[rank]);
   CU_RET_ERR(stat, "cudaMemcpy failed");

  stat = cudaMemcpyAsync(energygrp_elec[rank],energygrp_elec_d[rank],
        sizeof(real)*grppener->nener,cudaMemcpyDeviceToHost,streamBonded[rank]);
   CU_RET_ERR(stat, "cudaMemcpy failed");
   cudaStreamSynchronize(streamBonded[rank]);
  for (i = 0; i < grppener->nener; i++)
  {
      grppener->ener[egCOUL14][i] += energygrp_elec[rank][i];
      grppener->ener[egLJ14][i]  +=  energygrp_vdw[rank][i];
  }
  }
  cudaStreamSynchronize(streamBonded[rank]); // needed if we have no copies for bufferops sync
}
/*------------------------------------------------------------------------------*/

