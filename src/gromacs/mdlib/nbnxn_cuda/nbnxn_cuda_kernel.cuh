/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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
 *  \brief
 *  CUDA non-bonded kernel used through preprocessor-based code generation
 *  of multiple kernel flavors, see nbnxn_cuda_kernels.cuh.
 *
 *  NOTE: No include fence as it is meant to be included multiple times.
 *
 *  \author Szilárd Páll <pall.szilard@gmail.com>
 *  \ingroup module_mdlib
 */
#include "config.h"

#include "gromacs/math/utilities.h"
#include "gromacs/pbcutil/ishift.h"
/* Note that floating-point constants in CUDA code should be suffixed
 * with f (e.g. 0.5f), to stop the compiler producing intermediate
 * code that is in double precision.
 */

#if __CUDA_ARCH__ >= 300
/* Note: convenience macros, need to be undef-ed at the end of the file. */
#define REDUCE_SHUFFLE
/* On Kepler pre-loading i-atom types to shmem gives a few %,
   but on Fermi it does not */
#define IATYPE_SHMEM
#ifdef HAVE_CUDA_TEXOBJ_SUPPORT
#define USE_TEXOBJ
#endif
#endif

#if defined EL_EWALD_ANA || defined EL_EWALD_TAB
/* Note: convenience macro, needs to be undef-ed at the end of the file. */
#define EL_EWALD_ANY
#endif

#if defined EL_EWALD_ANY || defined EL_RF || defined LJ_EWALD || (defined EL_CUTOFF && defined CALC_ENERGIES)
/* Macro to control the calculation of exclusion forces in the kernel
 * We do that with Ewald (elec/vdw) and RF. Cut-off only has exclusion
 * energy terms.
 *
 * Note: convenience macro, needs to be undef-ed at the end of the file.
 */
#define EXCLUSION_FORCES
#endif

#if defined LJ_EWALD_COMB_GEOM || defined LJ_EWALD_COMB_LB
/* Note: convenience macro, needs to be undef-ed at the end of the file. */
#define LJ_EWALD
#endif


/*
   Kernel launch parameters:
    - #blocks   = #pair lists, blockId = pair list Id
    - #threads  = NTHREAD_Z * CL_SIZE^2
    - shmem     = see nbnxn_cuda.cu:calc_shmem_required()

    Each thread calculates an i force-component taking one pair of i-j atoms.
 */

/* Kernel launch bounds as function of NTHREAD_Z.
 * - CC 3.5/5.2: NTHREAD_Z=1, (64, 16) bounds
 * - CC 3.7:     NTHREAD_Z=2, (128, 16) bounds
 *
 * Note: convenience macros, need to be undef-ed at the end of the file.
 */
#if __CUDA_ARCH__ == 370
#define NTHREAD_Z           (2)
#define MIN_BLOCKS_PER_MP   (16)
#else
#define NTHREAD_Z           (1)
#define MIN_BLOCKS_PER_MP   (16)
#endif
#define THREADS_PER_BLOCK   (CL_SIZE*CL_SIZE*NTHREAD_Z)

#if __CUDA_ARCH__ >= 350
__launch_bounds__(THREADS_PER_BLOCK, MIN_BLOCKS_PER_MP)
#else
__launch_bounds__(THREADS_PER_BLOCK)
#endif
#ifdef PRUNE_NBL
#ifdef CALC_ENERGIES
__global__ void NB_KERNEL_FUNC_NAME(nbnxn_kernel, _VF_prune_cuda)
#else
__global__ void NB_KERNEL_FUNC_NAME(nbnxn_kernel, _F_prune_cuda)
#endif
#else
#ifdef CALC_ENERGIES
__global__ void NB_KERNEL_FUNC_NAME(nbnxn_kernel, _VF_cuda)
#else
__global__ void NB_KERNEL_FUNC_NAME(nbnxn_kernel, _F_cuda)
#endif
#endif
(const cu_atomdata_t atdat,
 const cu_nbparam_t nbparam,
 const cu_plist_t plist,
 bool bCalcFshift)
{
    /* convenience variables */
    const nbnxn_sci_t *pl_sci       = plist.sci;
#ifndef PRUNE_NBL
    const
#endif
    nbnxn_cj4_t        *pl_cj4      = plist.cj4;
    const nbnxn_excl_t *excl        = plist.excl;
    const int          *atom_types  = atdat.atom_types;
    int                 ntypes      = atdat.ntypes;
    const float4       *xq          = atdat.xq;
    float3             *f           = atdat.f;
    const float3       *shift_vec   = atdat.shift_vec;
    float               rcoulomb_sq = nbparam.rcoulomb_sq;
#ifdef VDW_CUTOFF_CHECK
    float               rvdw_sq     = nbparam.rvdw_sq;
    float               vdw_in_range;
#endif
#ifdef LJ_EWALD
    float               lje_coeff2, lje_coeff6_6;
#endif
#ifdef EL_RF
    float two_k_rf              = nbparam.two_k_rf;
#endif
#ifdef EL_EWALD_TAB
    float coulomb_tab_scale     = nbparam.coulomb_tab_scale;
#endif
#ifdef EL_EWALD_ANA
    float beta2                 = nbparam.ewald_beta*nbparam.ewald_beta;
    float beta3                 = nbparam.ewald_beta*nbparam.ewald_beta*nbparam.ewald_beta;
#endif
#ifdef PRUNE_NBL
    float rlist_sq              = nbparam.rlist_sq;
#endif

#ifdef CALC_ENERGIES
#ifdef EL_EWALD_ANY
    float  beta        = nbparam.ewald_beta;
    float  ewald_shift = nbparam.sh_ewald;
#else
    float  c_rf        = nbparam.c_rf;
#endif /* EL_EWALD_ANY */
    float *e_lj        = atdat.e_lj;
    float *e_el        = atdat.e_el;
#endif /* CALC_ENERGIES */

    /* thread/block/warp id-s */
    unsigned int tidxi  = threadIdx.x;
    unsigned int tidxj  = threadIdx.y;
    unsigned int tidx   = threadIdx.y * blockDim.x + threadIdx.x;
#if NTHREAD_Z == 1
    unsigned int tidxz  = 0;
#else
    unsigned int tidxz  = threadIdx.z;
#endif
    unsigned int bidx   = blockIdx.x;
    unsigned int widx   = tidx / WARP_SIZE; /* warp index */

    int          sci, ci, cj, ci_offset,
                 ai, aj,
                 cij4_start, cij4_end,
                 typei, typej,
                 i, jm, j4, wexcl_idx;
    float        qi, qj_f,
                 r2, inv_r, inv_r2, inv_r6,
                 c6, c12,
                 int_bit,
                 F_invr;
#ifdef CALC_ENERGIES
    float        E_lj, E_el;
#endif
#if defined CALC_ENERGIES || defined LJ_POT_SWITCH
    float        E_lj_p;
#endif
    unsigned int wexcl, imask, mask_ji;
    float4       xqbuf;
    float3       xi, xj, rv, f_ij, fcj_buf;
    float3       fci_buf[NCL_PER_SUPERCL]; /* i force buffer */
    nbnxn_sci_t  nb_sci;

    /* shmem buffer for i x+q pre-loading */
    extern __shared__  float4 xqib[];
    /* shmem buffer for cj, for each warp separately */
    int *cjs     = ((int *)(xqib + NCL_PER_SUPERCL * CL_SIZE)) + tidxz * 2 * NBNXN_GPU_JGROUP_SIZE;
#ifdef IATYPE_SHMEM
    /* shmem buffer for i atom-type pre-loading */
    int *atib    = ((int *)(xqib + NCL_PER_SUPERCL * CL_SIZE)) + NTHREAD_Z * 2 * NBNXN_GPU_JGROUP_SIZE;
#endif

#ifndef REDUCE_SHUFFLE
    /* shmem j force buffer */
#ifdef IATYPE_SHMEM
    float *f_buf = (float *)(atib + NCL_PER_SUPERCL * CL_SIZE);
#else
    float *f_buf = (float *)(cjs + NTHREAD_Z * 2 * NBNXN_GPU_JGROUP_SIZE);
#endif
#endif

    nb_sci      = pl_sci[bidx];         /* my i super-cluster's index = current bidx */
    sci         = nb_sci.sci;           /* super-cluster */
    cij4_start  = nb_sci.cj4_ind_start; /* first ...*/
    cij4_end    = nb_sci.cj4_ind_end;   /* and last index of j clusters */

    if (tidxz == 0)
    {
        /* Pre-load i-atom x and q into shared memory */
        ci = sci * NCL_PER_SUPERCL + tidxj;
        ai = ci * CL_SIZE + tidxi;
        xqib[tidxj * CL_SIZE + tidxi] = xq[ai] + shift_vec[nb_sci.shift];
#ifdef IATYPE_SHMEM
        /* Pre-load the i-atom types into shared memory */
        atib[tidxj * CL_SIZE + tidxi] = atom_types[ai];
#endif
    }
    __syncthreads();

    for (ci_offset = 0; ci_offset < NCL_PER_SUPERCL; ci_offset++)
    {
        fci_buf[ci_offset] = make_float3(0.0f);
    }

#ifdef LJ_EWALD
    /* TODO: we are trading registers with flops by keeping lje_coeff-s, try re-calculating it later */
    lje_coeff2   = nbparam.ewaldcoeff_lj*nbparam.ewaldcoeff_lj;
    lje_coeff6_6 = lje_coeff2*lje_coeff2*lje_coeff2*ONE_SIXTH_F;
#endif /* LJ_EWALD */


#ifdef CALC_ENERGIES
    E_lj = 0.0f;
    E_el = 0.0f;

#if defined EXCLUSION_FORCES /* Ewald or RF */
    if (nb_sci.shift == CENTRAL && pl_cj4[cij4_start].cj[0] == sci*NCL_PER_SUPERCL)
    {
        /* we have the diagonal: add the charge and LJ self interaction energy term */
        for (i = 0; i < NCL_PER_SUPERCL; i++)
        {
#if defined EL_EWALD_ANY || defined EL_RF || defined EL_CUTOFF
            qi    = xqib[i * CL_SIZE + tidxi].w;
            E_el += qi*qi;
#endif

#if defined LJ_EWALD
#ifdef USE_TEXOBJ
            E_lj += tex1Dfetch<float>(nbparam.nbfp_texobj, atom_types[(sci*NCL_PER_SUPERCL + i)*CL_SIZE + tidxi]*(ntypes + 1)*2);
#else
            E_lj += tex1Dfetch(nbfp_texref, atom_types[(sci*NCL_PER_SUPERCL + i)*CL_SIZE + tidxi]*(ntypes + 1)*2);
#endif /* USE_TEXOBJ */
#endif /* LJ_EWALD */

        }

        /* divide the self term(s) equally over the j-threads, then multiply with the coefficients. */
#ifdef LJ_EWALD
        E_lj /= CL_SIZE*NTHREAD_Z;
        E_lj *= 0.5f*ONE_SIXTH_F*lje_coeff6_6;
#endif  /* LJ_EWALD */

#if defined EL_EWALD_ANY || defined EL_RF || defined EL_CUTOFF
        E_el /= CL_SIZE*NTHREAD_Z;
#if defined EL_RF || defined EL_CUTOFF
        E_el *= -nbparam.epsfac*0.5f*c_rf;
#else
        E_el *= -nbparam.epsfac*beta*M_FLOAT_1_SQRTPI; /* last factor 1/sqrt(pi) */
#endif
#endif                                                 /* EL_EWALD_ANY || defined EL_RF || defined EL_CUTOFF */
    }
#endif                                                 /* EXCLUSION_FORCES */

#endif                                                 /* CALC_ENERGIES */

    /* skip central shifts when summing shift forces */
    if (nb_sci.shift == CENTRAL)
    {
        bCalcFshift = false;
    }

    /* loop over the j clusters = seen by any of the atoms in the current super-cluster */
    for (j4 = cij4_start + tidxz; j4 < cij4_end; j4 += NTHREAD_Z)
    {
        wexcl_idx   = pl_cj4[j4].imei[widx].excl_ind;
        imask       = pl_cj4[j4].imei[widx].imask;
        wexcl       = excl[wexcl_idx].pair[(tidx) & (WARP_SIZE - 1)];

#ifndef PRUNE_NBL
        if (imask)
#endif
        {
            /* Pre-load cj into shared memory on both warps separately */
            if ((tidxj == 0 || tidxj == 4) && tidxi < NBNXN_GPU_JGROUP_SIZE)
            {
                cjs[tidxi + tidxj * NBNXN_GPU_JGROUP_SIZE / 4] = pl_cj4[j4].cj[tidxi];
            }

            /* Unrolling this loop
               - with pruning leads to register spilling;
               - on Kepler is much slower;
               - doesn't work on CUDA <v4.1
               Tested with nvcc 3.2 - 5.0.7 */
#if !defined PRUNE_NBL && __CUDA_ARCH__ < 300 && GMX_CUDA_VERSION >= 4010
#pragma unroll 4
#endif
            for (jm = 0; jm < NBNXN_GPU_JGROUP_SIZE; jm++)
            {
                /* ((1U << NCL_PER_SUPERCL) - 1U) is the i-cluster interaction
                 * mask for a super-cluster with all NCL_PER_SUPERCL bits set.
                 */
                if (imask & (((1U << NCL_PER_SUPERCL) - 1U) << (jm * NCL_PER_SUPERCL)))
                {
                    mask_ji = (1U << (jm * NCL_PER_SUPERCL));

                    cj      = cjs[jm + (tidxj & 4) * NBNXN_GPU_JGROUP_SIZE / 4];
                    aj      = cj * CL_SIZE + tidxj;

                    /* load j atom data */
                    xqbuf   = xq[aj];
                    xj      = make_float3(xqbuf.x, xqbuf.y, xqbuf.z);
                    qj_f    = nbparam.epsfac * xqbuf.w;
                    typej   = atom_types[aj];

                    fcj_buf = make_float3(0.0f);

                    /* The PME and RF kernels don't unroll with CUDA <v4.1. */
#if !defined PRUNE_NBL && !(GMX_CUDA_VERSION < 4010 && defined EXCLUSION_FORCES)
#pragma unroll 8
#endif
                    for (i = 0; i < NCL_PER_SUPERCL; i++)
                    {
                        if (imask & mask_ji)
                        {
                            ci_offset   = i;                     /* i force buffer offset */

                            ci      = sci * NCL_PER_SUPERCL + i; /* i cluster index */
                            ai      = ci * CL_SIZE + tidxi;      /* i atom index */

                            /* all threads load an atom from i cluster ci into shmem! */
                            xqbuf   = xqib[i * CL_SIZE + tidxi];
                            xi      = make_float3(xqbuf.x, xqbuf.y, xqbuf.z);

                            /* distance between i and j atoms */
                            rv      = xi - xj;
                            r2      = norm2(rv);

#ifdef PRUNE_NBL
                            /* If _none_ of the atoms pairs are in cutoff range,
                               the bit corresponding to the current
                               cluster-pair in imask gets set to 0. */
                            if (!__any(r2 < rlist_sq))
                            {
                                imask &= ~mask_ji;
                            }
#endif

                            int_bit = (wexcl & mask_ji) ? 1.0f : 0.0f;

                            /* cutoff & exclusion check */
#ifdef EXCLUSION_FORCES
                            if (r2 < rcoulomb_sq *
                                (nb_sci.shift != CENTRAL || ci != cj || tidxj > tidxi))
#else
                            if (r2 < rcoulomb_sq * int_bit)
#endif
                            {
                                /* load the rest of the i-atom parameters */
                                qi      = xqbuf.w;
#ifdef IATYPE_SHMEM
                                typei   = atib[i * CL_SIZE + tidxi];
#else
                                typei   = atom_types[ai];
#endif

                                /* LJ 6*C6 and 12*C12 */
#ifdef USE_TEXOBJ
                                c6      = tex1Dfetch<float>(nbparam.nbfp_texobj, 2 * (ntypes * typei + typej));
                                c12     = tex1Dfetch<float>(nbparam.nbfp_texobj, 2 * (ntypes * typei + typej) + 1);
#else
                                c6      = tex1Dfetch(nbfp_texref, 2 * (ntypes * typei + typej));
                                c12     = tex1Dfetch(nbfp_texref, 2 * (ntypes * typei + typej) + 1);
#endif                          /* USE_TEXOBJ */


                                /* avoid NaN for excluded pairs at r=0 */
                                r2      += (1.0f - int_bit) * NBNXN_AVOID_SING_R2_INC;

                                inv_r   = rsqrt(r2);
                                inv_r2  = inv_r * inv_r;
                                inv_r6  = inv_r2 * inv_r2 * inv_r2;
#if defined EXCLUSION_FORCES
                                /* We could mask inv_r2, but with Ewald
                                 * masking both inv_r6 and F_invr is faster */
                                inv_r6  *= int_bit;
#endif                          /* EXCLUSION_FORCES */

                                F_invr  = inv_r6 * (c12 * inv_r6 - c6) * inv_r2;
#if defined CALC_ENERGIES || defined LJ_POT_SWITCH
                                E_lj_p  = int_bit * (c12 * (inv_r6 * inv_r6 + nbparam.repulsion_shift.cpot)*ONE_TWELVETH_F -
                                                     c6 * (inv_r6 + nbparam.dispersion_shift.cpot)*ONE_SIXTH_F);
#endif

#ifdef LJ_FORCE_SWITCH
#ifdef CALC_ENERGIES
                                calculate_force_switch_F_E(nbparam, c6, c12, inv_r, r2, &F_invr, &E_lj_p);
#else
                                calculate_force_switch_F(nbparam, c6, c12, inv_r, r2, &F_invr);
#endif /* CALC_ENERGIES */
#endif /* LJ_FORCE_SWITCH */


#ifdef LJ_EWALD
#ifdef LJ_EWALD_COMB_GEOM
#ifdef CALC_ENERGIES
                                calculate_lj_ewald_comb_geom_F_E(nbparam, typei, typej, r2, inv_r2, lje_coeff2, lje_coeff6_6, int_bit, &F_invr, &E_lj_p);
#else
                                calculate_lj_ewald_comb_geom_F(nbparam, typei, typej, r2, inv_r2, lje_coeff2, lje_coeff6_6, &F_invr);
#endif                          /* CALC_ENERGIES */
#elif defined LJ_EWALD_COMB_LB
                                calculate_lj_ewald_comb_LB_F_E(nbparam, typei, typej, r2, inv_r2, lje_coeff2, lje_coeff6_6,
#ifdef CALC_ENERGIES
                                                               int_bit, &F_invr, &E_lj_p
#else
                                                               0, &F_invr, NULL
#endif /* CALC_ENERGIES */
                                                               );
#endif /* LJ_EWALD_COMB_GEOM */
#endif /* LJ_EWALD */

#ifdef VDW_CUTOFF_CHECK
                                /* Separate VDW cut-off check to enable twin-range cut-offs
                                 * (rvdw < rcoulomb <= rlist)
                                 */
                                vdw_in_range  = (r2 < rvdw_sq) ? 1.0f : 0.0f;
                                F_invr       *= vdw_in_range;
#ifdef CALC_ENERGIES
                                E_lj_p       *= vdw_in_range;
#endif
#endif                          /* VDW_CUTOFF_CHECK */

#ifdef LJ_POT_SWITCH
#ifdef CALC_ENERGIES
                                calculate_potential_switch_F_E(nbparam, c6, c12, inv_r, r2, &F_invr, &E_lj_p);
#else
                                calculate_potential_switch_F(nbparam, c6, c12, inv_r, r2, &F_invr, &E_lj_p);
#endif /* CALC_ENERGIES */
#endif /* LJ_POT_SWITCH */

#ifdef CALC_ENERGIES
                                E_lj    += E_lj_p;
#endif


#ifdef EL_CUTOFF
#ifdef EXCLUSION_FORCES
                                F_invr  += qi * qj_f * int_bit * inv_r2 * inv_r;
#else
                                F_invr  += qi * qj_f * inv_r2 * inv_r;
#endif
#endif
#ifdef EL_RF
                                F_invr  += qi * qj_f * (int_bit*inv_r2 * inv_r - two_k_rf);
#endif
#if defined EL_EWALD_ANA
                                F_invr  += qi * qj_f * (int_bit*inv_r2*inv_r + pmecorrF(beta2*r2)*beta3);
#elif defined EL_EWALD_TAB
                                F_invr  += qi * qj_f * (int_bit*inv_r2 -
#ifdef USE_TEXOBJ
                                                        interpolate_coulomb_force_r(nbparam.coulomb_tab_texobj, r2 * inv_r, coulomb_tab_scale)
#else
                                                        interpolate_coulomb_force_r(r2 * inv_r, coulomb_tab_scale)
#endif /* USE_TEXOBJ */
                                                        ) * inv_r;
#endif /* EL_EWALD_ANA/TAB */

#ifdef CALC_ENERGIES
#ifdef EL_CUTOFF
                                E_el    += qi * qj_f * (int_bit*inv_r - c_rf);
#endif
#ifdef EL_RF
                                E_el    += qi * qj_f * (int_bit*inv_r + 0.5f * two_k_rf * r2 - c_rf);
#endif
#ifdef EL_EWALD_ANY
                                /* 1.0f - erff is faster than erfcf */
                                E_el    += qi * qj_f * (inv_r * (int_bit - erff(r2 * inv_r * beta)) - int_bit * ewald_shift);
#endif                          /* EL_EWALD_ANY */
#endif
                                f_ij    = rv * F_invr;

                                /* accumulate j forces in registers */
                                fcj_buf -= f_ij;

                                /* accumulate i forces in registers */
                                fci_buf[ci_offset] += f_ij;
                            }
                        }

                        /* shift the mask bit by 1 */
                        mask_ji += mask_ji;
                    }

                    /* reduce j forces */
#ifdef REDUCE_SHUFFLE
                    reduce_force_j_warp_shfl(fcj_buf, f, tidxi, aj);
#else
                    /* store j forces in shmem */
                    f_buf[                  tidx] = fcj_buf.x;
                    f_buf[    FBUF_STRIDE + tidx] = fcj_buf.y;
                    f_buf[2 * FBUF_STRIDE + tidx] = fcj_buf.z;

                    reduce_force_j_generic(f_buf, f, tidxi, tidxj, aj);
#endif
                }
            }
#ifdef PRUNE_NBL
            /* Update the imask with the new one which does not contain the
               out of range clusters anymore. */
            pl_cj4[j4].imei[widx].imask = imask;
#endif
        }
    }

    float fshift_buf = 0.0f;

    /* reduce i forces */
    for (ci_offset = 0; ci_offset < NCL_PER_SUPERCL; ci_offset++)
    {
        ai  = (sci * NCL_PER_SUPERCL + ci_offset) * CL_SIZE + tidxi;
#ifdef REDUCE_SHUFFLE
        reduce_force_i_warp_shfl(fci_buf[ci_offset], f,
                                 &fshift_buf, bCalcFshift,
                                 tidxj, ai);
#else
        f_buf[                  tidx] = fci_buf[ci_offset].x;
        f_buf[    FBUF_STRIDE + tidx] = fci_buf[ci_offset].y;
        f_buf[2 * FBUF_STRIDE + tidx] = fci_buf[ci_offset].z;
        __syncthreads();
        reduce_force_i(f_buf, f,
                       &fshift_buf, bCalcFshift,
                       tidxi, tidxj, ai);
        __syncthreads();
#endif
    }

    /* add up local shift forces into global mem, tidxj indexes x,y,z */
#ifdef REDUCE_SHUFFLE
    if (bCalcFshift && (tidxj & 3) < 3)
    {
        atomicAdd(&(atdat.fshift[nb_sci.shift].x) + (tidxj & ~4), fshift_buf);
    }
#else
    if (bCalcFshift && tidxj < 3)
    {
        atomicAdd(&(atdat.fshift[nb_sci.shift].x) + tidxj, fshift_buf);
    }
#endif

#ifdef CALC_ENERGIES
#ifdef REDUCE_SHUFFLE
    /* reduce the energies over warps and store into global memory */
    reduce_energy_warp_shfl(E_lj, E_el, e_lj, e_el, tidx);
#else
    /* flush the energies to shmem and reduce them */
    f_buf[              tidx] = E_lj;
    f_buf[FBUF_STRIDE + tidx] = E_el;
    reduce_energy_pow2(f_buf + (tidx & WARP_SIZE), e_lj, e_el, tidx & ~WARP_SIZE);
#endif
#endif
}

#undef REDUCE_SHUFFLE
#undef IATYPE_SHMEM
#undef USE_TEXOBJ

#undef NTHREAD_Z
#undef MIN_BLOCKS_PER_MP
#undef THREADS_PER_BLOCK

#undef EL_EWALD_ANY
#undef EXCLUSION_FORCES
#undef LJ_EWALD
