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

/*! \internal \file
 *  \brief
 *  CUDA non-bonded kernel used through preprocessor-based code generation
 *  of multiple kernel flavors for CC 2.x, see nbnxn_cuda_kernels.cuh.
 *
 *  NOTE: No include fence as it is meant to be included multiple times.
 *
 *  \author Szilárd Páll <pall.szilard@gmail.com>
 *  \author Berk Hess <hess@kth.se>
 *  \ingroup module_mdlib
 */

#include "gromacs/gpu_utils/cuda_arch_utils.cuh"
#include "gromacs/gpu_utils/cuda_kernel_utils.cuh"
#include "gromacs/math/utilities.h"
#include "gromacs/pbcutil/ishift.h"
/* Note that floating-point constants in CUDA code should be suffixed
 * with f (e.g. 0.5f), to stop the compiler producing intermediate
 * code that is in double precision.
 */

#if GMX_PTX_ARCH >= 300
#error "nbnxn_cuda_kernel_fermi.cuh included with GMX_PTX_ARCH >= 300"
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

#if defined LJ_COMB_GEOM || defined LJ_COMB_LB
#define LJ_COMB
#endif

/*
   Kernel launch parameters:
    - #blocks   = #pair lists, blockId = pair list Id
    - #threads  = c_clSize^2
    - shmem     = see nbnxn_cuda.cu:calc_shmem_required_nonbonded()

    Each thread calculates an i force-component taking one pair of i-j atoms.
 */

/**@{*/
/*! \brief Definition of kernel launch configuration parameters for CC 2.x.
 */

/* Kernel launch bounds, 16 blocks/multiprocessor can be kept in flight. */
#define THREADS_PER_BLOCK   (c_clSize*c_clSize)

__launch_bounds__(THREADS_PER_BLOCK)
#ifdef PRUNE_NBL
#ifdef CALC_ENERGIES
__global__ void NB_KERNEL_FUNC_NAME(nbnxn_kernel, _VF_prune_cuda)
#else
__global__ void NB_KERNEL_FUNC_NAME(nbnxn_kernel, _F_prune_cuda)
#endif /* CALC_ENERGIES */
#else
#ifdef CALC_ENERGIES
__global__ void NB_KERNEL_FUNC_NAME(nbnxn_kernel, _VF_cuda)
#else
__global__ void NB_KERNEL_FUNC_NAME(nbnxn_kernel, _F_cuda)
#endif /* CALC_ENERGIES */
#endif /* PRUNE_NBL */
(const cu_atomdata_t atdat,
 const cu_nbparam_t nbparam,
 const cu_plist_t plist,
 bool bCalcFshift)
#ifdef FUNCTION_DECLARATION_ONLY
;     /* Only do function declaration, omit the function body. */
#else
{
    /* convenience variables */
    const nbnxn_sci_t *pl_sci       = plist.sci;
#ifndef PRUNE_NBL
    const
#endif
    nbnxn_cj4_t        *pl_cj4      = plist.cj4;
    const nbnxn_excl_t *excl        = plist.excl;
#ifndef LJ_COMB
    const int          *atom_types  = atdat.atom_types;
    int                 ntypes      = atdat.ntypes;
#else
    const float2       *lj_comb     = atdat.lj_comb;
    float2              ljcp_i, ljcp_j;
#endif
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
#ifdef EL_EWALD_ANA
    float beta2                 = nbparam.ewald_beta*nbparam.ewald_beta;
    float beta3                 = nbparam.ewald_beta*nbparam.ewald_beta*nbparam.ewald_beta;
#endif
#ifdef PRUNE_NBL
    float rlist_sq              = nbparam.rlistOuter_sq;
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
    unsigned int bidx   = blockIdx.x;
    unsigned int widx   = tidx / warp_size; /* warp index */

    int          sci, ci, cj,
                 ai, aj,
                 cij4_start, cij4_end;
#ifndef LJ_COMB
    int          typei, typej;
#endif
    int          i, jm, j4, wexcl_idx;
    float        qi, qj_f,
                 r2, inv_r, inv_r2;
#if !defined LJ_COMB_LB || defined CALC_ENERGIES
    float        inv_r6, c6, c12;
#endif
#ifdef LJ_COMB_LB
    float        sigma, epsilon;
#endif
    float        int_bit,
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
    float3       fci_buf[c_numClPerSupercl]; /* i force buffer */
    nbnxn_sci_t  nb_sci;

    /*! i-cluster interaction mask for a super-cluster with all c_numClPerSupercl=8 bits set */
    const unsigned superClInteractionMask = ((1U << c_numClPerSupercl) - 1U);

    /*********************************************************************
     * Set up shared memory pointers.
     * sm_nextSlotPtr should always be updated to point to the "next slot",
     * that is past the last point where data has been stored.
     */
    extern __shared__  char sm_dynamicShmem[];
    char                   *sm_nextSlotPtr = sm_dynamicShmem;
    static_assert(sizeof(char) == 1, "The shared memory offset calculation assumes that char is 1 byte");

    /* shmem buffer for i x+q pre-loading */
    float4 *xqib    = (float4 *)sm_nextSlotPtr;
    sm_nextSlotPtr += (c_numClPerSupercl * c_clSize * sizeof(*xqib));

    /* shmem buffer for cj, for each warp separately */
    int *cjs        = (int *)(sm_nextSlotPtr);
    sm_nextSlotPtr += (c_nbnxnGpuClusterpairSplit * c_nbnxnGpuJgroupSize * sizeof(*cjs));

    /* shmem j force buffer */
    float *f_buf    = (float *)(sm_nextSlotPtr);
    sm_nextSlotPtr += (c_clSize * c_clSize * 3*sizeof(*f_buf));
    /*********************************************************************/

    nb_sci      = pl_sci[bidx];         /* my i super-cluster's index = current bidx */
    sci         = nb_sci.sci;           /* super-cluster */
    cij4_start  = nb_sci.cj4_ind_start; /* first ...*/
    cij4_end    = nb_sci.cj4_ind_end;   /* and last index of j clusters */

    {
        /* Pre-load i-atom x and q into shared memory */
        ci = sci * c_numClPerSupercl + tidxj;
        ai = ci * c_clSize + tidxi;

        xqbuf    = xq[ai] + shift_vec[nb_sci.shift];
        xqbuf.w *= nbparam.epsfac;
        xqib[tidxj * c_clSize + tidxi] = xqbuf;
    }
    __syncthreads();

    for (i = 0; i < c_numClPerSupercl; i++)
    {
        fci_buf[i] = make_float3(0.0f);
    }

#ifdef LJ_EWALD
    /* TODO: we are trading registers with flops by keeping lje_coeff-s, try re-calculating it later */
    lje_coeff2   = nbparam.ewaldcoeff_lj*nbparam.ewaldcoeff_lj;
    lje_coeff6_6 = lje_coeff2*lje_coeff2*lje_coeff2*c_oneSixth;
#endif


#ifdef CALC_ENERGIES
    E_lj = 0.0f;
    E_el = 0.0f;

#ifdef EXCLUSION_FORCES /* Ewald or RF */
    if (nb_sci.shift == CENTRAL && pl_cj4[cij4_start].cj[0] == sci*c_numClPerSupercl)
    {
        /* we have the diagonal: add the charge and LJ self interaction energy term */
        for (i = 0; i < c_numClPerSupercl; i++)
        {
#if defined EL_EWALD_ANY || defined EL_RF || defined EL_CUTOFF
            qi    = xqib[i * c_clSize + tidxi].w;
            E_el += qi*qi;
#endif

#ifdef LJ_EWALD
    #if DISABLE_CUDA_TEXTURES
            E_lj += LDG(&nbparam.nbfp[atom_types[(sci*c_numClPerSupercl + i)*c_clSize + tidxi]*(ntypes + 1)*2]);
    #else
            E_lj += tex1Dfetch(nbfp_texref, atom_types[(sci*c_numClPerSupercl + i)*c_clSize + tidxi]*(ntypes + 1)*2);
    #endif
#endif
        }

        /* divide the self term(s) equally over the j-threads, then multiply with the coefficients. */
#ifdef LJ_EWALD
        E_lj /= c_clSize;
        E_lj *= 0.5f*c_oneSixth*lje_coeff6_6;
#endif

#if defined EL_EWALD_ANY || defined EL_RF || defined EL_CUTOFF
        /* Correct for epsfac^2 due to adding qi^2 */
        E_el /= nbparam.epsfac*c_clSize;
#if defined EL_RF || defined EL_CUTOFF
        E_el *= -0.5f*c_rf;
#else
        E_el *= -beta*M_FLOAT_1_SQRTPI; /* last factor 1/sqrt(pi) */
#endif
#endif                                  /* EL_EWALD_ANY || defined EL_RF || defined EL_CUTOFF */
    }
#endif                                  /* EXCLUSION_FORCES */

#endif                                  /* CALC_ENERGIES */

#ifdef EXCLUSION_FORCES
    const int nonSelfInteraction = !(nb_sci.shift == CENTRAL & tidxj <= tidxi);
#endif

    /* loop over the j clusters = seen by any of the atoms in the current super-cluster */
    for (j4 = cij4_start; j4 < cij4_end; j4++)
    {
        wexcl_idx   = pl_cj4[j4].imei[widx].excl_ind;
        imask       = pl_cj4[j4].imei[widx].imask;
        wexcl       = excl[wexcl_idx].pair[(tidx) & (warp_size - 1)];

#ifndef PRUNE_NBL
        if (imask)
#endif
        {
            /* Pre-load cj into shared memory on both warps separately */
            if ((tidxj == 0 | tidxj == 4) & (tidxi < c_nbnxnGpuJgroupSize))
            {
                cjs[tidxi + tidxj * c_nbnxnGpuJgroupSize/c_splitClSize] = pl_cj4[j4].cj[tidxi];
            }

            /* Unrolling this loop with pruning leads to register spilling;
               Tested with up to nvcc 7.5 */
#if !defined PRUNE_NBL
#pragma unroll 4
#endif
            for (jm = 0; jm < c_nbnxnGpuJgroupSize; jm++)
            {
                if (imask & (superClInteractionMask << (jm * c_numClPerSupercl)))
                {
                    mask_ji = (1U << (jm * c_numClPerSupercl));

                    cj      = cjs[jm + (tidxj & 4) * c_nbnxnGpuJgroupSize/c_splitClSize];
                    aj      = cj * c_clSize + tidxj;

                    /* load j atom data */
                    xqbuf   = xq[aj];
                    xj      = make_float3(xqbuf.x, xqbuf.y, xqbuf.z);
                    qj_f    = xqbuf.w;
#ifndef LJ_COMB
                    typej   = atom_types[aj];
#else
                    ljcp_j  = lj_comb[aj];
#endif

                    fcj_buf = make_float3(0.0f);

#if !defined PRUNE_NBL
#pragma unroll 8
#endif
                    for (i = 0; i < c_numClPerSupercl; i++)
                    {
                        if (imask & mask_ji)
                        {
                            ci      = sci * c_numClPerSupercl + i; /* i cluster index */
                            ai      = ci * c_clSize + tidxi;       /* i atom index */

                            /* all threads load an atom from i cluster ci into shmem! */
                            xqbuf   = xqib[i * c_clSize + tidxi];
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
                            if ((r2 < rcoulomb_sq) * (nonSelfInteraction | (ci != cj)))
#else
                            if ((r2 < rcoulomb_sq) * int_bit)
#endif
                            {
                                /* load the rest of the i-atom parameters */
                                qi      = xqbuf.w;

#ifndef LJ_COMB
                                /* LJ 6*C6 and 12*C12 */
                                typei   = atom_types[ai];
                                fetch_nbfp_c6_c12(c6, c12, nbparam, ntypes * typei + typej);
#else
                                ljcp_i  = lj_comb[ai];
#ifdef LJ_COMB_GEOM
                                c6      = ljcp_i.x * ljcp_j.x;
                                c12     = ljcp_i.y * ljcp_j.y;
#else
                                /* LJ 2^(1/6)*sigma and 12*epsilon */
                                sigma   = ljcp_i.x + ljcp_j.x;
                                epsilon = ljcp_i.y * ljcp_j.y;
#if defined CALC_ENERGIES || defined LJ_FORCE_SWITCH || defined LJ_POT_SWITCH
                                convert_sigma_epsilon_to_c6_c12(sigma, epsilon, &c6, &c12);
#endif
#endif                          /* LJ_COMB_GEOM */
#endif                          /* LJ_COMB */

                                // Ensure distance do not become so small that r^-12 overflows
                                r2      = max(r2, NBNXN_MIN_RSQ);

                                inv_r   = rsqrt(r2);
                                inv_r2  = inv_r * inv_r;
#if !defined LJ_COMB_LB || defined CALC_ENERGIES
                                inv_r6  = inv_r2 * inv_r2 * inv_r2;
#ifdef EXCLUSION_FORCES
                                /* We could mask inv_r2, but with Ewald
                                 * masking both inv_r6 and F_invr is faster */
                                inv_r6  *= int_bit;
#endif                          /* EXCLUSION_FORCES */

                                F_invr  = inv_r6 * (c12 * inv_r6 - c6) * inv_r2;
#if defined CALC_ENERGIES || defined LJ_POT_SWITCH
                                E_lj_p  = int_bit * (c12 * (inv_r6 * inv_r6 + nbparam.repulsion_shift.cpot)*c_oneTwelveth -
                                                     c6 * (inv_r6 + nbparam.dispersion_shift.cpot)*c_oneSixth);
#endif
#else                           /* !LJ_COMB_LB || CALC_ENERGIES */
                                float sig_r  = sigma*inv_r;
                                float sig_r2 = sig_r*sig_r;
                                float sig_r6 = sig_r2*sig_r2*sig_r2;
#ifdef EXCLUSION_FORCES
                                sig_r6 *= int_bit;
#endif                          /* EXCLUSION_FORCES */

                                F_invr  = epsilon * sig_r6 * (sig_r6 - 1.0f) * inv_r2;
#endif                          /* !LJ_COMB_LB || CALC_ENERGIES */

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

#ifdef LJ_POT_SWITCH
#ifdef CALC_ENERGIES
                                calculate_potential_switch_F_E(nbparam, inv_r, r2, &F_invr, &E_lj_p);
#else
                                calculate_potential_switch_F(nbparam, inv_r, r2, &F_invr, &E_lj_p);
#endif /* CALC_ENERGIES */
#endif /* LJ_POT_SWITCH */

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
                                                        interpolate_coulomb_force_r(nbparam, r2 * inv_r)) * inv_r;
#endif                          /* EL_EWALD_ANA/TAB */

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
                                fci_buf[i] += f_ij;
                            }
                        }

                        /* shift the mask bit by 1 */
                        mask_ji += mask_ji;
                    }

                    /* reduce j forces */
                    /* store j forces in shmem */
                    f_buf[                   tidx] = fcj_buf.x;
                    f_buf[    c_fbufStride + tidx] = fcj_buf.y;
                    f_buf[2 * c_fbufStride + tidx] = fcj_buf.z;

                    reduce_force_j_generic(f_buf, f, tidxi, tidxj, aj);
                }
            }
#ifdef PRUNE_NBL
            /* Update the imask with the new one which does not contain the
               out of range clusters anymore. */
            pl_cj4[j4].imei[widx].imask = imask;
#endif
        }
    }

    /* skip central shifts when summing shift forces */
    if (nb_sci.shift == CENTRAL)
    {
        bCalcFshift = false;
    }

    float fshift_buf = 0.0f;

    /* reduce i forces */
    for (i = 0; i < c_numClPerSupercl; i++)
    {
        ai  = (sci * c_numClPerSupercl + i) * c_clSize + tidxi;
        f_buf[                   tidx] = fci_buf[i].x;
        f_buf[    c_fbufStride + tidx] = fci_buf[i].y;
        f_buf[2 * c_fbufStride + tidx] = fci_buf[i].z;
        __syncthreads();
        reduce_force_i(f_buf, f,
                       &fshift_buf, bCalcFshift,
                       tidxi, tidxj, ai);
        __syncthreads();
    }

    /* add up local shift forces into global mem, tidxj indexes x,y,z */
    if (bCalcFshift && tidxj < 3)
    {
        atomicAdd(&(atdat.fshift[nb_sci.shift].x) + tidxj, fshift_buf);
    }

#ifdef CALC_ENERGIES
    /* flush the energies to shmem and reduce them */
    f_buf[               tidx] = E_lj;
    f_buf[c_fbufStride + tidx] = E_el;
    reduce_energy_pow2(f_buf + (tidx & warp_size), e_lj, e_el, tidx & ~warp_size);
#endif
}
#endif /* FUNCTION_DECLARATION_ONLY */

#undef THREADS_PER_BLOCK

#undef EL_EWALD_ANY
#undef EXCLUSION_FORCES
#undef LJ_EWALD

#undef LJ_COMB
