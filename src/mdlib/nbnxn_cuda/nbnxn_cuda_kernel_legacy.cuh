/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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

#include "maths.h"
/* Note that floating-point constants in CUDA code should be suffixed
 * with f (e.g. 0.5f), to stop the compiler producing intermediate
 * code that is in double precision.
 */

/* Analytical Ewald is not implemented for the legacy kernels (as it is anyway
   slower than the tabulated kernel on Fermi). */
#ifdef EL_EWALD_ANA
#error Trying to generate Analytical Ewald legacy kernels which is not implemented in the legacy kernels!
#endif

/*
   Kernel launch parameters:
    - #blocks   = #pair lists, blockId = pair list Id
    - #threads  = CL_SIZE^2
    - shmem     = CL_SIZE^2 * sizeof(float)

    Each thread calculates an i force-component taking one pair of i-j atoms.
 */
#if __CUDA_ARCH__ >= 350
__launch_bounds__(64,16)
#endif
#ifdef PRUNE_NBL
#ifdef CALC_ENERGIES
__global__ void NB_KERNEL_FUNC_NAME(k_nbnxn, _ener_prune_legacy)
#else
__global__ void NB_KERNEL_FUNC_NAME(k_nbnxn, _prune_legacy)
#endif
#else
#ifdef CALC_ENERGIES
__global__ void NB_KERNEL_FUNC_NAME(k_nbnxn, _ener_legacy)
#else
__global__ void NB_KERNEL_FUNC_NAME(k_nbnxn, _legacy)
#endif
#endif
            (const cu_atomdata_t atdat,
             const cu_nbparam_t nbparam,
             const cu_plist_t plist,
             bool bCalcFshift)
{
    /* convenience variables */
    const nbnxn_sci_t *pl_sci   = plist.sci;
#ifndef PRUNE_NBL
    const
#endif
    nbnxn_cj4_t *pl_cj4         = plist.cj4;
    const nbnxn_excl_t *excl    = plist.excl;
    const int *atom_types       = atdat.atom_types;
    int ntypes                  = atdat.ntypes;
    const float4 *xq            = atdat.xq;
    float3 *f                   = atdat.f;
    const float3 *shift_vec     = atdat.shift_vec;
    float rcoulomb_sq           = nbparam.rcoulomb_sq;
#ifdef VDW_CUTOFF_CHECK
    float rvdw_sq               = nbparam.rvdw_sq;
    float vdw_in_range;
#endif
#ifdef EL_RF
    float two_k_rf              = nbparam.two_k_rf;
#endif
#ifdef EL_EWALD_TAB
    float coulomb_tab_scale     = nbparam.coulomb_tab_scale;
#endif
#ifdef PRUNE_NBL
    float rlist_sq              = nbparam.rlist_sq;
#endif

#ifdef CALC_ENERGIES
    float lj_shift    = nbparam.sh_invrc6;
#ifdef EL_EWALD_TAB
    float beta        = nbparam.ewald_beta;
    float ewald_shift = nbparam.sh_ewald;
#else
    float c_rf        = nbparam.c_rf;
#endif
    float *e_lj       = atdat.e_lj;
    float *e_el       = atdat.e_el;
#endif

    /* thread/block/warp id-s */
    unsigned int tidxi  = threadIdx.x;
    unsigned int tidxj  = threadIdx.y;
    unsigned int tidx   = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int bidx   = blockIdx.x;
    unsigned int widx   = tidx / WARP_SIZE; /* warp index */

    int sci, ci, cj, ci_offset,
        ai, aj,
        cij4_start, cij4_end,
        typei, typej,
        i, cii, jm, j4, nsubi, wexcl_idx;
    float qi, qj_f,
          r2, inv_r, inv_r2, inv_r6,
          c6, c12,
          int_bit,
#ifdef CALC_ENERGIES
          E_lj, E_el, E_lj_p,
#endif
          F_invr;
    unsigned int wexcl, imask, mask_ji;
    float4 xqbuf;
    float3 xi, xj, rv, f_ij, fcj_buf, fshift_buf;
    float3 fci_buf[NCL_PER_SUPERCL];    /* i force buffer */
    nbnxn_sci_t nb_sci;

    /* shmem buffer for i x+q pre-loading */
    extern __shared__  float4 xqib[];
    /* shmem j force buffer */
    float *f_buf = (float *)(xqib + NCL_PER_SUPERCL * CL_SIZE);

    nb_sci      = pl_sci[bidx];         /* my i super-cluster's index = current bidx */
    sci         = nb_sci.sci;           /* super-cluster */
    cij4_start  = nb_sci.cj4_ind_start; /* first ...*/
    cij4_end    = nb_sci.cj4_ind_end;   /* and last index of j clusters */

    /* Store the i-atom x and q in shared memory */
    /* Note: the thread indexing here is inverted with respect to the
       inner-loop as this results in slightly higher performance */
    ci = sci * NCL_PER_SUPERCL + tidxi;
    ai = ci * CL_SIZE + tidxj;
    xqib[tidxi * CL_SIZE + tidxj] = xq[ai] + shift_vec[nb_sci.shift];
    __syncthreads();

    for(ci_offset = 0; ci_offset < NCL_PER_SUPERCL; ci_offset++)
    {
        fci_buf[ci_offset] = make_float3(0.0f);
    }

#ifdef CALC_ENERGIES
    E_lj = 0.0f;
    E_el = 0.0f;

#if defined EL_EWALD_TAB || defined EL_RF
    if (nb_sci.shift == CENTRAL && pl_cj4[cij4_start].cj[0] == sci*NCL_PER_SUPERCL)
    {
        /* we have the diagonal: add the charge self interaction energy term */
        for (i = 0; i < NCL_PER_SUPERCL; i++)
        {
            qi    = xqib[i * CL_SIZE + tidxi].w;
            E_el += qi*qi;
        }
        /* divide the self term equally over the j-threads */
        E_el /= CL_SIZE;
#ifdef EL_RF
        E_el *= -nbparam.epsfac*0.5f*c_rf;
#else
        E_el *= -nbparam.epsfac*beta*M_FLOAT_1_SQRTPI; /* last factor 1/sqrt(pi) */
#endif
    }
#endif
#endif

    /* skip central shifts when summing shift forces */
    if (nb_sci.shift == CENTRAL)
    {
        bCalcFshift = false;
    }

    fshift_buf = make_float3(0.0f);

    /* loop over the j clusters = seen by any of the atoms in the current super-cluster */
    for (j4 = cij4_start; j4 < cij4_end; j4++)
    {
        wexcl_idx   = pl_cj4[j4].imei[widx].excl_ind;
        imask       = pl_cj4[j4].imei[widx].imask;
        wexcl       = excl[wexcl_idx].pair[(tidx) & (WARP_SIZE - 1)];

#ifndef PRUNE_NBL
        if (imask)
#endif
        {
            /* nvcc >v4.1 doesn't like this loop, it refuses to unroll it */
#if CUDA_VERSION >= 4010
            #pragma unroll 4
#endif
            for (jm = 0; jm < NBNXN_GPU_JGROUP_SIZE; jm++)
            {
                mask_ji = (imask >> (jm * CL_SIZE)) & supercl_interaction_mask;
                if (mask_ji)
                {
                    nsubi = __popc(mask_ji);

                    cj      = pl_cj4[j4].cj[jm];
                    aj      = cj * CL_SIZE + tidxj;

                    /* load j atom data */
                    xqbuf   = xq[aj];
                    xj      = make_float3(xqbuf.x, xqbuf.y, xqbuf.z);
                    qj_f    = nbparam.epsfac * xqbuf.w;
                    typej   = atom_types[aj];

                    fcj_buf = make_float3(0.0f);

                    /* loop over the i-clusters in sci */
                    /* #pragma unroll 8
                       -- nvcc doesn't like my code, it refuses to unroll it
                       which is a pity because here unrolling could help.  */
                    for (cii = 0; cii < nsubi; cii++)
                    {
                        i = __ffs(mask_ji) - 1;
                        mask_ji &= ~(1U << i);

                        ci_offset   = i;    /* i force buffer offset */

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
                            imask &= ~(1U << (jm * NCL_PER_SUPERCL + i));
                        }
#endif

                        int_bit = ((wexcl >> (jm * NCL_PER_SUPERCL + i)) & 1) ? 1.0f : 0.0f;

                        /* cutoff & exclusion check */
#if defined EL_EWALD_TAB || defined EL_RF
                        if (r2 < rcoulomb_sq *
                            (nb_sci.shift != CENTRAL || ci != cj || tidxj > tidxi))
#else
                        if (r2 < rcoulomb_sq * int_bit)
#endif
                        {
                            /* load the rest of the i-atom parameters */
                            qi      = xqbuf.w;
                            typei   = atom_types[ai];

                            /* LJ 6*C6 and 12*C12 */
                            c6      = tex1Dfetch(nbfp_texref, 2 * (ntypes * typei + typej));
                            c12     = tex1Dfetch(nbfp_texref, 2 * (ntypes * typei + typej) + 1);

                            /* avoid NaN for excluded pairs at r=0 */
                            r2      += (1.0f - int_bit) * NBNXN_AVOID_SING_R2_INC;

                            inv_r   = rsqrt(r2);
                            inv_r2  = inv_r * inv_r;
                            inv_r6  = inv_r2 * inv_r2 * inv_r2;
#if defined EL_EWALD_TAB || defined EL_RF
                            /* We could mask inv_r2, but with Ewald
                             * masking both inv_r6 and F_invr is faster */
                            inv_r6  *= int_bit;
#endif

                            F_invr  = inv_r6 * (c12 * inv_r6 - c6) * inv_r2;

#ifdef CALC_ENERGIES
                            E_lj_p  = int_bit * (c12 * (inv_r6 * inv_r6 - lj_shift * lj_shift) * 0.08333333f - c6 * (inv_r6 - lj_shift) * 0.16666667f);
#endif

#ifdef VDW_CUTOFF_CHECK
                            /* this enables twin-range cut-offs (rvdw < rcoulomb <= rlist) */
                            vdw_in_range = (r2 < rvdw_sq) ? 1.0f : 0.0f;
                            F_invr  *= vdw_in_range;
#ifdef CALC_ENERGIES
                            E_lj_p  *= vdw_in_range;
#endif
#endif
#ifdef CALC_ENERGIES
                            E_lj    += E_lj_p;
#endif


#ifdef EL_CUTOFF
                            F_invr  += qi * qj_f * inv_r2 * inv_r;
#endif
#ifdef EL_RF
                            F_invr  += qi * qj_f * (int_bit*inv_r2 * inv_r - two_k_rf);
#endif
#ifdef EL_EWALD_TAB
                            F_invr  += qi * qj_f * (int_bit*inv_r2 - interpolate_coulomb_force_r(r2 * inv_r, coulomb_tab_scale)) * inv_r;
#endif /* EL_EWALD_TAB */

#ifdef CALC_ENERGIES
#ifdef EL_CUTOFF
                            E_el    += qi * qj_f * (inv_r - c_rf);
#endif
#ifdef EL_RF
                            E_el    += qi * qj_f * (int_bit*inv_r + 0.5f * two_k_rf * r2 - c_rf);
#endif
#ifdef EL_EWALD_TAB
                            /* 1.0f - erff is faster than erfcf */
                            E_el    += qi * qj_f * (inv_r * (int_bit - erff(r2 * inv_r * beta)) - int_bit * ewald_shift);
#endif
#endif
                            f_ij    = rv * F_invr;

                            /* accumulate j forces in registers */
                            fcj_buf -= f_ij;

                            /* accumulate i forces in registers */
                            fci_buf[ci_offset] += f_ij;
                        }
                    }

                    /* store j forces in shmem */
                    f_buf[                  tidx] = fcj_buf.x;
                    f_buf[    FBUF_STRIDE + tidx] = fcj_buf.y;
                    f_buf[2 * FBUF_STRIDE + tidx] = fcj_buf.z;

                    /* reduce j forces */
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

    /* reduce i forces */
    for(ci_offset = 0; ci_offset < NCL_PER_SUPERCL; ci_offset++)
    {
        ai  = (sci * NCL_PER_SUPERCL + ci_offset) * CL_SIZE + tidxi;
        f_buf[                  tidx] = fci_buf[ci_offset].x;
        f_buf[    FBUF_STRIDE + tidx] = fci_buf[ci_offset].y;
        f_buf[2 * FBUF_STRIDE + tidx] = fci_buf[ci_offset].z;
        __syncthreads();
        reduce_force_i(f_buf, f,
                       &fshift_buf, bCalcFshift,
                       tidxi, tidxj, ai);
        __syncthreads();
    }

    /* add up local shift forces into global mem */
    if (bCalcFshift && tidxj == 0)
    {
        atomicAdd(&atdat.fshift[nb_sci.shift].x, fshift_buf.x);
        atomicAdd(&atdat.fshift[nb_sci.shift].y, fshift_buf.y);
        atomicAdd(&atdat.fshift[nb_sci.shift].z, fshift_buf.z);
    }

#ifdef CALC_ENERGIES
    /* flush the energies to shmem and reduce them */
    f_buf[              tidx] = E_lj;
    f_buf[FBUF_STRIDE + tidx] = E_el;
    reduce_energy_pow2(f_buf + (tidx & WARP_SIZE), e_lj, e_el, tidx & ~WARP_SIZE);
#endif
}
