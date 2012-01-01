/*  Launch parameterts:
    - #blocks   = #neighbor lists, blockId = neigbor_listId
    - #threads  = CLUSTER_SIZE^2
    - shmem     = CLUSTER_SIZE^2 * sizeof(float)
    - local mem = 4 bytes !!! 

    Each thread calculates an i force-component taking one pair of i-j atoms.
 */
#ifdef PRUNE_NBL
#ifdef CALC_ENERGIES
__global__ void FUNCTION_NAME(k_nbnxn, ener_prune_2)
#else
__global__ void FUNCTION_NAME(k_nbnxn, prune_2)
#endif
#else
#ifdef CALC_ENERGIES
__global__ void FUNCTION_NAME(k_nbnxn, ener_2)
#else
__global__ void FUNCTION_NAME(k_nbnxn, 2)
#endif
#endif
            (const cu_atomdata_t atomdata, 
            const cu_nb_params_t nb_params, 
            const cu_nblist_t nblist,
            gmx_bool calc_fshift)
{
    /* convenience variables */
    const nbnxn_sci_t *nbl_sci  = nblist.sci;
#ifndef PRUNE_NBL
    const
#endif
    nbnxn_cj4_t *nbl_cj4        = nblist.cj4;
    const nbnxn_excl_t    *excl = nblist.excl;
    const int *atom_types       = atomdata.atom_types;
    int ntypes                  = atomdata.ntypes;
    const float4 *xq            = atomdata.xq;
    float4 *f                   = atomdata.f;
    const float3 *shift_vec     = atomdata.shift_vec;
    float rcoulomb_sq           = nb_params.rcoulomb_sq;
#ifdef EL_EWALD
    float rvdw_sq               = nb_params.rvdw_sq;
#endif
#ifdef EL_RF
    float two_k_rf              = nb_params.two_k_rf;
#endif
#ifdef EL_EWALD
    float coulomb_tab_scale     = nb_params.coulomb_tab_scale;
#endif
#ifdef PRUNE_NBL
    float rlist_sq              = nb_params.rlist_sq;
#endif

#ifdef CALC_ENERGIES
    float lj_shift  = nb_params.lj_shift;
#ifdef EL_EWALD
    float beta      = nb_params.ewald_beta;
#endif
#ifdef EL_RF
    float c_rf      = nb_params.c_rf;
#endif
    float *e_lj     = atomdata.e_lj;
    float *e_el     = atomdata.e_el;
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
        i, cii, jm, j4,
        nsubi;
    float qi, qj_f,
          r2, inv_r, inv_r2, inv_r6,
          c6, c12,
#ifdef CALC_ENERGIES
          E_lj, E_el,
#endif
          F_invr; 
    float4  xqbuf;
    float3  xi, xj, rv;
    float3  shift;
    float3  f_ij, fcj_buf, fbuf_shift;
    nbnxn_sci_t nb_sci;
    unsigned int wexcl, int_bit;
    int wexcl_idx;
    unsigned imask, imask_j;
#ifdef PRUNE_NBL
    unsigned imask_prune;
#endif

    extern __shared__ float forcebuf[]; /* j force buffer */
    float3 fci_buf[NSUBCELL];           /* i force buffer */

    nb_sci      = nbl_sci[bidx];        /* cluster index */
    sci         = nb_sci.sci;           /* i cluster index = current block index */
    cij4_start  = nb_sci.cj4_ind_start; /* first ...*/
    cij4_end    = nb_sci.cj4_ind_end;   /* and last index of j clusters */

    shift       = shift_vec[nb_sci.shift];

    for(ci_offset = 0; ci_offset < NSUBCELL; ci_offset++)
    {
        fci_buf[ci_offset] = make_float3(0.0f);
    }

#ifdef CALC_ENERGIES
    E_lj = 0.0f;
    E_el = 0.0f;

#if defined EL_EWALD || defined EL_RF
    if (nb_sci.shift == CENTRAL && nbl_cj4[cij4_start].cj[0] == sci*NSUBCELL)
    {
        /* we have the diagonal: add the charge self interaction energy term */
        for (i = 0; i < NSUBCELL; i++)
        {
            ci    = sci * NSUBCELL + i;
            ai    = ci * CLUSTER_SIZE + tidxi;  /* i atom index */
            qi    = xq[ai].w;
            E_el += qi*qi;
        }
        /* divide the self term equally over the j-threads */
        E_el /= CLUSTER_SIZE;
#ifdef EL_RF
        E_el *= -nb_params.epsfac*0.5f*c_rf;
#else
        E_el *= -nb_params.epsfac*beta*0.56418958f; /* last factor 1/sqrt(pi) */
#endif
    }
#endif
#endif

    /* skip central shifts when summing shift forces */
    if (nb_sci.shift == CENTRAL)
    {
        calc_fshift = FALSE;
    }

    fbuf_shift = make_float3(0.0f);

    /* loop over the j clusters = seen by any of the atoms in the current super-cluster */
    for (j4 = cij4_start; j4 < cij4_end; j4++)
    {
        wexcl_idx   = nbl_cj4[j4].imei[widx].excl_ind;
        imask       = nbl_cj4[j4].imei[widx].imask;
        wexcl       = excl[wexcl_idx].pair[(tidx) & (WARP_SIZE - 1)];

#ifndef PRUNE_NBL
        if (imask)
#endif
        {
#ifdef PRUNE_NBL
            imask_prune = imask;
#endif

            /* #pragma unroll 4 
               -- nvcc doesn't like my code, it refuses to unroll it */
            for (jm = 0; jm < 4; jm++)
            {
                imask_j = (imask >> (jm * 8)) & 255U;
                if (imask_j)
                {
                    nsubi = __popc(imask_j);

                    cj      = nbl_cj4[j4].cj[jm];
                    aj      = cj * CLUSTER_SIZE + tidxj;

                    /* load j atom data */
                    xqbuf   = xq[aj];
                    xj      = make_float3(xqbuf.x, xqbuf.y, xqbuf.z);
                    qj_f    = nb_params.epsfac * xqbuf.w;
                    typej   = atom_types[aj];
                    xj      -= shift;

                    fcj_buf = make_float3(0.0f);

                    /* loop over the i-clusters in sci */
                    /* #pragma unroll 8 
                       -- nvcc doesn't like my code, it refuses to unroll it
                       which is a pity because here unrolling could help.  */
                    for (cii = 0; cii < nsubi; cii++)
                    {
                        i = __ffs(imask_j) - 1;
                        imask_j &= ~(1U << i);

                        ci_offset   = i;                       /* i force buffer offset */ 
                        ci          = sci * NSUBCELL + i;      /* i cluster index */
                        ai          = ci * CLUSTER_SIZE + tidxi;  /* i atom index */

                        /* all threads load an atom from i cluster ci into shmem! */
                        xqbuf   = xq[ai];
                        xi      = make_float3(xqbuf.x, xqbuf.y, xqbuf.z);

                        /* distance between i and j atoms */
                        rv      = xi - xj;
                        r2      = norm2(rv);

#ifdef PRUNE_NBL
                        /* If _none_ of the atoms pairs are in cutoff range,
                           the bit corresponding to the current cluster-pair 
                           in imask gets set to 0.
                         */
                        if (!__any(r2 < rlist_sq))
                        {
                            imask_prune &= ~(1U << (jm * NSUBCELL + i));
                        }
#endif

                        int_bit = ((wexcl >> (jm * NSUBCELL + i)) & 1);

                        /* cutoff & exclusion check */
#if defined EL_EWALD || defined EL_RF
                        /* small r2 check to avoid invr6 overflow */
                        if (r2 < rcoulomb_sq *
                            (nb_sci.shift != CENTRAL || ci != cj || tidxj > tidxi) *
                            (r2 > 1.0e-12f))
#else
                        if (r2 < rcoulomb_sq * int_bit)
#endif
                        {
                            /* load the rest of the i-atom parameters */
                            qi      = xqbuf.w;
                            typei   = atom_types[ai];

                            /* LJ 6*C6 and 12*C12 */
                            c6      = tex1Dfetch(tex_nbfp, 2 * (ntypes * typei + typej));
                            c12     = tex1Dfetch(tex_nbfp, 2 * (ntypes * typei + typej) + 1);

                            inv_r       = rsqrt(r2);
                            inv_r2      = inv_r * inv_r;
                            inv_r6      = inv_r2 * inv_r2 * inv_r2;
#if defined EL_EWALD || defined EL_RF
                            /* We could mask inv_r2, but with Ewald
                             * masking both inv_r6 and F_invr is faster */
                            inv_r6      *= int_bit;
#endif
#ifdef EL_EWALD
                            /* this enables twin-range cut-offs (rvdw < rcoulomb <= rlist) */
                            inv_r6      *= r2 < rvdw_sq;
#endif

                            F_invr      = inv_r6 * (c12 * inv_r6 - c6) * inv_r2;

#ifdef CALC_ENERGIES
                            E_lj        += (inv_r6 + lj_shift) * (0.08333333f * c12 * (inv_r6 + lj_shift) - 0.16666667f * c6);
#endif

#ifdef EL_CUTOFF
                            F_invr      += qi * qj_f * inv_r2 * inv_r;
#endif
#ifdef EL_RF
                            F_invr      += qi * qj_f * (int_bit*inv_r2 * inv_r - two_k_rf);
#endif
#ifdef EL_EWALD
                            F_invr      += qi * qj_f * (int_bit*inv_r2 - interpolate_coulomb_force_r(r2 * inv_r, coulomb_tab_scale)) * inv_r;
#endif

#ifdef CALC_ENERGIES
#ifdef EL_CUTOFF
                            E_el        += qi * qj_f * inv_r;
#endif
#ifdef EL_RF
                            E_el        += qi * qj_f * (int_bit*inv_r + 0.5f * two_k_rf * r2 - c_rf);
#endif
#ifdef EL_EWALD
                            /* 1.0f - erff is faster than erfcf */
                            E_el        += qi * qj_f * inv_r * (int_bit - erff(r2 * inv_r * beta));
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
                    forcebuf[                 tidx] = fcj_buf.x;
                    forcebuf[    STRIDE_DIM + tidx] = fcj_buf.y;
                    forcebuf[2 * STRIDE_DIM + tidx] = fcj_buf.z;

                    /* reduce j forces */
                    reduce_force_j_generic(forcebuf, f, tidxi, tidxj, aj);
                }
            }
#ifdef PRUNE_NBL
            /* Update the imask with the new one which does not contain the 
               out of range clusters anymore. Only the first thread in the 
               warp writes. ATM this gives a minor, but consistent improvement. */
            if (tidx & (WARP_SIZE - 1))
            {
                nbl_cj4[j4].imei[widx].imask = imask_prune;
            }
#endif
        }
    }

    /* reduce i forces */
    for(ci_offset = 0; ci_offset < NSUBCELL; ci_offset++)
    {
        ai  = (sci * NSUBCELL + ci_offset) * CLUSTER_SIZE + tidxi;
        forcebuf[                 tidx] = fci_buf[ci_offset].x;
        forcebuf[    STRIDE_DIM + tidx] = fci_buf[ci_offset].y;
        forcebuf[2 * STRIDE_DIM + tidx] = fci_buf[ci_offset].z;
        __syncthreads();
        reduce_force_i(forcebuf, f,
                       &fbuf_shift, calc_fshift,
                       tidxi, tidxj, ai);
        __syncthreads();
    }

    /* add up local shift forces */
    if (calc_fshift && tidxj == 0)
    {
        atomicAdd(&atomdata.f_shift[nb_sci.shift].x, fbuf_shift.x);
        atomicAdd(&atomdata.f_shift[nb_sci.shift].y, fbuf_shift.y);
        atomicAdd(&atomdata.f_shift[nb_sci.shift].z, fbuf_shift.z);
    }

#ifdef CALC_ENERGIES
    /* flush the partial energies to shmem and sum them up */
    forcebuf[             tidx] = E_lj;
    forcebuf[STRIDE_DIM + tidx] = E_el;
    reduce_energy_pow2(forcebuf, e_lj, e_el, tidx);
#endif
}
