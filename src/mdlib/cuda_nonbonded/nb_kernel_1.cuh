/*  Launch parameterts:
    - #blocks   = #neighbor lists, blockId = neigbor_listId
    - #threads  = CELL_SIZE^2
    - shmem     = (1 + NSUBCELL) * CELL_SIZE^2 * 3 * sizeof(float)

    Each thread calculates an i force-component taking one pair of i-j atoms.
 */
/*
TODO:
  - fix GPU_FACEL
  - make utility functions "static inline __device___"
 */
#ifdef PRUNE_NBL
#ifdef CALC_ENERGIES
__global__ void FUNCTION_NAME(k_calc_nb, forces_energies_prunenbl_1)
#else
__global__ void FUNCTION_NAME(k_calc_nb, forces_prunenbl_1)
#endif
#else
#ifdef CALC_ENERGIES
__global__ void FUNCTION_NAME(k_calc_nb, forces_energies_1)
#else
__global__ void FUNCTION_NAME(k_calc_nb, forces_1)
#endif
#endif
            (const cu_atomdata_t atomdata,
            const cu_nb_params_t nb_params,
            const cu_nblist_t nblist,
            gmx_bool calc_fshift)
{
    /* convenience variables */
    const gmx_nbl_ci_t *nbl_ci  = nblist.ci;
#ifndef PRUNE_NBL
    const
#endif
    gmx_nbl_sj4_t *nbl_sj4      = nblist.sj4;
    const gmx_nbl_excl_t *excl  = nblist.excl;
    const int *atom_types       = atomdata.atom_types;
    int ntypes                  = atomdata.ntypes;
    const float4 *xq            = atomdata.xq;
    float4 *f                   = atomdata.f;
    const float3 *shift_vec     = atomdata.shift_vec;
    float cutoff_sq             = nb_params.cutoff_sq;
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
#ifdef EL_EWALD
    float beta = nb_params.ewald_beta;
#endif
#ifdef EL_RF
    float c_rf = nb_params.c_rf;
#endif
    float *e_lj = atomdata.e_lj;
    float *e_el = atomdata.e_el;
#endif

    /* thread/block/warp id-s */
    unsigned int tidxi  = threadIdx.x;
    unsigned int tidxj  = threadIdx.y;
    unsigned int tidx   = threadIdx.y * blockDim.x + threadIdx.x;
    unsigned int bidx   = blockIdx.x;
    unsigned int widx   = tidx / WARP_SIZE; /* warp index */

    int ci, si, sj, si_offset,
        ai, aj,
        cij4_start, cij4_end,
        typei, typej,
        i, sii, jm, j4,
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
    float3  f_ij, fsj_buf, fbuf_shift;
    gmx_nbl_ci_t nb_ci;
    unsigned int wexcl, excl_bit;
    int wexcl_idx;
    unsigned imask, imask_j;
#ifdef PRUNE_NBL
    unsigned imask_prune;
#endif

    extern __shared__ float forcebuf[]; /* force buffer */

    nb_ci       = nbl_ci[bidx];         /* cell index */
    ci          = nb_ci.ci;             /* i cell index = current block index */
    cij4_start  = nb_ci.sj4_ind_start;  /* first ...*/
    cij4_end    = nb_ci.sj4_ind_end;    /* and last index of j cells */

    shift       = shift_vec[nb_ci.shift];

    for(si_offset = 0; si_offset < NSUBCELL; si_offset++)
    {
        forcebuf[                 (1 + si_offset) * STRIDE_SI + tidx] = 0.0f;
        forcebuf[    STRIDE_DIM + (1 + si_offset) * STRIDE_SI + tidx] = 0.0f;
        forcebuf[2 * STRIDE_DIM + (1 + si_offset) * STRIDE_SI + tidx] = 0.0f;
    }

#ifdef CALC_ENERGIES
    E_lj = 0.0f;
    E_el = 0.0f;
#endif

    /* skip central shifts when summing shift forces */
    if (nb_ci.shift == CENTRAL)
    {
        calc_fshift = FALSE;
    }

    fbuf_shift = make_float3(0.0f);

    /* loop over the j sub-cells = seen by any of the atoms in the current cell */
    for (j4 = cij4_start; j4 < cij4_end; j4++)
    {
        wexcl_idx   = nbl_sj4[j4].imei[widx].excl_ind;
        imask       = nbl_sj4[j4].imei[widx].imask;
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

                    sj      = nbl_sj4[j4].sj[jm];
                    aj      = sj * CELL_SIZE + tidxj;

                    /* load j atom data */
                    xqbuf   = xq[aj];
                    xj      = make_float3(xqbuf.x, xqbuf.y, xqbuf.z);
                    qj_f    = GPU_FACEL * xqbuf.w;
                    typej   = atom_types[aj];
                    xj      -= shift;

                    fsj_buf = make_float3(0.0f);

                    /* loop over i sub-cells in ci */
                    /* #pragma unroll 8 
                       -- nvcc doesn't like my code,  it refuses to unroll it
                       which is a pity because here unrolling could help.  */
                    for (sii = 0; sii < nsubi; sii++)
                    {
                        i = __ffs(imask_j) - 1;
                        imask_j &= ~(1U << i);

                        si_offset   = i;                       /* i force buffer offset */ 
                        si          = ci * NSUBCELL + i;       /* i subcell index */
                        ai          = si * CELL_SIZE + tidxi;  /* i atom index */

                        /* all threads load an atom from i cell si into shmem! */
                        xqbuf   = xq[ai];
                        xi      = make_float3(xqbuf.x, xqbuf.y, xqbuf.z);

                        /* distance between i and j atoms */
                        rv      = xi - xj;
                        r2      = norm2(rv);

#ifdef PRUNE_NBL
                        /* If _none_ of the atoms pairs are in cutoff range,
                           the bit corresponding to the current sub-cell-pair 
                           in imask gets set to 0.
                         */
                        if (!__any(r2 < rlist_sq))
                        {
                            imask_prune &= ~(1U << (jm * NSUBCELL + i));
                        }
#endif

                        excl_bit = ((wexcl >> (jm * NSUBCELL + i)) & 1);

                        /* cutoff & exclusion check */
                        if (r2 < cutoff_sq * excl_bit)
                        {
                            /* load the rest of the i-atom parameters */
                            qi      = xqbuf.w;
                            typei   = atom_types[ai];

                            /* LJ C6 and C12 */
                            c6      = tex1Dfetch(tex_nbfp, 2 * (ntypes * typei + typej));
                            c12     = tex1Dfetch(tex_nbfp, 2 * (ntypes * typei + typej) + 1);

                            inv_r       = rsqrt(r2);
                            inv_r2      = inv_r * inv_r;
                            inv_r6      = inv_r2 * inv_r2 * inv_r2;

                            F_invr      = inv_r6 * (12.0f * c12 * inv_r6 - 6.0f * c6) * inv_r2;
 
#ifdef CALC_ENERGIES
                            E_lj        += inv_r6 * (c12 * inv_r6 - c6);
#endif

#ifdef EL_CUTOFF
                            F_invr      += qi * qj_f * inv_r2 * inv_r;
#endif
#ifdef EL_RF
                            F_invr      += qi * qj_f * (inv_r2 * inv_r - two_k_rf);
#endif
#ifdef EL_EWALD
                            F_invr      += qi * qj_f * interpolate_coulomb_force_r(r2 * inv_r, coulomb_tab_scale) * inv_r;
#endif

#ifdef CALC_ENERGIES
#ifdef EL_CUTOFF
                            E_el        += qi * qj_f * inv_r;
#endif
#ifdef EL_RF
                            E_el        += qi * qj_f * (inv_r + 0.5f * two_k_rf * r2 - c_rf);
#endif
#ifdef EL_EWALD
                            E_el        += qi * qj_f * inv_r * erfcf(r2 * inv_r * beta);
#endif
#endif
                            f_ij    = rv * F_invr;

                            /* accumulate j forces in registers */
                            fsj_buf -= f_ij;

                            /* accumulate i forces in shmem */
                            forcebuf[                 (1 + si_offset) * STRIDE_SI + tidx] += f_ij.x;
                            forcebuf[    STRIDE_DIM + (1 + si_offset) * STRIDE_SI + tidx] += f_ij.y;
                            forcebuf[2 * STRIDE_DIM + (1 + si_offset) * STRIDE_SI + tidx] += f_ij.z;
                        }
                    }

                    /* store j forces in shmem */
                    forcebuf[                 tidx] = fsj_buf.x;
                    forcebuf[    STRIDE_DIM + tidx] = fsj_buf.y;
                    forcebuf[2 * STRIDE_DIM + tidx] = fsj_buf.z;

                    /* reduce j forces */
                    reduce_force_j_generic(forcebuf, f, tidxi, tidxj, aj);
                }
            }
#ifdef PRUNE_NBL
            /* Update the imask with the new one which does not contain the 
               out of range sub-cells anymore. Only the first thread in the 
               warp writes. ATM this gives a minor, but consistent improvement. */
            if (tidx & (WARP_SIZE - 1))
            {
                nbl_sj4[j4].imei[widx].imask = imask_prune;
            }
#endif
        }
    }
    __syncthreads();

    /* reduce i forces */
    for(si_offset = 0; si_offset < NSUBCELL; si_offset++)
    {
        ai  = (ci * NSUBCELL + si_offset) * CELL_SIZE + tidxi;  /* i atom index */
        reduce_force_i(forcebuf + (1 + si_offset) * STRIDE_SI, f,
                       &fbuf_shift, calc_fshift,
                       tidxi, tidxj, ai);
    }

    /* add up local shift forces */
    if (calc_fshift && tidxj == 0)
    {
        atomicAdd(&atomdata.f_shift[nb_ci.shift].x, fbuf_shift.x);
        atomicAdd(&atomdata.f_shift[nb_ci.shift].y, fbuf_shift.y);
        atomicAdd(&atomdata.f_shift[nb_ci.shift].z, fbuf_shift.z);
    }

#ifdef CALC_ENERGIES
    /* flush the partial energies to shmem and sum them up */
    forcebuf[             tidx] = E_lj;
    forcebuf[STRIDE_DIM + tidx] = E_el;
    reduce_energy_pow2(forcebuf, e_lj, e_el, tidx);
#endif
}
